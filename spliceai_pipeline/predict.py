"""The one SpliceAI prediction path shared by the reference score database
builder, its live fallback and the Sequence-to-cDNA tab (spec 5.1).

SpliceAI's own get_delta_scores() (spliceai/utils.py) predicts a 10,101-base
window around each variant: the five models' mean, with every base outside the
gene's interval replaced by N, the input reverse-complemented for a minus-strand
gene and the output flipped back to genomic order. predict() computes exactly
that for any stretch of sequence at once. The models take inputs of any length,
their convolutions reach exactly 5,000 bases to each side and their output is
cropped by 5,000 per side, so the prediction for a base depends only on the
input from 5,000 bases before it to 5,000 after it -- a gene predicted in long
chunks gives what get_delta_scores computes base by base (tests/test_predict.py
checks this against get_delta_scores itself).

The VCF tab keeps calling get_delta_scores unchanged; this module is for
everything that needs whole stretches of per-base probabilities.

The models are called through one tf.function each, traced once for any input
length. Keras' predict() would retrace for every new length -- every gene has
its own -- and is 2.5x slower here; the results differ from predict()'s by
float rounding only (7e-9 measured).

TensorFlow is only imported when a Predictor is made, so the constants and
verify_models() are cheap to import (refdb.py reads CONTEXT).
"""
import hashlib
import os

import numpy as np

from .transcripts import spliceai_package_dir

# Each output base depends on the input this far to either side.
CONTEXT = 5000

# Output bases per model call. Longer chunks spend less on the 2 x 5,000 bases
# of context each one needs; this one is ~10% context and uses ~0.5 GB of RAM.
DEFAULT_CHUNK = 100_000

# SHA-256 of the five model files in the pinned SpliceAI 1.3.1 wheel
# (spliceai_gui.spliceai_setup.SPLICEAI_WHEEL_SHA256).
MODEL_SHA256 = {
    "spliceai1.h5": "e1fd5adcef7489d604b10e79c40078ef790d51ef048c4ce3869c9119ac5de42b",
    "spliceai2.h5": "6ab042b82ab966b6d3582cb31b96f0859ea08a864f168d69e83aa14450a3b66e",
    "spliceai3.h5": "e2e790bde53dfdf410c6dc434a86122a7d12f3f38dc2ef45d85986e9ecf22fad",
    "spliceai4.h5": "ca88ac9e58e69ba6fdeed319b72f063f164c9abf7392eaccef903e94c1d99dd6",
    "spliceai5.h5": "791cd22c62a80a08d2ca674615a93ce8159d7b55bd157cfef2983b1bd6b41391",
}


class ModelMismatchError(RuntimeError):
    """The installed SpliceAI model files aren't the pinned ones."""


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_models():
    """{file name: SHA-256} of the installed models; raises ModelMismatchError
    unless all five are the pinned files."""
    folder = os.path.join(spliceai_package_dir(), "models")
    found = {}
    for name, expected in MODEL_SHA256.items():
        path = os.path.join(folder, name)
        found[name] = _sha256(path) if os.path.exists(path) else "missing"
        if found[name] != expected:
            raise ModelMismatchError(
                f"SpliceAI model {name} is not the published SpliceAI 1.3.1 file "
                f"(SHA-256 {found[name]}, expected {expected}). Reinstall SpliceAI."
            )
    return found


class Predictor:
    """The five models, the reference FASTA and SpliceAI's gene table, loaded
    once (through SpliceAI's own Annotator)."""

    def __init__(self, fasta_path, build):
        import tensorflow as tf

        from .score import BUILD_TO_ANNOTATION, Annotator, one_hot_encode

        self._one_hot = one_hot_encode
        self.annotator = Annotator(fasta_path, BUILD_TO_ANNOTATION[build])
        spec = [tf.TensorSpec([None, None, 4], tf.float32)]
        self._calls = [tf.function(lambda x, m=m: m(x, training=False), input_signature=spec)
                       for m in self.annotator.models]

    def _mean(self, x):
        import tensorflow as tf

        x = tf.constant(x, tf.float32)
        return np.mean([call(x).numpy() for call in self._calls], axis=0)

    def predict(self, fetch, seq_len, strand, mask_start, mask_end, start, end, chunk=DEFAULT_CHUNK,
                checkpoint=None):
        """Acceptor and donor probabilities, float32 array of shape (n, 2), for
        the 1-based positions start..end (inclusive) of a sequence, in genomic
        order.

        fetch(a, b) returns the bases a..b (1-based, inclusive) of a sequence
        seq_len long -- the reference chromosome, or later a patient haplotype.
        Input bases outside mask_start..mask_end, and beyond either end of the
        sequence, are N, as in get_delta_scores. checkpoint(), if given, is
        called between chunks (pause / cancel)."""
        out = np.empty((end - start + 1, 2), np.float32)
        for s in range(start, end + 1, chunk):
            e = min(s + chunk - 1, end)
            x = self._one_hot(masked_bases(fetch, seq_len, mask_start, mask_end, s - CONTEXT, e + CONTEXT))[None, :]
            if strand == "-":
                x = x[:, ::-1, ::-1]
            y = self._mean(np.ascontiguousarray(x))
            if strand == "-":
                y = y[:, ::-1]
            out[s - start:e - start + 1] = y[0, :, 1:3]
            if checkpoint is not None:
                checkpoint()
        return out


def masked_bases(fetch, seq_len, mask_start, mask_end, a, b):
    """Bases a..b (1-based, inclusive; may run past either end of the sequence)
    with N wherever they fall outside the mask or outside the sequence."""
    lo, hi = max(a, mask_start, 1), min(b, mask_end, seq_len)
    if lo > hi:
        return "N" * (b - a + 1)
    return "N" * (lo - a) + fetch(lo, hi) + "N" * (b - hi)
