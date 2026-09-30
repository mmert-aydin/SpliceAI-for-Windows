"""Live SpliceAI scoring (CPU, via TensorFlow/Keras).

spliceai.utils.one_hot_encode() calls np.fromstring() in binary mode, which
was hard-removed in NumPy 2.0 (raises ValueError: "The binary mode of
fromstring is removed, use frombuffer instead"). SpliceAI 1.3.1 (the latest
release on PyPI) has not been updated for this, so it is patched here with an
equivalent frombuffer-based implementation. Verified against the known
ANKRD26 chr10:27326999 T>C value (DS_AG=0.29) before relying on it.

Also note: spliceai's own declared dependency on pysam is sidestepped
entirely -- pysam has no Windows wheel either, but it is only imported by
spliceai's own __main__.py CLI, never by spliceai.utils (Annotator /
get_delta_scores), which is all this module uses.

One more compatibility issue: spliceai/__init__.py calls signal.signal(SIGINT,
...) at import time (to make Ctrl+C exit its own CLI silently). Python only
allows signal.signal() from the main thread -- fine for the CLI, but this
module is deliberately imported lazily (see cli.py), so in the GUI the first
import of `spliceai` happens inside a background QThread and raises
"ValueError: signal only works in main thread of the main interpreter"
without the guard below. Since this project never uses spliceai's own CLI,
that handler is irrelevant here regardless of which thread sets it up.
"""
import logging
import os
import re
import signal as _signal_module

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

import numpy as np
import keras

keras.config.disable_interactive_logging()

BUILD_TO_ANNOTATION = {"hg19": "grch37", "hg38": "grch38"}
DIST_VAR = 50


def _import_spliceai_utils_safely():
    real_signal = _signal_module.signal
    _signal_module.signal = lambda *a, **k: None
    try:
        import spliceai.utils as su
    finally:
        _signal_module.signal = real_signal
    return su


def _patch_one_hot_encode():
    su = _import_spliceai_utils_safely()

    def one_hot_encode(seq):
        base_map = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        seq = seq.upper().replace("A", "\x01").replace("C", "\x02")
        seq = seq.replace("G", "\x03").replace("T", "\x04").replace("N", "\x00")
        return base_map[np.frombuffer(seq.encode("latin-1"), np.int8) % 5]

    su.one_hot_encode = one_hot_encode


_patch_one_hot_encode()

# Imported only after _patch_one_hot_encode() (above) has run.
from spliceai.utils import Annotator, get_delta_scores

logging.getLogger().setLevel(logging.ERROR)

# SpliceAI's own "Skipping record (<reason>)" warnings, as shown in the output's
# source column.
SKIP_REASONS = {
    "ref issue": "skipped: REF doesn't match reference FASTA",
    "near chromosome end": "skipped: too close to chromosome end",
    "ref too long": "skipped: REF allele too long for SpliceAI",
    "fasta issue": "skipped: reference sequence unavailable",
    "bad input": "skipped: invalid record",
}


class _SkipReasonCollector(logging.Handler):
    """Catches the reason SpliceAI gives when it declines to score a variant.
    get_delta_scores() only logs a warning and returns nothing then -- with
    logging at ERROR (above) that's invisible, and the variant would come out
    with blank scores exactly like one that simply overlaps no gene."""

    def __init__(self):
        super().__init__(logging.WARNING)
        self.reason = None

    def emit(self, record):
        m = re.match(r"Skipping record \(([^)]+)\)", record.getMessage())
        if m and self.reason is None:
            self.reason = SKIP_REASONS.get(m.group(1), f"skipped: {m.group(1)}")


class _Record:
    __slots__ = ("chrom", "pos", "ref", "alts")

    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = [alt]


class LiveScorer:
    def __init__(self, fasta_path, build, mode):
        annotation = BUILD_TO_ANNOTATION[build]
        self._annotator = Annotator(fasta_path, annotation)
        self._mask = 1 if mode == "masked" else 0

    def score(self, chrom, pos, ref, alt):
        """Returns a list of 'ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'
        annotation strings (one per overlapping gene), possibly empty."""
        return self.score_with_reason(chrom, pos, ref, alt)[0]

    def score_with_reason(self, chrom, pos, ref, alt):
        """(annotations, skip_reason): like score(), plus -- when SpliceAI
        declined to score the variant -- why (one of SKIP_REASONS' values);
        None otherwise, including for a variant that just overlaps no gene."""
        record = _Record(chrom, pos, ref, alt)
        root = logging.getLogger()
        collector = _SkipReasonCollector()
        level = root.level
        root.addHandler(collector)
        root.setLevel(logging.WARNING)
        try:
            annotations = get_delta_scores(record, self._annotator, DIST_VAR, self._mask)
        finally:
            root.removeHandler(collector)
            root.setLevel(level)
        return annotations, (None if annotations else collector.reason)
