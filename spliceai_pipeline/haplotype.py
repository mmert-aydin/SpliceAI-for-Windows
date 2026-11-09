"""The patient haplotype and its coordinate map (spec 7.1).

The haplotype is a reference region with the chosen edits applied, in
genomic orientation (predict.py turns a minus-strand gene around itself). An
edit replaces reference bases start..end (1-based, inclusive; end = start - 1
for an insertion, which goes in front of start) with what the patient has:
an inversion's bases are the reverse complement of the segment (or the
patient's own bases there), a tandem duplication is an inserted copy.

CoordMap lists the haplotype as ordered segments: haplotype start and end,
the reference position they come from, orientation (+1 as in the reference,
-1 inverted, 0 inserted: no reference counterpart) and a label (INS1, DUP1,
INV1, ...). It converts in both directions; a base inside inserted sequence
is named after it, e.g. INS1+23 (its 23rd base). The viewer's shared columns
come from it.

predict_patient() is the patient's prediction (spec 7.2): SpliceAI's output
at a base depends only on the 5,000 bases on each side of it, so only the
haplotype within 5,000 bases of an edit is predicted again (predict.py, the
gene's mask lifted onto the haplotype); everywhere else the patient's values
are the reference's, taken from the database through the coordinate map.
ref_aligned() puts them back on reference positions the way get_delta_scores
does for a VCF variant.
"""
import bisect
import time
from dataclasses import dataclass, field

import numpy as np

from .align import reverse_complement
from .predict import CONTEXT

# Model input per second (each predicted stretch plus 5,000 bases each side),
# measured with predict.py on hg38 chr21 (Phase 1, Ryzen 5 5600); for the
# estimate shown before a run.
INPUT_BASES_PER_SECOND = 53_000

_PREFIX = {"insertion": "INS", "duplication": "DUP", "inversion": "INV", "complex": "CPX",
           "SNV": "SNV", "deletion": "DEL"}


@dataclass(frozen=True)
class Edit:
    start: int
    end: int            # start - 1 for an insertion
    alt: str            # the patient's bases ("" for a deletion)
    kind: str           # SNV, insertion, deletion, duplication, inversion, complex
    label: str = ""

    @classmethod
    def from_event(cls, event):
        """An align.Event as an edit."""
        return cls(event.start, event.end, event.inserted, event.kind)


@dataclass(frozen=True)
class Segment:
    hap_start: int      # 1-based, inclusive
    hap_end: int
    ref_start: int      # orientation +1: hap_start's reference base; -1: hap_end's; 0: None
    ref_end: int
    orientation: int
    label: str = ""


class CoordMap:
    def __init__(self, segments, region_start, region_end):
        self.segments = segments
        self.region_start, self.region_end = region_start, region_end
        self.length = segments[-1].hap_end if segments else 0
        self._starts = [s.hap_start for s in segments]
        # Segments with a reference counterpart are in reference order too (edits don't overlap).
        self._placed = [s for s in segments if s.orientation]
        self._ref_starts = [s.ref_start for s in self._placed]

    def segment_at(self, h):
        if not 1 <= h <= self.length:
            raise IndexError(f"haplotype position {h} outside 1..{self.length}")
        return self.segments[bisect.bisect_right(self._starts, h) - 1]

    def to_ref(self, h):
        """(reference position or None, label) of haplotype position h: e.g.
        (1234, "") for an unchanged base, (1250, "INV1") inside an inversion,
        (1240, "SNV1") at a substituted base, (None, "INS1+23") inside an insertion."""
        s = self.segment_at(h)
        if s.orientation == 1:
            return s.ref_start + (h - s.hap_start), s.label
        if s.orientation == -1:
            return s.ref_end - (h - s.hap_start), s.label
        return None, f"{s.label}+{h - s.hap_start + 1}"

    def to_hap(self, r):
        """Haplotype position of reference position r, or None if it was
        deleted (or lies outside the region). An inverted base maps to where it
        ended up; a duplicated base to its original copy."""
        i = bisect.bisect_right(self._ref_starts, r) - 1
        if i < 0 or r > self._placed[i].ref_end:
            return None
        s = self._placed[i]
        return s.hap_start + (r - s.ref_start if s.orientation == 1 else s.ref_end - r)

    def name(self, h, chrom=""):
        """How the viewer names haplotype position h: "chr7:1,234", "chr7:1,250 (INV1)" or "INS1+23"."""
        r, label = self.to_ref(h)
        if r is None:
            return label
        where = f"{chrom + ':' if chrom else ''}{r:,}"
        return f"{where} ({label})" if label else where

    def lift(self, r, inward=1):
        """to_hap(r), or for a deleted r the nearest base kept in direction
        `inward` (+1: to the right, -1: to the left); None if there is none."""
        while self.region_start <= r <= self.region_end:
            h = self.to_hap(r)
            if h is not None:
                return h
            r += inward
        return None

    def edited_spans(self):
        """Haplotype stretches (start, end) that differ from the reference:
        every edited segment, and where bases were deleted, the two bases on
        either side of the gap."""
        spans = [(s.hap_start, s.hap_end) for s in self.segments if s.label]
        expected = self.region_start
        for i, s in enumerate(self.segments):
            if s.orientation != 1 or s.label:
                expected = None
                continue
            if expected is not None and s.ref_start != expected:
                spans.append((max(1, s.hap_start - 1), s.hap_start))
            expected = s.ref_end + 1
        last = self.segments[-1] if self.segments else None
        if last and last.orientation == 1 and not last.label and last.ref_end < self.region_end:
            spans.append((last.hap_end, last.hap_end))
        return sorted(spans)


def build(ref_seq, region_start, edits):
    """(haplotype sequence, CoordMap) of the reference region ref_seq (first
    base at region_start) with `edits` applied. Raises ValueError for edits
    that overlap or reach outside the region."""
    region_end = region_start + len(ref_seq) - 1
    edits = sorted(edits, key=lambda e: (e.start, e.end))
    counts, labelled = {}, []
    for e in edits:
        if not e.label:
            prefix = _PREFIX.get(e.kind, "EDIT")
            counts[prefix] = counts.get(prefix, 0) + 1
            e = Edit(e.start, e.end, e.alt, e.kind, f"{prefix}{counts[prefix]}")
        labelled.append(e)
    for a, b in zip(labelled, labelled[1:]):
        if b.start <= a.end or (a.start == b.start and a.end < a.start and b.end < b.start):
            raise ValueError(f"{a.label} and {b.label} overlap")
    pieces, segments = [], []
    hap = 0     # haplotype bases so far
    ref = region_start  # next reference base to copy

    def add(seq, ref_start, ref_end, orientation, label=""):
        nonlocal hap
        if seq:
            segments.append(Segment(hap + 1, hap + len(seq), ref_start, ref_end, orientation, label))
            pieces.append(seq)
            hap += len(seq)

    for e in labelled:
        if e.start < region_start or e.end > region_end or e.start > region_end + 1:
            raise ValueError(f"{e.label} ({e.start:,}..{e.end:,}) reaches outside the region "
                             f"{region_start:,}..{region_end:,}")
        if e.start > ref:
            add(ref_seq[ref - region_start:e.start - region_start], ref, e.start - 1, 1)
        n = e.end - e.start + 1
        if e.kind == "inversion":
            alt = e.alt or reverse_complement(ref_seq[e.start - region_start:e.end - region_start + 1])
            if len(alt) != n:
                raise ValueError(f"{e.label}: an inversion keeps its length ({n:,}, not {len(alt):,})")
            add(alt, e.start, e.end, -1, e.label)
        elif n > 0 and len(e.alt) == n:
            add(e.alt, e.start, e.end, 1, e.label)      # substitution: base for base
        else:
            add(e.alt, None, None, 0, e.label)           # inserted (a deletion adds nothing)
        ref = e.end + 1
    if ref <= region_end:
        add(ref_seq[ref - region_start:], ref, region_end, 1)
    return "".join(pieces), CoordMap(segments, region_start, region_end)


# Column kinds for the viewer's shared axis (columns()).
SAME, SUBSTITUTED, INSERTED, DELETED, INVERTED = range(5)


def columns(cmap, start, end):
    """The viewer's shared horizontal axis (spec 8.2), like a multiple
    alignment of the reference and the haplotype over reference positions
    start..end: (ref, hap, kind) arrays, one entry per column. ref / hap are
    the 1-based positions (0 where that side has a gap), kind one of SAME,
    SUBSTITUTED, INSERTED (a gap in the reference rows), DELETED (a gap in the
    patient rows) or INVERTED (reference order; the patient's base is the
    complement, read the other way). Every reference position in start..end
    is one column, in order; inserted bases sit right after the reference
    base they follow."""
    ref, hap, kind = [], [], []
    expected = cmap.region_start

    def deleted(upto):
        if upto > expected:
            ref.append(np.arange(expected, upto))
            hap.append(np.zeros(upto - expected, np.int64))
            kind.append(np.full(upto - expected, DELETED, np.int8))

    for s in cmap.segments:
        n = s.hap_end - s.hap_start + 1
        if s.orientation == 0:
            ref.append(np.zeros(n, np.int64))
            hap.append(np.arange(s.hap_start, s.hap_end + 1))
            kind.append(np.full(n, INSERTED, np.int8))
            continue
        deleted(s.ref_start)
        ref.append(np.arange(s.ref_start, s.ref_end + 1))
        if s.orientation == 1:
            hap.append(np.arange(s.hap_start, s.hap_end + 1))
            kind.append(np.full(n, SUBSTITUTED if s.label else SAME, np.int8))
        else:
            hap.append(np.arange(s.hap_end, s.hap_start - 1, -1))
            kind.append(np.full(n, INVERTED, np.int8))
        expected = s.ref_end + 1
    deleted(cmap.region_end + 1)
    ref, hap, kind = (np.concatenate(x) if x else np.zeros(0, np.int64) for x in (ref, hap, kind))
    # Keep start..end, and inserted columns between the first and last kept reference column.
    placed = np.nonzero((ref >= start) & (ref <= end))[0]
    if not len(placed):
        return ref[:0], hap[:0], kind[:0]
    keep = slice(placed[0], placed[-1] + 1)
    return ref[keep], hap[keep], kind[keep].astype(np.int8)


def apply_events(ref_seq, region_start, events):
    """build() with the align.Event objects chosen (apply=True)."""
    return build(ref_seq, region_start, [Edit.from_event(e) for e in events if e.apply])


# --- the patient's prediction (spec 7.2) ---------------------------------------------

@dataclass
class PatientPrediction:
    """SpliceAI's acceptor and donor probabilities along the patient
    haplotype, over the gene's stored region lifted onto it."""
    start: int                  # haplotype positions start..end (1-based, inclusive)
    end: int
    values: np.ndarray          # float32 (end - start + 1, 2): acceptor, donor
    mask: tuple                 # the gene's mask, lifted (haplotype positions)
    windows: list = field(default_factory=list)     # (start, end) stretches predicted again
    seconds: float = 0.0

    def at(self, h):
        return self.values[h - self.start]


def patient_region(cmap, plan):
    """(start, end, mask) on the haplotype for the gene's stored region and mask (plan_gene())."""
    start, end = cmap.lift(plan["start"], 1), cmap.lift(plan["end"], -1)
    a, b = cmap.lift(plan["mask"][0], 1), cmap.lift(plan["mask"][1], -1)
    if None in (start, end, a, b) or start > end:
        raise ValueError("the chosen edits delete the whole gene")
    return start, end, (min(a, b), max(a, b))


def recompute_windows(cmap, start, end, context=CONTEXT):
    """The haplotype stretches within start..end whose prediction can differ
    from the reference's: within `context` bases of an edit, merged."""
    windows = []
    for a, b in cmap.edited_spans():
        a, b = max(start, a - context), min(end, b + context)
        if a > b:
            continue
        if windows and a <= windows[-1][1] + 1:
            windows[-1] = (windows[-1][0], max(windows[-1][1], b))
        else:
            windows.append((a, b))
    return windows


def estimate_seconds(windows):
    """Rough time to predict these windows (each needs its 5,000 bases of context on both sides)."""
    return sum(b - a + 1 + 2 * CONTEXT for a, b in windows) / INPUT_BASES_PER_SECOND


def predict_patient(predictor, hap, cmap, plan, wt_values, checkpoint=None, on_progress=None):
    """The patient's prediction for one gene: wt_values are the reference's
    (plan["start"]..plan["end"], as refdb.ReferenceScores.gene() gives them),
    predictor a predict.Predictor of that genome build. Only recompute_windows()
    is predicted; checkpoint() is called between model calls (pause, cancel),
    on_progress(bases done, bases to do) after each window."""
    t0 = time.time()
    start, end, mask = patient_region(cmap, plan)
    values = np.empty((end - start + 1, 2), np.float32)
    # Unchanged bases: the reference's values, through the coordinate map.
    for s in cmap.segments:
        a, b = max(s.hap_start, start), min(s.hap_end, end)
        if a > b or s.orientation != 1 or s.label:
            continue
        r = s.ref_start + (a - s.hap_start)
        values[a - start:b - start + 1] = wt_values[r - plan["start"]:r - plan["start"] + b - a + 1]
    windows = recompute_windows(cmap, start, end)
    total, done = sum(b - a + 1 for a, b in windows), 0
    for a, b in windows:
        values[a - start:b - start + 1] = predictor.predict(lambda x, y: hap[x - 1:y], len(hap), plan["strand"],
                                                            mask[0], mask[1], a, b, checkpoint=checkpoint)
        done += b - a + 1
        if on_progress:
            on_progress(done, total)
    return PatientPrediction(start, end, values, mask, windows, time.time() - t0)


def ref_aligned(prediction, cmap, start, end):
    """The patient's values at reference positions start..end, the way
    get_delta_scores lines them up for a VCF variant: a deleted base gets 0,
    and bases inserted after a reference base are pooled into it (the
    maximum), as an insertion's are into its anchor base."""
    out = np.zeros((end - start + 1, 2), np.float32)
    for s in cmap.segments:
        if s.orientation == 0:
            continue
        r0, r1 = max(s.ref_start, start), min(s.ref_end, end)
        if r0 > r1:
            continue
        if s.orientation == 1:
            # Haplotype positions of r0..r1, clipped to the prediction.
            h0, h1 = s.hap_start + (r0 - s.ref_start), s.hap_start + (r1 - s.ref_start)
            lo, hi = max(h0, prediction.start), min(h1, prediction.end)
            if lo <= hi:
                out[r0 + (lo - h0) - start:r0 + (hi - h0) - start + 1] =                     prediction.values[lo - prediction.start:hi - prediction.start + 1]
        else:
            for r in range(r0, r1 + 1):
                h = s.hap_start + (s.ref_end - r)
                if prediction.start <= h <= prediction.end:
                    out[r - start] = prediction.at(h)
    for s in cmap.segments:
        if s.orientation != 0 or s.hap_start == 1:
            continue
        anchor, _ = cmap.to_ref(s.hap_start - 1)
        if anchor is None or not start <= anchor <= end:
            continue
        a, b = max(s.hap_start, prediction.start), min(s.hap_end, prediction.end)
        if a <= b:
            pooled = prediction.values[a - prediction.start:b - prediction.start + 1].max(axis=0)
            out[anchor - start] = np.maximum(out[anchor - start], pooled)
    return out
