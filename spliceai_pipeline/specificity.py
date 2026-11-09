"""Sequence-specificity control (spec 14 C3): is the predicted effect of an
event specific to its sequence, or would any sequence of that kind there do
the same?

For each applied event, the local prediction is run again with N random
replacements, the other events kept as they are:
- an insertion, a duplication, an inversion or a complex change: a random
  sequence of the same length and the same dinucleotide composition
  (Altschul & Erickson 1985), at the same position;
- a deletion: a deletion of the same size at a random position in the same
  intron (inside the transcript, if it touches an exon).
Substitutions of a base or two have no such replacements and are left out.

Each event is measured against the sample without it (the other events
kept; the wild type when it is the only one). Reported per event: how often
the top-1 transcript differs from that one, and the maximum delta (max
|difference| per base within 5,000 bases of the change, lined up as
get_delta_scores lines a variant up), next to the observed values. An
empirical comparison, not a p-value.
"""
import random
from dataclasses import dataclass, field

import numpy as np

from . import haplotype, splice_paths
from .cdna import to_reference_exons
from .predict import CONTEXT
from .texts import tx

DEFAULT_N = 20


def dinucleotide_shuffle(seq, rng):
    """A random sequence with exactly seq's dinucleotide counts, first and
    last base: a random Eulerian path through seq's dinucleotide graph whose
    last exits form a random tree towards the last base."""
    if len(seq) < 3:
        return seq
    edges = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    last = seq[-1]
    while True:
        # Each base but the last leaves for the last time along a random edge;
        # those last exits must lead every base to `last` (no cycle).
        exit_of = {v: rng.choice(out) for v, out in edges.items() if v != last}
        ok = True
        for v in exit_of:
            seen, x = set(), v
            while x != last:
                if x in seen:
                    ok = False
                    break
                seen.add(x)
                x = exit_of[x]
            if not ok:
                break
        if ok:
            break
    order = {}
    for v, out in edges.items():
        rest = list(out)
        if v in exit_of:
            rest.remove(exit_of[v])
        rng.shuffle(rest)
        order[v] = rest + ([exit_of[v]] if v in exit_of else [])
    out, x = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        x = order[x].pop(0)
        out.append(x)
    return "".join(out)


def counts(seq):
    return sorted(a + b for a, b in zip(seq, seq[1:]))


@dataclass
class EventControl:
    label: str                  # the event, as the differences table names it
    kind: str
    observed_changed: bool      # the observed top-1 differs from the WT top-1
    observed_delta: float
    changed: int = 0            # of the replicates
    deltas: list = field(default_factory=list)
    note: str = ""

    @property
    def n(self):
        return len(self.deltas)

    def text(self):
        if self.note:
            return tx("{label}: {note}", label=self.label, note=self.note)
        median = float(np.median(self.deltas)) if self.deltas else 0.0
        return tx("{label}: observed: top-1 {top}, max delta {delta:.2f}; random: {changed}/{n} changed, median max "
                  "delta {median:.2f} (range {low:.2f}-{high:.2f})", label=self.label,
                  top=tx("changed") if self.observed_changed else tx("unchanged"), delta=self.observed_delta,
                  changed=self.changed, n=self.n, median=median, low=min(self.deltas), high=max(self.deltas))


def _top(result, values, start, cmap):
    from .cdna_analysis import anchors

    c, s = result.gene.canonical, result.settings
    a = anchors(result.gene)
    plus = c.strand == "+"
    tss, tes = cmap.lift(a["transcript start"], 1 if plus else -1), cmap.lift(a["transcript end"], -1 if plus else 1)
    path = splice_paths.top_paths(values[:, 0], values[:, 1], start, c.strand, tss, tes, k=1, threshold=s.threshold,
                                  min_exon=s.min_exon, min_intron=s.min_intron)[0]
    return _key(to_reference_exons(path.exons, cmap))


def _key(exons):
    return [tuple(x) if x else None for x in exons]


def aligned(result, pred, cmap):
    """A prediction on the reference's positions (the gene's stored region),
    lined up as get_delta_scores lines a variant up."""
    return haplotype.ref_aligned(pred, cmap, result.plan["start"], result.plan["end"])


def max_delta(result, values, base, edit):
    """max |values - base| within CONTEXT bases of `edit` (reference positions, the gene's stored region)."""
    p0 = result.plan["start"]
    lo = max(min(edit.start, edit.end) - CONTEXT - p0, 0)
    hi = min(max(edit.start, edit.end) + CONTEXT - p0 + 1, len(base))
    return float(np.max(np.abs(values[lo:hi] - base[lo:hi]))) if hi > lo else 0.0


def _intron(result, e):
    """(lo, hi) of the canonical intron holding edit e, else the transcript's span."""
    exons = sorted(result.gene.canonical.exons)
    for (a, b), (c, d) in zip(exons, exons[1:]):
        if b < e.start and e.end < c:
            return b + 1, c - 1
    return exons[0][0], exons[-1][1]


def replacement(result, e, rng):
    """A random stand-in for edit e, or None if it has none."""
    ref = result.compat.ref_seq
    r0 = result.compat.region[0]
    n = e.end - e.start + 1
    if e.kind == "deletion":
        lo, hi = _intron(result, e)
        if hi - lo + 1 <= n:
            return None
        start = rng.randint(lo, hi - n + 1)
        return haplotype.Edit(start, start + n - 1, "", "deletion")
    if e.kind == "inversion":
        seq = e.alt or haplotype.reverse_complement(ref[e.start - r0:e.end - r0 + 1])
        return haplotype.Edit(e.start, e.end, dinucleotide_shuffle(seq, rng), "complex")
    if len(e.alt) >= 3:                 # insertion, duplication, complex
        return haplotype.Edit(e.start, e.end, dinucleotide_shuffle(e.alt, rng), e.kind)
    return None


def run(result, scores, n=DEFAULT_N, seed=1, checkpoint=None, on_progress=None):
    """[EventControl] for result's applied events (cdna_analysis.Analysis);
    scores: the build's refdb.ReferenceScores (its predictor). Each event is
    measured against the sample without it (the other events kept; the wild
    type when it is the only one): whether the top-1 transcript changes, and
    the maximum delta within CONTEXT bases of it."""
    rng = random.Random(seed)
    by = {(t.kind, t.rank): t for t in result.transcripts}
    region0, ref_seq = result.compat.region[0], result.compat.ref_seq

    def predict(edits):
        hap, cmap = haplotype.build(ref_seq, region0, edits)
        pred = haplotype.predict_patient(scores.predictor(), hap, cmap, result.plan, result.wt_values,
                                         checkpoint=checkpoint)
        return _top(result, pred.values, pred.start, cmap), aligned(result, pred, cmap)

    observed_top = _key(by[("Patient", 1)].ref_exons)
    observed = aligned(result, result.patient, result.cmap)
    out = []
    jobs = len(result.edits) * n
    done = 0
    for i, e in enumerate(result.edits):
        size = max(e.end - e.start + 1, len(e.alt))
        label = tx("The {kind} at {pos:,} ({size:,} bp)", kind=e.kind, pos=e.start, size=size)
        others = [x for j, x in enumerate(result.edits) if j != i]
        if replacement(result, e, random.Random(0)) is None:
            out.append(EventControl(label, e.kind, False, 0.0,
                                    note=tx("no random replacement for this kind or size of change; left out")))
            done += n
            continue
        if others:
            base_top, base = predict(others)
        else:
            base_top, base = _key(by[("WT predicted", 1)].exons), result.wt_values
        ctl = EventControl(label, e.kind, observed_top != base_top, max_delta(result, observed, base, e))
        out.append(ctl)
        for _ in range(n):
            if checkpoint:
                checkpoint()
            stand_in = replacement(result, e, rng)
            try:
                top, values = predict(others + [stand_in])
            except ValueError:          # a random deletion overlapping another edit: left out
                done += 1
                continue
            ctl.changed += top != base_top
            ctl.deltas.append(max_delta(result, values, base, stand_in))
            done += 1
            if on_progress:
                on_progress(done, jobs)
        if not ctl.deltas:              # every stand-in overlapped another edit: nothing to compare with
            ctl.note = tx("every random replacement overlapped another change; left out")
    return out


# --- breakpoint sensitivity (spec 14 C6) ----------------------------------------------------

SHIFTS = (-50, -10, 10, 50)


@dataclass
class Breakpoints:
    """The breakpoint sensitivity of one event: whether the sample's top-1
    transcript stays the same with either breakpoint shifted."""
    label: str
    tried: list = field(default_factory=list)       # [side ("start" / "end"), shift, changed]
    note: str = ""

    def text(self):
        if self.note:
            return tx("{label}: {note}", label=self.label, note=self.note)
        changed = [f"{side} {shift:+d}" for side, shift, c in self.tried if c]
        same = len(self.tried) - len(changed)
        return (tx("{label}: the top transcript stays the same for {same} of {n} shifted breakpoints; it changes "
                   "with {changed} bp", label=self.label, same=same, n=len(self.tried), changed=", ".join(changed))
                if changed else
                tx("{label}: the top transcript stays the same for {same} of {n} shifted breakpoints",
                   label=self.label, same=same, n=len(self.tried)))


def shifted(e, side, shift, ref_seq, region0):
    """Edit e with its start or end breakpoint moved by `shift` bases; None if
    that leaves nothing or the edit has no such breakpoint."""
    r = lambda a, b: ref_seq[a - region0:b - region0 + 1]
    if e.end < e.start:                         # inserted before e.start
        if e.kind == "duplication":
            a, b = e.start - len(e.alt), e.start - 1
            direct = e.alt == r(a, b)
            if direct or e.alt == haplotype.reverse_complement(r(a, b)):
                a, b = (a + shift, b) if side == "start" else (a, b + shift)
                if a > b:
                    return None
                copy = r(a, b)
                return haplotype.Edit(b + 1, b, copy if direct else haplotype.reverse_complement(copy), e.kind)
        if side == "end":
            return None                         # a plain insertion has one breakpoint
        return haplotype.Edit(e.start + shift, e.start + shift - 1, e.alt, e.kind)
    if e.start == e.end and len(e.alt) == 1:
        return None                             # a substitution
    a, b = (e.start + shift, e.end) if side == "start" else (e.start, e.end + shift)
    if a > b:
        return None
    if e.kind == "inversion":
        return haplotype.Edit(a, b, "", "inversion")
    return haplotype.Edit(a, b, "" if e.kind == "deletion" else e.alt, e.kind)


def breakpoint_sensitivity(result, scores, shifts=SHIFTS, checkpoint=None, on_progress=None):
    """[Breakpoints] for result's applied events: each breakpoint moved by
    each of `shifts` bases (the other events kept), the top-1 transcript
    compared with the sample's own."""
    by = {(t.kind, t.rank): t for t in result.transcripts}
    region0, ref_seq = result.compat.region[0], result.compat.ref_seq
    observed = _key(by[("Patient", 1)].ref_exons)
    jobs = [(i, side, s) for i in range(len(result.edits)) for side in ("start", "end") for s in shifts]
    out = [Breakpoints(f"The {e.kind} at {e.start:,} ({max(e.end - e.start + 1, len(e.alt)):,} bp)")
           for e in result.edits]
    for done, (i, side, s) in enumerate(jobs, 1):
        if checkpoint:
            checkpoint()
        stand_in = shifted(result.edits[i], side, s, ref_seq, region0)
        if stand_in is not None:
            edits = [stand_in if j == i else x for j, x in enumerate(result.edits)]
            try:
                hap, cmap = haplotype.build(ref_seq, region0, edits)
            except ValueError:                  # shifted out of the region or onto another event
                stand_in = None
        if stand_in is not None:
            pred = haplotype.predict_patient(scores.predictor(), hap, cmap, result.plan, result.wt_values,
                                             checkpoint=checkpoint)
            out[i].tried.append([side, s, _top(result, pred.values, pred.start, cmap) != observed])
        if on_progress:
            on_progress(done, len(jobs))
    for b in out:
        if not b.tried:
            b.note = tx("no breakpoint to shift (a substitution)")
    return out
