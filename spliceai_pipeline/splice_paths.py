"""Transcript reconstruction from SpliceAI's per-base probabilities (spec
7.3); the same algorithm for the wild type and the patient.

A transcript runs from its start anchor to its end anchor (from the canonical
transcript, lifted onto the sequence; nothing else is taken from the
annotation), alternating a donor (the last base of an exon) and an acceptor
(the first base of the next one) on the transcript's strand -- where
SpliceAI puts them, checked against the canonical exons of genes on both
strands. Candidate sites are bases whose acceptor or donor probability is at
least the site threshold.

Every candidate is a use / no-use event: a path scores the sum of log(p)
over the sites it uses plus the sum of log(1 - p) over the candidates it
doesn't, with p clipped away from 0 and 1. Scoring the used sites alone
would always favour skipping exons; this way skipping a strong exon or
keeping an intron pays for every strong site it ignores. Since the score is
a constant plus the sum of logit(p) over the used sites, an exact k-best
dynamic programme over the candidates in transcript order finds the top
paths in O(n k log k). The relative score is a softmax across the paths
reported: a ranking, not an isoform frequency.
"""
import heapq
import math
from dataclasses import dataclass, field

import numpy as np

SITE_THRESHOLD = 0.05
MIN_EXON = 3        # internal exons (spec 17: 3 nt keeps micro-exons)
MIN_INTRON = 60
TOP = 5
# p is clipped to [CLIP, 1 - CLIP]: float16 storage rounds p >= 0.99976 to 1.0.
CLIP = 1e-4


@dataclass
class Path:
    exons: list                 # [(start, end)] genomic, inclusive, in transcript order
    score: float                # log-likelihood as above
    sites: list = field(default_factory=list)   # [(kind, position, p)] used, in transcript order
    relative: float = 0.0       # softmax across the paths reported

    @property
    def introns(self):
        return len(self.exons) - 1


def _logit(p):
    p = min(max(p, CLIP), 1 - CLIP)
    return math.log(p) - math.log(1 - p)


def top_paths(acceptor, donor, offset, strand, tx_start, tx_end, k=TOP, threshold=SITE_THRESHOLD,
              min_exon=MIN_EXON, min_intron=MIN_INTRON):
    """The k best transcripts. acceptor, donor: probabilities at positions
    offset, offset + 1, ... (genomic order, as predict.py gives them);
    tx_start, tx_end: the transcript's first and last base (genomic; on the
    minus strand tx_start is the larger). Returns [Path], best first."""
    sign = 1 if strand == "+" else -1
    length = (tx_end - tx_start) * sign      # transcript coordinate t = (x - tx_start) * sign, 0..length
    if length < 0:
        raise ValueError("the transcript's end anchor lies before its start anchor")
    to_index = lambda t: tx_start + sign * t - offset
    lo, hi = sorted((to_index(0), to_index(length)))
    if lo < 0 or hi >= len(acceptor):
        raise ValueError("the transcript's anchors lie outside the predicted stretch")
    acc = np.asarray(acceptor[lo:hi + 1], np.float64)
    don = np.asarray(donor[lo:hi + 1], np.float64)
    if sign < 0:
        acc, don = acc[::-1], don[::-1]
    # Candidates, in transcript order: (t, kind, p).
    cands = sorted([(int(t), "D", float(don[t])) for t in np.nonzero(don >= threshold)[0] if t < length]
                   + [(int(t), "A", float(acc[t])) for t in np.nonzero(acc >= threshold)[0] if t > 0])
    constant = sum(math.log(1 - min(max(p, CLIP), 1 - CLIP)) for _, _, p in cands)

    # Per node, its k best partial paths: (score, predecessor node, predecessor's rank).
    START = -1
    best = {START: [(0.0, None, None)]}
    donors_ready, acceptors_ready = [], []          # running top-k over predecessors now in reach
    donors_sorted = [(t, i) for i, (t, kind, _) in enumerate(cands) if kind == "D"]
    acceptors_sorted = [(t, i) for i, (t, kind, _) in enumerate(cands) if kind == "A"]
    di = ai = 0

    def merge(running, node):
        entries = [(s, node, rank) for rank, (s, _, _) in enumerate(best[node])]
        return heapq.nlargest(k, running + entries, key=lambda e: e[0])

    acceptor_pool = []
    for i, (t, kind, p) in enumerate(cands):
        w = _logit(p)
        if kind == "A":
            while di < len(donors_sorted) and donors_sorted[di][0] <= t - min_intron - 1:
                donors_ready = merge(donors_ready, donors_sorted[di][1])
                di += 1
            best[i] = [(s + w, node, rank) for s, node, rank in donors_ready]
        else:
            while ai < len(acceptors_sorted) and acceptors_sorted[ai][0] <= t - min_exon + 1:
                acceptor_pool = merge(acceptor_pool, acceptors_sorted[ai][1])
                ai += 1
            pool = heapq.nlargest(k, acceptor_pool + [(0.0, START, 0)], key=lambda e: e[0])
            best[i] = [(s + w, node, rank) for s, node, rank in pool]
    # The end anchor: the last exon from any acceptor, or one exon from start to end.
    ends = [(0.0, START, 0)] + [(s, i, rank) for _, i in acceptors_sorted for rank, (s, _, _) in enumerate(best[i])]
    ends = heapq.nlargest(k, ends, key=lambda e: e[0])

    paths = []
    for score, node, rank in ends:
        sites = []
        while node != START:
            t, kind, p = cands[node]
            sites.append((kind, t, p))
            _, node, rank = best[node][rank]
        sites.reverse()
        bounds = [0] + [t for _, t, _ in sites] + [length]
        exons_t = list(zip(bounds[0::2], bounds[1::2]))
        exons = [tuple(sorted((tx_start + sign * a, tx_start + sign * b))) for a, b in exons_t]
        paths.append(Path(exons, constant + score,
                          [("acceptor" if kind == "A" else "donor", tx_start + sign * t, p) for kind, t, p in sites]))
    top = max(p.score for p in paths)
    weights = [math.exp(p.score - top) for p in paths]
    for p, w in zip(paths, weights):
        p.relative = w / sum(weights)
    return paths


def path_score(exons, acceptor, donor, offset, strand, threshold=SITE_THRESHOLD):
    """The score top_paths() gives a path with these exons (for checks): its
    sites must be candidates; the remaining candidates count as unused."""
    tx_start = exons[0][0] if strand == "+" else exons[0][1]
    tx_end = exons[-1][1] if strand == "+" else exons[-1][0]
    used = set()
    for n, (a, b) in enumerate(exons):
        first, last = (a, b) if strand == "+" else (b, a)
        if n:
            used.add(("A", first))
        if n < len(exons) - 1:
            used.add(("D", last))
    lo, hi = sorted((tx_start, tx_end))
    score = 0.0
    for x in range(lo, hi + 1):
        for kind, arr in (("A", acceptor), ("D", donor)):
            p = float(arr[x - offset])
            at_edge = (kind == "A" and x == tx_start) or (kind == "D" and x == tx_end)
            if p < threshold or at_edge:
                if (kind, x) in used:
                    raise ValueError(f"{kind} at {x} isn't a candidate site")
                continue
            pc = min(max(p, CLIP), 1 - CLIP)
            score += math.log(pc) if (kind, x) in used else math.log(1 - pc)
    return score
