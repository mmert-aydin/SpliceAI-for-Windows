"""Aligning a pasted patient sequence to the reference (spec 6.1, 6.2).

parse_sequence() reads the paste (raw or FASTA). align() places it on a
reference region and lists where it differs:

- 500 bases and longer: 21-mers that occur once in the region and once in the
  paste are seeds, in both orientations. Runs of seeds on one diagonal become
  blocks, extended base by base while the bases still match, and the forward
  blocks are chained by collinearity (an exact best chain in O(n log n)).
  Between two chained blocks, a gap where one side is empty is a plain
  insertion or deletion; other gaps are aligned with parasail (global,
  affine gaps). Reverse blocks inside a gap are an inversion. The paste's two
  ends are aligned semi-globally; an end that doesn't match (under 90%
  identity) is left out and reported.
- 50 to 499 bases: exact dynamic programming against the whole region (the
  paste aligned end to end, the region's ends free), in both orientations.

Indels are left-aligned with normalize.normalize_variant, and an insertion
equal to the reference right after it is a tandem duplication. Coordinates
are 1-based, inclusive; an event's pos/ref/alt are VCF-style (indels keep
their anchor base), its start/end the reference bases it replaces (start =
end + 1 for an insertion: nothing replaced, inserted before start).
compatibility() is the "Check compatibility" report, validate_events() the
check of the events chosen to build the patient haplotype (haplotype.py).
"""
import re
from dataclasses import dataclass, field

import numpy as np

from .normalize import normalize_variant
from .texts import tx

K = 21
MIN_LENGTH = 50
SEED_FROM = 500
# parasail's DNA matrix (match +5, mismatch -4) with affine gaps: a long
# insertion or deletion is one gap, not many.
GAP_OPEN, GAP_EXTEND = 10, 1
# Gaps with both sides longer than this many cells aren't aligned base by
# base: they become one complex replacement (a rearrangement, not a variant).
MAX_GAP_CELLS = 40_000_000
# A gap whose base-by-base alignment matches fewer of its aligned columns than
# this is one replacement (a delins or a rearrangement), not scattered changes.
LOW_IDENTITY = 0.75
END_MIN_IDENTITY = 0.90
APPLY_MIN_SIZE = 3

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_CODE = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def parse_sequence(text):
    """The bases of a pasted sequence: raw or one FASTA record; case,
    whitespace and digits (position numbers) don't matter. Raises ValueError,
    worded for the user, for other characters, several records, or fewer
    than MIN_LENGTH bases."""
    lines = text.splitlines()
    headers = [i for i, line in enumerate(lines) if line.lstrip().startswith(">")]
    if len(headers) > 1:
        raise ValueError(tx("The text holds {n} FASTA sequences; paste one.", n=len(headers)))
    if headers and any(line.strip() for line in lines[:headers[0]]):
        raise ValueError(tx("There is text before the FASTA header line (>...)."))
    body = "".join(lines[headers[0] + 1:] if headers else lines)
    seq = re.sub(r"[\s\d]+", "", body).upper()
    bad = sorted(set(seq) - set("ACGTN"))
    if bad:
        shown = ", ".join(repr(c) for c in bad[:8])
        raise ValueError(tx("The sequence contains characters other than A, C, G, T and N: {shown}.", shown=shown))
    if len(seq) < MIN_LENGTH:
        raise ValueError(tx("The sequence has {n} bases; at least {least} are needed to place it reliably on the "
                            "reference.", n=len(seq), least=MIN_LENGTH))
    return seq


def parse_region(text):
    """"chr7:117,480,025-117,668,665" (commas, spaces and "chr" optional) ->
    ("7", 117480025, 117668665). Raises ValueError."""
    m = re.fullmatch(r"\s*(?:chr)?([0-9]{1,2}|X|Y|M|MT)\s*:\s*([\d,]+)\s*-\s*([\d,]+)\s*", text, re.IGNORECASE)
    if not m:
        raise ValueError(tx("{text} isn't a region like chr7:117,480,025-117,668,665.", text=repr(text)))
    chrom, start, end = m.group(1).upper(), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if not 1 <= start <= end:
        raise ValueError(tx("{text}: the start has to be at least 1 and not after the end.", text=repr(text)))
    return chrom, start, end


def parse_coordinates(text):
    """A region (parse_region()) or a single position, "chr7:117,559,590" ->
    (chrom, start, end, single), start == end for a position. 1-based,
    inclusive. Raises ValueError."""
    m = re.fullmatch(r"\s*(?:chr)?([0-9]{1,2}|X|Y|M|MT)\s*:\s*([\d,]+)\s*", text, re.IGNORECASE)
    if m:
        pos = int(m.group(2).replace(",", ""))
        if pos < 1:
            raise ValueError(tx("{text}: a position is at least 1.", text=repr(text)))
        return m.group(1).upper(), pos, pos, True
    try:
        return (*parse_region(text), False)
    except ValueError:
        raise ValueError(tx("{text} is neither a region like chr7:117,480,025-117,668,665 nor a position like "
                            "chr7:117,559,590.", text=repr(text))) from None


# A single position: the paste is looked for within this many bases of it,
# plus twice its own length, before the whole gene is searched (spec 14 C5).
SEARCH_FLANK = 10_000
FLANK = 5000


def analysed_region(gene_start, gene_end, chrom_length, window=None, flank=FLANK):
    """(start, end) the paste is aligned to: the coordinate window if one was
    given (spec 6.1), else the gene (its mask) with `flank` bases on each
    side, within the chromosome. Raises ValueError for a window off the
    chromosome."""
    if window is not None:
        start, end = window
        if end > chrom_length:
            raise ValueError(tx("The window ends at {end:,}, past the chromosome's end ({length:,}).", end=end,
                                length=chrom_length))
        return start, end
    return max(1, gene_start - flank), min(chrom_length, gene_end + flank)


@dataclass
class Event:
    """One difference between the paste and the reference."""
    kind: str           # SNV, insertion, deletion, duplication, inversion, complex
    pos: int            # VCF-style
    ref: str
    alt: str
    start: int          # reference bases replaced: start..end (end = start - 1 for an insertion)
    end: int
    inserted: str       # what the patient has instead of start..end
    size: int
    query_start: int = 0    # where it is in the paste (1-based), for display
    apply: bool = False
    note: str = ""
    # Spec 14 C8: what it looks like -- only a label; the haplotype comes from kind/start/end/inserted.
    label: str = ""
    microhomology: str = ""             # bases the breakpoints can slide over (left-aligned here)
    breakpoints: list = field(default_factory=list)     # [[first, last] possible place of each breakpoint]
    blocks: list = field(default_factory=list)          # [[q0, q1, ref_start, ref_end, +1/-1]]: inserted bases' sources


@dataclass
class Alignment:
    """Where a paste lies on a reference region, and how it differs."""
    query_length: int
    reverse: bool                   # the paste aligns as its reverse complement
    ref_start: int = 0              # covered reference interval (1-based), 0 if nothing aligned
    ref_end: int = 0
    query_start: int = 0            # aligned part of the paste (1-based, of the paste as aligned)
    query_end: int = 0
    matches: int = 0
    columns: int = 0
    score: int = 0
    events: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    method: str = ""

    @property
    def identity(self):
        return self.matches / self.columns if self.columns else 0.0


# --- seeds and blocks --------------------------------------------------------------

def _kmers(seq):
    """(values, valid) of every K-mer of seq: 2 bits per base; a K-mer with N is invalid."""
    codes = _CODE[np.frombuffer(seq.encode("ascii"), np.uint8)]
    n = len(codes) - K + 1
    if n <= 0:
        return np.zeros(0, np.int64), np.zeros(0, bool)
    bad = codes == 4
    c = np.where(bad, 0, codes).astype(np.int64)
    values = np.zeros(n, np.int64)
    for j in range(K):
        values = (values << 2) | c[j:j + n]
    bad_run = np.convolve(bad.astype(np.int32), np.ones(K, np.int32), "valid") > 0
    return values, ~bad_run


def _unique(values, valid):
    """{value: position} of the K-mers that occur exactly once, as sorted arrays."""
    pos = np.nonzero(valid)[0]
    vals = values[pos]
    u, first, counts = np.unique(vals, return_index=True, return_counts=True)
    keep = counts == 1
    return u[keep], pos[first[keep]]


def _seed_pairs(query, ref_unique):
    """(query positions, reference positions) of shared unique K-mers."""
    qv, qp = _unique(*_kmers(query))
    rv, rp = ref_unique
    _, qi, ri = np.intersect1d(qv, rv, assume_unique=True, return_indices=True)
    return qp[qi], rp[ri]


def _blocks(q, r):
    """Runs of seeds on one diagonal and consecutive in the paste, as blocks
    [q0, q1) x [r0, r1) (0-based, half-open), sorted by q0."""
    if len(q) == 0:
        return []
    order = np.lexsort((q, r - q))
    q, r = q[order], r[order]
    d = r - q
    breaks = np.nonzero((np.diff(d) != 0) | (np.diff(q) != 1))[0] + 1
    starts = np.r_[0, breaks]
    ends = np.r_[breaks, len(q)] - 1
    blocks = [[int(q[s]), int(q[e]) + K, int(r[s]), int(r[e]) + K] for s, e in zip(starts, ends)]
    return sorted(blocks)


def _chain(blocks):
    """The heaviest chain (by block length) whose blocks increase in both the
    paste and the reference: exact, with a Fenwick tree over reference starts."""
    if not blocks:
        return []
    ranks = {v: i + 1 for i, v in enumerate(sorted({b[2] for b in blocks}))}
    tree_val = [0] * (len(ranks) + 1)
    tree_idx = [-1] * (len(ranks) + 1)
    best, prev = [0] * len(blocks), [-1] * len(blocks)

    def query(i):  # max over ranks 1..i
        v, idx = 0, -1
        while i > 0:
            if tree_val[i] > v:
                v, idx = tree_val[i], tree_idx[i]
            i -= i & -i
        return v, idx

    def update(i, v, idx):
        while i < len(tree_val):
            if v > tree_val[i]:
                tree_val[i], tree_idx[i] = v, idx
            i += i & -i

    pending = []    # blocks with the same q0 enter the tree together, so none chains onto another
    for n, b in enumerate(blocks):
        if pending and blocks[pending[0]][0] != b[0]:
            for m in pending:
                update(ranks[blocks[m][2]], best[m], m)
            pending = []
        v, idx = query(ranks[b[2]] - 1)
        best[n], prev[n] = v + (b[1] - b[0]), idx
        pending.append(n)
    n = max(range(len(blocks)), key=best.__getitem__)
    chain = []
    while n >= 0:
        chain.append(list(blocks[n]))
        n = prev[n]
    chain.reverse()
    # Consecutive blocks may overlap (a repeat at a boundary): trim the later one.
    out = []
    for b in chain:
        if out:
            over = max(out[-1][1] - b[0], out[-1][3] - b[2], 0)
            b = [b[0] + over, b[1], b[2] + over, b[3]]
            if b[1] - b[0] <= 0 or b[3] - b[2] <= 0:
                continue
        out.append(b)
    return out


def _extend(chain, query, ref):
    """Grows each block base by base into the gaps around it while the bases match."""
    for i, b in enumerate(chain):
        q_lo = chain[i - 1][1] if i else 0
        r_lo = chain[i - 1][3] if i else 0
        while b[0] > q_lo and b[2] > r_lo and query[b[0] - 1] == ref[b[2] - 1]:
            b[0] -= 1
            b[2] -= 1
        q_hi = chain[i + 1][0] if i + 1 < len(chain) else len(query)
        r_hi = chain[i + 1][2] if i + 1 < len(chain) else len(ref)
        while b[1] < q_hi and b[3] < r_hi and query[b[1]] == ref[b[3]]:
            b[1] += 1
            b[3] += 1
    return chain


# --- aligning gaps ---------------------------------------------------------------------

def _parasail():
    import parasail

    return parasail


def _ops(cigar):
    """[(op, length)] from parasail's decoded CIGAR (=, X, I: paste only, D: reference only)."""
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([=XIDM])", cigar)]


def _walk(ops, query, ref, q0, r0, diffs, counts):
    """Follows an alignment from paste position q0 and reference position r0
    (0-based), adding each run of non-matching columns to diffs as
    (q, r, ref bases, paste bases) and counting matches and columns."""
    q, r = q0, r0
    run = None
    for op, n in ops:
        if op == "=":
            if run:
                diffs.append(run)
                run = None
            q, r = q + n, r + n
            counts[0] += n
        else:
            run = run or [q, r, "", ""]
            if op == "X":
                run[2] += ref[r:r + n]
                run[3] += query[q:q + n]
                q, r = q + n, r + n
            elif op == "I":
                run[3] += query[q:q + n]
                q += n
            else:
                run[2] += ref[r:r + n]
                r += n
        counts[1] += n
    if run:
        diffs.append(run)
    return q, r


def _same_length(a, b, q0, r0, diffs, counts):
    """paste a against reference b of the same length, base for base, when
    they differ at few positions (the common case: a repeat between two seeds
    on one diagonal). False if they differ too much for that."""
    x = np.frombuffer(a.encode("ascii"), np.uint8) != np.frombuffer(b.encode("ascii"), np.uint8)
    where = np.nonzero(x)[0]
    if len(where) > max(3, len(a) // 100):
        return False
    runs = np.split(where, np.nonzero(np.diff(where) != 1)[0] + 1) if len(where) else []
    for run in runs:
        i, j = int(run[0]), int(run[-1]) + 1
        diffs.append([q0 + i, r0 + i, b[i:j], a[i:j]])
    counts[0] += len(a) - len(where)
    counts[1] += len(a)
    return True


def _gap(query, ref, q0, q1, r0, r1, inverted, diffs, counts):
    """What lies between two chained blocks: the same sequence, a few
    substitutions, an inversion (reverse blocks inside the gap on both
    sides), or anything else, aligned with parasail."""
    a, b = query[q0:q1], ref[r0:r1]
    if len(a) == len(b) and _same_length(a, b, q0, r0, diffs, counts):
        return
    inside = [x for x in inverted if q0 <= x[0] and x[1] <= q1 and r0 <= x[2] and x[3] <= r1
              and x[1] - x[0] >= 2 * K]
    if inside:
        iq0, iq1 = min(x[0] for x in inside), max(x[1] for x in inside)
        ir0, ir1 = min(x[2] for x in inside), max(x[3] for x in inside)
        iq0, iq1, ir0, ir1 = _grow_inversion(query, ref, q0, q1, r0, r1, iq0, iq1, ir0, ir1)
        _global(query, ref, q0, iq0, r0, ir0, diffs, counts)
        diffs.append([iq0, ir0, ref[ir0:ir1], query[iq0:iq1], True])
        counts[1] += max(iq1 - iq0, ir1 - ir0)
        _global(query, ref, iq1, q1, ir1, r1, diffs, counts)
    else:
        _global(query, ref, q0, q1, r0, r1, diffs, counts)


def _global(query, ref, q0, q1, r0, r1, diffs, counts):
    """Aligns paste[q0:q1] with ref[r0:r1] end to end."""
    a, b = query[q0:q1], ref[r0:r1]
    if not a and not b:
        return
    if a == b:
        counts[0] += len(a)
        counts[1] += len(a)
        return
    if not a or not b or len(a) * len(b) > MAX_GAP_CELLS:
        diffs.append([q0, r0, b, a])
        counts[1] += max(len(a), len(b))
        return
    result = _parasail().nw_trace_scan_sat(a, b, GAP_OPEN, GAP_EXTEND, _parasail().dnafull)
    ops = _ops(_decoded(result))
    # Two unrelated stretches still align, with gaps scattered to catch chance matches: many separate
    # changes and few matching columns (gap columns counted) mean one replacement.
    runs = sum(1 for i, (op, _) in enumerate(ops) if op != "=" and (i == 0 or ops[i - 1][0] == "="))
    columns = sum(n for _, n in ops)
    if min(len(a), len(b)) >= 10 and runs >= 4 and sum(n for op, n in ops if op == "=") < LOW_IDENTITY * columns:
        diffs.append([q0, r0, b, a])
        counts[1] += max(len(a), len(b))
        return
    _walk(ops, query, ref, q0, r0, diffs, counts)


def _strip_end_gaps(ops, start=True, end=True):
    """Reference-only columns at the free end(s) of a semi-global alignment:
    (reference bases skipped at the start, the remaining ops)."""
    skip = 0
    while start and ops and ops[0][0] == "D":
        skip += ops.pop(0)[1]
    while end and ops and ops[-1][0] == "D":
        ops.pop()
    return skip, ops


def _end(query, ref, q0, q1, r_anchor, before, diffs, counts):
    """Aligns an end of the paste (paste[q0:q1]) next to a block that starts
    (before=True) or ends at reference position r_anchor, the reference free
    on the far side. Returns the reference position it reaches, or None if
    it doesn't match well enough (it is then left out)."""
    n = q1 - q0
    if n == 0:
        return r_anchor
    margin = n + 50
    if before:
        lo = max(0, r_anchor - margin)
        window = ref[lo:r_anchor]
        fn = _parasail().sg_db_trace_scan_sat
    else:
        lo = r_anchor
        window = ref[r_anchor:r_anchor + margin]
        fn = _parasail().sg_de_trace_scan_sat
    if not window:
        return None
    result = fn(query[q0:q1], window, GAP_OPEN, GAP_EXTEND, _parasail().dnafull)
    # Only the far side is free; a deletion right next to the block is real.
    skip, ops = _strip_end_gaps(_ops(_decoded(result)), start=before, end=not before)
    local_diffs, local_counts = [], [0, 0]
    _, r_end = _walk(ops, query, ref, q0, lo + skip, local_diffs, local_counts)
    if not local_counts[1] or local_counts[0] / local_counts[1] < END_MIN_IDENTITY:
        return None
    diffs.extend(local_diffs)
    counts[0] += local_counts[0]
    counts[1] += local_counts[1]
    return lo + skip if before else r_end


def _decoded(result):
    cigar = result.cigar.decode
    return cigar.decode() if isinstance(cigar, bytes) else cigar


# --- events --------------------------------------------------------------------------------

def _overlap(a, b):
    """Events a and b (a first) overlap -- two insertions at one place too."""
    return b.start <= a.end or (a.start == b.start and a.end < a.start and b.end < b.start)


def _events(fetch, chrom, ref_offset, diffs):
    """The differences as events in reference order. Left-aligning an indel
    can move it onto its neighbour (a deletion next to an inversion, say);
    such an indel stays where the alignment put it."""
    make = lambda d, left=True: _event(fetch, chrom, ref_offset, d[0], d[1], d[2], d[3], inverted=len(d) > 4,
                                       left_align=left)
    found = [(d, make(d)) for d in diffs]
    while True:
        found.sort(key=lambda x: (x[1].start, x[1].end))
        clash = next((i for i in range(len(found) - 1) if _overlap(found[i][1], found[i + 1][1])), None)
        if clash is None:
            return [e for _, e in found]
        for j in (clash + 1, clash):
            d, e = found[j]
            plain = make(d, left=False)
            if (plain.start, plain.end) != (e.start, e.end):
                found[j] = (d, plain)
                break
        else:
            return [e for _, e in found]


def _event(fetch, chrom, ref_offset, q, r, ref_bases, alt_bases, inverted=False, left_align=True):
    """An Event from a difference at 0-based paste position q and region position r."""
    start = ref_offset + r
    if inverted:
        if len(ref_bases) != len(alt_bases):        # an indel inside the inverted stretch: one replacement
            return Event("complex", start, ref_bases, alt_bases, start, start + len(ref_bases) - 1, alt_bases,
                         max(len(ref_bases), len(alt_bases)), query_start=q + 1, apply=True,
                         label="inversion with changes", note="inverted, with an indel inside: applied as one "
                                                              "replacement")
        return Event("inversion", start, ref_bases, alt_bases, start, start + len(ref_bases) - 1, alt_bases,
                     len(ref_bases), query_start=q + 1, apply=True)
    if len(ref_bases) == len(alt_bases) == 1:
        return Event("SNV", start, ref_bases, alt_bases, start, start, alt_bases, 1, query_start=q + 1)
    # VCF style: an indel gets the base before it as anchor, then is left-aligned.
    pos, ref, alt = start, ref_bases, alt_bases
    if not ref or not alt:
        pos -= 1
        anchor = fetch(chrom, pos, pos)
        ref, alt = anchor + ref, anchor + alt
    if left_align:
        pos, ref, alt = normalize_variant(fetch, chrom, pos, ref, alt)
    shared = 0
    while shared < min(len(ref), len(alt)) and ref[shared] == alt[shared]:
        shared += 1
    rep_start = pos + shared
    replaced, inserted = ref[shared:], alt[shared:]
    if not replaced:
        kind, size = "insertion", len(inserted)
        after = fetch(chrom, rep_start, rep_start + len(inserted) - 1)
        if after == inserted:
            kind = "duplication"
    elif not inserted:
        kind, size = "deletion", len(replaced)
    else:
        kind, size = "complex", max(len(replaced), len(inserted))
    event = Event(kind, pos, ref, alt, rep_start, rep_start + len(replaced) - 1, inserted, size, query_start=q + 1)
    if kind == "duplication":
        event.note = f"tandem duplication of {rep_start:,}..{rep_start + size - 1:,}"
    return event


def set_default_apply(events, min_size=APPLY_MIN_SIZE):
    """Spec 6.2: events of min_size or more and every inversion and duplication.
    With none of those, the SNVs are ticked too: a lone SNV is the variant, not
    noise around a larger change. Small events left unticked are reported by
    cdna_analysis.check_notices()."""
    large = [e.kind in ("inversion", "duplication") or e.size >= min_size for e in events]
    for e, big in zip(events, large):
        e.apply = big or (not any(large) and e.kind == "SNV")


# --- what the events look like (spec 14 C8) --------------------------------------------------------
# Labels only: an event's kind, start, end and inserted bases build the
# haplotype, so the paste is always reproduced exactly, whatever the label.

SOURCE_MIN = 20         # inserted bases needed to look for where they come from
SOURCE_SHARE = 0.9      # one source covering this much of them explains them
ADJACENT = 10           # a source within this many bases of the breakpoint is next to it


def _grow(seq, ref, b):
    """A seed block [q0, q1, r0, r1] grown while the bases still match."""
    q0, q1, r0, r1 = b
    while q0 > 0 and r0 > 0 and seq[q0 - 1] == ref[r0 - 1]:
        q0, r0 = q0 - 1, r0 - 1
    while q1 < len(seq) and r1 < len(ref) and seq[q1] == ref[r1]:
        q1, r1 = q1 + 1, r1 + 1
    return q0, q1, r0, r1


def sources(seq, ref, ref_unique):
    """Where stretches of `seq` come from in the region `ref`: [(q0, q1, r0,
    r1, orientation)] (0-based, half-open; r on ref, the stretch reversed and
    complemented for orientation -1), longest first among those that don't
    overlap in seq, then in seq order."""
    found = []
    for orientation, s in ((1, seq), (-1, reverse_complement(seq))):
        for b in _blocks(*_seed_pairs(s, ref_unique)):
            q0, q1, r0, r1 = _grow(s, ref, b)
            if orientation < 0:
                q0, q1 = len(seq) - q1, len(seq) - q0
            found.append((q0, q1, r0, r1, orientation))
    chosen = []
    for q0, q1, r0, r1, o in sorted(found, key=lambda p: p[0] - p[1]):
        # Growing base by base may run a few bases into the next piece by chance: trim, don't drop.
        for c0, c1, *_ in chosen:
            if c0 <= q0 < c1:
                cut = c1 - q0
                q0, r0, r1 = q0 + cut, r0 + cut * (o > 0), r1 - cut * (o < 0)
            if c0 < q1 <= c1:
                cut = q1 - c0
                q1, r1, r0 = q1 - cut, r1 - cut * (o > 0), r0 + cut * (o < 0)
        if q1 - q0 > 0 and all(q1 <= c[0] or c[1] <= q0 for c in chosen):
            chosen.append((q0, q1, r0, r1, o))
    return sorted(chosen)


def _covered(ref_unique_positions, i):
    """A unique K-mer of the region covers its base i."""
    import bisect

    k = bisect.bisect_right(ref_unique_positions, i)
    return k > 0 and ref_unique_positions[k - 1] > i - K


def _slide(e, ref, offset, positions):
    """Microhomology (a deletion) or repeat context (an insertion): how far the
    left-aligned event could move right and stay the same change."""
    s = e.start - offset
    if e.kind == "deletion":
        n = e.end - e.start + 1
        k = 0
        while k < n and s + n + k < len(ref) and ref[s + k] == ref[s + n + k]:
            k += 1
        bases = ref[s:s + k]
    else:
        seq, k = e.inserted, 0
        while seq and k < 10 * len(seq) and s + k < len(ref) and seq[k % len(seq)] == ref[s + k]:
            k += 1
        bases = (seq * (k // len(seq) + 1))[:k] if seq else ""
    notes = []
    if k:
        e.microhomology = bases
        e.breakpoints = [[e.start, e.start + k], [e.end, e.end + k]]
        notes.append(f"{'microhomology' if e.kind == 'deletion' else 'repeat'} {bases[:20]}"
                     f"{'...' if k > 20 else ''} ({k} bp): the breakpoints can lie anywhere over {k} bases; shown "
                     "left-aligned")
    if e.kind == "deletion" and not (_covered(positions, s) and _covered(positions, e.end - offset)):
        notes.append("a breakpoint lies in repeated sequence (no unique 21-mer covers it): its place is less certain")
    e.note = "; ".join([x for x in [e.note] + notes if x])


def _blocks_text(e, offset):
    parts, at = [], 0
    for q0, q1, r0, r1, o in e.blocks:
        if q0 > at:
            parts.append(f"{at + 1}-{q0} novel")
        a, b = (r0, r1) if o > 0 else (r1, r0)
        parts.append(f"{q0 + 1}-{q1} = {a:,}-{b:,} ({'+' if o > 0 else '-'})")
        at = q1
    if at < len(e.inserted):
        parts.append(f"{at + 1}-{len(e.inserted)} novel")
    return "blocks of the inserted bases: " + "; ".join(parts)


def _label_inserted(e, ref, offset, ref_unique):
    """The label of the bases an insertion or a replacement puts in, from where they come from."""
    seq = e.inserted
    s, t = e.start - offset, e.end - offset          # replaced ref[s..t] (none for an insertion)
    replaces = t >= s
    if replaces and len(seq) == t - s + 1:
        inv = reverse_complement(ref[s:t + 1])
        if sum(x != y for x, y in zip(seq, inv)) <= 0.05 * len(seq):
            return "inversion", ""
    found = sources(seq, ref, ref_unique) if len(seq) >= SOURCE_MIN else []
    found = [p for p in found if p[1] - p[0] >= SOURCE_MIN]
    e.blocks = [[q0, q1, offset + r0, offset + r1 - 1, o] for q0, q1, r0, r1, o in found]
    best = max(found, key=lambda p: p[1] - p[0]) if found else None
    if not found:
        return ("delins (novel sequence)" if replaces else
                "tandem duplication" if e.kind == "duplication" else "insertion (novel sequence)"), ""
    if best[1] - best[0] < SOURCE_SHARE * len(seq):
        return "complex replacement", _blocks_text(e, offset)
    q0, q1, r0, r1, o = best
    source = f"{offset + r0:,}..{offset + r1 - 1:,}"
    before, after = (s, t + 1) if replaces else (s, s)
    adjacent = abs(r1 - before) <= ADJACENT or abs(r0 - after) <= ADJACENT
    if not replaces:
        if o > 0:
            return ("tandem duplication", "") if adjacent else ("templated insertion", f"a copy of {source}")
        return (("inverted duplication", f"an inverted copy of {source}, right next to it") if adjacent else
                ("inverted templated insertion", f"an inverted copy of {source}"))
    if o > 0:
        return "deletion with templated insertion", f"the inserted bases are a copy of {source}"
    if r0 < t + 1 and s < r1:
        return "deletion-inversion-deletion", f"the inserted bases are {source}, inverted; the rest deleted"
    if adjacent:
        return "deletion with inverted duplication", f"the inserted bases are an inverted copy of {source}, next to it"
    return "deletion with inverted templated insertion", f"the inserted bases are an inverted copy of {source}"


def classify(events, ref, offset):
    """Labels the events of a paste aligned to the region `ref` (first base
    at `offset`), in place: SNV, deletion (with its microhomology), insertion
    (novel; a tandem, inverted or templated copy), inversion, the deletion
    and inversion around one breakpoint pair, a delins or a templated
    replacement, else a complex replacement with the blocks its bases come
    from. Labels never change what is applied."""
    ref_unique = _unique(*_kmers(ref))
    positions = sorted(int(p) for p in ref_unique[1])
    for e in events:
        if e.label:
            continue
        if e.kind in ("SNV", "inversion"):
            e.label = e.kind
        elif e.kind == "deletion":
            e.label = "deletion"
            _slide(e, ref, offset, positions)
        else:
            e.label, detail = _label_inserted(e, ref, offset, ref_unique)
            # In repeats the alignment finds no inverted seeds and gives an exact inversion as a replacement.
            if e.label == "inversion" and e.inserted == reverse_complement(ref[e.start - offset:e.end - offset + 1]):
                e.kind = "inversion"
            if e.kind in ("insertion", "duplication"):
                _slide(e, ref, offset, positions)
            e.note = "; ".join(x for x in (e.note, detail) if x)
    # A deletion right next to an inversion is one rearrangement with it.
    events = sorted(events, key=lambda e: (e.start, e.end))
    for i, e in enumerate(events):
        if e.kind != "inversion":
            continue
        left = i > 0 and events[i - 1].kind == "deletion" and e.start - events[i - 1].end <= 3
        right = i + 1 < len(events) and events[i + 1].kind == "deletion" and events[i + 1].start - e.end <= 3
        if left or right:
            name = "deletion-inversion-deletion" if left and right else "deletion-inversion" if left else \
                "inversion-deletion"
            for x in [events[i - 1]] * left + [e] + [events[i + 1]] * right:
                x.label = name
                x.note = "; ".join(y for y in (x.note, f"one {name} with its neighbour(s)") if y)


# --- align ----------------------------------------------------------------------------------

def align(query, ref, ref_offset, fetch, chrom, apply_min_size=APPLY_MIN_SIZE):
    """Aligns a paste (parse_sequence()) to the reference region `ref`, whose
    first base is at ref_offset on `chrom`. fetch(chrom, a, b) reads the
    reference (for left-aligning indels, which may move outside the region).
    Returns an Alignment; a paste that aligns as its reverse complement is
    aligned that way and marked reverse."""
    if len(query) < SEED_FROM:
        return _align_short(query, ref, ref_offset, fetch, chrom, apply_min_size)
    ref_unique = _unique(*_kmers(ref))
    forward = _chain(_blocks(*_seed_pairs(query, ref_unique)))
    rc = reverse_complement(query)
    backward = _chain(_blocks(*_seed_pairs(rc, ref_unique)))
    covered = lambda chain: sum(b[1] - b[0] for b in chain)
    reverse = covered(backward) > 2 * covered(forward) and covered(backward) > 0.5 * len(query)
    seq, chain = (rc, backward) if reverse else (query, forward)
    result = Alignment(len(query), reverse, method="seeds (unique 21-mers) and parasail")
    if not chain:
        result.notes.append("No 21-base stretch of the sequence occurs once in the region: it doesn't align here.")
        return result
    chain = _extend(chain, seq, ref)
    # Reverse blocks of the other orientation, for inversions inside gaps.
    other = reverse_complement(seq)
    inverted = [[len(seq) - b[1], len(seq) - b[0], b[2], b[3]]
                for b in _blocks(*_seed_pairs(other, ref_unique))]
    diffs, counts = [], [0, 0]
    events = []
    start_word, end_word = ("last", "first") if reverse else ("first", "last")      # of the paste as given
    first_r = _end(seq, ref, 0, chain[0][0], chain[0][2], True, diffs, counts)
    q_start = 0 if first_r is not None else chain[0][0]
    if first_r is None:
        first_r = chain[0][2]
        result.notes.append(f"The {start_word} {chain[0][0]:,} bases of the sequence don't match the reference next to "
                            "the aligned part; they are left out.")
    for i, b in enumerate(chain):
        counts[0] += b[1] - b[0]
        counts[1] += b[1] - b[0]
        if i + 1 == len(chain):
            break
        nb = chain[i + 1]
        _gap(seq, ref, b[1], nb[0], b[3], nb[2], inverted, diffs, counts)
    last = chain[-1]
    last_r = _end(seq, ref, last[1], len(seq), last[3], False, diffs, counts)
    q_end = len(seq) if last_r is not None else last[1]
    if last_r is None:
        last_r = last[3]
        result.notes.append(f"The {end_word} {len(seq) - last[1]:,} bases of the sequence don't match the reference next "
                            "to the aligned part; they are left out.")
    events += _events(fetch, chrom, ref_offset, diffs)
    set_default_apply(events, apply_min_size)
    classify(events, ref, ref_offset)
    result.events = sorted(events, key=lambda e: (e.start, e.end))
    result.ref_start, result.ref_end = ref_offset + first_r, ref_offset + last_r - 1
    result.query_start, result.query_end = q_start + 1, q_end
    result.matches, result.columns = counts
    return result


def _grow_inversion(seq, ref, q0, q1, r0, r1, iq0, iq1, ir0, ir1):
    """Moves an inversion's ends outward while the paste still matches the
    reverse complement of the reference there (within the gap)."""
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    # Paste position iq0 - 1 pairs with reference position ir1 (inverted), and so on.
    while iq0 > q0 and ir1 < r1 and comp.get(seq[iq0 - 1]) == ref[ir1]:
        iq0, ir1 = iq0 - 1, ir1 + 1
    while iq1 < q1 and ir0 > r0 and comp.get(seq[iq1]) == ref[ir0 - 1]:
        iq1, ir0 = iq1 + 1, ir0 - 1
    return iq0, iq1, ir0, ir1


def _align_short(query, ref, ref_offset, fetch, chrom, apply_min_size):
    """50-499 bases: the best end-to-end placement of the paste in the region,
    found exactly (score pass over the whole region, then a traceback in a
    window around where it ends), in both orientations."""
    p = _parasail()
    best = None
    for reverse, seq in ((False, query), (True, reverse_complement(query))):
        scan = p.sg_dx_scan_sat(seq, ref, GAP_OPEN, GAP_EXTEND, p.dnafull)
        if best is None or scan.score > best[0]:
            best = (scan.score, reverse, seq, scan.end_ref)
    score, reverse, seq, end_ref = best
    lo = max(0, end_ref - 2 * len(seq) - 50)
    window = ref[lo:end_ref + 1]
    traced = p.sg_dx_trace_scan_sat(seq, window, GAP_OPEN, GAP_EXTEND, p.dnafull)
    skip, ops = _strip_end_gaps(_ops(_decoded(traced)))
    diffs, counts = [], [0, 0]
    _, r_end = _walk(ops, seq, ref, 0, lo + skip, diffs, counts)
    result = Alignment(len(query), reverse, ref_offset + lo + skip, ref_offset + r_end - 1, 1, len(seq),
                       counts[0], counts[1], score, method="dynamic programming (parasail) over the whole region")
    result.events = _events(fetch, chrom, ref_offset, diffs)
    set_default_apply(result.events, apply_min_size)
    classify(result.events, ref, ref_offset)
    if result.identity < END_MIN_IDENTITY:
        result.notes.append(f"Only {100 * result.identity:.1f}% identity at its best place: the sequence may not "
                            "belong to this region.")
    return result


# --- compatibility report and chosen events ----------------------------------------------

def compatibility(alignment, gene_start, gene_end, spliceai=None, canonical=None, other=None):
    """The "Check compatibility" report (spec 6.2) as lines of text and the
    problems that block a run: identity, the covered reference interval, the
    parts of the gene the paste doesn't cover, orientation, a better fit to
    the other genome build (other: (build, Alignment) or None), and whether
    SpliceAI's gene interval and the canonical transcript differ
    (spliceai, canonical: (start, end) or None)."""
    lines, problems, warnings = [], [], []
    a = alignment
    if not a.columns:
        problems.append(tx("The sequence doesn't align to the region."))
        return {"lines": lines + a.notes, "problems": problems, "warnings": warnings}
    lines.append(tx("Identity {pct:.2f}% over {n:,} aligned columns ({method}).", pct=100 * a.identity, n=a.columns,
                    method=a.method))
    q = (a.query_length - a.query_end + 1, a.query_length - a.query_start + 1) if a.reverse else (a.query_start, a.query_end)
    lines.append(tx("Covered reference: {start:,}..{end:,} ({n:,} bases; bases {q0:,}..{q1:,} of the sequence).",
                    start=a.ref_start, end=a.ref_end, n=a.ref_end - a.ref_start + 1, q0=q[0], q1=q[1]))
    uncovered = [(s, e) for s, e in ((gene_start, min(gene_end, a.ref_start - 1)),
                                     (max(gene_start, a.ref_end + 1), gene_end)) if s <= e]
    if uncovered:
        warnings.append(tx("The sequence doesn't cover the whole gene: {parts} come from the reference and are "
                           "treated as normal.", parts=tx("{a} and {b}", a=f"{uncovered[0][0]:,}..{uncovered[0][1]:,}",
                                                     b=f"{uncovered[1][0]:,}..{uncovered[1][1]:,}")
                           if len(uncovered) == 2 else f"{uncovered[0][0]:,}..{uncovered[0][1]:,}"))
    if other is not None:
        build, oa = other
        if oa.columns and (oa.score > a.score * 1.02 or oa.identity > a.identity + 0.005):
            warnings.append(tx("The sequence fits {build} better ({pct:.2f}% identity there): check the genome build.",
                               build=build, pct=100 * oa.identity))
    if spliceai and canonical and tuple(spliceai) != tuple(canonical):
        lines.append(tx("SpliceAI's gene interval ({s0:,}..{s1:,}) and the canonical transcript ({c0:,}..{c1:,}) "
                        "differ; the analysis covers both.", s0=spliceai[0], s1=spliceai[1], c0=canonical[0],
                        c1=canonical[1]))
    lines.append(tx("{n} difference(s) found.", n=len(a.events)))
    return {"lines": lines + a.notes, "problems": problems, "warnings": warnings}


def validate_events(events, region_start, region_end, anchors=()):
    """Problems and warnings for the events chosen to build the haplotype:
    overlapping events and breakpoints outside the analysed region block it;
    an event touching an anchor (anchors: [(name, position)], e.g. the
    transcript start and end and the ATG) is a warning."""
    problems, warnings = [], []
    chosen = sorted((e for e in events if e.apply), key=lambda e: (e.start, e.end))
    for a, b in zip(chosen, chosen[1:]):
        # An insertion (end = start - 1) right at another event's edge still overlaps nothing.
        if b.start <= a.end or (b.start == a.start and a.start > a.end and b.start > b.end):
            problems.append(tx("{kind_a} at {a:,} and {kind_b} at {b:,} overlap; choose one.", kind_a=a.kind,
                               a=a.start, kind_b=b.kind, b=b.start))
    for e in chosen:
        if e.start < region_start or max(e.end, e.start - 1) > region_end:
            problems.append(tx("The {kind} at {pos:,} reaches outside the analysed region ({start:,}..{end:,}).",
                               kind=e.kind, pos=e.start, start=region_start, end=region_end))
        named = set()
        for name, pos in anchors:
            if e.start - 1 <= pos <= e.end + 1 and name not in named:
                named.add(name)
                warnings.append(tx("The {kind} at {start:,} touches the {name} ({pos:,}).", kind=e.kind, start=e.start,
                                   name=name, pos=pos))
    return {"problems": problems, "warnings": warnings}
