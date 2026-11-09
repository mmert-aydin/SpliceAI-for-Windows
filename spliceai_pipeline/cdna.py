"""cDNA, protein and comparisons of transcripts (spec 7.4).

A transcript is a list of exons (genomic, inclusive, in transcript order) on
a sequence -- the reference genome for the reference and wild-type
transcripts, the patient haplotype for the patient's. build() splices it
(reverse-complemented on the minus strand) and translates from the anchored
ATG to the first stop: no ATG at the anchor means "ATG lost" (no start is
invented), no stop reached means "no stop codon". summarize() compares it
with the reference transcript: frame, premature stop, NMD by the 50-55 nt
rule, protein length. concordance() is the exon table (identical, start or
end shifted by n, missing, extra), hgvs_r() and hgvs_p() approximate
HGVS-style descriptions -- predicted, approximate, never a diagnosis --
and pairwise() an alignment of two cDNAs or proteins for display.
"""
from dataclasses import dataclass, field

from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from .align import reverse_complement
from .texts import tx, tx_join

STOPS = ("TAA", "TAG", "TGA")
# A stop more than this many bases before the last exon-exon junction makes
# the mRNA a target of nonsense-mediated decay (the 50-55 nt rule; 55 is its
# conservative end); one in the first START_PROXIMAL bases of the coding
# sequence may escape it by re-initiation.
NMD_DISTANCE = 55
START_PROXIMAL = 150
LONG_INSERT = 20
FRAME_PROBE = 30        # bases of the CDS looked for past a second change


@dataclass
class CDNA:
    exons: list                 # [(start, end)] genomic, transcript order
    strand: str
    sequence: str               # transcript orientation
    exon_starts: list           # 0-based index in `sequence` where each exon begins
    cds_start: int = None       # index of the A of the ATG; None: ATG lost
    cds_end: int = None         # index after the stop codon; None: no stop codon reached
    protein: str = ""           # without the stop
    notes: list = field(default_factory=list)
    c1: int = None              # the reference's c. numbering, from the annotated CDS: index of c.1 ...
    c_end: int = None           # ... and after its last base. Unset: cds_start, cds_end

    @property
    def stop_start(self):
        return self.cds_end - 3 if self.cds_end is not None else None

    @property
    def numbering(self):
        """(index of c.1, index after the last c. position); HGVS takes both from
        the annotated CDS, which translation may not reach (a CTG start, a TGA read
        through)."""
        return (self.cds_start if self.c1 is None else self.c1, self.cds_end if self.c_end is None else self.c_end)

    def index(self, x):
        """cDNA index (0-based) of genomic position x, or None if it isn't in an exon."""
        for (a, b), o in zip(self.exons, self.exon_starts):
            if a <= x <= b:
                return o + (x - a if self.strand == "+" else b - x)
        return None

    def regions(self):
        """[(name, start, end)] of 5'UTR, CDS and 3'UTR in the cDNA (0-based, half-open)."""
        if self.cds_start is None:
            return [("UTR or non-coding", 0, len(self.sequence))]
        end = self.cds_end if self.cds_end is not None else len(self.sequence)
        parts = [("5'UTR", 0, self.cds_start), ("CDS", self.cds_start, end), ("3'UTR", end, len(self.sequence))]
        return [p for p in parts if p[2] > p[1]]


def build(fetch, strand, exons, atg):
    """The CDNA of `exons` on a sequence read with fetch(a, b) (1-based,
    inclusive), translated from the ATG anchored at genomic position atg (its
    A; on the minus strand the highest of its three bases)."""
    pieces, starts, o = [], [], 0
    for a, b in exons:
        s = fetch(a, b).upper()
        pieces.append(s if strand == "+" else reverse_complement(s))
        starts.append(o)
        o += b - a + 1
    cdna = CDNA(list(exons), strand, "".join(pieces), starts)
    i = cdna.index(atg) if atg is not None else None
    if i is None:
        cdna.notes.append("ATG lost: the anchored start codon isn't in this transcript")
        return cdna
    return translate(cdna, i)


def translate(cdna, i):
    """Translates cdna from index i (the anchored ATG's A) to the first stop."""
    codon = cdna.sequence[i:i + 3]
    if codon != "ATG":
        cdna.notes.append(f"ATG lost: the anchored start codon reads {codon or 'nothing'} here")
        return cdna
    cdna.cds_start = i
    for j in range(i, len(cdna.sequence) - 2, 3):
        if cdna.sequence[j:j + 3] in STOPS:
            cdna.cds_end = j + 3
            break
    end = cdna.cds_end - 3 if cdna.cds_end is not None else i + (len(cdna.sequence) - i) // 3 * 3
    cdna.protein = str(Seq(cdna.sequence[i:end]).translate())
    if cdna.cds_end is None:
        cdna.notes.append("no stop codon: translation runs to the end of the transcript")
    return cdna


def reference(fetch, transcript):
    """The CDNA of an annotated transcript (transcripts.Transcript)."""
    exons = list(transcript.exons) if transcript.strand == "+" else list(reversed(transcript.exons))
    atg = None
    if transcript.cds:
        atg = transcript.cds[0] if transcript.strand == "+" else transcript.cds[1]
    out = build(fetch, transcript.strand, exons, atg)
    if transcript.cds:
        last = transcript.cds[1] if transcript.strand == "+" else transcript.cds[0]
        out.c1, out.c_end = out.index(atg), out.index(last) + 1
    return out


@dataclass
class Summary:
    exon_count: int
    length: int
    protein_length: int
    frame: str                  # "no change", "in frame", "frameshift", "ATG lost", "non-coding"
    ptc: bool = False           # a stop before the reference's
    ptc_codon: int = None       # its codon number
    nmd: str = ""               # "expected", "escapes: ..." or ""
    notes: list = field(default_factory=list)


def _frameshift(cdna, ref):
    """Does cdna read out of the reference's frame at the end of its CDS? The
    net change inside the CDS counts: past the stop nothing can restore the
    frame. A second change past the CDS shows as a suffix too short, and then
    the CDS's last stretch that cdna still has gives the shift."""
    a, b = ref.sequence[ref.cds_start:], cdna.sequence[cdna.cds_start:]
    i = _common_prefix(a, b)
    cds = (ref.cds_end or len(ref.sequence)) - ref.cds_start
    if len(a) - _common_suffix(a, b, min(len(a), len(b)) - i) > cds:
        for end in range(cds, i + FRAME_PROBE - 1, -FRAME_PROBE):
            at = b.find(a[end - FRAME_PROBE:end], i)
            if at >= 0:
                return (at - (end - FRAME_PROBE)) % 3 != 0
    return (len(b) - len(a)) % 3 != 0


def summarize(cdna, ref):
    """Exon count, length, frame relative to the reference CDS, premature
    stop and NMD (50-55 nt rule, with last-exon and start-proximal escape),
    protein length."""
    s = Summary(len(cdna.exons), len(cdna.sequence), len(cdna.protein), "")
    s.notes = list(cdna.notes)
    if ref.cds_start is None:
        s.frame = "non-coding"
        return s
    if cdna.cds_start is None:
        s.frame = "ATG lost"
        return s
    a, b = ref.sequence[ref.cds_start:], cdna.sequence[cdna.cds_start:]
    i = _common_prefix(a, b)
    ref_cds_length = (ref.cds_end or len(ref.sequence)) - ref.cds_start
    if a == b or i >= ref_cds_length:
        s.frame = "no change"
    else:
        s.frame = "frameshift" if _frameshift(cdna, ref) else "in frame"
    if cdna.cds_end is not None and ref.cds_end is not None and s.frame != "no change":
        # The stop is the reference's if it sits as far from the transcript's
        # end and everything from it on is the same (an in-frame deletion
        # keeps it); it is premature if it comes earlier.
        ref_tail, tail = len(ref.sequence) - ref.stop_start, len(cdna.sequence) - cdna.stop_start
        same_stop = tail == ref_tail and _common_suffix(ref.sequence, cdna.sequence,
                                                        min(len(ref.sequence), len(cdna.sequence))) >= ref_tail
        if not same_stop and tail > ref_tail:
            s.ptc, s.ptc_codon = True, len(cdna.protein) + 1
    if cdna.cds_end is not None and (s.ptc or s.frame == "frameshift"):
        s.nmd = nmd(cdna)
    return s


def nmd(cdna):
    """NMD by the 50-55 nt rule for the transcript's stop codon."""
    if cdna.cds_end is None:
        return "not assessed: no stop codon"
    if len(cdna.exons) < 2:
        return "escapes: single exon, no exon-exon junction"
    last_junction = cdna.exon_starts[-1]
    stop = cdna.stop_start
    if stop >= last_junction:
        return "escapes: stop codon in the last exon"
    distance = last_junction - cdna.cds_end
    if distance <= NMD_DISTANCE:
        return f"escapes: stop codon {distance} nt before the last exon-exon junction (50-55 nt rule)"
    if stop - cdna.cds_start < START_PROXIMAL:
        return (f"expected by the 50-55 nt rule (stop {distance} nt before the last junction), but the stop is "
                f"start-proximal (codon {(stop - cdna.cds_start) // 3 + 1}) and may escape it")
    return f"expected: stop codon {distance} nt before the last exon-exon junction (50-55 nt rule)"


def _common_prefix(a, b):
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a, b, limit):
    j = 0
    while j < limit and a[-1 - j] == b[-1 - j]:
        j += 1
    return j


# --- exon concordance -----------------------------------------------------------------

def concordance(ref_exons, other_exons, strand):
    """Per exon, objectively: rows of {"ref": n or None, "other": m or None,
    "status": ...} with status "identical", "start shifted n", "end shifted
    n" (both may apply; + is further 3'), "missing" or "extra (novel)", and a
    summary line. Exons are genomic (start, end) in transcript order; an
    exon given as None (it lies in inserted sequence) counts as extra."""
    sign = 1 if strand == "+" else -1
    first = lambda e: e[0] if strand == "+" else e[1]
    last = lambda e: e[1] if strand == "+" else e[0]
    overlaps = lambda e, f: e is not None and f is not None and e[0] <= f[1] and f[0] <= e[1]
    rows, used = [], set()
    for n, r in enumerate(ref_exons, 1):
        hits = [m for m, o in enumerate(other_exons, 1) if overlaps(r, o)]
        if not hits:
            rows.append({"ref": n, "other": None, "status": tx("missing")})
            continue
        m = hits[0]
        used.update(hits)
        o = other_exons[m - 1]
        parts = []
        ds, de = (first(o) - first(r)) * sign, (last(o) - last(r)) * sign
        if ds:
            parts.append(tx("start shifted {n:+d}", n=ds))
        if de:
            parts.append(tx("end shifted {n:+d}", n=de))
        status = tx_join(", ", parts) if parts else tx("identical")
        shared = [k for k, x in enumerate(ref_exons, 1) if k != n and overlaps(x, o)]
        if shared:
            status = tx("{status} (one exon with reference exon {exons})", status=status,
                        exons=", ".join(map(str, shared)))
        rows.append({"ref": n, "other": m, "status": status, "start_shift": ds, "end_shift": de})
    for m, o in enumerate(other_exons, 1):
        if m not in used:
            rows.append({"ref": None, "other": m, "status": tx("extra (novel)")})
    counts = {}
    for row in rows:
        key = row["status"]
        key = "shifted" if "shifted" in key else ("identical" if key.startswith("identical") else key)
        counts[key] = counts.get(key, 0) + 1
    order = ("identical", "shifted", "missing", "extra (novel)")
    summary = ", ".join(f"{counts[k]} {k}" for k in order if k in counts)
    return rows, summary


def to_reference_exons(exons, cmap):
    """Patient exons (haplotype positions) as reference exons for the
    concordance table: both ends through the coordinate map; an exon with an
    end in inserted sequence keeps the nearest mapped bases, and one wholly
    inside inserted sequence becomes None (extra)."""
    out = []
    for a, b in exons:
        lo, hi = sorted((a, b))
        ra = next((cmap.to_ref(h)[0] for h in range(lo, hi + 1) if cmap.to_ref(h)[0] is not None), None)
        rb = next((cmap.to_ref(h)[0] for h in range(hi, lo - 1, -1) if cmap.to_ref(h)[0] is not None), None)
        out.append(None if ra is None else tuple(sorted((ra, rb))))
    return out


# --- approximate HGVS -------------------------------------------------------------------

def c_position(ref, i):
    """HGVS c. numbering of index i of the reference cDNA."""
    first, end = ref.numbering
    if first is None:
        return f"{i + 1}"           # non-coding: n. numbering
    if i < first:
        return f"-{first - i}"
    end = end if end is not None else len(ref.sequence)
    if i < end:
        return f"{i - first + 1}"
    return f"*{i - end + 1}"


def genomic_of_index(ref, i):
    """Genomic position of cDNA index i (0-based)."""
    k = max(j for j, o in enumerate(ref.exon_starts) if o <= i)
    a, b = ref.exons[k]
    o = i - ref.exon_starts[k]
    if not 0 <= o <= b - a:
        raise ValueError(f"index {i} is past the transcript's end")
    return a + o if ref.strand == "+" else b - o


def c_to_genomic(ref, text):
    """Genomic position of an HGVS c. position on the reference transcript:
    "123", "-45" (5'UTR), "*67" (3'UTR), each with an intronic offset like
    "3718-2477" or "123+5" (+ from an exon's last base, - before its first,
    in transcript direction). Raises ValueError for anything else."""
    import re

    m = re.fullmatch(r"\s*(?:c\.)?([-*]?)(\d+)(?:([+-])(\d+))?\s*", text)
    first, end = ref.numbering
    if not m or first is None:
        raise ValueError(f"{text!r} isn't a c. position this program can resolve")
    prefix, n, sign, offset = m.group(1), int(m.group(2)), m.group(3), m.group(4)
    if n == 0:
        raise ValueError(f"{text!r}: HGVS has no position 0")
    if prefix == "-":
        i = first - n
    elif prefix == "*":
        if end is None:
            raise ValueError("the reference transcript has no stop codon, so c.* positions don't exist")
        i = end - 1 + n
    else:
        i = first + n - 1
    if not 0 <= i < len(ref.sequence):
        raise ValueError(f"c.{text} lies outside the reference transcript")
    x = genomic_of_index(ref, i)
    if sign:
        k = max(j for j, o in enumerate(ref.exon_starts) if o <= i)
        last = ref.exon_starts[k + 1] - 1 if k + 1 < len(ref.exons) else len(ref.sequence) - 1
        if (sign == "+" and i != last) or (sign == "-" and i != ref.exon_starts[k]):
            raise ValueError(f"c.{text}: an intronic offset counts from an exon's {'last' if sign == '+' else 'first'} "
                             "base, and this one isn't")
        d = int(offset) * (1 if sign == "+" else -1)
        x += d if ref.strand == "+" else -d
    return x


def genomic_to_c(ref, x):
    """HGVS c. position of genomic position x on the reference transcript
    ("123", "-45", "*67", or intronic "3718-2477"); None outside the transcript."""
    i = ref.index(x)
    if i is not None:
        return c_position(ref, i)
    sign = 1 if ref.strand == "+" else -1
    # In transcript direction: the last base of the exon before x, the first of the one after.
    t = lambda p: p * sign
    before = [(e, o) for e, o in zip(ref.exons, ref.exon_starts) if t(e[1] if sign > 0 else e[0]) < t(x)]
    after = [(e, o) for e, o in zip(ref.exons, ref.exon_starts) if t(e[0] if sign > 0 else e[1]) > t(x)]
    if not before or not after:
        return None
    (e1, o1), (e2, o2) = before[-1], after[0]
    last_base = e1[1] if sign > 0 else e1[0]
    first_base = e2[0] if sign > 0 else e2[1]
    up, down = abs(x - last_base), abs(first_base - x)
    if up <= down:
        return f"{c_position(ref, o1 + (e1[1] - e1[0]))}+{up}"
    return f"{c_position(ref, o2)}-{down}"


def _rna(seq):
    return seq.lower().replace("t", "u") if len(seq) <= LONG_INSERT else f"[{len(seq)} nt]"


def hgvs_r(ref, cdna, number=None):
    """r. description of how a transcript differs from `ref`, 3'-most, as one
    change: r.(=), >, del, ins, delins. Positions are ref's own c. numbering,
    or number(i) for its index i -- e.g. the reference transcript's c.
    numbering for a wild-type prediction (number_on()), so that every change
    reads in the numbering of the literature."""
    number = number or (lambda k: c_position(ref, k))
    a, b = ref.sequence, cdna.sequence
    i = _common_prefix(a, b)
    j = _common_suffix(a, b, min(len(a), len(b)) - i)
    gone, new = a[i:len(a) - j], b[i:len(b) - j]
    span = lambda: number(i) + ("" if len(gone) == 1 else "_" + number(i + len(gone) - 1))
    if not gone and not new:
        return "r.(=)"
    if not new:
        return f"r.({span()}del)"
    if not gone:
        return f"r.({number(i - 1)}_{number(i)}ins{_rna(new)})"
    if len(gone) == len(new) == 1:
        return f"r.({number(i)}{_rna(gone)}>{_rna(new)})"
    return f"r.({span()}delins{_rna(new)})"


def number_on(reference, other):
    """number(i) for hgvs_r(): index i of `other` (a transcript on the reference
    genome) in `reference`'s c. numbering, through its genomic position
    (intronic positions as c.123+45)."""
    def number(i):
        if not 0 <= i < len(other.sequence):
            return "?"
        return genomic_to_c(reference, genomic_of_index(other, i)) or "?"
    return number


def hgvs_p(ref, cdna):
    """p. description of the protein change (three-letter codes): p.(=),
    substitutions, deletions, insertions, frameshifts with the new stop,
    stops gained and lost (extensions), and p.(Met1?) when the ATG is lost."""
    if ref.cds_start is None:
        return "p.?"
    if cdna.cds_start is None:
        return "p.(Met1?)"
    r, p = ref.protein, cdna.protein
    if r == p and (cdna.cds_end is not None) == (ref.cds_end is not None):
        return "p.(=)"
    i = _common_prefix(r, p)
    aa = lambda s, k: seq3(s[k]) if k < len(s) else "Ter"
    if i >= len(p) and cdna.cds_end is not None and i < len(r):
        return f"p.({aa(r, i)}{i + 1}Ter)"       # the stop comes where the reference's protein went on
    if i >= len(r) and len(p) > len(r) and ref.cds_end is not None:
        tail = f"Ter{len(p) - len(r) + 1}" if cdna.cds_end is not None else "Ter?"
        return f"p.(Ter{len(r) + 1}{aa(p, len(r))}ext{tail})"    # the stop is lost and the protein goes on
    if _frameshift(cdna, ref):
        new_stop = f"Ter{len(p) - i + 1}" if cdna.cds_end is not None else "Ter?"
        return f"p.({aa(r, i)}{i + 1}{aa(p, i)}fs{new_stop})"
    j = _common_suffix(r, p, min(len(r), len(p)) - i)
    gone, new = r[i:len(r) - j], p[i:len(p) - j]
    if cdna.cds_end is not None and len(p) < len(r) and not new and j == 0:
        return f"p.({aa(r, i)}{i + 1}Ter)"
    three = lambda s: "".join(seq3(c) for c in s)
    span = lambda: f"{aa(r, i)}{i + 1}" + ("" if len(gone) == 1 else f"_{aa(r, i + len(gone) - 1)}{i + len(gone)}")
    if not gone:
        return f"p.({aa(r, i - 1)}{i}_{aa(r, i)}{i + 1}ins{three(new)})"
    if not new:
        return f"p.({span()}del)"
    if len(gone) == len(new) == 1:
        return f"p.({aa(r, i)}{i + 1}{aa(p, i)})"
    return f"p.({span()}delins{three(new)})"


# --- pairwise alignment for display -------------------------------------------------------------

def pairwise(a, b, protein=False):
    """(a with gaps, markup, b with gaps): a global alignment for display
    (parasail; DNA +5/-4 or BLOSUM62, gap open 10, extend 1)."""
    import parasail

    matrix = parasail.blosum62 if protein else parasail.dnafull
    if not a or not b:
        return a or "-" * len(b), " " * max(len(a), len(b)), b or "-" * len(a)
    # Two transcripts of one gene differ in a few places: only the stretch between their common start and end
    # is aligned (a 14 kb cDNA pair in milliseconds instead of seconds; spec 15 item 2).
    i = _common_prefix(a, b)
    j = _common_suffix(a, b, min(len(a), len(b)) - i)
    ma, mb = a[i:len(a) - j], b[i:len(b) - j]
    if ma and mb:
        tb = parasail.nw_trace_scan_sat(ma, mb, 10, 1, matrix).traceback
        top, mid, bottom = tb.query, tb.comp, tb.ref
    else:
        top, mid, bottom = ma or "-" * len(mb), " " * max(len(ma), len(mb)), mb or "-" * len(ma)
    return (a[:i] + top + a[len(a) - j:], "|" * i + mid + "|" * j, b[:i] + bottom + b[len(b) - j:])
