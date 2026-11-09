"""Events given as coordinates (spec 14 C6): the "From coordinates" input of
the Sequence-to-cDNA tab. Each line is one event, 1-based and inclusive as in
VCF files and genome browsers:

    deletion chr17:17228330-17228462            (optionally followed by the deleted bases)
    tandem duplication chr17:17228330-17228462  (a copy right after the stretch)
    inverted tandem duplication chr17:17228330-17228462
    inversion chr17:17228330-17228462
    insertion chr17:17228330 ACGTACGT           (inserted after the position, as in VCF)

or HGVS g. on the chromosome's RefSeq accession, which fixes the build:

    NC_000017.11:g.17228330_17228462del         (also ...delACGT, dup, inv)
    NC_000017.11:g.17228330_17228331insACGT
    NC_000017.11:g.17228330_17228335delinsACG
    NC_000017.11:g.17228330A>G

Lines starting with # and empty lines are skipped. The events are known
exactly, so they become the differences table directly, with no alignment;
the reference is still checked: every position on the chromosome and in the
region analysed for the gene, no N, and any reference bases given (deleted
bases, a substitution's reference) equal to the FASTA's.
"""
import re
from dataclasses import dataclass

from .align import Event, reverse_complement
from .texts import tx

# The chromosomes' RefSeq accessions, from UCSC's chromAlias tables
# (hgdownload.soe.ucsc.edu/goldenPath/<build>/bigZips/<build>.chromAlias.txt,
# 2026-10-28). hg19's chrM is NC_001807.4, not the rCRS of hg38.
_GRCH37 = {"1": 10, "2": 11, "3": 11, "4": 11, "5": 9, "6": 11, "7": 13, "8": 10, "9": 11, "10": 10, "11": 9,
           "12": 11, "13": 10, "14": 8, "15": 9, "16": 9, "17": 10, "18": 9, "19": 9, "20": 10, "21": 8, "22": 10,
           "X": 10, "Y": 9}
_GRCH38 = {"1": 11, "2": 12, "3": 12, "4": 12, "5": 10, "6": 12, "7": 14, "8": 11, "9": 12, "10": 11, "11": 10,
           "12": 12, "13": 11, "14": 9, "15": 10, "16": 10, "17": 11, "18": 10, "19": 10, "20": 11, "21": 9, "22": 11,
           "X": 11, "Y": 10}


def _accession(chrom, version):
    number = {"X": 23, "Y": 24}.get(chrom) or int(chrom)
    return f"NC_{number:06d}.{version}"


ACCESSIONS = {acc: (build, chrom) for build, table in (("hg19", _GRCH37), ("hg38", _GRCH38))
              for chrom, version in table.items() for acc in [_accession(chrom, version)]}
ACCESSIONS["NC_001807.4"] = ("hg19", "M")
ACCESSIONS["NC_012920.1"] = ("hg38", "M")

KINDS = {"deletion": "deletion", "del": "deletion", "tandem duplication": "duplication", "duplication": "duplication",
         "dup": "duplication", "inverted tandem duplication": "inverted duplication",
         "inverted duplication": "inverted duplication", "inversion": "inversion", "inv": "inversion",
         "insertion": "insertion", "ins": "insertion"}
_PLAIN = re.compile(r"(?P<kind>[a-z ]+?)\s+(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|M|MT)\s*:\s*(?P<a>[\d,]+)"
                    r"(?:\s*-\s*(?P<b>[\d,]+))?(?:\s+(?P<seq>[ACGTN]+))?(?P<inv>\s+inverted)?", re.IGNORECASE)
_HGVS = re.compile(r"(?P<acc>NC_\d+\.\d+)\s*:\s*g\.(?P<a>\d+)(?:_(?P<b>\d+))?"
                   r"(?:(?P<op>delins|del|dup|inv|ins)(?P<seq>[ACGTN]*)|(?P<sub_ref>[ACGT])>(?P<sub_alt>[ACGT]))",
                   re.IGNORECASE)


@dataclass
class CoordEvent:
    line: int
    text: str
    kind: str           # deletion, duplication, inverted duplication, inversion, insertion, delins, SNV
    chrom: str          # "17"
    start: int          # the reference stretch (an insertion: the base it follows, start == end)
    end: int
    sequence: str = ""  # inserted bases (insertion, delins, SNV's alternative)
    ref: str = ""       # reference bases given, to check


def parse(text, build):
    """[CoordEvent] of the "From coordinates" box. Raises ValueError naming the
    line for anything it can't read, or an HGVS accession of the other build."""
    out = []
    for n, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        m = _HGVS.fullmatch(line)
        if m:
            out.append(_hgvs(n, line, m, build))
            continue
        m = _PLAIN.fullmatch(line)
        if not m or m.group("kind").strip().lower() not in KINDS:
            raise ValueError(tx("Line {n} ({line}): not an event like \"deletion chr17:17228330-17228462\", "
                                "\"insertion chr17:17228330 ACGT\" or NC_000017.11:g.17228330_17228462del.", n=n,
                                line=repr(line)))
        kind = KINDS[m.group("kind").strip().lower()]
        if m.group("inv") and kind == "duplication":
            kind = "inverted duplication"
        a = int(m.group("a").replace(",", ""))
        b = int(m.group("b").replace(",", "")) if m.group("b") else None
        seq = (m.group("seq") or "").upper()
        chrom = "M" if m.group("chrom").upper() == "MT" else m.group("chrom").upper()
        if kind == "insertion":
            if b is not None or not seq:
                raise ValueError(tx("Line {n}: an insertion is a position and the inserted bases, e.g. "
                                    "\"insertion chr17:17228330 ACGT\" (inserted after that position).", n=n))
            out.append(CoordEvent(n, line, kind, chrom, a, a, sequence=seq))
            continue
        if b is None:
            b = a
        if kind != "deletion" and seq:
            raise ValueError(tx("Line {n}: only a deletion takes bases after its coordinates (the deleted ones).", n=n))
        out.append(CoordEvent(n, line, kind, chrom, a, b, ref=seq))
    for e in out:
        if e.start < 1 or e.end < e.start:
            raise ValueError(tx("Line {n}: the start has to be at least 1 and not after the end.", n=e.line))
    return out


def _hgvs(n, line, m, build):
    acc = m.group("acc").upper()
    if acc not in ACCESSIONS:
        raise ValueError(tx("Line {n}: {acc} isn't a GRCh37 or GRCh38 chromosome accession.", n=n, acc=acc))
    acc_build, chrom = ACCESSIONS[acc]
    if acc_build != build:
        raise ValueError(tx("Line {n}: {acc} is a {assembly} chromosome, but {build} is chosen.", n=n, acc=acc,
                            assembly="GRCh38 (hg38)" if acc_build == "hg38" else "GRCh37 (hg19)", build=build))
    a = int(m.group("a"))
    b = int(m.group("b")) if m.group("b") else a
    if m.group("sub_ref"):
        if m.group("b"):
            raise ValueError(tx("Line {n}: a substitution is at one position.", n=n))
        return CoordEvent(n, line, "SNV", chrom, a, a, sequence=m.group("sub_alt").upper(),
                          ref=m.group("sub_ref").upper())
    op, seq = m.group("op").lower(), m.group("seq").upper()
    if op == "ins":
        if b != a + 1 or not seq:
            raise ValueError(tx("Line {n}: an HGVS insertion names the two bases it lies between (g.100_101insACGT) "
                                "and the inserted bases.", n=n))
        return CoordEvent(n, line, "insertion", chrom, a, a, sequence=seq)
    if op == "delins":
        if not seq:
            raise ValueError(tx("Line {n}: delins needs the inserted bases.", n=n))
        return CoordEvent(n, line, "delins", chrom, a, b, sequence=seq)
    kind = {"del": "deletion", "dup": "duplication", "inv": "inversion"}[op]
    if seq and len(seq) != b - a + 1:
        raise ValueError(tx("Line {n}: {bases} bases given for {positions:,} positions.", n=n, bases=len(seq),
                            positions=b - a + 1))
    return CoordEvent(n, line, kind, chrom, a, b, ref=seq)


def to_events(parsed, chrom, region, fetch, chrom_length):
    """align.Event objects (apply=True) of parsed events on `chrom` (the
    gene's), checked against the reference: fetch(a, b) reads the
    chromosome; region (start, end): the region analysed for the gene, which
    every event has to lie in. Raises ValueError (naming the line) when a
    check fails."""
    events = []
    for e in parsed:
        n = e.line
        if e.chrom != chrom:
            raise ValueError(tx("Line {n}: chr{chrom}, but the gene is on chr{gene_chrom}.", n=n, chrom=e.chrom,
                                gene_chrom=chrom))
        last = e.end + (1 if e.kind == "insertion" else 0)
        if e.end > chrom_length:
            raise ValueError(tx("Line {n}: {end:,} is past the chromosome's end ({length:,}).", n=n, end=e.end,
                                length=chrom_length))
        if e.start < region[0] or last > region[1]:
            raise ValueError(tx("Line {n}: {start:,}..{end:,} is outside the region analysed for this gene "
                                "({r0:,}..{r1:,}).", n=n, start=e.start, end=e.end, r0=region[0], r1=region[1]))
        ref = fetch(e.start, e.end).upper()
        if "N" in ref:
            raise ValueError(tx("Line {n}: the reference has N (unknown bases) in {start:,}..{end:,}.", n=n,
                                start=e.start, end=e.end))
        if e.ref and e.ref != ref:
            cut = lambda s: s[:30] + ("..." if len(s) > 30 else "")
            raise ValueError(tx("Line {n}: the reference at {where} is {ref}, not {given}.", n=n,
                                where=f"chr{chrom}:{e.start:,}" + (f"-{e.end:,}" if e.end != e.start else ""),
                                ref=cut(ref), given=cut(e.ref)))
        size = e.end - e.start + 1
        note = "given as coordinates"
        if e.kind == "deletion":
            anchor = fetch(e.start - 1, e.start - 1).upper()
            ev = Event("deletion", e.start - 1, anchor + ref, anchor, e.start, e.end, "", size)
        elif e.kind in ("duplication", "inverted duplication"):
            copy = ref if e.kind == "duplication" else reverse_complement(ref)
            ev = Event("duplication", e.end, ref[-1], ref[-1] + copy, e.end + 1, e.end, copy, size)
            note = (f"tandem duplication of {e.start:,}..{e.end:,}" if e.kind == "duplication" else
                    f"inverted tandem duplication of {e.start:,}..{e.end:,} (its reverse complement right after "
                    "it)") + ", given as coordinates"
        elif e.kind == "inversion":
            ev = Event("inversion", e.start, ref, reverse_complement(ref), e.start, e.end, reverse_complement(ref),
                       size)
        elif e.kind == "insertion":
            ev = Event("insertion", e.start, ref, ref + e.sequence, e.start + 1, e.start, e.sequence,
                       len(e.sequence))
        elif e.kind == "delins":
            ev = Event("complex", e.start, ref, e.sequence, e.start, e.end, e.sequence, max(size, len(e.sequence)))
        else:
            ev = Event("SNV", e.start, ref, e.sequence, e.start, e.start, e.sequence, 1)
        ev.apply, ev.note = True, note
        events.append(ev)
    events.sort(key=lambda x: (x.start, x.end))
    for a, b in zip(events, events[1:]):
        if b.start <= a.end or (a.start == b.start and a.end < a.start and b.end < b.start):
            raise ValueError(tx("The {kind_a} at {a:,} and the {kind_b} at {b:,} overlap.", kind_a=a.kind, a=a.start,
                                kind_b=b.kind, b=b.start))
    return events
