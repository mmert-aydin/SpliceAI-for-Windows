"""What the viewer shows, computed once, without Qt (spec 8.2): the shared
columns (haplotype.columns), the bases and per-base tracks on them, every
transcript's exons as column runs coloured by the reference exon they match,
the differences, the columns without scores, and the hover text. The
viewer's widgets only draw it; the tests check it.

A ViewData comes from a finished analysis (cdna_analysis.Analysis, or one
reopened from a results bundle), or -- per-base view only -- from a CSV the
viewer exported ("Export visible region as CSV").

Spec 15: a transcript row finds its exon blocks through position-to-column
maps and works out which column is in which exon only for the columns asked
for, so a 2.3 Mb gene opens at once; flipped() is the same view in transcript
orientation (columns reversed, bases complemented) for a minus-strand gene.

Spec 16: event_info describes each applied event (cdna_analysis.event_summary);
details() is the side panel's sections, with the transcripts that agree on
a column merged into one line ("WT predicted 1-5, Sample 1-4: intron 13").
"""
import copy
import csv
import dataclasses
import threading
from dataclasses import dataclass, field

import numpy as np

from . import cdna, cdna_analysis, haplotype
from .texts import tx

TRACKS = ("ref_acc", "ref_don", "pat_acc", "pat_don", "delta_acc", "delta_don")
TRACK_NAMES = {"ref_acc": "reference acceptor", "ref_don": "reference donor", "pat_acc": "sample acceptor",
               "pat_don": "sample donor", "delta_acc": "delta acceptor", "delta_don": "delta donor"}


@dataclass
class ExonBlock:
    first: int              # column indices, inclusive
    last: int
    number: int             # in transcript order
    identity: int           # the reference exon it matches (1-based), 0 for none (novel)
    start: int = 0          # its genomic (or haplotype) coordinates
    end: int = 0
    status: str = ""        # its concordance with the reference
    duplicate: bool = False     # in inserted bases, a copy of reference exon `identity`


@dataclass
class TranscriptRow:
    name: str
    group: str              # "reference", "alternative", "wt", "patient"
    rank: int
    given: np.ndarray = None    # per-column membership read from a CSV; else worked out when asked for
    blocks: list = field(default_factory=list)
    relative: float = None
    on_haplotype: bool = False
    weak: bool = False      # a sample candidate with weak evidence (cdna_analysis.evidence)
    exons: list = field(default_factory=list)       # [(start, end)] in transcript order, on `positions`
    positions: np.ndarray = None    # per column: the reference (or haplotype) position the row is placed by
    cds: tuple = None       # (first, last) column of the coding part, for an annotated transcript; UTRs are thinner

    @property
    def membership(self):
        """Per column: exon number (transcript order), -1 intron, 0 outside."""
        if self.given is None:
            self.given = self.membership_range(0, len(self.positions) - 1)
        return self.given

    def membership_range(self, c0, c1):
        """membership[c0:c1 + 1], worked out for those columns only."""
        if self.given is not None:
            return self.given[c0:c1 + 1]
        pos = self.positions[c0:c1 + 1]
        out = np.zeros(len(pos), np.int16)
        if not self.exons:
            return out
        order = np.argsort([e[0] for e in self.exons])
        starts = np.array([self.exons[i][0] for i in order], np.int64)
        ends = np.array([self.exons[i][1] for i in order], np.int64)
        out[(pos >= starts[0]) & (pos <= ends.max())] = -1
        idx = np.searchsorted(starts, pos, side="right") - 1
        inside = idx >= 0
        inside[inside] &= pos[inside] <= ends[idx[inside]]
        out[inside] = (order + 1)[idx[inside]]
        return out

    def flipped(self, n):
        f = copy.copy(self)
        f.positions = None if self.positions is None else self.positions[::-1]
        f.given = None if self.given is None else self.given[::-1]
        f.blocks = [dataclasses.replace(b, first=n - 1 - b.last, last=n - 1 - b.first) for b in self.blocks]
        f.cds = None if self.cds is None else (n - 1 - self.cds[1], n - 1 - self.cds[0])
        return f


_COMPLEMENT = np.arange(128, dtype=np.uint32)
for _a, _b in ("AT", "TA", "CG", "GC", "at", "ta", "cg", "gc"):
    _COMPLEMENT[ord(_a)] = ord(_b)


def _complement(bases):
    """A '<U1' array of bases, complemented."""
    codes = np.ascontiguousarray(bases).view(np.uint32)
    return _COMPLEMENT[np.minimum(codes, 127)].view("<U1")


VIEW_FLANK = 3000      # bases of reference shown beyond the gene on each side
FLANK_WHY = "outside the gene: reference only, not scored"


class ViewData:
    def __init__(self):
        self.n = 0
        self.limited = ""           # set when built from a CSV: what the viewer can't show then
        self._base = None           # flipped(): the view it was flipped from
        self.reversed = False       # columns in transcript orientation, 5'->3', for a minus-strand gene
        self._lock = threading.Lock()   # the lazy caches, filled from the viewer's drawing thread too

    # --- from an analysis ------------------------------------------------------------------
    @classmethod
    def from_analysis(cls, result):
        self = cls()
        self.result = result
        g, plan, cmap, pred = result.gene, result.plan, result.cmap, result.patient
        self.chrom, self.strand, self.build, self.gene = "chr" + g.chrom, g.strand, result.build, g.name
        self.threshold = result.settings.threshold
        self.source = result.wt_source
        region0, region_seq = result.compat.region[0], result.compat.ref_seq
        # Reference beyond the gene on both sides (not scored, hatched), so an event near its start or end isn't
        # squeezed against the edge of the view.
        lo = max(region0, cmap.region_start, plan["start"] - VIEW_FLANK)
        hi = min(region0 + len(region_seq) - 1, cmap.region_end, plan["end"] + VIEW_FLANK)
        self.ref, self.hap, self.kind = haplotype.columns(cmap, lo, hi)
        self.n = len(self.ref)
        has_ref, has_hap = self.ref > 0, self.hap > 0
        self.ref_base = np.full(self.n, "-", "<U1")
        self.pat_base = np.full(self.n, "-", "<U1")
        ref_arr = np.frombuffer(region_seq.encode("ascii"), "S1").astype("<U1")
        hap_arr = np.frombuffer(result.hap.encode("ascii"), "S1").astype("<U1")
        self.ref_base[has_ref] = ref_arr[self.ref[has_ref] - region0]
        self.pat_base[has_hap] = hap_arr[self.hap[has_hap] - 1]
        # Bases the patient's sequence covered (the rest is the reference, shown muted).
        a = result.compat.alignment
        self.provided = (self.kind == haplotype.INSERTED) | ((self.ref >= a.ref_start) & (self.ref <= a.ref_end))
        self.provided &= has_hap
        t = {name: np.full(self.n, np.nan, np.float32) for name in TRACKS}
        wt = result.wt_values
        self.gene_range = (plan["start"], plan["end"])      # beyond it: the VIEW_FLANK, reference only
        scored = has_ref & (self.ref >= plan["start"]) & (self.ref <= plan["end"])
        t["ref_acc"][scored] = wt[self.ref[scored] - plan["start"], 0]
        t["ref_don"][scored] = wt[self.ref[scored] - plan["start"], 1]
        inside = has_hap & (self.hap >= pred.start) & (self.hap <= pred.end)
        t["pat_acc"][inside] = pred.values[self.hap[inside] - pred.start, 0]
        t["pat_don"][inside] = pred.values[self.hap[inside] - pred.start, 1]
        t["delta_acc"] = t["pat_acc"] - t["ref_acc"]
        t["delta_don"] = t["pat_don"] - t["ref_don"]
        self.tracks = t
        self._ref_cols = np.nonzero(has_ref)[0]
        self._refs = self.ref[self._ref_cols]      # the reference positions shown, ascending
        self._col_of_hap = np.full(int(self.hap.max()) + 2, -1, np.int64)
        self._col_of_hap[self.hap[has_hap]] = np.nonzero(has_hap)[0]
        self.rows = self._rows(result)
        self.events = []
        self.event_info = cdna_analysis.event_summary(result) if result.transcripts else []
        for e in a.events:
            if e.apply:
                first, last = self.columns_of(e.start, e.end)
                if e.end < e.start:     # an insertion: its inserted columns
                    first, last = self.inserted_after(e.start - 1)
                self.events.append((first, last, f"{e.kind} at {e.start:,} ({e.size:,} bp)"))
        if len(self.event_info) != len(self.events):
            self.event_info = []
        self.chrom_end_padded = [tuple(x) for x in plan.get("chrom_end_padded") or []]
        return self

    def _rows(self, result):
        g = result.gene
        rows = []
        by_kind = {(t.kind, t.rank): t for t in result.transcripts}
        ref_t = by_kind.get(("Reference", 0))
        canonical = ref_t.exons if ref_t else []
        # The canonical's exons are identities 1..n; another transcript's exon takes
        # the identity of the canonical exon it overlaps, else 0 (novel).
        identity = lambda e: next((k for k, c in enumerate(canonical, 1) if e and e[0] <= c[1] and c[0] <= e[1]), 0)
        for t in result.transcripts:
            group = {"Reference": "reference", "WT predicted": "wt", "Patient": "patient"}[t.kind]
            name = t.kind if not t.rank else f"{'Sample' if group == 'patient' else t.kind} {t.rank}"
            if group == "reference":
                name = f"Reference ({g.canonical.accession})"
            status = {row["other"]: row["status"] for row in t.vs_ref[1] if row["other"]}
            identities = [identity(e) for e in t.ref_exons]
            copies = cdna_analysis.in_insertion(result, t) if group == "patient" else {}
            for number, (copy, _at, _n) in copies.items():
                identities[number - 1] = copy
            rows.append(self._row(name, group, t.rank, t.exons, identities, on_hap=group == "patient",
                                  relative=t.relative, status=status, copies=copies))
            if group == "reference" and g.canonical.cds:
                rows[-1].cds = self.columns_of(*g.canonical.cds)
            rows[-1].weak = "weak evidence" in cdna_analysis.evidence(t, result.settings)
        for alt in g.transcripts[1:]:
            exons = list(alt.exons) if alt.strand == "+" else list(reversed(alt.exons))
            rows.append(self._row(alt.accession, "alternative", 0, exons, [identity(e) for e in exons]))
            if alt.cds:
                rows[-1].cds = self.columns_of(*alt.cds)
        return rows

    def _row(self, name, group, rank, exons, identities, on_hap=False, relative=None, status=None, copies=()):
        row = TranscriptRow(name, group, rank, relative=relative, on_haplotype=on_hap, exons=list(exons),
                            positions=self.hap if on_hap else self.ref)
        for number, (e, ident) in enumerate(zip(exons, identities), 1):
            for first, last in self._spans(e[0], e[1], on_hap):
                row.blocks.append(ExonBlock(first, last, number, ident, e[0], e[1],
                                            (status or {}).get(number, ""), bool(ident) and number in copies))
        return row

    def _spans(self, a, b, on_hap):
        """[(first, last)] columns of positions a..b of the haplotype or the reference, one pair for each run
        of adjacent columns (an inversion or a deletion can part a sample exon's bases); none if none is shown."""
        if on_hap:
            cols = self._col_of_hap[max(a, 0):max(min(b + 1, len(self._col_of_hap)), 0)]
            cols = np.sort(cols[cols >= 0])
            if not len(cols):
                return []
            cuts = np.nonzero(np.diff(cols) > 1)[0]
            return [(int(x), int(y)) for x, y in zip(np.r_[cols[0], cols[cuts + 1]], np.r_[cols[cuts], cols[-1]])]
        i, j = int(np.searchsorted(self._refs, a)), int(np.searchsorted(self._refs, b, side="right")) - 1
        return [(int(self._ref_cols[i]), int(self._ref_cols[j]))] if i <= j else []

    def flipped(self):
        """This view in the other orientation (spec 15 item 1): columns reversed and
        bases complemented, e.g. a minus-strand gene in transcript orientation
        (5' to 3', exon 1 on the left). Positions, scores and exon numbers stay
        what they are, so labels and tooltips read the same."""
        f = copy.copy(self)
        f._base, f.reversed, n = self, not self.reversed, self.n
        r = lambda a: a[::-1]
        f.ref, f.hap, f.kind, f.provided = r(self.ref), r(self.hap), r(self.kind), r(self.provided)
        f.ref_base, f.pat_base = _complement(r(self.ref_base)), _complement(r(self.pat_base))
        f.tracks = {k: r(v) for k, v in self.tracks.items()}
        for cache in ("_pyramids", "_unscored"):
            f.__dict__.pop(cache, None)
        f.rows = [row.flipped(n) for row in self.rows]
        f.events = [(n - 1 - b, n - 1 - a, text) for a, b, text in self.events]      # event_info: the same
        return f

    # --- coordinates --------------------------------------------------------------------------------
    def columns_of(self, start, end):
        """(first, last) columns of reference positions start..end (clipped)."""
        if self._base is not None:
            a, b = self._base.columns_of(start, end)
            return self.n - 1 - b, self.n - 1 - a
        refs = self._refs
        a = int(np.searchsorted(refs, start))
        b = int(np.searchsorted(refs, end, side="right")) - 1
        a, b = min(max(a, 0), len(refs) - 1), min(max(b, 0), len(refs) - 1)
        return int(self._ref_cols[a]), int(self._ref_cols[max(a, b)])

    def column_of_ref(self, pos):
        return self.columns_of(pos, pos)[0]

    def inserted_after(self, pos):
        """(first, last) inserted columns right after reference position pos (or pos's column twice)."""
        if self._base is not None:
            a, b = self._base.inserted_after(pos)
            return self.n - 1 - b, self.n - 1 - a
        c = self.column_of_ref(pos)
        last = c
        while last + 1 < self.n and self.kind[last + 1] == haplotype.INSERTED:
            last += 1
        if last > c:
            return c + 1, last
        # A deleted pos: haplotype.columns() puts inserted bases before the deleted ones.
        last = c
        while last >= 0 and self.kind[last] == haplotype.DELETED:
            last -= 1
        first = last
        while first >= 0 and self.kind[first] == haplotype.INSERTED:
            first -= 1
        return (first + 1, last) if last > first and last < c else (c, c)

    def label(self, c):
        """The ruler's name for column c: a genomic position, or INS1+23 inside inserted sequence."""
        if self.ref[c]:
            return f"{self.ref[c]:,}"
        if hasattr(self, "result"):
            return self.result.cmap.to_ref(int(self.hap[c]))[1]
        return getattr(self, "hap_labels", {}).get(c, "inserted")

    def position_text(self, c):
        """Column c's place, short: "chr21:14,533,713 · c.57+12345", or the inserted base's label."""
        if not self.ref[c]:
            return tx("inserted: {label}", label=self.label(c))
        text = f"{self.chrom}:{self.ref[c]:,}"
        if hasattr(self, "result") and self.result.transcripts:
            c_pos = cdna.genomic_to_c(self.result.transcripts[0].cdna, int(self.ref[c]))
            if c_pos:
                return tx("{place} · c.{c}", place=text, c=c_pos)
            return tx("{place} · outside the reference transcript", place=text)
        return text

    def canonical_region(self, c):
        """Column c on the canonical transcript, by its exon blocks: "exon 3", "intron 13", "" outside."""
        ref = next((r for r in self.rows if r.group == "reference"), None)
        if ref is None:
            return ""
        blocks = sorted(ref.blocks, key=lambda b: b.number)
        blk = next((b for b in blocks if b.first <= c <= b.last), None)
        if blk is not None:
            return tx("exon {n}", n=blk.number)
        k = next((a.number for a, b in zip(blocks, blocks[1:])
                  if min(a.last, b.last) < c < max(a.first, b.first) and b.number == a.number + 1), None)
        return tx("intron {n}", n=k) if k else ""

    def details(self, c):
        """The side panel for column c (spec 16 item 3): [(section, [(key, value)])] -- Position, Bases,
        Scores, and Transcripts, where the rows that agree on c share one line."""
        position = []
        if self.ref[c]:
            position.append((tx("Genomic"), f"{self.chrom}:{self.ref[c]:,}"))
            if hasattr(self, "result") and self.result.transcripts:
                c_pos = cdna.genomic_to_c(self.result.transcripts[0].cdna, int(self.ref[c]))
                position.append(("c.", f"c.{c_pos}" if c_pos else tx("outside the reference transcript")))
        else:
            position.append((tx("Inserted"), self.label(c)))
        region = self.canonical_region(c)
        if region:
            position.append((tx("Canonical"), region))
        status = {haplotype.INVERTED: tx("inverted"), haplotype.SUBSTITUTED: tx("changed"),
                  haplotype.INSERTED: tx("inserted"), haplotype.DELETED: tx("deleted")}.get(int(self.kind[c]), tx("same"))
        if self.hap[c] and not self.provided[c]:
            status = tx("{status} (from the reference)", status=status)
        bases = [(tx("Reference"), self.ref_base[c]), (tx("Sample"), self.pat_base[c]), (tx("Status"), status)]
        f = lambda v: "--" if np.isnan(v) else f"{v:.3f}"
        t = self.tracks
        scores = [(tx("Ref. acceptor"), f(t["ref_acc"][c])), (tx("Ref. donor"), f(t["ref_don"][c])),
                  (tx("Sample acceptor"), f(t["pat_acc"][c])), (tx("Sample donor"), f(t["pat_don"][c])),
                  (tx("Delta acceptor"), f(t["delta_acc"][c])), (tx("Delta donor"), f(t["delta_don"][c]))]
        return [(tx("Position"), position), (tx("Bases"), bases), (tx("Scores"), scores),
                (tx("Transcripts"), [(what, rows) for rows, what in self.memberships(c, region)])]

    def memberships(self, c, region=None):
        """[(rows, what)]: which transcript rows have column c in which exon or intron, rows that agree merged:
        ("WT predicted 1-5, Sample 1-4", "intron 13"). Exons are named after the reference exon they match."""
        region = self.canonical_region(c) if region is None else region
        groups = {}
        for row in self.rows:
            m = int(row.membership_range(c, c)[0])
            if not m:
                continue
            if m > 0:
                blk = next((b for b in row.blocks if b.first <= c <= b.last), None)
                what = tx("novel exon (its exon {n})", n=m)
                if blk is not None and blk.identity:
                    what = (tx("duplicate of exon {n}", n=blk.identity) if blk.duplicate else
                            tx("exon {n}", n=blk.identity))
            else:
                what = region if region.startswith("intron") else \
                    tx("intron (in reference {region})", region=region) if region else tx("intron")
            groups.setdefault(what, []).append(row)
        out = []
        for what, rows in groups.items():
            names, ranks = [], {}
            for row in rows:
                if row.group in ("wt", "patient") and row.rank:
                    ranks.setdefault(row.name.rsplit(" ", 1)[0], []).append(row.rank)
                else:
                    names.append(row.name)
            names += [f"{prefix} {cdna_analysis._numbers(r)}" for prefix, r in ranks.items()]
            out.append((", ".join(names), what))
        return out

    def hover(self, c):
        """details(c) as plain text, one "key: value" per line."""
        return details_text(self.details(c))

    def max_abs_delta(self, a, b):
        """The largest |delta| (acceptor or donor) over columns a..b, from the pyramids (O(log n))."""
        best = 0.0
        for track in ("delta_acc", "delta_don"):
            levels, lo, hi, k = self.pyramid(track), max(a, 0), min(b, self.n - 1), 0
            while lo <= hi and k < len(levels):
                level = levels[k]
                if lo % 2:
                    best = max(best, abs(float(level[lo])))
                    lo += 1
                if hi % 2 == 0 and lo <= hi:
                    best = max(best, abs(float(level[hi])))
                    hi -= 1
                lo, hi, k = lo // 2, (hi - 1) // 2 if hi >= 0 else -1, k + 1
        return best

    def unscored(self):
        """[(first, last, reason)]: column runs without a reference score, or
        scored with N in place of the sequence past a chromosome end (5,000
        bases of context missing). The viewer hatches them."""
        if self._base is not None:
            return [(self.n - 1 - b, self.n - 1 - a, why) for a, b, why in self._base.unscored()]
        out = []
        missing = np.isnan(self.tracks["ref_acc"]) & (self.ref > 0)
        gene = getattr(self, "gene_range", None)
        if gene is not None:        # the flank shown past the gene is never scored: not a gap in the scores
            beyond = (self.ref > 0) & ((self.ref < gene[0]) | (self.ref > gene[1]))
            missing &= ~beyond
            out += _runs(beyond, FLANK_WHY)
        for a, b, why in _runs(missing, "no score"):
            out.append((a, b, why))
        for lo, hi in getattr(self, "chrom_end_padded", []):
            a, b = self.columns_of(lo, hi)
            out.append((a, b, "within 5 kb of the chromosome's end: scored with N as context"))
        return out

    def unscored_transcribed(self):
        """unscored() runs that fall on an exon or intron of the canonical or a
        curated alternative transcript (the completeness rule of spec 8.2)."""
        rows = [r for r in self.rows if r.group in ("reference", "alternative")]
        bad = []
        for a, b, why in self.unscored():
            if why != FLANK_WHY and any((r.membership_range(a, b) != 0).any() for r in rows):
                bad.append((a, b, why))
        return bad

    # --- per-pixel maxima for zoomed-out drawing ---------------------------------------------------
    def pyramid(self, track):
        """[level 0, 1, ...]: per 2^k columns the largest value (for deltas,
        the largest change either way, sign kept), built once per track."""
        with self._lock:
            return self._pyramid(track)

    def _pyramid(self, track):
        cache = self.__dict__.setdefault("_pyramids", {})
        if track not in cache:
            v = np.nan_to_num(self.tracks[track], nan=0.0).astype(np.float32)
            levels = [v]
            while len(levels[-1]) > 1:
                x = levels[-1]
                if len(x) % 2:
                    x = np.append(x, 0.0)
                a, b = x[0::2], x[1::2]
                levels.append(np.where(np.abs(a) >= np.abs(b), a, b))
            cache[track] = levels
        return cache[track]

    # --- CSV of a stretch of columns, and back ------------------------------------------------------
    def csv_header(self):
        return (["column", "genomic position", "haplotype position", "reference base", "sample base",
                 "reference acceptor", "reference donor", "sample acceptor", "sample donor", "delta acceptor",
                 "delta donor"] + [r.name for r in self.rows])

    def write_csv(self, path, first=0, last=None):
        """One row per column (spec 8.2's "Export visible region as CSV")."""
        last = self.n - 1 if last is None else last
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fh.write(f"# {self.gene} {self.build} {self.chrom} {self.strand}, threshold {self.threshold}; "
                     f"SpliceAI Sequence-to-cDNA viewer export. {cdna_analysis.RESEARCH_NOTE}\n")
            w = csv.writer(fh)
            w.writerow(self.csv_header())
            num = lambda v: "" if np.isnan(v) else f"{v:.6f}"
            members = [r.membership_range(first, last) for r in self.rows]
            for c in range(first, last + 1):
                member = [("exon " + str(m) if m > 0 else "intron" if m < 0 else "")
                          for m in (int(x[c - first]) for x in members)]
                w.writerow([c, self.ref[c] or "", self.hap[c] or "", self.ref_base[c], self.pat_base[c]]
                           + [num(self.tracks[k][c]) for k in TRACKS] + member)

    @classmethod
    def from_csv(cls, path):
        """The per-base view of a CSV written by write_csv(): bases, tracks and
        exon membership, but no transcripts' sequences or comparisons."""
        self = cls()
        with open(path, newline="", encoding="utf-8") as fh:
            first = fh.readline()
            if not first.startswith("#") or "viewer export" not in first:
                raise ValueError(f"{path} isn't a CSV exported by the Sequence-to-cDNA viewer.")
            parts = first[2:].split()
            self.gene, self.build, self.chrom, self.strand = parts[0], parts[1], parts[2], parts[3].rstrip(",")
            self.threshold = float(first.split("threshold ")[1].split(";")[0])
            reader = csv.reader(fh)
            header = next(reader)
            body = list(reader)
        # Exports from before spec 14 C7 said "patient" where they now say "sample".
        if [h.replace("patient", "sample") for h in header[:11]] != cls().csv_header_base() or not body:
            raise ValueError(f"{path} doesn't have the viewer's columns.")
        self.n = len(body)
        col = lambda i, kind: np.array([kind(r[i]) if r[i] != "" else 0 for r in body], np.int64)
        self.ref, self.hap = col(1, int), col(2, int)
        self.kind = np.where(self.ref == 0, haplotype.INSERTED, np.where(self.hap == 0, haplotype.DELETED,
                                                                            haplotype.SAME)).astype(np.int8)
        self.ref_base = np.array([r[3] for r in body], "<U1")
        self.pat_base = np.array([r[4] for r in body], "<U1")
        self.provided = self.hap > 0
        self.tracks = {k: np.array([float(r[5 + i]) if r[5 + i] else np.nan for r in body], np.float32)
                       for i, k in enumerate(TRACKS)}
        self.source = "(from a CSV export)"
        self.rows = []
        for j, name in enumerate(header[11:], 11):
            m = np.array([int(r[j].split()[1]) if r[j].startswith("exon") else -1 if r[j] == "intron" else 0
                          for r in body], np.int16)
            group = ("reference" if name.startswith("Reference") else "wt" if name.startswith("WT")
                     else "patient" if name.startswith(("Sample", "Patient")) else "alternative")
            row = TranscriptRow(name, group, 0, m)
            for a, b, number in _exon_runs(m):
                row.blocks.append(ExonBlock(a, b, number, number if group == "reference" else 0))
            self.rows.append(row)
        self._ref_cols = np.nonzero(self.ref > 0)[0]
        self._refs = self.ref[self._ref_cols]
        self.events, self.event_info, self.chrom_end_padded = [], [], []
        self._col_of_hap = np.full(int(self.hap.max()) + 2, -1, np.int64)
        self._col_of_hap[self.hap[self.hap > 0]] = np.nonzero(self.hap > 0)[0]
        self.limited = ("Loaded from a CSV export: the per-base view only. Transcript sequences, comparisons and "
                        "the exon schematic need the results file (.spliceai-cdna).")
        return self

    @staticmethod
    def csv_header_base():
        return ["column", "genomic position", "haplotype position", "reference base", "sample base",
                "reference acceptor", "reference donor", "sample acceptor", "sample donor", "delta acceptor",
                "delta donor"]


def details_text(sections):
    """Side-panel sections as plain text."""
    return "\n".join(line for title, pairs in sections
                     for line in [title] + [f"  {k}: {v}" for k, v in pairs])


def _runs(mask, why):
    idx = np.nonzero(mask)[0]
    if not len(idx):
        return []
    breaks = np.nonzero(np.diff(idx) != 1)[0]
    starts, ends = np.r_[idx[0], idx[breaks + 1]], np.r_[idx[breaks], idx[-1]]
    return [(int(a), int(b), why) for a, b in zip(starts, ends)]


def _exon_runs(membership):
    out, c = [], 0
    while c < len(membership):
        m = int(membership[c])
        if m > 0:
            d = c
            while d + 1 < len(membership) and membership[d + 1] == m:
                d += 1
            out.append((c, d, m))
            c = d + 1
        else:
            c += 1
    return out
