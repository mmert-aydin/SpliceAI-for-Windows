"""The Sequence-to-cDNA analysis end to end, without Qt (spec 6, 7, 9): what
the tab's "Check compatibility" and "Run" do, so it can be tested and saved
(bundle.py) on its own.

check_compatibility() reads the paste, aligns it to the gene's region (or the
coordinate window) and returns the report, the differences and the anchors
to check them against. run_analysis() then takes the reference's values
(database, or computed live), the wild type's top transcripts, the patient
haplotype with the chosen edits, its prediction near the edits, the
patient's top transcripts, and every transcript's cDNA, protein and
comparison with the reference and the wild type, plus the warnings of spec 9.

Genes without their MANE Select transcript (spec 5.3, decided after Phase
7): where the build's table lacks the MANE Select transcript's locus but has
another NM_ there (hg19: RHBG and 4 more), that NM_ is the reference -- its
cDNA, its transcript start, end and ATG as anchors, and what everything is
compared with -- and the tab and the results bundle say so. With no NM_ at the
locus there is nothing to anchor a transcript to: only the explanation is
shown, and nothing is predicted.

Evidence (decided after Phase 8): a sample candidate whose every splice site
that differs from the WT prediction changed by less than Settings.weak_delta
(SpliceAI's probability at that base, sample minus reference) is labelled
"weak evidence"; the label never changes the ranking. lead() is the
relative-score gap between the sample's ranks 1 and 2.
"""
import re
import time
from dataclasses import dataclass, field
from types import SimpleNamespace

from . import align, cdna, haplotype, splice_paths, transcripts
from .texts import tx, tx_join

# One line for exported files, which travel without the startup pop-up (spec 14 C7).
RESEARCH_NOTE = ("Research use only: predictions from SpliceAI, not measurements; RNA (cDNA) confirmation is required "
                 "for any clinical conclusion.")

@dataclass
class Settings:
    threshold: float = splice_paths.SITE_THRESHOLD
    min_exon: int = splice_paths.MIN_EXON
    min_intron: int = splice_paths.MIN_INTRON
    top: int = splice_paths.TOP
    apply_min_size: int = align.APPLY_MIN_SIZE
    low_concordance: float = 0.9
    weak_delta: float = 0.5


class GeneTable:
    """SpliceAI's genes of one build with their transcripts, found by
    SpliceAI's symbol, the canonical transcript's or any curated one's."""

    def __init__(self, build, refseq_dir, mane_path):
        from .mane import load_mane_select

        self.build = build
        self.genes = transcripts.load_genes(build, refseq_dir, mane_path)
        self.mane_by_symbol = {symbol: acc for acc, symbol in load_mane_select(mane_path).values()}
        self._by = {}
        for g in self.genes:
            names = {g.name} | {t.symbol for t in g.transcripts}
            for n in names:
                self._by.setdefault(n.upper(), []).append(g)

    def find(self, symbol):
        """Genes a symbol names (usually one; several at a read-through or nested locus)."""
        hits = self._by.get(symbol.strip().upper(), [])
        exact = [g for g in hits if g.name.upper() == symbol.strip().upper()]
        return exact or list({id(g): g for g in hits}.values())

    def at(self, chrom, start, end=None):
        """Genes whose mask overlaps chrom:start..end ("7" or "chr7"): MANE Select ones first, then those
        with another canonical transcript, then the rest."""
        end = start if end is None else end
        chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
        hits = [g for g in self.genes if g.chrom == chrom.upper() and g.mask[0] <= end and start <= g.mask[1]]
        return sorted(hits, key=lambda g: (g.canonical is None, g.rule != transcripts.RULE_MANE, g.name))

    def reference_note(self, gene):
        """(kind, text) for a gene without its MANE Select transcript (spec
        5.3), else (None, None). kind "none": no NM_ at the locus, nothing can
        be predicted; "fallback": another NM_ is the reference, and text says
        which."""
        mane = self.mane_by_symbol.get(gene.name) or (gene.canonical and self.mane_by_symbol.get(gene.canonical.symbol))
        if gene.canonical is None:
            text = tx("There is no curated RefSeq transcript (NM_) at {gene}'s locus in the {build} table, so no "
                      "reference cDNA can be built, and there is no transcript start, end and ATG to anchor predicted "
                      "transcripts. Nothing can be predicted for this gene.", gene=gene.name, build=self.build)
            if mane:
                text = tx("{text} Its MANE Select transcript {mane} exists, but the UCSC {build} table lacks it here.",
                          text=text, mane=mane, build=self.build)
            return "none", text
        if gene.rule == transcripts.RULE_LONGEST_CDS and mane:
            return "fallback", tx("No MANE Select transcript is available for this gene on {build}. The reference "
                                  "cDNA and the anchors (transcript start, end and ATG) are taken from {accession}.",
                                  build=self.build, accession=gene.canonical.accession)
        return None, None


def anchors(gene):
    """{name: genomic position} of the canonical transcript: transcript start
    and end (in transcript orientation), the ATG's A and the stop codon's last base."""
    c = gene.canonical
    plus = c.strand == "+"
    out = {"transcript start": c.start if plus else c.end, "transcript end": c.end if plus else c.start}
    if c.cds:
        out["ATG"] = c.cds[0] if plus else c.cds[1]
        out["stop codon"] = c.cds[1] if plus else c.cds[0]
    return out


@dataclass
class Compatibility:
    sequence: str
    region: tuple               # (start, end): the analysed region, where the haplotype is built
    window: tuple               # (start, end) the paste was aligned to
    ref_seq: str                # the analysed region's reference bases
    alignment: object           # align.Alignment
    report: dict                # align.compatibility()
    anchors: dict
    found: tuple = None         # a single position given: (start, end) of the reference the paste was found on
    search: str = ""            # ...and the report's line on where it was looked for


def found(alignment, length):
    """The paste was found where it was aligned: at least half of it aligned, at 90% identity or more."""
    a = alignment
    return bool(a.columns) and a.query_end - a.query_start + 1 >= length / 2 and a.identity >= 0.9


def check_compatibility(text, gene, fetch, chrom_length, window=None, other=None, position=None):
    """"Check compatibility": the paste read, aligned to the gene's region
    (or the coordinate window), and reported on. fetch(a, b) reads the
    reference chromosome; other: (build, align.Alignment) of the same paste on
    the other genome build, or None. position (spec 14 C5): a single
    position to look around -- within twice the paste's length plus
    align.SEARCH_FLANK bases first, the whole gene if it isn't found there;
    the range found is in .found and the report. Raises ValueError for a
    paste that can't be read."""
    sequence = align.parse_sequence(text)
    if position is not None:
        span = 2 * len(sequence) + align.SEARCH_FLANK
        near = (max(1, position - span), min(chrom_length, position + span))
        compat = check_compatibility(text, gene, fetch, chrom_length, window=near, other=other)
        near_found = True
        if not found(compat.alignment, len(sequence)):
            compat = check_compatibility(text, gene, fetch, chrom_length, other=other)
            near_found = False
        a = compat.alignment
        where = f"chr{gene.chrom}:{a.ref_start:,}-{a.ref_end:,}"
        if found(a, len(sequence)):
            compat.found = (a.ref_start, a.ref_end)
            line = (tx("Found within {span:,} bases of {position:,}: {where}.", span=span, position=position,
                       where=where) if near_found else
                    tx("Not found within {span:,} bases of {position:,}; found in the whole gene: {where}.", span=span,
                       position=position, where=where))
        else:
            line = tx("Not found within {span:,} bases of {position:,}, nor in the whole gene.", span=span,
                      position=position)
        compat.search = line
        compat.report["lines"].insert(0, line)
        return compat
    region = align.analysed_region(gene.mask[0], gene.mask[1], chrom_length)
    window = align.analysed_region(gene.mask[0], gene.mask[1], chrom_length, window=window)
    ref_seq = fetch(*region).upper()
    window_seq = ref_seq[window[0] - region[0]:window[1] - region[0] + 1] if region[0] <= window[0] and window[1] <= region[1] \
        else fetch(*window).upper()
    chrom = "chr" + gene.chrom
    a = align.align(sequence, window_seq, window[0], lambda c, x, y: fetch(x, y).upper(), chrom)
    c = gene.canonical
    report = align.compatibility(a, gene.mask[0], gene.mask[1], spliceai=(gene.start, gene.end),
                                 canonical=(c.start, c.end) if c else None, other=other)
    return Compatibility(sequence, region, window, ref_seq, a, report, anchors(gene) if c else {})


def from_coordinates(text, gene, build, fetch, chrom_length):
    """Spec 14 C6: the sample built from events given as coordinates
    (coord_events.py), as a Compatibility: the analysed region of the
    reference with the events applied, the events as the differences, no
    alignment. Raises ValueError for an event that can't be read or doesn't
    match the reference."""
    from . import coord_events

    parsed = coord_events.parse(text, build)
    if not parsed:
        raise ValueError(tx("No event given: one per line, e.g. \"deletion chr17:17228330-17228462\"."))
    region = align.analysed_region(gene.mask[0], gene.mask[1], chrom_length)
    ref_seq = fetch(*region).upper()
    events = coord_events.to_events(parsed, gene.chrom, region, lambda a, b: fetch(a, b), chrom_length)
    align.classify(events, ref_seq, region[0])      # labels and where inserted bases come from; never what is applied
    hap, _ = haplotype.apply_events(ref_seq, region[0], events)
    changed = sum(max(e.end - e.start + 1, len(e.inserted)) for e in events)
    columns = max(len(hap), len(ref_seq))
    a = align.Alignment(len(hap), False, region[0], region[1], 1, len(hap), columns - changed, columns,
                        events=events, method="from coordinates (no alignment)")
    report = {"lines": [tx("Built from {n} event(s) given as coordinates, with no alignment: the sample is the "
                           "reference {region} with them applied, and the reference bases were checked.", n=len(events),
                           region=f"chr{gene.chrom}:{region[0]:,}-{region[1]:,}"),
                        tx("{n} difference(s).", n=len(events))],
              "problems": [], "warnings": []}
    return Compatibility(hap, region, region, ref_seq, a, report, anchors(gene) if gene.canonical else {})


def chosen_check(compat):
    """align.validate_events() for the events chosen now, against the analysed
    region and the anchors. The ATG and the stop are single bases there; their
    codons' middle bases are added, so that its reach of one base takes in all three."""
    a = compat.anchors
    points = list(a.items())
    if "ATG" in a and "stop codon" in a:
        step = 1 if a["stop codon"] > a["ATG"] else -1
        points += [("ATG", a["ATG"] + step), ("stop codon", a["stop codon"] - step)]
    return align.validate_events(compat.alignment.events, *compat.region, anchors=points)


@dataclass
class TranscriptResult:
    kind: str                   # "Reference", "WT predicted", "Patient"
    rank: int                   # 0 for the reference
    relative: float = None
    score: float = None
    exons: list = field(default_factory=list)       # on its own sequence (the patient's: haplotype positions)
    ref_exons: list = field(default_factory=list)   # on the reference (None for an exon inside inserted sequence)
    cdna: object = None
    summary: object = None
    vs_ref: tuple = ("", [], "", "")                # (summary line, rows, r., p.)
    vs_wt: tuple = ("", [], "", "")
    change_vs_ref: str = ""
    change_vs_wt: str = ""
    max_delta: float = None     # sample candidates: max_delta(); None = no site differs from the WT prediction


@dataclass
class Analysis:
    build: str
    gene: object                # transcripts.Gene
    settings: Settings
    compat: Compatibility
    edits: list                 # haplotype.Edit, the ones applied
    hap: str
    cmap: object
    plan: dict
    wt_values: object
    wt_source: str
    patient: object             # haplotype.PatientPrediction
    transcripts: list = field(default_factory=list)
    reference_note: tuple = (None, None)
    warnings: list = field(default_factory=list)
    seconds: float = 0.0
    specificity: list = field(default_factory=list)     # specificity.EventControl, once that control has run
    breakpoints: list = field(default_factory=list)     # specificity.Breakpoints, once that check has run


def exon_labels(base_exons, canonical):
    """Names for the exons of a transcript compared with: "Exon n" after the
    canonical exon each one matches (so a change always reads in the
    canonical transcript's numbering), "its extra exon (start-end)" otherwise
    -- also for a second exon within the same canonical one, so two exons
    never share a name."""
    out, used = [], set()
    for e in base_exons:
        n = next((k for k, c in enumerate(canonical, 1) if e and e[0] <= c[1] and c[0] <= e[1]), None)
        n = None if n in used else n
        used.add(n)
        out.append(f"Exon {n}" if n else f"its extra exon ({e[0]:,}-{e[1]:,})" if e else "its exon in inserted sequence")
    return out


def change_text(rows, other_exons, strand, summary=None, r="", p="", base_exons=(), labels=None):
    """The short structured "Change vs ..." text, e.g. "Exon 8 skipped |
    r.(1043_1170del) | in frame | p.(Val348_Gln390del)". labels: names of
    the compared transcript's exons (exon_labels(); default "Exon n"). An
    extra exon within the compared transcript (base_exons) is a pseudoexon."""
    name = lambda n: labels[n - 1] if labels else f"Exon {n}"
    parts, retained = [], set()
    lo = min((min(e) for e in base_exons if e), default=None)
    hi = max((max(e) for e in base_exons if e), default=None)
    for row in rows:
        status, n = row["status"], row["ref"]
        if status == "missing":
            parts.append(f"{name(n)} skipped")
        elif status == "extra (novel)":
            e = other_exons[row["other"] - 1]
            if not e:
                parts.append("Exon in inserted sequence")
            else:
                inside = lo is not None and lo <= e[0] and e[1] <= hi
                parts.append(f"{'Pseudoexon' if inside else 'Novel exon'} {e[1] - e[0] + 1} nt")
        elif status != "identical":
            if "one exon with reference exon" in status:
                partner = int(status.rsplit(" ", 1)[1].rstrip(")").split(",")[0])
                retained.add(min(n, partner))
                continue
            for piece in status.split(", "):
                what, _, shift = piece.rpartition(" ")
                site = "cryptic acceptor" if what.startswith("start") else "cryptic donor"
                parts.append(f"{name(n)} {'start' if what.startswith('start') else 'end'} {shift} ({site})")
    for k in sorted(retained):
        a, b = name(k), name(k + 1)
        if a.startswith("Exon ") and b == f"Exon {int(a.split()[1]) + 1}":
            parts.append(f"Intron {a.split()[1]} retained")
        else:
            parts.append(f"Intron between {a[0].lower() + a[1:]} and {b[0].lower() + b[1:]} retained")
    if not parts and (not r or r == "r.(=)"):
        return "No change"
    if summary is not None and summary.frame not in ("", "no change"):
        parts.append(summary.frame)
    return " | ".join(parts + [x for x in (r, p) if x and x not in ("r.(=)", "p.(=)")])


def run_analysis(build, gene, compat, scores, fetch, table, settings=None, checkpoint=None, on_progress=None):
    """The whole analysis (see module docstring). scores: a
    refdb.ReferenceScores of the build; fetch(a, b) reads the reference
    chromosome; table: the build's GeneTable (for the reference note).
    on_progress(fraction 0..1, step text). Uses the events chosen
    (apply=True) in compat.alignment.events."""
    t0 = time.time()
    settings = settings or Settings()
    note = table.reference_note(gene)
    if note[0] == "none":
        raise ValueError(note[1])
    say = on_progress or (lambda fraction, text: None)
    say(0.02, tx("Reference values (database or computed live)"))
    wt, plan, source = scores.gene(gene.name, checkpoint=checkpoint,
                                   on_progress=lambda f: say(0.02 + 0.2 * f, tx("Computing the reference live")))
    edits = [haplotype.Edit.from_event(e) for e in compat.alignment.events if e.apply]
    hap, cmap = haplotype.build(compat.ref_seq, compat.region[0], edits)
    say(0.25, tx("Sample haplotype: predicting near the edits"))
    pred = haplotype.predict_patient(scores.predictor(), hap, cmap, plan, wt, checkpoint=checkpoint,
                                     on_progress=lambda d, n: say(0.25 + 0.6 * d / max(n, 1), tx("Sample haplotype")))
    result = Analysis(build, gene, settings, compat, edits, hap, cmap, plan, wt, source, pred, reference_note=note)
    say(0.88, tx("Transcripts, cDNA and protein"))
    _transcripts(result, fetch, settings)
    result.warnings = warnings(result)
    result.seconds = time.time() - t0
    say(1.0, tx("Done"))
    return result


def _transcripts(result, fetch, settings):
    gene, cmap, pred, wt, plan = result.gene, result.cmap, result.patient, result.wt_values, result.plan
    c = gene.canonical
    plus = c.strand == "+"
    anchor = anchors(gene)
    tss, tes = anchor["transcript start"], anchor["transcript end"]
    kw = dict(k=settings.top, threshold=settings.threshold, min_exon=settings.min_exon, min_intron=settings.min_intron)
    wild = splice_paths.top_paths(wt[:, 0], wt[:, 1], plan["start"], c.strand, tss, tes, **kw)
    h_tss, h_tes = cmap.lift(tss, 1 if plus else -1), cmap.lift(tes, -1 if plus else 1)
    patient = splice_paths.top_paths(pred.values[:, 0], pred.values[:, 1], pred.start, c.strand, h_tss, h_tes, **kw)
    describe(result, lambda a, b: fetch(a, b).upper(), wild, patient)


def describe(result, ref_fetch, wild, patient):
    """Fills result.transcripts from the wild type's and the patient's paths
    (splice_paths.Path, or anything with exons, relative and score): cDNA,
    protein, summaries and comparisons. Needs only the reference region and
    the haplotype, so a saved analysis (bundle.py) can be rebuilt."""
    gene, cmap = result.gene, result.cmap
    c = gene.canonical
    plus = c.strand == "+"
    atg = anchors(gene).get("ATG")
    hap_fetch = lambda a, b: result.hap[a - 1:b]
    reference = cdna.reference(ref_fetch, c)
    ref_exons = list(c.exons) if plus else list(reversed(c.exons))
    h_atg = cmap.to_hap(atg) if atg is not None else None

    def compare(t, base, base_exons, number=None):
        rows, line = cdna.concordance(base_exons, t.ref_exons, c.strand)
        return line, rows, cdna.hgvs_r(base[0], t.cdna, number), cdna.hgvs_p(base[0], t.cdna)

    out = []
    ref_t = TranscriptResult("Reference", 0, exons=ref_exons, ref_exons=ref_exons, cdna=reference)
    out.append(ref_t)
    for rank, path in enumerate(wild, 1):
        out.append(TranscriptResult("WT predicted", rank, path.relative, path.score, path.exons, path.exons,
                                    cdna.build(ref_fetch, c.strand, path.exons, atg)))
    for rank, path in enumerate(patient, 1):
        out.append(TranscriptResult("Patient", rank, path.relative, path.score, path.exons,
                                    cdna.to_reference_exons(path.exons, cmap),
                                    cdna.build(hap_fetch, c.strand, path.exons, h_atg)))
    wt_top = out[1]
    for t in out:
        t.summary = cdna.summarize(t.cdna, reference)
        t.vs_ref = compare(t, (reference,), ref_exons)
        t.change_vs_ref = change_text(t.vs_ref[1], t.ref_exons, c.strand, t.summary, t.vs_ref[2], t.vs_ref[3],
                                      ref_exons)
        if t.kind != "Reference":
            # Against the WT prediction, but named and numbered as the canonical transcript is.
            t.vs_wt = compare(t, (wt_top.cdna,), wt_top.ref_exons, cdna.number_on(reference, wt_top.cdna))
            t.change_vs_wt = change_text(t.vs_wt[1], t.ref_exons, c.strand, cdna.summarize(t.cdna, wt_top.cdna),
                                         t.vs_wt[2], t.vs_wt[3], [e for e in wt_top.ref_exons if e],
                                         exon_labels(wt_top.ref_exons, ref_exons))
        if t.kind == "Patient":
            t.max_delta = max_delta(result, t, wt_top)
    result.transcripts = out


def _sites(exons, strand):
    """{(0 acceptor / 1 donor, position)} a transcript splices at: each exon's
    first base (in transcript orientation) an acceptor but the first exon's,
    its last base a donor but the last exon's."""
    out = set()
    for k, (a, b) in enumerate(exons):
        first, last = (a, b) if strand == "+" else (b, a)
        if k:
            out.add((0, first))
        if k < len(exons) - 1:
            out.add((1, last))
    return out


def max_delta(result, t, wt_top):
    """The largest change in SpliceAI's probability (|sample - reference| at
    the same base) over the splice sites where sample candidate t differs from
    the WT prediction wt_top: the sites only t uses and those only the WT
    uses. None if no site differs. A site in inserted sequence has no
    reference base (reference 0); a reference site deleted in the sample has
    sample 0."""
    strand, cmap = result.gene.canonical.strand, result.cmap
    wt, w0 = result.wt_values, result.plan["start"]
    pv, p0 = result.patient.values, result.patient.start
    at = lambda values, start, kind, x: float(values[x - start, kind]) if x is not None and 0 <= x - start < len(
        values) else 0.0
    wt_sites = _sites(wt_top.exons, strand)
    mine = {(kind, h, cmap.to_ref(h)[0]) for kind, h in _sites(t.exons, strand)}
    deltas = [abs(at(pv, p0, kind, h) - at(wt, w0, kind, r)) for kind, h, r in mine if (kind, r) not in wt_sites]
    used = {(kind, r) for kind, _h, r in mine}
    deltas += [abs(at(wt, w0, kind, r) - at(pv, p0, kind, cmap.to_hap(r))) for kind, r in wt_sites
               if (kind, r) not in used]
    return max(deltas) if deltas else None


def evidence(t, settings):
    """The results table's "Evidence vs WT" for transcript t."""
    if t.kind != "Patient":
        return ""
    if t.max_delta is None:
        return "same as the WT prediction"
    text = f"max delta {t.max_delta:.2f}"
    if t.max_delta < settings.weak_delta:
        text = f"weak evidence (delta < {settings.weak_delta:g}); " + text
    return text


def lead(result):
    """(gap, relative of rank 1, of rank 2) of the sample's top two candidates; None with fewer than two."""
    top = [t for t in result.transcripts if t.kind == "Patient"][:2]
    if len(top) < 2:
        return None
    return top[0].relative - top[1].relative, top[0].relative, top[1].relative


def lead_text(result):
    gap = lead(result)
    if gap is None:
        return ""
    return tx("Sample rank 1 leads rank 2 by {gap:.3f} in relative score ({a:.3f} vs {b:.3f}).", gap=gap[0], a=gap[1],
              b=gap[2])


def rna_check(result, t):
    """Spec 15 item 10, for sample candidate t when it differs from the WT
    prediction: where RT-PCR primers could go -- the nearest canonical exons
    on both sides of the change that the reference transcript and t both
    have unchanged -- and the expected product sizes, from the cDNAs. The
    sizes run from the start of the first of those exons to the end of the
    second; primers placed inside them shorten both products by the same
    amount, so the difference stays. None when t is the WT prediction's own
    transcript. No primers are designed."""
    if t.kind != "Patient" or t.max_delta is None:
        return None
    ref = _canonical(result)
    by = {(x.kind, x.rank): x for x in result.transcripts}
    wt, reference = by.get(("WT predicted", 1)), by[("Reference", 0)]
    changed = []
    for r in t.vs_wt[1]:
        if r["status"] == "identical":
            continue
        if r["ref"] is not None and wt is not None:
            changed.append(wt.ref_exons[r["ref"] - 1])
        if r["other"] is not None:
            e = t.ref_exons[r["other"] - 1]
            if e is None:           # an exon in inserted sequence: between the exons around it
                k = r["other"] - 1
                before = next((x for x in reversed(t.ref_exons[:k]) if x), None)
                after = next((x for x in t.ref_exons[k + 1:] if x), None)
                if before and after:
                    e = (min(max(before), max(after)) + 1, max(min(before), min(after)) - 1)
            if e:
                changed.append(e)
    if not changed:
        return None
    lo, hi = min(min(e) for e in changed), max(max(e) for e in changed)
    same = {r["ref"]: r["other"] for r in t.vs_ref[1] if r["ref"] is not None and r["status"] == "identical"}
    plus = result.gene.canonical.strand == "+"
    upstream = [k for k in same if (ref[k - 1][1] < lo if plus else ref[k - 1][0] > hi)]
    downstream = [k for k in same if (ref[k - 1][0] > hi if plus else ref[k - 1][1] < lo)]
    if not upstream or not downstream:
        side = "upstream" if not upstream else "downstream"
        return {"up": max(upstream) if upstream else None, "down": min(downstream) if downstream else None,
                "reference": None, "predicted": None,
                "text": f"no unchanged exon {side} of the change for a primer"}
    up, down = max(upstream), min(downstream)

    def size(cdna_, first, last):
        a, b = cdna_.exons[last]
        return cdna_.exon_starts[last] + (b - a + 1) - cdna_.exon_starts[first]

    size_ref = size(reference.cdna, up - 1, down - 1)
    size_pred = size(t.cdna, same[up] - 1, same[down] - 1)
    return {"up": up, "down": down, "reference": size_ref, "predicted": size_pred,
            "text": f"exon {up} -> exon {down}: reference {size_ref:,} bp, predicted {size_pred:,} bp"}


def _numbers(numbers):
    """'3-6' for consecutive numbers, else '3, 5 and 6'."""
    numbers = sorted(numbers)
    if len(numbers) > 1 and numbers[-1] - numbers[0] == len(numbers) - 1:
        return f"{numbers[0]}-{numbers[-1]}"
    words = [str(n) for n in numbers]
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def _canonical(result):
    c = result.gene.canonical
    return list(c.exons) if c.strand == "+" else list(reversed(c.exons))


def in_insertion(result, t):
    """Sample candidate t's exons wholly inside inserted bases: {exon number: (canonical exon it is a copy of or 0,
    reference position of the insertion, length)}. A copy is an insertion whose bases stand in the reference
    (the nearest place, if several), with the exon's bases on a canonical exon there."""
    out = {}
    cmap, ref_seq, r0 = result.cmap, result.compat.ref_seq, result.compat.region[0]
    canonical = _canonical(result)
    for number, (exon, on_ref) in enumerate(zip(t.exons, t.ref_exons), 1):
        if on_ref is not None:
            continue
        lo, hi = sorted(exon)
        s = cmap.segment_at(lo)
        i = cmap.segments.index(s)
        before = cmap.segments[i - 1] if i else None
        at = before.ref_end + 1 if before and before.orientation else cmap.segments[i + 1].ref_start
        bases = result.hap[s.hap_start - 1:s.hap_end]
        hits, k = [], ref_seq.find(bases) if len(bases) >= align.SOURCE_MIN else -1
        while k >= 0:
            hits.append(r0 + k)
            k = ref_seq.find(bases, k + 1)
        copy = 0
        if hits:
            source = min(hits, key=lambda x: abs(x - at)) + lo - s.hap_start
            copy = next((k for k, (a, b) in enumerate(canonical, 1) if source <= b and a <= source + hi - lo), 0)
        out[number] = (copy, at, hi - lo + 1)
    return out


def intron_of(ref, lo, hi):
    """k if lo..hi lies in intron k (between exons k and k + 1 of ref, in transcript order), else None."""
    return next((k for k in range(1, len(ref)) if min(ref[k - 1][1], ref[k][1]) < lo
                 and hi < max(ref[k - 1][0], ref[k][0])), None)


def wt_regions(result):
    """Where the WT prediction disagrees with the reference, per region of the
    canonical transcript (spec 15 item 9; C2 reads the same): {("exon", k):
    status} for an exon not predicted identically, {("intron", k): what} for
    intron k (between exons k and k + 1) when the WT prediction has an exon
    there that the reference doesn't, or keeps it in one exon."""
    ref = _canonical(result)
    wt = next((t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1), None)
    out = {}
    for r in (wt.vs_ref[1] if wt else []):
        if r["ref"] is not None and r["status"] != "identical":
            out[("exon", r["ref"])] = r["status"]
            m = re.search(r"one exon with reference exon ([\d, ]+)", r["status"])
            for other in (int(x) for x in m.group(1).split(",")) if m else ():
                out[("intron", min(other, r["ref"]))] = tx("kept in one exon by the WT prediction (retained)")
        elif r["ref"] is None and r["other"] is not None:
            e = wt.exons[r["other"] - 1]
            k = intron_of(ref, min(e), max(e))
            if k:
                out[("intron", k)] = tx("an exon of {n:,} nt in the WT prediction only", n=max(e) - min(e) + 1)
    return out


def region_sites(result, kind, k):
    """SpliceAI's WT probability at the reference's own splice sites bounding a
    region of wt_regions(): exon k's acceptor and donor, or intron k's donor
    (exon k) and acceptor (exon k + 1). The transcript's start and end are not
    splice sites. [(name, probability)], e.g. [("acceptor of exon 7", 0.31)]."""
    ref, strand = _canonical(result), result.gene.canonical.strand
    wt, w0 = result.wt_values, result.plan["start"]
    wanted = [("acceptor", k), ("donor", k)] if kind == "exon" else [("donor", k), ("acceptor", k + 1)]
    out = []
    for site, j in wanted:
        if not 1 <= j <= len(ref) or (site, j) in (("acceptor", 1), ("donor", len(ref))):
            continue
        a, b = ref[j - 1]
        x = (a if strand == "+" else b) if site == "acceptor" else (b if strand == "+" else a)
        if 0 <= x - w0 < len(wt):
            name = tx("acceptor of exon {n}", n=j) if site == "acceptor" else tx("donor of exon {n}", n=j)
            out.append((name, float(wt[x - w0, 0 if site == "acceptor" else 1])))
    return out


def locate(result, lo, hi):
    """(where, first exon, last exon) of the canonical exons around reference
    positions lo..hi (C2): the two exons on both sides of the intron it lies
    in plus one more on each side; in an exon, that exon, its neighbours and
    one more on each side; before or after the transcript, its first or last
    two exons."""
    c, ref = result.gene.canonical, _canonical(result)
    n = len(ref)
    hit = [k for k, (a, b) in enumerate(ref, 1) if a <= hi and lo <= b]
    intron = intron_of(ref, lo, hi)
    if hit:
        where, first, last = tx("exon {n}", n=_numbers(hit)), hit[0] - 2, hit[-1] + 2
    elif intron:
        where, first, last = tx("intron {n}", n=intron), intron - 1, intron + 2
    else:
        before = (hi < c.start) if c.strand == "+" else (lo > c.end)
        where = tx("before the transcript") if before else tx("after the transcript")
        first, last = (1, 2) if before else (n - 1, n)
    return where, max(1, first), min(n, last)


def local_concordance(result):
    """Spec 14 C2: per applied edit, whether the WT prediction matches the
    reference for the canonical exons around it (locate()) and the introns
    between them (wt_regions()). [dict(edit, where, exons, differing: {exon:
    status}, introns: {intron: what}, text)]."""
    regions = wt_regions(result)
    out = []
    for e in result.edits:
        where, first, last = locate(result, min(e.start, e.end), max(e.start, e.end))
        exons = list(range(first, last + 1))
        differing = {k: regions[("exon", k)] for k in exons if ("exon", k) in regions}
        introns = {k: regions[("intron", k)] for k in range(first, last) if ("intron", k) in regions}
        size = max(e.end - e.start + 1, len(e.alt))
        if not differing and not introns:
            text = tx("Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction matches the reference "
                      "for exons {exons} ({n} of {n}).", kind=e.kind, pos=e.start, size=size, where=where,
                      exons=_numbers(exons), n=len(exons))
        else:
            same = [k for k in exons if k not in differing]
            which = tx_join("; ", [tx("exon {k} ({status})", k=k, status=s) for k, s in differing.items()]
                            + [tx("intron {k} ({status})", k=k, status=s) for k, s in introns.items()])
            text = (tx("Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction differs from the "
                       "reference for {which}, and matches it for exons {same} ({n_same} of {n} exons identical).",
                       kind=e.kind, pos=e.start, size=size, where=where, which=which, same=_numbers(same),
                       n_same=len(same), n=len(exons)) if same else
                    tx("Around the {kind} at {pos:,} ({size:,} bp, {where}): the WT prediction differs from the "
                       "reference for {which} (0 of {n} exons identical).", kind=e.kind, pos=e.start, size=size,
                       where=where, which=which, n=len(exons)))
        out.append({"edit": e, "where": where, "exons": exons, "differing": differing, "introns": introns,
                    "text": text})
    return out


COLUMNS = ["Type", "Rank", "Relative score", "Evidence vs WT", "Exons", "cDNA length", "Frame", "PTC/NMD",
           "Protein length", "RNA check (RT-PCR)", "Change vs Reference", "Change vs WT predicted"]


def table_rows(result):
    """The results table (spec 8.1), one dict per transcript, as the tab, the
    CSV export and the bundle's summary show it."""
    rows = []
    kind, note = result.reference_note
    accession = result.gene.canonical.accession if result.gene.canonical else ""
    for t in result.transcripts:
        s = t.summary
        ptc = ""
        if s.ptc:
            ptc = f"PTC codon {s.ptc_codon}" + (f"; NMD {s.nmd.split(':')[0]}" if s.nmd else "")
        elif s.nmd:
            ptc = f"NMD {s.nmd.split(':')[0]}"
        row = {"Type": "Sample" if t.kind == "Patient" else t.kind, "Rank": t.rank or "",
               "Relative score": "" if t.relative is None else f"{t.relative:.3f}",
               "Evidence vs WT": evidence(t, result.settings), "Exons": len(t.exons), "cDNA length": len(t.cdna.sequence), "Frame": s.frame, "PTC/NMD": ptc,
               "Protein length": len(t.cdna.protein) if t.cdna.cds_start is not None else "",
               "RNA check (RT-PCR)": (rna_check(result, t) or {}).get("text", ""),
               "Change vs Reference": t.change_vs_ref, "Change vs WT predicted": t.change_vs_wt}
        if t.kind == "Reference" and kind == "fallback":
            row["Type"] = f"Reference ({accession})"
        rows.append(row)
    return rows


# Spec 15 item 5: every notice is critical (it can stop Run), a caution (shown as "?", never stops Run) or
# general information that is always true -- the last is no longer a warning: it is in the tab's introduction
# and help (spelled out in spliceai_gui/cdna_help.py).
CRITICAL, CAUTION = "critical", "caution"


def check_notices(compat):
    """What the check found about this input, by severity: [(level, text,
    confirmable)]. Critical problems stop Run; a confirmable one (a change
    touching the ATG or the transcript's start or end) only until "Run
    anyway" is ticked. Cautions never stop it."""
    checked = chosen_check(compat)
    out = [(CRITICAL, p, False) for p in compat.report["problems"] + checked["problems"]]
    out += [(CRITICAL, w, True) for w in checked["warnings"]]
    out += [(CAUTION, w, False) for w in compat.report["warnings"]]
    return out


def result_notices(result):
    """[(level, text, rows)] for the results: rows are the (kind, rank) of the
    table rows each notice belongs next to."""
    samples = [(t.kind, t.rank) for t in result.transcripts if t.kind == "Patient"]
    out = [(level, text, samples) for level, text, _ in check_notices(result.compat)]
    if len(result.transcripts) > 1:
        rows = result.transcripts[1].vs_ref[1]
        ref_rows = [r for r in rows if r["ref"] is not None]
        same = sum(r["status"] == "identical" for r in ref_rows)
        if ref_rows and same / len(ref_rows) < result.settings.low_concordance:
            out.append((CAUTION, tx("Low wild-type concordance: only {same} of {n} reference exons are predicted "
                                    "identically for the reference sequence. Predictions for the sample are less "
                                    "reliable for this gene.", same=same, n=len(ref_rows)), [("WT predicted", 1)]))
    for d in local_concordance(result):
        if d["differing"] or d["introns"]:
            out.append((CAUTION, tx("{text} Interpret the sample prediction there with caution.", text=d["text"]),
                        samples))
    for t in result.transcripts:
        if t.kind == "Patient" and "weak evidence" in evidence(t, result.settings):
            out.append((CAUTION, tx("Weak evidence: every splice site where sample candidate {rank} differs from the "
                                    "WT prediction changed by less than {delta:g}.", rank=t.rank,
                                    delta=result.settings.weak_delta),
                        [(t.kind, t.rank)]))
    return out


def warnings(result):
    """The notices' texts, kept with the results (bundle.py)."""
    return [text for _, text, _ in result_notices(result)]

# --- the applied events, described (spec 16 items 1 and 2) -------------------------------------------------

ABBREVIATION = {"deletion": "DEL", "insertion": "INS", "duplication": "DUP", "inversion": "INV", "complex": "DELINS",
                "SNV": "SNV"}


def _region_at(result, x):
    """(kind, k, utr) of reference position x on the canonical transcript: ("exon", k, "5'UTR" / "3'UTR" /
    ""), ("intron", k, ""), ("upstream" / "downstream", 0, "")."""
    c, ref = result.gene.canonical, _canonical(result)
    k = next((k for k, (a, b) in enumerate(ref, 1) if a <= x <= b), None)
    if k:
        c_pos = cdna.genomic_to_c(result.transcripts[0].cdna, x) or ""
        return "exon", k, "5'UTR" if c_pos.startswith("-") else "3'UTR" if c_pos.startswith("*") else ""
    k = intron_of(ref, x, x)
    if k:
        return "intron", k, ""
    return ("upstream" if (x < c.start if c.strand == "+" else x > c.end) else "downstream"), 0, ""


def _region_name(r):
    kind, k, utr = r
    if kind == "exon":
        return tx("{utr} in exon {n}", utr=utr, n=k) if utr else tx("exon {n}", n=k)
    if kind == "intron":
        return tx("intron {n}", n=k)
    return tx("upstream of the transcript") if kind == "upstream" else tx("downstream of the transcript")


def event_location(result, e):
    """(where, c. range, nearest sites) of an applied event on the canonical transcript, e.g. ("intron 13",
    "c.2000+50_c.2000+182", "donor of exon 13: 50 bp; acceptor of exon 14: 1,180 bp")."""
    c, ref = result.gene.canonical, _canonical(result)
    lo, hi = (e.end, e.start) if e.end < e.start else (e.start, e.end)       # an insertion: the two bases around it
    first, last = (lo, hi) if c.strand == "+" else (hi, lo)                 # in transcript order
    a, b = _region_at(result, first), _region_at(result, last)
    hit = [k for k, (x, y) in enumerate(ref, 1) if x <= hi and lo <= y]
    if a == b or (a[:2] == b[:2] and a[0] == "exon"):
        where = _region_name(a) if a == b else tx("exon {n}", n=a[1])
    elif e.end < e.start:
        where = tx("between {a} and {b}", a=_region_name(a), b=_region_name(b))
    elif a[0] == b[0] == "exon":
        where = tx("spanning exons {a} to {b}", a=a[1], b=b[1])
    elif not hit:
        where = tx("spanning {a} to {b}", a=_region_name(a), b=_region_name(b))
    else:
        where = (tx("spanning {a} to {b} (exons {n})", a=_region_name(a), b=_region_name(b), n=_numbers(hit))
                 if len(hit) > 1 else
                 tx("spanning {a} to {b} (exon {n})", a=_region_name(a), b=_region_name(b), n=_numbers(hit)))
    ref_cdna = result.transcripts[0].cdna
    c1, c2 = cdna.genomic_to_c(ref_cdna, first), cdna.genomic_to_c(ref_cdna, last)
    c_range = "" if c1 is None or c2 is None else f"c.{c1}" if c1 == c2 else f"c.{c1}_c.{c2}"
    distance = lambda s: 0 if lo <= s <= hi else min(abs(s - lo), abs(s - hi))
    sites = []
    for what, ks in (("donor", range(1, len(ref))), ("acceptor", range(2, len(ref) + 1))):
        # In transcript order a donor is an exon's last base, an acceptor its first.
        at = lambda k: (ref[k - 1][1] if c.strand == "+" else ref[k - 1][0]) if what == "donor" else \
            (ref[k - 1][0] if c.strand == "+" else ref[k - 1][1])
        k = min(ks, key=lambda k: distance(at(k)), default=None)
        if k is not None:
            d = distance(at(k))
            site = tx("donor of exon {n}", n=k) if what == "donor" else tx("acceptor of exon {n}", n=k)
            sites.append(tx("{site}: at it", site=site) if d == 0 else tx("{site}: {d:,} bp", site=site, d=d))
    return where, c_range, tx_join("; ", sites)


def event_origin(e, chrom, ref_seq="", offset=0):
    """Spec 16 item 2: where the bases an insertion, duplication, inversion or replacement puts in come from,
    from the aligner's blocks (align.classify): (text, [(start, end, +1 / -1)] the source regions). ref_seq
    (first base at offset): the region the blocks were looked for in, to recognise a short tandem copy."""
    if e.kind in ("SNV", "deletion"):
        return "", []
    if e.kind == "inversion" or e.label.startswith("inversion"):
        return (tx("inverted copy of {region} (in place)", region=f"{chrom}:{e.start:,}-{e.end:,}"),
                [(e.start, e.end, -1)])
    n = len(e.inserted)
    if e.blocks:
        parts, at, sources = [], 0, []
        for q0, q1, r0, r1, o in e.blocks:
            if q0 > at:
                parts.append((at, q0, tx("novel")))
            adjacent = abs(r1 + 1 - e.start) <= align.ADJACENT or abs(r0 - (e.end + 1)) <= align.ADJACENT
            region = f"{chrom}:{r0:,}-{r1:,}"
            parts.append((q0, q1, tx("inverted copy of {region}", region=region) if o < 0 else
                          tx("tandem copy of {region}", region=region) if adjacent else
                          tx("templated from {region}", region=region)))
            sources.append((r0, r1, o))
            at = q1
        if at < n:
            parts.append((at, n, tx("novel")))
        if len(parts) == 1:
            return parts[0][2], sources
        return tx_join("; ", [tx("bases {a}-{b}: {what}", a=a + 1, b=b, what=what) for a, b, what in parts]), sources
    if e.end < e.start and ref_seq:
        # Too short (or too repetitive) for blocks: a copy right before or after the insertion point.
        s = e.start - offset
        for r0, o in ((s - n, 1), (s, 1), (s - n, -1), (s, -1)):
            piece = ref_seq[max(r0, 0):r0 + n] if r0 >= 0 else ""
            if len(piece) == n and (piece if o > 0 else align.reverse_complement(piece)) == e.inserted:
                start = offset + r0
                region = f"{chrom}:{start:,}-{start + n - 1:,}"
                return (tx("tandem copy of {region}", region=region) if o > 0 else
                        tx("inverted copy of {region}", region=region), [(start, start + n - 1, o)])
    if n < align.SOURCE_MIN:
        return tx("too short to trace (under {n} bases)", n=align.SOURCE_MIN), []
    return tx("novel (not found in the gene region)"), []


def short_sequence(seq, head=6, tail=7):
    """Spec 17 item 3: "CGATCG…ACGATCG, 193 bp"; a short sequence whole."""
    return seq if len(seq) <= head + tail + 3 else f"{seq[:head]}…{seq[-tail:]}, {len(seq):,} bp"


def difference_places(compat, gene):
    """(location, origin) of each of compat's differences, before any run: "intron 4, c.375+120" and
    event_origin()'s text; the differences table shows them."""
    c = gene.canonical
    if c is None:
        return [("", "") for _ in compat.alignment.events]
    r0, ref_seq = compat.region[0], compat.ref_seq
    reference = cdna.reference(lambda a, b: ref_seq[a - r0:b - r0 + 1], c)
    stand_in = SimpleNamespace(gene=gene, transcripts=[SimpleNamespace(cdna=reference)])
    out = []
    for e in compat.alignment.events:
        where, c_range, _sites = event_location(stand_in, e)
        out.append((f"{where}, {c_range}" if c_range else where, event_origin(e, "chr" + gene.chrom, ref_seq, r0)[0]))
    return out


def compat_card(compat, gene):
    """Spec 17 item 3: the compatibility summary, [(label, text)]: identity, the covered region with what is
    missing of the transcript, the orientation (information, never a warning) and the number of differences.
    Empty when the sequence doesn't align; the technical details stay in compat.report["lines"]."""
    a = compat.alignment
    if not a.columns:
        return []
    chrom, c = "chr" + gene.chrom, gene.canonical
    span = (c.start, c.end) if c else tuple(gene.mask)
    gaps = [(s, e) for s, e in ((span[0], min(span[1], a.ref_start - 1)), (max(span[0], a.ref_end + 1), span[1]))
            if s <= e]
    region = f"{chrom}:{a.ref_start:,}-{a.ref_end:,}"
    if not gaps:
        covered = tx("{region}; whole transcript covered", region=region)
    else:
        ref = (list(c.exons) if c.strand == "+" else list(reversed(c.exons))) if c else []
        parts = []
        for s, e in gaps:
            hit = [k for k, (x, y) in enumerate(ref, 1) if x <= e and s <= y]
            gap = f"{chrom}:{s:,}-{e:,}"
            parts.append(gap if not hit else tx("{gap} (exons {n})", gap=gap, n=_numbers(hit)) if len(hit) > 1 else
                         tx("{gap} (exon {n})", gap=gap, n=_numbers(hit)))
        covered = (tx("{region}; missing {a} and {b}", region=region, a=parts[0], b=parts[1]) if len(parts) == 2 else
                   tx("{region}; missing {a}", region=region, a=parts[0]))
    strand = "+" if gene.strand == "+" else "−"
    if a.method.startswith("from coordinates"):
        how = tx("the sample was built from the reference with the events given.")
    elif a.reverse == (gene.strand == "-"):
        how = tx("the sequence was given in transcript orientation (5'→3'), as in Ensembl/NCBI gene exports, and is "
                 "used as such.")
    elif gene.strand == "-":
        how = tx("the sequence was given as the reference (+) strand, as in genome browsers, and is used as such.")
    else:
        how = tx("the sequence was given as its reverse complement (3'→5') and is used turned around.")
    orientation = tx("{gene} is on the {strand} strand; {how}", gene=gene.name, strand=strand, how=how)
    return [(tx("Identity"), f"{100 * a.identity:.2f}%"), (tx("Covered region"), covered),
            (tx("Orientation"), orientation), (tx("Differences"), f"{len(a.events):,}")]


def flanking_exons(result, lo, hi):
    """(first, last) canonical exons around reference positions lo..hi: the exons it touches and their
    neighbours, the two around the intron it lies in, or the first or last exon before or after the
    transcript."""
    c, ref = result.gene.canonical, _canonical(result)
    hit = [k for k, (a, b) in enumerate(ref, 1) if a <= hi and lo <= b]
    k = intron_of(ref, lo, hi)
    if hit:
        return max(1, hit[0] - 1), min(len(ref), hit[-1] + 1)
    if k:
        return k, k + 1
    return (1, 1) if (hi < c.start if c.strand == "+" else lo > c.end) else (len(ref), len(ref))


def _local_change(result, t, lo, hi):
    """What transcript t does around reference positions lo..hi (flanking_exons()), as a "Change vs
    Reference" text: "exon 14 skipped", "no change"."""
    first, last = flanking_exons(result, lo, hi)
    ref = _canonical(result)
    g_lo, g_hi = min(min(e) for e in ref[first - 1:last]), max(max(e) for e in ref[first - 1:last])
    near = lambda e: not e or (g_lo <= min(e) and max(e) <= g_hi)
    rows = [r for r in t.vs_ref[1] if (r["ref"] is not None and first <= r["ref"] <= last and r["status"] != "identical")
            or (r["ref"] is None and r["other"] is not None and near(t.ref_exons[r["other"] - 1]))]
    text = change_text(rows, t.ref_exons, result.gene.canonical.strand, base_exons=ref)
    return text[0].lower() + text[1:]


def event_summary(result):
    """One dict per applied event, in reference order: kind, columns-free description for the viewer's
    labels and summary card, the tab, the CSV and the bundle -- "label" ("DEL 133 bp · chr19:... · intron 13
    (c....) · donor of exon 13: 50 bp; ..."), "origin" and its "sources", "effect" (what sample 1 does
    there, and the WT prediction if it differs from the reference there too) and "line" (all of it)."""
    chrom = "chr" + result.gene.chrom
    by_kind = {(t.kind, t.rank): t for t in result.transcripts}
    sample, wt = by_kind.get(("Patient", 1)), by_kind.get(("WT predicted", 1))
    out = []
    for e in (e for e in result.compat.alignment.events if e.apply):
        lo, hi = (e.end, e.start) if e.end < e.start else (e.start, e.end)
        size = max(e.end - e.start + 1, len(e.inserted))
        coords = f"{chrom}:{lo:,}_{hi:,}" if e.end < e.start else (
            f"{chrom}:{lo:,}" if lo == hi else f"{chrom}:{lo:,}-{hi:,}")
        where, c_range, sites = event_location(result, e)
        label = " · ".join(x for x in (f"{ABBREVIATION.get(e.kind, e.kind.upper())} {size:,} bp", coords,
                                       where + (f" ({c_range})" if c_range else ""), sites) if x)
        origin, sources = event_origin(e, chrom, result.compat.ref_seq, result.compat.region[0])
        effect = ""
        if sample is not None:
            mine = _local_change(result, sample, lo, hi)
            there = _local_change(result, wt, lo, hi) if wt is not None else "no change"
            effect = (f"sample 1: {mine}" if there == "no change" else
                      f"sample 1: same as WT predicted 1 there, {mine} vs the reference" if mine == there else
                      f"sample 1: {mine} (WT predicted 1 there: {there})")
        line = " · ".join(x for x in (label, f"origin: {origin}" if origin else "", effect) if x)
        out.append({"event": e, "kind": e.kind, "label": label, "where": where, "c": c_range, "sites": sites,
                    "origin": origin, "sources": sources, "effect": effect, "line": line})
    return out

