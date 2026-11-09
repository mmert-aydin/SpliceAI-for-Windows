"""From a variant scored in the VCF tab to the Sequence-to-cDNA tab (spec 16
item 4), without Qt.

- event_line(): the variant as a "From coordinates" line, HGVS g. on the
  chromosome's RefSeq accession (the VCF tab's normalized, left-aligned
  variant, its anchor base dropped).
- genotype(): the variant's genotype in the VCF text (first sample), and
  what it means for the other allele.
- engine_scores(): DS_AG / DS_AL / DS_DG / DS_DL and DP as SpliceAI's
  get_delta_scores would give them, from the cDNA engine's reference values
  and the sample's prediction: the same 101-base window, a deleted base
  scored 0 in the sample, an insertion's bases folded into its anchor's
  maximum, and (masked mode) SpliceAI's own masking at its annotated exon
  boundaries (refdb_validate.scores_from).
- compare(): those against the VCF tab's DS / DP under the 5.6 tolerance
  (the VCF tab prints 2 decimals: half a step plus 1e-5), and why they can
  differ where they do.
"""
import re
from dataclasses import dataclass, field

import numpy as np

from . import coord_events, refdb_validate
from .texts import tx

DIST = 50                       # SpliceAI's window, as the VCF tab scores (score.DIST_VAR)
VCF_TOLERANCE = 0.005 + 1e-5    # spec 5.6: half a rounding step of the 2-decimal VCF output, plus 1e-5
EDGE = 5000 + DIST              # get_delta_scores puts N past the gene's ends within this many bases
NAMES = ("AG", "AL", "DG", "DL")
MAX_QUEUE = 10


def _bare(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def accession(build, chrom):
    return next(acc for acc, (b, c) in coord_events.ACCESSIONS.items() if b == build and c == _bare(chrom))


def event_line(build, chrom, pos, ref, alt):
    """The VCF variant as one "From coordinates" line: SNV, insertion, deletion or delins in HGVS g."""
    ref, alt = ref.upper(), alt.upper()
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    r, a, start = ref[k:], alt[k:], pos + k
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    g = f"{accession(build, chrom)}:g."
    if not r and not a:
        raise ValueError(f"{chrom}:{pos} {ref}>{alt} changes nothing.")
    if not r:
        return f"{g}{start - 1}_{start}ins{a}"
    span = f"{start}" if len(r) == 1 else f"{start}_{start + len(r) - 1}"
    if not a:
        return f"{g}{span}del"
    if len(r) == len(a) == 1:
        return f"{g}{start}{r}>{a}"
    return f"{g}{span}delins{a}"


@dataclass
class Genotype:
    gt: str = ""                # as in the VCF, e.g. "0/1"
    sample: str = ""
    samples: int = 0
    meaning: str = ""           # "heterozygous", "homozygous", ...
    other_allele: str = ""      # what the other allele is expected to give

    @property
    def heterozygous(self):
        return self.meaning == "heterozygous"

    def text(self):
        if not self.gt:
            return self.meaning or tx("no genotype in the VCF")
        who = tx("{sample}, the first of {n} samples", sample=self.sample, n=self.samples) if self.samples > 1 \
            else self.sample
        text = (tx("{gt} ({meaning}; {who})", gt=self.gt, meaning=self.meaning, who=who) if who else
                tx("{gt} ({meaning})", gt=self.gt, meaning=self.meaning))
        return tx("{text}. {other}", text=text, other=self.other_allele) if self.other_allele else text


def genotype(vcf_text, chrom, pos, ref, alt, normalize=None):
    """The variant's genotype in the first sample of the VCF text. Matched on chromosome, position, REF and the
    ALT allele; failing that, records within 100 bases are normalized with normalize(chrom, pos, ref, alt) ->
    (pos, ref, alt), as the VCF tab does, and compared again."""
    samples, near = [], []
    for line in vcf_text.splitlines():
        if line.startswith("#CHROM"):
            samples = line.rstrip("\r\n").split("\t")[9:]
            continue
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) < 5 or _bare(cols[0]) != _bare(chrom):
            continue
        try:
            p = int(cols[1])
        except ValueError:
            continue
        alts = cols[4].split(",")
        if p == pos and cols[3].upper() == ref.upper() and alt.upper() in [x.upper() for x in alts]:
            return _genotype(cols, [x.upper() for x in alts].index(alt.upper()) + 1, samples)
        if abs(p - pos) <= 100:
            near.append((p, cols, alts))
    for p, cols, alts in near if normalize else []:
        for i, x in enumerate(alts, 1):
            try:
                if tuple(normalize(cols[0], p, cols[3], x)) == (pos, ref, alt):
                    return _genotype(cols, i, samples)
            except (ValueError, KeyError):
                continue
    return Genotype(meaning=tx("the record wasn't found in the VCF text"))


def _genotype(cols, allele, samples):
    if len(cols) < 10 or "GT" not in cols[8].split(":"):
        return Genotype(meaning=tx("no genotype (the VCF has no GT field)"))
    gt = cols[9].split(":")[cols[8].split(":").index("GT")]
    calls = re.split(r"[/|]", gt)
    mine = calls.count(str(allele))
    others = [c for c in calls if c not in (str(allele), "0", ".")]
    g = Genotype(gt, samples[0] if samples else "", len(samples))
    if all(c == "." for c in calls):
        g.meaning = tx("genotype missing")
    elif len(calls) == 1:
        g.meaning = tx("hemizygous") if mine else tx("not carried")
    elif mine >= 2:
        g.meaning = tx("homozygous")
        g.other_allele = tx("Both alleles carry it: no WT transcript is expected from this locus.")
    elif mine == 1 and not others:
        g.meaning = tx("heterozygous")
        g.other_allele = tx("The other allele is expected to give the WT transcript.")
    elif mine == 1:
        g.meaning = tx("heterozygous, another ALT on the other allele")
        g.other_allele = tx("The other allele carries another ALT of this record (allele {allele}), so it isn't "
                            "expected to give the WT transcript.", allele=others[0])
    else:
        g.meaning = tx("not carried by this sample")
    return g


def engine_scores(result, pos, ref, alt, masked=False, exons=None):
    """(DS [AG, AL, DG, DL], DP, the four delta series) for the VCF variant pos/ref/alt from the analysis
    `result` (its reference values and the sample's prediction), as get_delta_scores computes them.
    exons: SpliceAI's exons of the gene (1-based) for masked mode; default the gene's own."""
    plan, wt, pred, cmap = result.plan, result.wt_values, result.patient, result.cmap
    n = 2 * DIST + 1
    y_ref, y_alt = np.zeros((1, n, 3)), np.zeros((1, n, 3))
    value = lambda h: pred.values[h - pred.start] if h is not None and pred.start <= h <= pred.end else (0.0, 0.0)
    for i, x in enumerate(range(pos - DIST, pos + DIST + 1)):
        if plan["start"] <= x <= plan["end"]:
            y_ref[0, i, 1:] = wt[x - plan["start"]]
        h = cmap.to_hap(x)
        if x == pos and h is None and len(ref) != len(alt) and 1 in (len(ref), len(alt)) and cmap.to_hap(pos - 1):
            h = cmap.to_hap(pos - 1) + 1        # a replaced base: get_delta_scores has the first alternative base here
        if x == pos and len(ref) == 1 and len(alt) > 1 and h is not None:
            # get_delta_scores: an insertion's bases fold into its anchor, as their maximum
            y_alt[0, i, 1:] = np.max([value(h + j) for j in range(len(alt))], axis=0)
        else:
            y_alt[0, i, 1:] = value(h)          # a deleted base: 0
    bounds = np.array(sorted({b for e in (exons or result.gene.exons) for b in e}))
    dist_exon = int(min(bounds - pos, key=abs)) if len(bounds) else 10 ** 9
    return refdb_validate.scores_from(y_ref, y_alt, dist_exon, masked)


@dataclass
class Comparison:
    vcf: tuple                  # (DS [4], DP [4]) from the VCF tab
    engine: tuple               # (DS [4], DP [4]) from the engine
    problems: list = field(default_factory=list)    # where they differ beyond the tolerance
    reasons: list = field(default_factory=list)     # why they can

    @property
    def agree(self):
        return not self.problems

    def lines(self, mode=None):
        """Header text: [VCF tab line, engine line] -- 'DS_AG 0.29 (DP 12)' ...; mode: the VCF tab's, named."""
        out = []
        vcf = tx("VCF tab ({mode})", mode=mode) if mode else tx("VCF tab")
        for name, (ds, dp) in ((vcf, self.vcf), (tx("cDNA engine"), self.engine)):
            out.append(tx("{name}: {scores}", name=name, scores=", ".join(
                f"DS_{k} {'.' if d is None else f'{d:.2f}'} (DP {'.' if p is None else p})"
                for k, d, p in zip(NAMES, ds, dp))))
        return out


def compare(row, engine, mode="raw", gene=None, others_applied=False):
    """The VCF tab's row (merge.ScoreRow) against engine_scores(): problems under the 5.6 tolerance and the
    reasons they can differ (spec 16 item 4)."""
    ds_v = [row.DS_AG, row.DS_AL, row.DS_DG, row.DS_DL]
    dp_v = [row.DP_AG, row.DP_AL, row.DP_DG, row.DP_DL]
    ds_e, dp_e, series = engine
    c = Comparison((ds_v, dp_v), (ds_e, dp_e))
    if any(v is None for v in ds_v):
        c.problems.append(tx("the VCF tab has no score for this variant"))
        c.reasons.append(tx("SpliceAI scores no variant that replaces several bases by several others (a delins), and "
                            "none outside its genes or past its reference; the engine scores it all the same"))
        return c
    c.problems = refdb_validate.compare(ds_e, dp_e, series, ds_v, dp_v, VCF_TOLERANCE)[0]
    if not c.problems:
        return c
    if gene is not None and (row.pos - gene.start < EDGE or gene.end - row.pos < EDGE):
        c.reasons.append(tx("near the gene's edge: within {n:,} bases of SpliceAI's gene start or end, the VCF tab "
                            "(get_delta_scores) replaces the sequence beyond the gene with N, while the engine uses the "
                            "real sequence there (the database's extended mask)", n=EDGE))
    if mode == "masked":
        c.reasons.append(tx("the VCF tab scored in masked mode: a gain at an annotated splice site and a loss "
                            "elsewhere count 0 there; the engine's comparison follows SpliceAI's annotation, which can "
                            "place exon boundaries differently from the canonical transcript shown here"))
    if row.source and "precomputed" in row.source.lower():
        c.reasons.append(tx("the VCF tab's value comes from Illumina's precomputed files (their own annotation and "
                            "rounding), not from a live computation"))
    if others_applied:
        c.reasons.append(tx("other variants were applied on the same haplotype, so the engine's deltas include their "
                            "effect"))
    if not c.reasons:
        c.reasons.append(tx("no known reason (the conventions agree here): please report it"))
    return c


@dataclass
class Request:
    """One variant sent from the VCF tab: what the Sequence-to-cDNA tab builds and runs."""
    build: str
    gene: object                # transcripts.Gene
    row: object                 # merge.ScoreRow of the variant
    genotype: Genotype
    mode: str = "raw"           # the VCF tab's scoring mode for this row
    others: list = field(default_factory=list)      # ScoreRows applied too, together on one haplotype

    @property
    def lines(self):
        return [event_line(self.build, r.chrom, r.pos, r.ref, r.alt) for r in [self.row] + self.others]

    @property
    def name(self):
        r = self.row
        return f"{r.chrom}:{r.pos:,} {r.ref}>{r.alt}"

    def compare(self, result):
        r = self.row
        return compare(r, engine_scores(result, r.pos, r.ref, r.alt, masked=self.mode == "masked"), self.mode,
                       self.gene, bool(self.others))


def summary(request, comparison):
    """What the results bundle keeps of where the analysis came from (its manifest's "from_vcf")."""
    return {"variant": request.name, "gene": request.gene.name, "mode": request.mode,
            "genotype": request.genotype.text(), "events": request.lines,
            "other_variants_applied": [f"{r.chrom}:{r.pos} {r.ref}>{r.alt}" for r in request.others],
            "vcf_tab": comparison.lines()[0], "cdna_engine": comparison.lines()[1],
            "within_tolerance": comparison.agree, "differences": comparison.problems, "reasons": comparison.reasons}
