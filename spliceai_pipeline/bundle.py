"""The results bundle (spec 8.3): one file (a zip with its own extension)
holding a finished Sequence-to-cDNA analysis, readable by people and by this
program, which restores the table and the viewer from it without SpliceAI,
TensorFlow, a FASTA, the database or the annotation.

Inside, for people: summary.tsv (the results table), concordance.tsv,
transcripts.tsv (every transcript's exons), edits.tsv (the differences found,
and which were applied), reference_scores.tsv and sample_scores.tsv (per
base), cdna.fasta, protein.fasta, alignments.txt. For the program:
analysis.json, region.fa and haplotype.fa (the sequences everything is
rebuilt from), the per-base values exactly as .npy files, and manifest.json
(build, gene, transcripts, parameters, versions, and every file's SHA-256).
"""
import dataclasses
import hashlib
import io
import json
import os
import time
import zipfile

import numpy as np

from . import __version__, align, cdna, cdna_analysis, haplotype, specificity, transcripts

EXTENSION = ".spliceai-cdna"
FORMAT = 1
KIND = "SpliceAI Sequence-to-cDNA results"


class BundleError(ValueError):
    """A file that isn't a results bundle this program can read; worded for the user."""


def _tsv(header, rows):
    out = io.StringIO()
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join("" if v is None else str(v) for v in row) + "\n")
    return out.getvalue()


def _fasta(records):
    return "".join(f">{name}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n"
                   for name, seq in records)


def _npy(array):
    buf = io.BytesIO()
    np.save(buf, np.asarray(array), allow_pickle=False)
    return buf.getvalue()


def _gene_json(gene):
    return dataclasses.asdict(gene)


def _gene_from(d):
    make = lambda t: transcripts.Transcript(t["accession"], t["symbol"], t["chrom"], t["strand"], t["start"], t["end"],
                                          tuple(t["cds"]) if t["cds"] else None, tuple(tuple(e) for e in t["exons"]))
    g = transcripts.Gene(d["name"], d["chrom"], d["strand"], d["start"], d["end"], tuple(tuple(e) for e in d["exons"]))
    g.canonical = make(d["canonical"]) if d["canonical"] else None
    g.rule = d["rule"]
    g.transcripts = [make(t) for t in d["transcripts"]]
    g.excluded = [tuple(x) for x in d["excluded"]]
    g.other_mane = list(d["other_mane"])
    return g


def save(result, path, versions=None):
    """Writes the analysis to `path` (EXTENSION added if missing). versions:
    {"spliceai": ..., "models": {...}, "database": {...}} as known to the caller."""
    if not path.endswith(EXTENSION):
        path += EXTENSION
    g, a, cm = result.gene, result.compat.alignment, result.cmap
    region0 = result.compat.region[0]
    chrom = "chr" + g.chrom
    files = {}
    rows = cdna_analysis.table_rows(result)
    files["summary.tsv"] = _tsv(cdna_analysis.COLUMNS, [[r[c] for c in cdna_analysis.COLUMNS] for r in rows])
    conc = []
    for t in result.transcripts:
        for against, (line, crows, r, p) in (("reference", t.vs_ref), ("WT predicted 1", t.vs_wt)):
            for row in crows:
                conc.append([t.kind, t.rank, against, row["ref"], row["other"], row["status"]])
    files["concordance.tsv"] = _tsv(["type", "rank", "compared with", "reference exon", "this exon", "status"], conc)
    tx_rows = []
    for t in result.transcripts:
        for n, (e, r) in enumerate(zip(t.exons, t.ref_exons), 1):
            tx_rows.append([t.kind, t.rank, "" if t.relative is None else repr(t.relative),
                            "" if t.score is None else repr(t.score), n, e[0], e[1],
                            r[0] if r else "", r[1] if r else ""])
    files["transcripts.tsv"] = _tsv(["type", "rank", "relative score", "score", "exon", "start", "end",
                                     "reference start", "reference end"], tx_rows)
    files["edits.tsv"] = _tsv(["kind", "pos", "ref", "alt", "start", "end", "size", "inserted", "applied", "label",
                               "origin of the inserted bases", "note"],
                              [[e.kind, e.pos, e.ref, e.alt, e.start, e.end, e.size, e.inserted, e.apply, e.label,
                                cdna_analysis.event_origin(e, chrom, result.compat.ref_seq, region0)[0], e.note]
                               for e in a.events])
    plan = result.plan
    wt = result.wt_values
    ref_pos = np.arange(plan["start"], plan["end"] + 1)
    ref_bases = result.compat.ref_seq[plan["start"] - region0:plan["end"] - region0 + 1]
    files["reference_scores.tsv"] = _tsv(["chrom", "position", "base", "acceptor", "donor"],
                                         ([chrom, int(p), b, f"{x:.6f}", f"{y:.6f}"]
                                          for p, b, (x, y) in zip(ref_pos, ref_bases, wt)))
    pv = result.patient.values
    files["sample_scores.tsv"] = _tsv(["haplotype position", "reference position", "base", "acceptor", "donor"],
                                       ([h, cm.name(h, chrom), result.hap[h - 1], f"{x:.6f}", f"{y:.6f}"]
                                        for h, (x, y) in zip(range(result.patient.start, result.patient.end + 1), pv)))
    fasta, protein = [], []
    for t in result.transcripts:
        name = f"{t.kind.replace(' ', '_')}_{t.rank} {g.name}"
        fasta.append((name, t.cdna.sequence))
        if t.cdna.protein:
            protein.append((name, t.cdna.protein))
    files["cdna.fasta"] = _fasta(fasta)
    files["protein.fasta"] = _fasta(protein)
    lines = []
    by_kind = {(t.kind, t.rank): t for t in result.transcripts}
    for x, y in ((("Reference", 0), ("WT predicted", 1)), (("WT predicted", 1), ("Patient", 1))):
        if x in by_kind and y in by_kind:
            for what, sa, sb, protein_ in (("cDNA", by_kind[x].cdna.sequence, by_kind[y].cdna.sequence, False),
                                           ("protein", by_kind[x].cdna.protein, by_kind[y].cdna.protein, True)):
                top, mid, bottom = cdna.pairwise(sa, sb, protein=protein_)
                lines.append(f"# {what}: {x[0]} {x[1] or ''} vs {y[0]} {y[1]}\n{top}\n{mid}\n{bottom}\n")
    files["alignments.txt"] = "\n".join(lines)
    files["region.fa"] = _fasta([(f"{chrom}:{region0}-{result.compat.region[1]} reference", result.compat.ref_seq)])
    files["haplotype.fa"] = _fasta([("sample haplotype", result.hap)])
    files["wt_values.npy"] = _npy(wt)
    files["patient_values.npy"] = _npy(pv)
    state = {
        "build": result.build, "gene": _gene_json(g), "settings": dataclasses.asdict(result.settings),
        "compat": {"sequence": result.compat.sequence, "region": list(result.compat.region),
                   "window": list(result.compat.window), "report": result.compat.report,
                   "anchors": result.compat.anchors,
                   "alignment": {k: v for k, v in dataclasses.asdict(a).items() if k != "events"},
                   "events": [dataclasses.asdict(e) for e in a.events]},
        "edits": [dataclasses.asdict(e) for e in result.edits],
        "segments": [dataclasses.asdict(s) for s in cm.segments],
        "plan": plan, "wt_source": result.wt_source,
        "patient": {"start": result.patient.start, "end": result.patient.end, "mask": list(result.patient.mask),
                    "windows": [list(w) for w in result.patient.windows], "seconds": result.patient.seconds},
        "transcripts": [{"kind": t.kind, "rank": t.rank, "relative": t.relative, "score": t.score,
                         "exons": [list(e) for e in t.exons]} for t in result.transcripts],
        "reference_note": list(result.reference_note), "warnings": result.warnings, "seconds": result.seconds,
        "specificity": [dataclasses.asdict(c) for c in result.specificity],
        "breakpoints": [dataclasses.asdict(b) for b in result.breakpoints],
    }
    files["analysis.json"] = json.dumps(state, indent=1)
    data = {name: (body.encode("utf-8") if isinstance(body, str) else body) for name, body in files.items()}
    c = g.canonical
    manifest = {
        "kind": KIND, "format": FORMAT, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "build": result.build, "gene": g.name, "strand": g.strand, "canonical": c.accession if c else None,
        "transcripts": [t.accession for t in g.transcripts], "parameters": dataclasses.asdict(result.settings),
        # The note of spec 5.3, e.g. which NM_ stands in for a missing MANE Select transcript.
        "reference_note": result.reference_note[1],
        # Decided after Phase 8: each sample candidate's evidence, and how clearly rank 1 leads rank 2.
        "research_use": cdna_analysis.RESEARCH_NOTE,
        "evidence": {f"Sample {t.rank}": {"max_delta": t.max_delta,
                                             "label": cdna_analysis.evidence(t, result.settings)}
                     for t in result.transcripts if t.kind == "Patient"},
        "sample_lead": cdna_analysis.lead_text(result),
        # Spec 15 item 10: where RT-PCR primers could go, and the products' sizes.
        "rna_check": {f"Sample {t.rank}": cdna_analysis.rna_check(result, t)
                      for t in result.transcripts if t.kind == "Patient"},
        "local_concordance": [d["text"] for d in cdna_analysis.local_concordance(result)],
        # Spec 16: each applied event -- type, size, place, nearest sites, origin, what sample 1 does there.
        "events": [d["line"] for d in cdna_analysis.event_summary(result)],
        "sequence_specificity": [c.text() for c in result.specificity],
        "breakpoint_sensitivity": [b.text() for b in result.breakpoints],
        "software_version": __version__, **(versions or {}),
        "files": {name: hashlib.sha256(body).hexdigest() for name, body in data.items()},
    }
    with zipfile.ZipFile(path + ".tmp", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("manifest.json", json.dumps(manifest, indent=1))
        for name, body in data.items():
            z.writestr(name, body)
    os.replace(path + ".tmp", path)
    return path


def load(path):
    """(Analysis, manifest, [warnings]) from a bundle. Checks it is one, that
    every file matches the manifest's SHA-256, and says when it was made by
    another version. Raises BundleError otherwise."""
    try:
        z = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise BundleError(f"{path} isn't a results bundle (not a readable zip file: {exc}).") from exc
    with z:
        try:
            manifest = json.loads(z.read("manifest.json"))
        except (KeyError, ValueError) as exc:
            raise BundleError(f"{path} isn't a results bundle of this program (no readable manifest.json).") from exc
        if manifest.get("kind") != KIND:
            raise BundleError(f"{path} isn't a Sequence-to-cDNA results bundle.")
        if manifest.get("format") != FORMAT:
            raise BundleError(f"{path} has bundle format {manifest.get('format')}; this program reads format {FORMAT}.")
        data = {}
        for name, digest in manifest["files"].items():
            try:
                body = z.read(name)
            except KeyError as exc:
                raise BundleError(f"{path} is incomplete: {name} is missing.") from exc
            if hashlib.sha256(body).hexdigest() != digest:
                raise BundleError(f"{path} is damaged: {name} doesn't match its checksum.")
            data[name] = body
    notes = []
    if manifest.get("software_version") != __version__:
        notes.append(f"Made by version {manifest.get('software_version')} of this program; this is {__version__}.")
    state = json.loads(data["analysis.json"])
    read_fasta = lambda b: "".join(line for line in b.decode().splitlines() if not line.startswith(">"))
    gene = _gene_from(state["gene"])
    c = state["compat"]
    events = [align.Event(**e) for e in c["events"]]
    alignment = align.Alignment(**c["alignment"], events=events)
    ref_seq = read_fasta(data["region.fa"])
    compat = cdna_analysis.Compatibility(c["sequence"], tuple(c["region"]), tuple(c["window"]), ref_seq, alignment,
                                         c["report"], c["anchors"])
    segments = [haplotype.Segment(**s) for s in state["segments"]]
    cmap = haplotype.CoordMap(segments, c["region"][0], c["region"][1])
    p = state["patient"]
    patient = haplotype.PatientPrediction(p["start"], p["end"], np.load(io.BytesIO(data["patient_values.npy"])),
                                          tuple(p["mask"]), [tuple(w) for w in p["windows"]], p["seconds"])
    result = cdna_analysis.Analysis(
        state["build"], gene, cdna_analysis.Settings(**state["settings"]), compat,
        [haplotype.Edit(**e) for e in state["edits"]], read_fasta(data["haplotype.fa"]), cmap, state["plan"],
        np.load(io.BytesIO(data["wt_values.npy"])), state["wt_source"], patient,
        reference_note=tuple(state["reference_note"]), warnings=state["warnings"], seconds=state["seconds"],
        specificity=[specificity.EventControl(**c) for c in state.get("specificity", [])],
        breakpoints=[specificity.Breakpoints(**b) for b in state.get("breakpoints", [])])
    paths = {"WT predicted": [], "Patient": []}
    for t in state["transcripts"]:
        if t["kind"] in paths:
            paths[t["kind"]].append(_SavedPath([tuple(e) for e in t["exons"]], t["relative"], t["score"]))
    if gene.canonical is not None and state["transcripts"]:
        region0 = compat.region[0]
        cdna_analysis.describe(result, lambda a, b: ref_seq[a - region0:b - region0 + 1], paths["WT predicted"],
                               paths["Patient"])
    return result, manifest, notes


@dataclasses.dataclass
class _SavedPath:
    exons: list
    relative: float
    score: float
