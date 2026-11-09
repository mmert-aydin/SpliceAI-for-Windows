# Checks

Plain Python scripts, no test framework: each prints `PASS` / `FAIL` / `SKIP`
lines and exits non-zero if anything failed.

```
powershell -ExecutionPolicy Bypass -File tests\run-tests.ps1     # all of them
.venv\Scripts\python tests\test_pipeline.py                      # or one
```

| Script | Checks | Time |
|---|---|---|
| `test_pipeline.py` | three single regression variants (build, position, REF, ALT and their unrounded masked scores, checked as numbers), normalization (incl. complex indels), blank lines, chromosome names (MT/chrM, unknown, RefSeq-named FASTA), wrong-build stop, SpliceAI skip reasons, SnpEff per gene, ANKRD26 chr10:27326999 T>C = 0.29, the ClinVar demo VCF against its frozen 1.1.1 output (`tests/data`, masked and raw) | ~3 min |
| `test_normalize.py` | two synthetic repeat indels deep in ANKRD26 (hg19: chr10:27298042 ACAGG>A, chr10:27297926 G>GC) in every equivalent form: same sequence, same left-aligned variant, masked max 0.00 | ~1 min |
| `test_predict.py` | the reference score database's prediction (`predict.py`) equals the reference prediction inside SpliceAI's own `get_delta_scores`: both strands, windows across gene ends, overlapping genes, ANKRD26; the extended mask where it matches the VCF tab; chunk size changes nothing | ~1 min |
| `test_refdb_build.py` | transcript matching and masks (synthetic and real), 0-based/1-based conversions, what each gene stores, storage precision, a real build killed mid-way and resumed, adding a gene to a packed chromosome, changed settings refused, keep-awake, the validation's comparison rule and its delta-score check on built genes, model files checked | ~2 min |
| `test_refdb_reader.py` | a build steered the way the window steers it (pause, resume, stop, a gene asked for first on another chromosome, progress.json, genes from a file), the reader (packed and not-yet-packed genes, stretches, by position) and its refusals (other build, other model files, other format, failed verification, a damaged or missing file), the live fallback and its cache, verification against a validation set (and the real hg38 one, if built here), import, the speed benchmark | ~2 min |
| `test_refdb_gui.py` | the database in the window: the wizard (choice, check with benchmark, start), a background build paused, the program closed mid-build and started again, the build continued and verified (against the Phase 2 values, if built here), "Build this gene now", import (also by dropping the folder on the window), top-bar status, cache | ~4 min |
| `test_align_haplotype.py` | reading a pasted sequence and a region; the patient haplotype and its coordinate map (every base both ways, each kind of edit, region edges, 200 random edit sets, the ATG of a minus-strand gene lifted); synthetic patients from the real reference aligned back -- 7 known edits found as such (left-aligned, tandem duplication, exact inversion breakpoints), reverse-complemented, partial and junk-ended pastes, 25 random patients rebuilt exactly, short pastes by dynamic programming, the other build, the report, the chosen events' checks, a 2.2 Mb paste and the slowest short case | ~15 s |
| `test_align_complex.py` | complex rearrangements: 20 synthetic samples per type from the real hg38 reference -- delins, deletion-inversion-deletion, inverted duplication, inverted templated insertion, deletion with a templated insertion, deletion with an inverted duplication, a mosaic of two copies (complex replacement with its blocks), deletions with 1-10 bp microhomology (reported, left-aligned) and with both breakpoints in Alu repeats (the true breakpoint in the range reported); every sample rebuilt exactly (100% required), the share labelled as made reported (90% required per type) | ~25 s |
| `test_patient.py` | the patient's prediction on top of the database: predicting only within 5,000 bases of 7 edits equals predicting the whole haplotype (a 150 kb plus-strand and a 174 kb minus-strand gene; exact reference values and the database's float16 ones), progress, estimate, cancelling, a deleted gene refused; the engine's raw DS/DP at a variant's +/-50 window equal `get_delta_scores`' (40 SNVs, insertions, deletions, random and at splice sites, spec 5.6 rule) | ~1.5 min |
| `test_splice_cdna.py` | transcripts from SpliceAI's probabilities: synthetic canonical, forced skip, forced retention, pseudoexon, no skipping bias, minimum lengths, minus strand, the k best against brute force (60 random cases); real reference cDNAs on both strands (BRCA2 3,418 aa, APP 770 aa), exon skipping in frame and frameshift (r., p., premature stop), ATG lost, no stop codon, NMD (50-55 nt rule, last exon, start-proximal), exon concordance, approximate HGVS (a change against a WT prediction in the reference's numbering), pairwise alignment; the wild type's top path from the database is the canonical transcript (BRCA2, SAMSN1, MRPL39) | ~5 s |
| `test_cdna_analysis.py` | the engine under the Sequence-to-cDNA tab: genes by either symbol, the hg19 genes without a usable reference transcript (BOP1, RHBG), exons named in the canonical numbering, a synthetic patient checked and run (differences, defaults, the results table, warnings), the results bundle saved and loaded (identical table, coordinates and scores; damaged files refused), bedGraph written and read back through tabix | ~20 s |
| `test_wt_concordance.py` | the genome-wide WT concordance: each category of non-matching exon on synthetic genes (first or last exon edge, canonical site below the threshold, alternative exon -- a moved site or a skip another NM_ has -- other), the per-gene table written and read back, the shipped tables of both builds and their place in the `.spec`, hg38 chr21 from the database against the shipped rows | ~10 s |
| `test_coord_events.py` | events given as coordinates: every plain form and HGVS g. (del, dup, inv, ins, delins, substitution) read, the accession's build checked (UCSC chromAlias accessions), the reference checked (deleted bases, a substitution's base, outside the region, another chromosome, overlaps), the sample built equal to the haplotype of the same edits and rebuilt exactly by aligning it back, the breakpoint sensitivity run and its shifted breakpoints | ~30 s |
| `test_cdna_tab.py` | the tab in the real window: two tabs, the "no database" panel, build and gene, a synthetic patient checked and run, the differences table, results, warnings, status, CSV / bundle / bedGraph exports, a bundle dropped onto the window, a gene without an NM_, Cancel, help texts; spec 15: the introduction shown once and reopened from '?', its texts the help texts, the short Gene placeholder, 'Open viewer' in its accent colour (contrast checked), the strand on the gene line, notices by severity ('?' cautions that never block Run, a red '!' for critical problems, 'Run anyway' for a change touching the ATG), icons next to exactly the result rows a notice concerns, the live computation in the status line | ~5 min |
| `test_layout.py` | the Sequence-to-cDNA tab fits 1366 and 1920 px screens at 100%, 125% and 150% display scaling (no sideways scrolling, every button and help icon inside, Advanced open), the genome build on its own row, the Gene placeholder not cut off; laid out by the real Windows platform for its font metrics, never shown on screen | ~10 s |
| `test_viewer.py` | the viewer: the shared columns from the coordinate map (insertion, deletion, inversion, 300 random edit sets), what it shows (tracks, exons, the completeness rule, hover with c. positions, CSV out and back), the window offscreen, a results file opened with no SpliceAI, FASTA or database, and a 2.3 Mb gene opened in under 1 s with panning under 50 ms a frame; spec 16: each applied change described (type, size, place in the canonical transcript with c., nearest sites, origin of inserted bases, what sample 1 does there), the summary card, the stacked view opening on the first change with its flanking exons, the labels above the ruler, Focus, the delta row scaled to the view with its marks, exon numbers where there is room, the orientation above the view, the side panel as sections with merged transcript rows in a splitter whose width is kept, a tandem duplication's origin with Go to origin and its tooltip, the toolbar at 1366 px with 100/125/150% scaling (Windows platform) | ~1.5 min |
| `test_vcf_bridge.py` | from the VCF tab to the Sequence-to-cDNA tab: the variant as an HGVS g. event, its genotype from the VCF (also after normalization), the confirmation's contents, the prefilled tab and the run, the engine's DS / DP against the VCF tab's within the 5.6 tolerance (the three public regression variants and ClinVar demo variants, masked), the header with both and the genotype, a variant in two genes, other variants of the gene applied together, a queue of two saved as one results file each, the right-click menu; what the transcripts show is printed, not checked | ~3 min |
| `test_defaults.py` | on a fresh install SnpEff annotation and MANE Select are on; a choice saved earlier is kept; SnpEff switched off automatically because it wasn't found isn't saved as the user's choice | ~10 s |
| `test_spliceai_source.py` | from source SpliceAI is found in `python-packages` as in the installed program, with the same setup dialog when it's missing, and isn't in `.venv`; the pinned package installed from its file into an empty home and imported from there; `build.ps1`'s message | ~30 s |
| `test_screens.py` | the whole program and the viewer at 1280x720, 1366x768, 1600x900 and 1920x1080 at 100, 125, 150 and 175% (real Windows font metrics, never shown on screen): nothing cut off or scrolling sideways, Advanced open too, the schematic's labels inside it and off each other, tiles at the screen's ratio, text as sharp as a direct drawing; a screenshot of each view to `docs/phase_reports/phase_13_screens/` (`--no-screens` to skip them) | ~8 min |
| `test_build_selection.py` | automatic hg19/hg38 selection from the VCF header, per-build FASTA rows, settings migration | ~15 s |
| `test_downloads.py` | downloads verify certificates through Windows (truststore); a certificate failure explains what to do; live check that NCBI and UCSC answer | ~10 s |
| `test_reference_locator.py` | the reference genomes are found by themselves (which folders, in which order, .fai required); the FASTA rows stay folded away when they are; they open when one is missing | ~15 s |
| `test_drag_drop.py` | dragging files in from Explorer: a VCF anywhere on the window loads, a FASTA goes to its build's row, a MANE file / folder to its field, anything else is explained | ~10 s |
| `test_gui_features.py` | Pause / Resume / End (also during SnpEff), closing during a run, CPU display, MANE folder, SpliceAI check, scrolling | ~2 min |
| `test_review_phase14.py` | the independent review's findings, one section each: c. numbering from the annotated CDS, peaks in the heat strip, the whole ATG and stop codon in the warnings, no c.0, frame inside the CDS, stop-loss extensions, delins scores, exon blocks, event columns, a reversed paste's notes, the schematic's room | ~30 s |
| `test_benchmark.py` | the three levels of `benchmark.py`'s comparison on made-up exon chains, both strands, and the exact intervals | ~5 s |

`benchmark.py` and `benchmark_report.py` are not checks either: they run the
RNA-confirmed cases of `..\benchmark` (never committed) through the engine and
write the workbook.

`scenarios.py` is not one of the checks: it builds the synthetic corpus in
`..\test_scenarios` (never committed) and runs it through the engine, the
sequence-specificity control and the viewer's exports for the report
(`generate`, `check`, `run`, `report`; about half an hour in all).

**Needs**

- the build venv `.venv` (`installer\build.ps1` creates it);
- `hg19.fa` / `hg38.fa` with their `.fai` files in `~\SpliceAI_reference_data`,
  or set `SPLICEAI_REF_DIR`;
- internet once: the pinned SpliceAI wheel is downloaded into `tests\.cache`
  (git-ignored) and used from there -- never installed into `.venv`, which
  `build.ps1` refuses.
- SnpEff checks use `installer\thirdparty` (`installer\fetch-thirdparty.ps1`)
  and are skipped without it; the reference score database checks need its
  `mane` folder, and the RefSeq tables, which are in the repository's `refseq`
  folder (and in `installer\thirdparty\refseq`).

The GUI checks run the real window invisibly (Qt offscreen) with a throw-away
home folder, so your own settings are never touched. All test variants are
generated from the reference genome -- no patient data is needed.

**Your own data (optional, never committed).** Point these at a VCF (plain or
gzipped, whatever its name says) and its known-good output from an earlier run;
`test_pipeline.py` then scores the VCF's records at the positions the
known-good file lists -- keep that file small, it picks the variants out of a
whole exome -- and compares them, masked, as numbers: the unrounded scores
(read from inside SpliceAI's `get_delta_scores`) must be within the reference
database tolerance of the known-good values (2^-12 + 1e-5; plus half a print
step, 0.005, when the known-good file holds printed 2-decimal values), the DP
the same where the score is at least 0.05 and clearly peaked, and the printed
table the unrounded values rounded. A score next to a rounding boundary
(0.995332 prints as 1.00) can't flip the check either way.

```
$env:SPLICEAI_TEST_VCF = "C:\...\sample.vcf"
$env:SPLICEAI_TEST_EXPECTED = "C:\...\sample-known-good.tsv"
$env:SPLICEAI_TEST_BUILD = "hg38"
```

**The Sequence-to-cDNA validation** is a report, not
a check, and isn't in `run-tests.ps1`: `validate_cdna.py` writes the WT
concordance on 30 MANE genes, the literature and the synthetic controls
(default thresholds, nothing tuned) to the file given. It needs the hg38
database in `~\SpliceAI_reference_data\spliceai-refdb`. About 1 minute.

```
.venv\Scripts\python tests\validate_cdna.py validation.md
```

With `SPLICEAI_PATIENT_TSV` (gene, chrom, pos, ref, alt, build per line) and
`SPLICEAI_PATIENT_REPORT` set, it also appends those variants' predicted
transcripts to the report file named -- keep both outside the tracked files.

**The viewer's speed** is a report too: `perf_viewer.py`
opens the viewer on DMD (2.2 Mb) and MRPL39 (22 kb) and measures the time to
the first view and the frame times while panning and zooming, with Windows'
own platform and the window never shown (it needs the hg38 database).

```
.venv\Scripts\python tests\perf_viewer.py now
```

The clean-machine tests of Setup.exe are separate: `installer\test\run-sandbox.ps1`.
