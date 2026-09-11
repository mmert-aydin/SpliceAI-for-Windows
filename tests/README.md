# Checks

Plain Python scripts, no test framework: each prints `PASS` / `FAIL` / `SKIP`
lines and exits non-zero if anything failed.

```
powershell -ExecutionPolicy Bypass -File tests\run-tests.ps1     # all of them
.venv\Scripts\python tests\test_pipeline.py                      # or one
```

| Script | Checks | Time |
|---|---|---|
| `test_pipeline.py` | normalization (incl. complex indels), blank lines, chromosome names (MT/chrM, unknown, RefSeq-named FASTA), wrong-build stop, SpliceAI skip reasons, SnpEff per gene, ANKRD26 chr10:27326999 T>C = 0.29 | ~2 min |
| `test_build_selection.py` | automatic hg19/hg38 selection from the VCF header, per-build FASTA rows, settings migration | ~15 s |
| `test_gui_features.py` | Pause / Resume / End (also during SnpEff), closing during a run, CPU display, MANE folder, SpliceAI check, scrolling | ~2 min |

**Needs**

- the build venv `.venv` (`installer\build.ps1` creates it);
- `hg19.fa` / `hg38.fa` with their `.fai` files in `~\SpliceAI_reference_data`,
  or set `SPLICEAI_REF_DIR`;
- internet once: the pinned SpliceAI wheel is downloaded into `tests\.cache`
  (git-ignored) and used from there -- never installed into `.venv`, which
  `build.ps1` refuses.
- SnpEff checks use `installer\thirdparty` (`installer\fetch-thirdparty.ps1`)
  and are skipped without it.

The GUI checks run the real window invisibly (Qt offscreen) with a throw-away
home folder, so your own settings are never touched. All test variants are
generated from the reference genome -- no patient data is needed.

**Your own data (optional, never committed).** Point these at a VCF and its
known-good output from an earlier run; `test_pipeline.py` then scores the first
20 records and compares them:

```
$env:SPLICEAI_TEST_VCF = "C:\...\sample.vcf"
$env:SPLICEAI_TEST_EXPECTED = "C:\...\sample-known-good.tsv"
$env:SPLICEAI_TEST_BUILD = "hg38"
```

The clean-machine tests of Setup.exe are separate: `installer\test\run-sandbox.ps1`.
