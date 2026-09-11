# Third-Party Licenses

This program is an interface to, and a pipeline around, several external
tools and datasets, each licensed separately from this project and from
each other. None of them are stored in this repository. **SpliceAI is never
bundled or redistributed** — not in the repository, the app, or its Windows
installer: it is downloaded on each user's own PC, from PyPI, during Setup
(or later by the app, if Setup couldn't). The Windows installer
(`installer/`) does include SnpEff, the Temurin Java runtime, the NCBI MANE
summary and the Microsoft Visual C++ runtime DLLs, whose licenses permit
redistribution — see each section below.
This file exists to give the fuller licensing picture beyond the short
notice shown in the app's first-launch dialog; it is not legal advice, and
compliance with each license below is the user's own responsibility.

---

## SpliceAI

**What it is:** the deep-learning splice-effect prediction model this
project is built around — both its Python source (`spliceai` on PyPI,
installed via `pip install --no-deps spliceai`, used for live scoring) and,
if the user downloads them, its precomputed genome-wide score files and
trained model weights.

**License:** SpliceAI's source code and its trained models are licensed
**separately**, under different terms:

- **Source code:** [PolyForm Strict License 1.0.0](https://github.com/Illumina/SpliceAI/blob/master/LICENSE).
  This is a noncommercial license — it permits use for personal, academic,
  research, and similar noncommercial purposes, but does **not** permit
  redistributing the software or creating derivative works from it, and
  does not permit commercial use. It also includes a patent-defense
  termination clause (asserting the software infringes a patent you hold
  terminates your license).
- **Trained models / precomputed scores:** [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
  — academic and non-commercial use only. Illumina's own README states this
  directly: *"These annotations are free for academic and not-for-profit
  use; other use requires a commercial license from Illumina, Inc."*

**A note on a discrepancy worth knowing about:** the `spliceai` package
published on PyPI (version 1.3.1, the latest release, installed by this
project's own setup instructions) still declares itself "GPLv3" in its own
package metadata and README badge. That appears to be stale — Illumina's
GitHub repository's actual, current `LICENSE` file states PolyForm Strict
1.0.0 for the source and CC BY-NC 4.0 for the models, with no mention of
GPLv3. This project defers to the license file in Illumina's repository as
authoritative, since it's the more specific and more recently maintained
statement, but users should be aware the PyPI metadata for the exact
version they install doesn't reflect it.

**Citation:** Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et
al. "Predicting Splicing from Primary Sequence with Deep Learning." *Cell*.
2019;176(3):535-548.e24.

**Source:** [github.com/Illumina/SpliceAI](https://github.com/Illumina/SpliceAI)
· Precomputed scores: [basespace.illumina.com/s/otSPW8hnhaZR](https://basespace.illumina.com/s/otSPW8hnhaZR)

**What this project does:** installs the unmodified `spliceai` 1.3.1
package on the user's own PC. The Windows installer downloads it from PyPI
during Setup, checks it against the published SHA-256 and unpacks it into
the app's folder; if that isn't possible, the app offers the same download
later. It is never included in the installer or the app itself. The app
reads the precomputed score files if the user separately downloads them.
Neither the source, the model weights, nor
the precomputed score data are redistributed here — the user obtains each
directly from Illumina, under Illumina's terms above, and is solely
responsible for using them within academic/non-commercial bounds (or
obtaining a commercial license from Illumina if that applies to their use).

---

## This project's relationship to SpliceAI

### In plain terms

This project is one person's work to make SpliceAI easier to install, run,
and use locally on Windows — a bridge around the setup friction (WSL-free
packaging, a GUI, indexed lookup against the precomputed score files), not
a replacement, modification, or new version of SpliceAI itself. It isn't
monetized in any way — no advertising, no paid tier, no sale of the
software or of anything it produces. For the authoritative source on what
SpliceAI actually is and how it may be used, read
[Illumina's SpliceAI repository](https://github.com/Illumina/SpliceAI)
directly rather than taking this project's word for it.

### Formal statement

The following is quoted verbatim from Illumina's own SpliceAI repository —
not paraphrased or summarized:

> "SpliceAI source code is provided under the [PolyForm Strict License
> 1.0.0](https://github.com/Illumina/SpliceAI/blob/master/LICENSE). SpliceAI
> includes several third party packages provided under other open source
> licenses, please see [NOTICE](https://github.com/Illumina/SpliceAI/blob/master/NOTICE)
> for additional details. The trained models used by SpliceAI (located in
> this package at spliceai/models) are provided under the
> [CC BY NC 4.0](https://github.com/Illumina/SpliceAI/blob/master/LICENSE)
> license for academic and non-commercial use; other use requires a
> commercial license from Illumina, Inc. Purchase of AI scores and models
> for commercial use is available at AI_licensing@illumina.com."
>
> — Illumina, Inc., [SpliceAI README](https://github.com/Illumina/SpliceAI/blob/master/README.md)
> (this framing sentence is corroborated by the opening of Illumina's own
> [LICENSE file](https://github.com/Illumina/SpliceAI/blob/master/LICENSE),
> which states the same thing in shorter form)

> "This package annotates genetic variants with their predicted effect on
> splicing, as described in Jaganathan et al, Cell 2019 in press. The
> annotations for all possible substitutions, 1 base insertions, and 1-4
> base deletions within genes are available
> [here](https://basespace.illumina.com/s/otSPW8hnhaZR) for download. These
> annotations are free for academic and not-for-profit use; other use
> requires a commercial license from Illumina, Inc."
>
> — Illumina, Inc., [SpliceAI README](https://github.com/Illumina/SpliceAI/blob/master/README.md)

All rights to SpliceAI — its source code, trained models, and precomputed
annotations — belong to Illumina, Inc. This project claims no ownership
over any part of SpliceAI and holds no license to it beyond what any
individual user obtains directly from Illumina under the terms quoted
above.

What someone does with this tool is their own decision, but it is not a
decision this project can authorize — it is governed entirely by
Illumina's terms quoted above, which this project has no power to expand,
waive, or interpret on a user's behalf. This project's role is to state
those terms clearly, not to enforce or adjudicate them.

---

## SnpEff

**What it is:** a standalone Java tool, used here (optionally) to annotate
each variant with its affected transcript and coding-sequence (HGVS.c)
position.

**License:** [MIT License](https://github.com/pcingola/SnpEff/blob/master/LICENSE.md).
Copyright 2021, Pablo Cingolani.

> Permission is hereby granted, free of charge, to any person obtaining a
> copy of this software and associated documentation files (the
> "Software"), to deal in the Software without restriction, including
> without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to permit
> persons to whom the Software is furnished to do so, subject to the
> following conditions:
>
> The above copyright notice and this permission notice shall be included
> in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
> OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
> NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
> DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
> OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
> USE OR OTHER DEALINGS IN THE SOFTWARE.

**Citation:** Cingolani P, Platts A, Wang le L, et al. "A program for
annotating and predicting the effects of single nucleotide polymorphisms,
SnpEff." *Fly (Austin)*. 2012;6(2):80-92.

**Source:** [github.com/pcingola/SnpEff](https://github.com/pcingola/SnpEff)

**What this project does:** runs SnpEff as an unmodified standalone Java
tool (never vendored in this repository or listed as a Python dependency)
via its own `snpEff.jar` and RefSeq databases (`GRCh37.p13` / `GRCh38.p14`).
The Windows installer includes the unmodified SnpEff distribution (with its
`LICENSE.md`) and both databases, so annotation works on offline PCs; the
MIT license permits this. "Download SnpEff..." can still fetch it
separately.

---

## Eclipse Temurin (Java runtime)

**What it is:** a free, open-source build of OpenJDK from the Eclipse
Adoptium project. SnpEff is a Java program and needs Java 21 or newer.

**License:** [GNU General Public License, version 2, with the Classpath
Exception](https://openjdk.org/legal/gplv2+ce.html) (GPLv2+CE).

**Source:** [adoptium.net](https://adoptium.net) ·
[github.com/adoptium](https://github.com/adoptium)

**What this project does:** the Windows installer includes the official,
unmodified Temurin 21.0.12.1+1 JRE (with its own `legal/` license and
notice files) in a `java` folder inside the app, used solely to run SnpEff.
Corresponding source: [github.com/adoptium/jdk21u](https://github.com/adoptium/jdk21u),
tag `jdk-21.0.12.1+1`. For copies of the app without it, "Download
SnpEff..." downloads the same kind of JRE from Adoptium (checked against
Adoptium's published SHA-256) when the PC has no Java 21+. It is never
installed system-wide and doesn't change PATH.

---

## NCBI MANE (Matched Annotation from NCBI and EBI)

**What it is:** a joint NCBI/EBI reference dataset — one expert-curated,
representative RefSeq/Ensembl transcript pair ("MANE Select") per
protein-coding gene. Used here, if downloaded, to prefer that transcript
when a variant overlaps several, rather than an arbitrary or
computationally-predicted one.

**License:** standard public NCBI reference data. NCBI's own data policy
states: *"NCBI itself places no restrictions on the use or distribution of
the data"* it produces (this applies to MANE's summary file, hosted and
distributed via NCBI's FTP site). No MANE-specific license file was found
beyond this general NCBI policy; MANE is a joint NCBI/EBI product, and EBI
data is generally distributed under similarly open terms, but this project
did not locate a MANE-specific EBI statement to cite separately.

**Source:** [ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/)
· [ncbi.nlm.nih.gov/refseq/MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)

**What this project does:** downloads the current `MANE.GRCh38.vX.X.summary.txt.gz`
file (discovered dynamically from NCBI's own directory listing, not a
hardcoded version, so it doesn't silently go stale) into a local, gitignored
`mane/` folder. Not stored in this repository; the Windows installer
includes `MANE.GRCh38.v1.5.summary.txt.gz` unmodified so it works offline.

---

## Microsoft Visual C++ runtime (Windows installer only)

**What it is:** the Visual C++ 2015–2022 x64 runtime DLLs
(`msvcp140*.dll`, `vcruntime140*.dll`, `concrt140.dll`, `vcomp140.dll`,
`vccorlib140.dll`), needed by TensorFlow, Qt and the bundled Python.

**License:** Microsoft Visual Studio license terms. These files are on
Microsoft's "Distributable Code" list and may be redistributed with an
application ("app-local" deployment).

**What this project does:** the Windows installer copies version 14.44 of
these DLLs next to the app's exe, so the app runs on PCs without the Visual
C++ Redistributable and without admin rights. Nothing is installed
system-wide.

---

## UCSC Genome Browser reference sequence (hg19 / hg38)

**What it is:** the human reference genome FASTA files (`hg19.fa.gz` /
`hg38.fa.gz`), used as the reference sequence for normalization and live
scoring. Downloaded from UCSC's own goldenPath mirror rather than
NCBI/Ensembl directly, since it's a single, simple, pre-built per-build
FASTA download with no separate indexing/assembly step required.

**License:** UCSC's own stated position, from their data-use conditions
page, is that Genome Browser sequence and annotation data are *"freely
available for any use"*, subject to crediting the data contributor when
used in a publication, and to any additional restrictions noted for
specific species/tracks on their Credits page. This is a broad usage grant
in UCSC's own words, not a formal "public domain" declaration — for the
human reference genome specifically, the underlying sequence data itself
(produced by the public Genome Reference Consortium) has long been treated
as unrestricted, and UCSC's redistribution of it carries no additional
restriction beyond the attribution expectation above.

**Source:** [hgdownload.soe.ucsc.edu/goldenPath](https://hgdownload.soe.ucsc.edu/goldenPath/)
· [genome.ucsc.edu/conditions.html](https://genome.ucsc.edu/conditions.html)

**What this project does:** downloads the plain FASTA for the build the
user selects, decompresses it, and builds a `.fai` index (via `pyfaidx`, so
no `samtools` dependency) into a local, gitignored `data/` folder. Not
redistributed in this repository.
