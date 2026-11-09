## What the tool can and cannot show

**What these results prove, and what they don't.** The truth here is the sequence: each file was made with known changes, so "exact match" proves that the tool finds and places them. There is no truth for splicing in synthetic samples. Nobody has measured what these made-up changes do to RNA, so a predicted transcript here shows what the model says, never that it is right.

**It can show:**

- Where a sequence differs from the reference: the type, place, size and bases of each change, on either strand, in repeats too, and several changes in one sample. An unchanged file gives no differences.
- What SpliceAI predicts for the sequence as given: a lost splice site can give a skipped exon or a nearby cryptic site, an inverted exon can drop out, and a duplicated exon can be read twice.
- How firm the ground is: whether the prediction for the unchanged sequence matches the annotated exons around the change. Where it doesn't, the result says so.
- Whether an effect belongs to the sequence: the specificity control repeats the prediction with random stand-ins of the same length and composition.
- Where to look on RNA: the exons for an RT-PCR around the change, with the expected product sizes.

**It cannot show:**

- That a transcript exists. These are ranked predictions of one model, not measurements. Tissue, expression level and how much of each transcript is made are outside it, and RNA (cDNA) confirmation is required for any clinical conclusion.
- A change the unchanged prediction already misses. Where SpliceAI's top transcript leaves out an annotated exon, damaging that exon can read as "No change": SCN1A exon 5 is an example here, whose unchanged prediction already skips exons 2, 3, 5 and 6.
- Every micro-exon. With the default minimum internal exon of 3 nt these exons can be predicted, and TP53 exon 3 (22 nt) and the 3 nt exons of MYBPC3, RWDD4 and CDH23 tested here all are. But SpliceAI scores some micro-exons too low: across hg38, 597 of the 765 canonical internal exons under 25 nt are predicted. The low minimum can also add a few false short exons.
- That "No change" means no effect. It means the top transcript stays the same. All 24 intronic inversions read "No change", and most deep insertions do, but the model only looks at sequence near the change. It can't see distant effects, changes to RNA structure or gene regulation.
- That an effect is caused by the exact bases. In MLH1 intron 15 the real insertion gives a 114 nt pseudoexon, but 18 of 20 random stand-ins also changed the top transcript, with score changes up to the same size (0.69), so the effect looks tied to putting bases at that place, not to those bases. In TP53 intron 1, 12 of 20 stand-ins changed the top transcript and the real insertion didn't.
- A clear winner every time. Where the first two candidates are close, read them together: PTEN's exon 8 deletion leads by 0.002 in relative score.
