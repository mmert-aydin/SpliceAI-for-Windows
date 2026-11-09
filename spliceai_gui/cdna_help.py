"""Every explanation of the Sequence-to-cDNA tab, written once (spec 15 item 6):
the first-time introduction is made of these sections, and the tab's "?"
help dialogs show the same texts.

Sections: what the tab does (WHAT_IT_DOES), every parameter (PARAMETERS),
every result term (RESULT_TERMS), the general notes that are always true and
so are no longer shown as warnings (GENERAL_NOTES, spec 15 item 5), how the
tab classes its notices (NOTICES), and how to read the result
(HOW_TO_READ).
"""
import html
from .i18n import N_, t

WHAT_IT_DOES = (
    t('In short: the tab shows, one below the other, the reference exon boundaries, the transcript SpliceAI'
      ' predicts for the normal (WT) sequence and the transcripts it predicts for the sample, so that the e'
      "ffect of a change can be read from their differences.\n\nThis tab predicts how a sample's DNA sequence"
      ' for a single gene may be spliced, and the cDNA and protein that would result. The sample sequence i'
      's entered by pasting it, by loading it from a FASTA file, or by building it from events given as coo'
      'rdinates; the tab then shows three kinds of transcript, one below the other.\n\nThe reference transcri'
      "pt is the gene's annotated canonical transcript (MANE Select on hg38). The WT predicted transcript i"
      's predicted from the reference genome sequence alone: SpliceAI assigns each base a probability of be'
      'ing a splice acceptor and a splice donor, and the strongest sites are chained into transcripts. Only'
      ' the transcript start, the transcript end and the ATG are taken from the annotation; all exon bounda'
      'ries in between are predicted. The transcripts predicted for the sample are obtained by applying the'
      ' same method to the sample sequence; up to five candidates are shown, ranked from the most likely.\n\n'
      'Comparing the WT prediction with the reference shows how well the method reproduces the known transc'
      'ript of this gene. Where the two agree, a difference between the sample prediction and the WT predic'
      "tion is more likely to reflect the sample's change than a limitation of the method. Where the WT pre"
      'diction already differs from the reference, a predicted change in that region should be interpreted '
      'with more caution.')
)

# (name, default, what it means, what raising and lowering it does)
PARAMETERS = [
    (t("Site threshold"), "0.05",
     t("The lowest SpliceAI probability at which a base can be used as a splice site when transcripts are built."),
     t("A higher value keeps only strong sites: fewer and clearer candidates, at the risk of missing a weak but "
       "real exon (if one of its sites falls below the threshold, the exon can only be predicted as skipped). A "
       "lower value admits weaker sites: more candidates, with more noise.")),
    (t("Minimum internal exon"), "3 nt",
     t("The shortest internal exon a predicted transcript may contain. Micro-exons (internal exons of 3 to about 30 "
       "nucleotides) are common in brain and muscle transcripts; exon 3 of TP53, for example, is 22 nucleotides "
       "long. An exon shorter than this value cannot be predicted, even when SpliceAI scores both of its sites "
       "strongly. The default of 3 nt keeps micro-exons within reach: compared with 25 nt, it predicts 686 more "
       "exons correctly across all hg38 genes while adding 88 incorrect ones. SpliceAI tends to score micro-exons "
       "lower than longer exons, so not all of them are found. The first and last exons are not subject to this "
       "limit."),
     t("A higher value removes very short exons, including real micro-exons and a few incorrect ones. Values below "
       "3 nt are of no use, as no annotated hg38 exon is that short.")),
    (t("Minimum intron"), "60 nt",
     t("The shortest intron a predicted transcript may have."),
     t("A higher value excludes short introns and, with them, some cryptic sites close to an exon. A lower value "
       "allows very short stretches to be spliced out, which rarely occurs in reality.")),
    (t("Apply by default from"), "10 bp",
     t("Differences of this size or larger found by Check compatibility, as well as all inversions and "
       "duplications, are selected under Apply. Smaller differences are not selected, because they often represent "
       "sequencing noise or known polymorphisms (these are scored in the VCF tab)."),
     t("A higher value selects fewer differences by default; a lower value selects more. Each difference can still "
       "be selected or cleared individually.")),
    (t("Transcripts per sequence (Advanced)"), "5",
     t("The number of candidate transcripts shown for the WT prediction and for the sample."),
     t("A higher value also shows less likely alternatives; the relative scores are then shared among more "
       "transcripts, so each becomes smaller. A lower value shows only the most likely few.")),
    (t("Weak evidence below delta (Advanced)"), "0.5",
     t("A sample candidate is marked \"weak evidence\" when every splice site at which it differs from the WT "
       "prediction changed by less than this value. Jaganathan et al. (Cell, 2019) used delta-score cut-offs of 0.2 "
       "(high recall: most splice-altering variants pass, with more false positives), 0.5 (their recommended "
       "balance) and 0.8 (high precision: fewer false positives, more variants missed)."),
     t("A higher value marks more candidates as weak evidence; a lower value marks fewer. The ranking is not "
       "affected.")),
    (t("Random stand-ins (Advanced)"), "20",
     t("The number of random replacements the sequence-specificity control uses for each change."),
     t("More stand-ins give a more stable comparison but take longer, since each requires one more prediction "
       "(about a third of a second for a typical gene). Fewer stand-ins are faster but less stable.")),
]

# (term, explanation)
RESULT_TERMS = [
    (t("Relative score"),
     t("Each transcript's score is the log-likelihood of its splice sites under SpliceAI's probabilities: each "
       "candidate site the transcript uses contributes log(p), and each candidate site it ignores contributes log(1 "
       "- p). Skipping a strong exon or retaining an intron is therefore penalised by the strong sites it ignores. "
       "The relative score converts the scores of the transcripts shown into shares that sum to 1; it is a ranking, "
       "not an estimate of how often each transcript is produced. When all candidates carry the same change, they "
       "differ only elsewhere in the gene, and the relative score then reflects those other sites rather than the "
       "certainty of the change. The gap between the sample's rank 1 and rank 2, shown below the results table, "
       "indicates how clearly the top candidate leads.")),
    (t("Evidence vs WT and weak evidence"),
     t("For each sample candidate, the largest change in SpliceAI's probability (sample minus reference at the same "
       "base) at the splice sites where it differs from the WT prediction: sites it uses that the WT prediction "
       "does not, and sites it no longer uses. A candidate is marked \"weak evidence\" when every differing site "
       "changed by less than the weak-evidence setting; its difference from the WT prediction then arises from "
       "small shifts in probability, often at sites the sequence change barely affects, rather than from a strongly "
       "gained or lost site.")),
    (t("Gene WT concordance"),
     t("How many of this gene's canonical exons the WT prediction reproduces exactly, from a genome-wide run with "
       "the default parameters; it is shown on the gene line once the gene is selected. For genes with low "
       "concordance, the method's predictions are less reliable.")),
    (t("Local concordance"),
     t("For each change, whether the WT prediction matches the reference in the surrounding exons: the exons on "
       "both sides of the affected intron, one further exon on each side, and the introns between them. The exon "
       "schematic shades every exon and intron where they differ.")),
    (t("Sequence-specificity control"),
     t("An optional check that runs only when its button is clicked. For each change, the prediction is repeated "
       "with random stand-ins in place of the change: an insertion, duplication or inversion is replaced by a "
       "random sequence of the same length and dinucleotide composition at the same position, and a deletion by a "
       "deletion of the same size at a random position in the same intron. Each change is compared with the sample "
       "without it (the other changes kept; the wild type when it is the only change). The tab reports how often "
       "the top transcript changes with the stand-ins and their largest delta within 5,000 bases of the change, "
       "alongside the values for the real change. If random sequences alter the transcript about as often and as "
       "strongly as the real change, the predicted effect is probably related to the position (for example a weak "
       "site nearby that any insertion would shift) rather than to the particular sequence. If the real change "
       "stands out, its effect is more likely to be specific to its sequence. This is an empirical comparison based "
       "on a small number of random draws; it is not a statistical test and does not give a p-value.")),
    (t("Breakpoint sensitivity"),
     t("An optional check after Run: each breakpoint of each change is shifted by 10 and 50 bases in both "
       "directions, the prediction is repeated, and the result shows whether the top transcript remains the same. A "
       "prediction that changes with such small shifts depends on base-exact breakpoints.")),
    (t("RNA check (RT-PCR)"),
     t("For each sample candidate that differs from the WT prediction: the nearest exons on both sides of the "
       "change that are unchanged in both the reference transcript and the candidate, where RT-PCR primers could be "
       "placed (for example \"exon 7 -> exon 9\"), and the expected product sizes for the reference and the "
       "predicted transcript. Sizes are measured on the cDNAs from the start of the first of these exons to the end "
       "of the second; primers placed within the exons shorten both products by the same amount, so the difference "
       "between them is preserved. The program does not design primers.")),
]

GENERAL_NOTES = [
    t("These are predictions for research and hypothesis generation, not measurements. A predicted cDNA may differ "
      "from the actual one, and any change must be confirmed with RNA (cDNA) from the sample before a clinical "
      "conclusion is drawn."),
    t("SpliceAI was trained and evaluated mainly on substitutions and small insertions and deletions. Mid-size "
      "insertions and deletions, inversions and duplications fall outside its typical use, so predictions for them "
      "are exploratory."),
    t("CNV and SV calls from short-read sequencing often have imprecise breakpoints, whereas SpliceAI requires "
      "base-exact breakpoints: a few bases can create or abolish a splice site. Where possible, use breakpoints "
      "confirmed by long-read or Sanger sequencing, or run the breakpoint sensitivity check."),
    t("SpliceAI can miss some documented deep intronic variants (for example CFTR c.3718-2477C>T, whose 84 bp "
      "cryptic exon it does not predict); the absence of a predicted event therefore does not exclude an effect."),
    t("SpliceAI does not model tissue specificity, the choice of transcription start or polyadenylation site, "
      "regulatory factors or isoform usage. Many genes also undergo physiological alternative splicing, which the "
      "tab does not label."),
]

NOTICES = (
    t("The tab shows a notice only when it concerns the current input and affects how the result should be read. A "
      "\"?\" icon lists cautions (for example a sequence that does not cover the whole gene, or a WT prediction "
      "that differs from the reference near a change); hover over or click the icon to read them. A red \"!\" "
      "marks a critical problem (for example a sequence that does not align, a reference base that does not match, "
      "or a change affecting the ATG or the transcript's start or end); only these can prevent Run. In the "
      "results, the same icons appear next to the rows concerned. Information that always applies is stated once, "
      "here and in the introduction, rather than repeated as a warning each time.")
)

HOW_TO_READ = (
    t("A change in the sample is a stronger hypothesis when four conditions hold together: the WT prediction "
      "matches the reference around the change (local concordance); the delta at the changed sites is at least 0.5 "
      "(it is not marked as weak evidence); the sample's rank 1 clearly leads rank 2; and the sequence-specificity "
      "control shows that random sequences at the same position do not produce the same change.\n"
      "\n"
      "It is a weaker hypothesis when the WT prediction is discordant near the change (for example a cryptic exon "
      "in that intron present only in the WT prediction), when it is marked as weak evidence, or when random "
      "sequences produce the same change.\n"
      "\n"
      "In all cases, confirm the change with RNA (cDNA) from the sample; the RNA check column suggests where "
      "RT-PCR primers could be placed and which product sizes to expect.")
)

SEQUENCE_NOTE = (t("Ideally, the sequence covers the whole gene. Shorter sequences can also be used: regions not covered are "
                   "taken from the reference and treated as normal, and every position always receives 5,000 bases of reference "
                   "context on each side."))


def _parameters_text():
    e = html.escape
    return "".join(f"<p><b>{e(name)}</b> ({t('default')} {e(default)}): {e(meaning)} {e(effect)}</p>"
                   for name, default, meaning, effect in PARAMETERS)


def _terms_html(terms):
    """Terms in bold, each its own paragraph (rich text: QMessageBox shows it as such)."""
    e = html.escape
    return "".join(f"<p><b>{e(term)}</b>: {e(text)}</p>" for term, text in terms)


# Reading the results, in three dialogs: the scores, the optional checks, and how to read the result.
SCORE_TERMS = RESULT_TERMS[:4]
CHECK_TERMS = RESULT_TERMS[4:]


def _notes_text():
    return "\n\n".join(GENERAL_NOTES)


# The "?" help dialogs (plain text), from the same sections.
TAB_HELP_TEXT = WHAT_IT_DOES + "\n\n" + t("Always keep in mind:") + "\n\n" + _notes_text()
PARAMETERS_HELP_TEXT = _parameters_text()
SCORES_HELP_TEXT = _terms_html(SCORE_TERMS)
CHECKS_HELP_TEXT = _terms_html(CHECK_TERMS)
HOW_TO_READ_HELP_TEXT = "".join(f"<p>{html.escape(p)}</p>" for p in HOW_TO_READ.split("\n\n"))
CHECK_HELP_TEXT = (
    t("\"Check compatibility\" places the sequence on the reference around the selected gene (or within the "
      "coordinate window, if one was given) before any prediction is made.\n"
      "\n"
      "The summary shows the identity, the part of the reference covered by the sequence, the part of the "
      "transcript that is missing (taken from the reference and treated as normal), the orientation in which the "
      "sequence was given, and the number of differences. A sequence given in transcript orientation, as exported "
      "by Ensembl and NCBI, is the usual case for a minus-strand gene and requires no action. Details contains the "
      "technical information: the alignment method, the aligned columns, the positions in the sequence, and where "
      "SpliceAI's gene interval differs from the canonical transcript. A notice is shown when the sequence clearly "
      "fits the other genome build better.\n"
      "\n"
      "All differences are listed: substitutions, insertions, deletions, duplications (insertions repeating the "
      "adjacent sequence), inversions and complex rearrangements. Insertions and deletions are shifted to their "
      "leftmost equivalent position, and each difference is labelled by type. Location gives its position in the "
      "canonical transcript, Inserted the inserted bases (shortened in the middle; hover to see all of them, "
      "right-click to copy them), and Origin the source of those bases. Only the differences selected under Apply "
      "are used to build the sample's sequence: by default, those of 10 bp or more and all inversions and "
      "duplications, because smaller differences often represent sequencing noise or known polymorphisms. The "
      "selection can be changed as needed.\n"
      "\n") + NOTICES
)
COORDS_HELP_TEXT = (
    t("Instead of pasting a sequence, the events themselves can be entered, one per line, and the tab builds the "
      "sample's sequence from the reference. Supported events: a deletion (chr17:17228330-17228462, optionally "
      "followed by the deleted bases, which are then checked), a tandem duplication (a copy immediately after the "
      "stretch; write \"inverted tandem duplication\" or add \"inverted\" for a reverse-complemented copy), an "
      "inversion, or an insertion (a position and the inserted bases, inserted after that position, as in VCF). "
      "HGVS g. notation on the chromosome's RefSeq accession is also accepted, e.g. "
      "NC_000017.11:g.17228330_17228462del, ...dup, ...inv, g.100_101insACGT, g.100_105delinsACG or g.100A>G; the "
      "accession must belong to the selected genome build.\n"
      "\n"
      "Coordinates are 1-based and inclusive, as in VCF files and genome browsers. \"Generate sequence\" builds "
      "the sample's sequence for the gene region from the reference and these events and places it in the paste "
      "box (\"Save as FASTA...\" saves it). As the events are known exactly, they go directly into the differences "
      "table; the reference is still checked: every event must lie within the gene region, and any reference bases "
      "given must match the reference genome.\n"
      "\n") + GENERAL_NOTES[2]
)
COORDINATES_HELP_TEXT = (
    t("Optional. A range, such as chr7:117,480,025-117,668,665, restricts the alignment to that part of the "
      "reference, which is faster for a short sequence. A single position, such as chr7:117,559,590, makes the "
      "check search for the sequence around it: first within twice the sequence length plus 10,000 bases on each "
      "side, then across the whole gene if it is not found there. The range found is then entered here.\n"
      "\n"
      "Coordinates are 1-based and inclusive, as in VCF files and genome browsers. If no gene has been entered, "
      "entering coordinates suggests the gene or genes located there.")
)


# Shown first in the introduction while there is no reference score database -- an install from GitHub,
# whose Setup has no database beside it to copy (the USB drive's Setup copies it).
FIRST_STEPS = (
    t("This tab requires two resources that the VCF tab does not: SpliceAI's predictions for the reference genome, "
      "stored in a reference score database (about 1.1 GB per genome build), and the reference genome itself. No "
      "database is installed on this computer yet; it is provided on the program's USB drive and is not "
      "distributed through GitHub.\n"
      "\n"
      "1. The reference genome (hg19.fa, hg38.fa): VCF tab > Settings shows whether it is present; if not, "
      "Advanced settings > Download... obtains it from UCSC (about 1 GB per build, once).\n"
      "\n"
      "2. The database, in one of the following ways:\n"
      "• With the USB drive: click \"From a folder (USB drive)...\" in this tab and select spliceai-refdb\\hg19, "
      "then spliceai-refdb\\hg38 on the drive. Every file is checked before use. (Running Setup from the drive "
      "does the same.)\n"
      "• Without the USB drive: \"Compute it on this computer...\" computes it with SpliceAI. This takes many "
      "hours per genome build and runs in the background; the program remains usable, and the computation resumes "
      "after a restart.\n"
      "• Without a database: with \"Continue without it\", each gene is computed when it is analysed (a minute or "
      "two, kept for later use).\n"
      "\n"
      "Purpose of the database: each analysis compares the sample's sequence with the reference, and predicting "
      "the reference for a whole gene is the slowest step. Once this has been done, each analysis only needs to "
      "predict the sample.")
)


def intro_html(first_steps=False):
    """The first-time introduction (spec 15 item 6), as HTML for a text browser; first_steps puts
    FIRST_STEPS (getting the database) first."""
    e = html.escape
    paras = lambda text: "".join(f"<p>{e(p)}</p>" for p in text.split("\n\n"))
    h = lambda level, text: f"<h{level}>{e(t(text))}</h{level}>"
    out = [h(2, N_("The Sequence-to-cDNA Tab"))]
    if first_steps:
        out += [h(3, N_("Before the First Analysis")), paras(FIRST_STEPS).replace("\n", "<br>")]
    out += [h(3, N_("What It Does")), paras(WHAT_IT_DOES), h(3, N_("The Parameters"))]
    for name, default, meaning, effect in PARAMETERS:
        out.append(f"<p><b>{e(name)}</b> ({t('default')} {e(default)}). {e(meaning)} {e(effect)}</p>")
    out.append(h(3, N_("The Results")))
    for term, text in RESULT_TERMS:
        out.append(f"<p><b>{e(term)}</b>. {e(text)}</p>")
    out += [h(3, N_("Always Keep in Mind")), "".join(f"<p>{e(n)}</p>" for n in GENERAL_NOTES),
            h(3, N_("Notices")), paras(NOTICES), h(3, N_("How to Read the Results")), paras(HOW_TO_READ)]
    return "".join(out)
