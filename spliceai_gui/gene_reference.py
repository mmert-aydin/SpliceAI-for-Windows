"""Local, offline reference of known gene symbols, for validating pasted gene lists.

Source: the gene-annotation files bundled with the spliceai package itself
(spliceai/annotations/grch37.txt and grch38.txt), rather than a separately
embedded HGNC dump. Two reasons:

1. They're already installed as part of this app's own dependencies -- no
   extra file to ship, no network fetch at runtime.
2. They ARE the exact set of genes this pipeline can ever produce a score
   for (Annotator.get_name_and_strand() only looks up genes present in one
   of these two files). A gene symbol that isn't in either file could never
   appear in a results table anyway, so validating against this list is more
   directly useful here than validating against the ~43,000-symbol HGNC
   master table, which includes many withdrawn symbols, non-coding loci,
   and other entries with no splice annotation at all.

The tradeoff: this is the gene set as of whatever GENCODE/RefSeq snapshot
spliceai's authors built their annotation from, not a continuously updated
HGNC feed -- a very recently renamed gene symbol might not be recognized.
"""
import functools
import importlib.resources
import re

import pandas as pd

_SEPARATOR_RE = re.compile(r"[,;\s]+")


@functools.lru_cache(maxsize=1)
def known_gene_symbols():
    """Union of gene symbols from spliceai's grch37 + grch38 annotation files."""
    symbols = set()
    for build in ("grch37", "grch38"):
        path = importlib.resources.files("spliceai") / "annotations" / f"{build}.txt"
        with importlib.resources.as_file(path) as p:
            df = pd.read_csv(p, sep="\t", usecols=["#NAME"])
        symbols.update(df["#NAME"].astype(str))
    return symbols


def parse_gene_tokens(text):
    """Split pasted text into gene-symbol tokens, tolerating newline/comma/
    semicolon/whitespace separators in any mix."""
    return [t for t in _SEPARATOR_RE.split(text.strip()) if t]


class GeneValidationResult:
    def __init__(self, valid, invalid):
        self.valid = valid
        self.invalid = invalid


def validate_genes(tokens):
    reference = known_gene_symbols()
    upper_to_canonical = {s.upper(): s for s in reference}

    seen_valid, seen_invalid = set(), set()
    valid, invalid = [], []
    for token in tokens:
        canonical = upper_to_canonical.get(token.upper())
        if canonical is not None:
            if canonical not in seen_valid:
                seen_valid.add(canonical)
                valid.append(canonical)
        elif token not in seen_invalid:
            seen_invalid.add(token)
            invalid.append(token)

    return GeneValidationResult(valid=valid, invalid=invalid)
