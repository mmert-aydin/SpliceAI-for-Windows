"""Results table model: displays merge.ScoreRow objects with proper numeric
sorting and combined gene-search + score-threshold filtering.
"""
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt

from . import theme

SORT_ROLE = Qt.UserRole + 1

# Shown in a monospaced font: every column whose value is read as digits or as
# an identifier, so decimal points and accessions line up down the column.
MONO_COLUMNS = {
    "variant_id", "snpeff_transcript", "snpeff_hgvs_c", "chrom", "pos", "alteration",
    "ref", "alt", "allele_fraction", "max_score",
    "DS_AG", "DP_AG", "DS_AL", "DP_AL", "DS_DG", "DP_DG", "DS_DL", "DP_DL",
}

# (header label, ScoreRow attribute, kind) for each results-table column, in
# display order; kind ("str"/"int"/"float") picks formatting and sort handling.

# Order: what identifies the variant (gene, transcript, where it sits), then
# the answer (max score), then the eight scores it came from, then where the
# answer came from. Columns can still be dragged into any other order, and
# where they are left is remembered (config "column_order").
COLUMNS = [
    ("Gene", "gene", "str"),
    ("Transcript", "snpeff_transcript", "str"),
    ("Variant Id", "variant_id", "str"),
    ("Region", "snpeff_region", "str"),
    ("Hgvs", "snpeff_hgvs_c", "str"),
    ("Max Score", "max_score", "float"),
    ("DS AG", "DS_AG", "float"),
    ("DP AG", "DP_AG", "int"),
    ("DS AL", "DS_AL", "float"),
    ("DP AL", "DP_AL", "int"),
    ("DS DG", "DS_DG", "float"),
    ("DP DG", "DP_DG", "int"),
    ("DS DL", "DS_DL", "float"),
    ("DP DL", "DP_DL", "int"),
    ("Source", "source", "str"),
    ("Crom", "chrom", "str"),
    ("Position", "pos", "int"),
    ("Alteration", "alteration", "str"),
    ("Ref", "ref", "str"),
    ("Alt", "alt", "str"),
    ("Allele Fraction", "allele_fraction", "float"),
]


def _column_index(attr):
    return next(i for i, (_, a, _) in enumerate(COLUMNS) if a == attr)


# Positions (in COLUMNS) of the columns the GUI refers to by index, looked up
# by attribute name.
VARIANT_ID_COLUMN = _column_index("variant_id")
GENE_COLUMN = _column_index("gene")
MAX_SCORE_COLUMN = _column_index("max_score")




# Header tooltips, keyed by ScoreRow attribute name; columns not listed here
# get no tooltip.
COLUMN_TOOLTIPS = {
    "snpeff_region": (
        "Where in the gene the variant falls, from SnpEff, for the same transcript shown in "
        "the Transcript column: \"Intron 5/26\" is the 5th of 26 introns, \"Exon 3/27\" the 3rd "
        "of 27 exons. \"Upstream\" is within 5 kb before the transcript starts (the promoter "
        "region); \"Downstream\", \"5' UTR\", \"3' UTR\" and \"Intergenic\" as named. Blank "
        "without SnpEff."
    ),
    "allele_fraction": (
        "Allele fraction, read from the input VCF's AF field (INFO or FORMAT/sample) if "
        "present; otherwise computed from AD (allelic depth) as alt / (ref + alt) if "
        "that's present instead. Blank if the VCF has neither."
    ),
}


def sort_key_value(row, attr, kind):
    """The same None-handling ScoreTableModel.data()'s SORT_ROLE branch uses.
    Public (not module-private) because main_window uses it directly as a
    sort key for a plain Python sorted() call on the full row list -- see
    ScoreFilterProxyModel's docstring for why sorting happens there rather
    than through QSortFilterProxyModel's own sort()/lessThan()."""
    value = getattr(row, attr)
    if kind in ("float", "int"):
        return value if value is not None else float("-inf")
    return value or ""


class ScoreTableModel(QAbstractTableModel):
    def __init__(self, rows=None, parent=None):
        super().__init__(parent)
        self._rows = rows or []

    def set_rows(self, rows):
        self.beginResetModel()
        self._rows = rows
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal:
            if role == Qt.DisplayRole:
                return COLUMNS[section][0]
            if role == Qt.ToolTipRole:
                return COLUMN_TOOLTIPS.get(COLUMNS[section][1])
        return None

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        row = self._rows[index.row()]
        _, attr, kind = COLUMNS[index.column()]
        value = getattr(row, attr)

        if role == Qt.DisplayRole:
            if value is None:
                return ""
            if kind == "float":
                return f"{value:.2f}"
            return str(value)

        if role == SORT_ROLE:
            return sort_key_value(row, attr, kind)

        if role == Qt.FontRole and attr in MONO_COLUMNS:
            return theme.mono_font()

        if role == Qt.TextAlignmentRole and kind in ("float", "int"):
            return int(Qt.AlignRight | Qt.AlignVCenter)

        return None

    def row_at(self, row_index):
        return self._rows[row_index]


class ScoreFilterProxyModel(QSortFilterProxyModel):
    """Combines a gene-name substring search with a max_score >= threshold filter.

    Deliberately filter-only -- sorting is NOT done through this proxy's own
    sort()/lessThan() mechanism. QSortFilterProxyModel's built-in sort requires
    one Python-level lessThan() callback per comparison (~n*log(n) of them),
    and with dynamicSortFilter (Qt's default), every later filter change
    re-triggers that full re-sort too. Profiled at ~2 million lessThan() calls
    and >30s for a single filter-clear on a ~121k-row result set with a sort
    active -- consistently too slow for interactive use, no matter how cheap
    each callback body is, purely from the number of Python/C++ round-trips.

    Instead, main_window pre-sorts the full row list directly in Python (see
    table_model.sort_key_value + main_window._on_sort_indicator_changed) --
    a single native Timsort pass, no per-comparison callback overhead -- and
    this proxy is kept permanently at sort column -1 (see
    main_window._build_results_group's `self.proxy_model.sort(-1)`), so it
    only ever filters; filtering an already-sorted source preserves that
    order for free, with no lessThan() calls at all.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._quick_genes = []
        self._threshold = None
        self._gene_set = None
        self.setSortRole(SORT_ROLE)

    def set_quick_gene_filters(self, genes):
        """genes: iterable of gene-name substrings the user added one at a time via the
        "Quick gene filter" add-box (OR'd together -- a row matches if its gene contains
        any of them), or empty/None to clear. Independent of set_gene_set()'s bulk exact-match
        import list; the two combine with AND, consistent with how threshold combines with both."""
        self._quick_genes = [g.strip().lower() for g in genes if g and g.strip()] if genes else []
        self.invalidateFilter()

    def set_threshold(self, threshold):
        """threshold: float or None to clear."""
        self._threshold = threshold
        self.invalidateFilter()

    def set_gene_set(self, genes):
        """genes: iterable of exact gene symbols (matches gene_reference's canonical
        casing, which is also the pipeline's own annotation casing), or None to clear."""
        self._gene_set = set(genes) if genes else None
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        if self._quick_genes or self._gene_set:
            row = self.sourceModel().row_at(source_row)
            gene_value = row.gene or ""
            gene_value_lower = gene_value.lower()

            if self._quick_genes and not any(g in gene_value_lower for g in self._quick_genes):
                return False
            if self._gene_set and gene_value not in self._gene_set:
                return False

        if self._threshold is not None:
            row = self.sourceModel().row_at(source_row)
            score_value = row.max_score
            if score_value is None or score_value < self._threshold:
                return False

        return True
