"""Results table model: displays merge.ScoreRow objects with proper numeric
sorting and combined gene-search + score-threshold filtering.
"""
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt

SORT_ROLE = Qt.UserRole + 1

# (header label, ScoreRow attribute, kind) for each results-table column, in
# display order; kind ("str"/"int"/"float") picks formatting and sort handling.
#
# (Original comments weren't recoverable from the compiled program; blank and
# comment lines like these keep line numbers where they were, which Python 3.13
# bakes into class definitions.)

COLUMNS = [
    ("Variant Id", "variant_id", "str"),
    ("Gene", "gene", "str"),
    ("Transcript", "snpeff_transcript", "str"),
    ("Hgvs", "snpeff_hgvs_c", "str"),
    ("Crom", "chrom", "str"),
    ("Position", "pos", "int"),
    ("Alteration", "alteration", "str"),
    ("Allele Fraction", "allele_fraction", "float"),
    ("Max Score", "max_score", "float"),
    ("Ref", "ref", "str"),
    ("Alt", "alt", "str"),
    ("DS AG", "DS_AG", "float"),
    ("DP AG", "DP_AG", "int"),
    ("DS AL", "DS_AL", "float"),
    ("DP AL", "DP_AL", "int"),
    ("DS DG", "DS_DG", "float"),
    ("DP DG", "DP_DG", "int"),
    ("DS DL", "DS_DL", "float"),
    ("DP DL", "DP_DL", "int"),
    ("Source", "source", "str"),
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
