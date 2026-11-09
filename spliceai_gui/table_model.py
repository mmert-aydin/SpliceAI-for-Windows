"""Results table model: displays merge.ScoreRow objects with proper numeric
sorting and combined gene-search + score-threshold filtering.
"""
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRectF, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPen, QPixmap

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
    ("Zygosity", "zygosity", "str"),
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
    "zygosity": (
        "The genotype of the VCF's first sample for this ALT allele, from its GT field: "
        "Het (one copy), Hom (both copies: no WT transcript expected from this locus) or "
        "Hemi (a single-copy call carrying it, e.g. chrX in a male). Blank if the VCF has no GT."
    ),
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


# vcfio.zygosities()' codes -> (label, filled dots, dots drawn, meaning)
ZYGOSITY = {
    "het": ("Het", 1, 2, "Heterozygous: one copy of this ALT"),
    "hom": ("Hom", 2, 2, "Homozygous: both copies carry this ALT"),
    "hemi": ("Hemi", 1, 1, "Hemizygous: a single-copy call carrying this ALT"),
}
_ICONS = {}


def zygosity_icon(code):
    """Two dots for the two copies of the locus, filled where the ALT is: Hom
    both, Het one, Hemi a single dot. Drawn, not a file, so it's sharp at any
    display scaling."""
    if code not in ZYGOSITY:
        return None
    if code not in _ICONS:
        _, filled, drawn, _ = ZYGOSITY[code]
        icon = QIcon()
        for scale in (1, 2):
            pm = QPixmap(22 * scale, 12 * scale)
            pm.fill(Qt.transparent)
            p = QPainter(pm)
            p.setRenderHint(QPainter.Antialiasing, True)
            p.scale(scale, scale)
            color = QColor(theme.TEAL)
            p.setPen(QPen(color, 1.4))
            for i in range(drawn):
                p.setBrush(color if i < filled else Qt.NoBrush)
                p.drawEllipse(QRectF(1.5 + i * 10, 1.5, 9, 9))
            p.end()
            pm.setDevicePixelRatio(scale)
            icon.addPixmap(pm)
        _ICONS[code] = icon
    return _ICONS[code]


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

        if attr == "zygosity":
            if role == Qt.DisplayRole:
                return ZYGOSITY[value][0] if value in ZYGOSITY else (value or "")
            if role == Qt.DecorationRole:
                return zygosity_icon(value)
            if role == Qt.ToolTipRole and value in ZYGOSITY:
                return ZYGOSITY[value][3]

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


class ScoreFilterProxyModel:
    """Combines a gene-name substring search with a max_score >= threshold filter
    over a ScoreTableModel's rows, keeping the source order.

    Plain Python, not a QSortFilterProxyModel: that one calls a Python
    filterAcceptsRow() once per source row and, after invalidateFilter(),
    emits a rowsRemoved/rowsInserted pair per changed range -- on a WGS result
    (millions of rows) adding one gene froze the window for minutes. Here a
    filter change is one list comprehension, and the gene test is decided once
    per distinct gene name (tens of thousands) rather than once per row.
    on_changed is called exactly once per change.

    Sorting is done by main_window on the full row list (see
    main_window._on_sort_indicator_changed); filtering keeps that order.
    """

    def __init__(self, source_model, on_changed=None):
        self._source = source_model
        self._on_changed = on_changed
        self._quick_genes = []
        self._threshold = None
        self._gene_set = None
        self._genes = set()
        self._rows = []
        source_model.modelReset.connect(self._on_source_reset)

    def _on_source_reset(self):
        self._genes = {r.gene or "" for r in self._source._rows}
        self.refilter()

    def set_quick_gene_filters(self, genes):
        """genes: iterable of gene-name substrings the user added one at a time via the
        "Quick gene filter" add-box (OR'd together -- a row matches if its gene contains
        any of them), or empty/None to clear. Independent of set_gene_set()'s bulk exact-match
        import list; the two combine with AND, consistent with how threshold combines with both."""
        self._quick_genes = [g.strip().lower() for g in genes if g and g.strip()] if genes else []
        self.refilter()

    def set_threshold(self, threshold):
        """threshold: float or None to clear."""
        self._threshold = threshold
        self.refilter()

    def set_gene_set(self, genes):
        """genes: iterable of exact gene symbols (matches gene_reference's canonical
        casing, which is also the pipeline's own annotation casing), or None to clear."""
        self._gene_set = set(genes) if genes else None
        self.refilter()

    def refilter(self):
        rows = self._source._rows
        if self._quick_genes or self._gene_set:
            allowed = {
                g for g in self._genes
                if (not self._quick_genes or any(q in g.lower() for q in self._quick_genes))
                and (not self._gene_set or g in self._gene_set)
            }
            rows = [r for r in rows if (r.gene or "") in allowed]
        if self._threshold is not None:
            t = self._threshold
            rows = [r for r in rows if r.max_score is not None and r.max_score >= t]
        self._rows = rows
        if self._on_changed:
            self._on_changed()

    def rowCount(self):
        return len(self._rows)

    def row_at(self, row_index):
        return self._rows[row_index]
