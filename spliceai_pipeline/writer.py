"""Writing the merged score table to CSV/TSV."""
import csv

from .merge import ScoreRow

FIELDNAMES = [
    "chrom", "pos", "ref", "alt", "gene",
    "DS_AG", "DP_AG", "DS_AL", "DP_AL", "DS_DG", "DP_DG", "DS_DL", "DP_DL",
    "max_score", "allele_fraction", "source", "snpeff_transcript", "snpeff_hgvs_c",
]

# (Original comments weren't recoverable from the compiled program; blank and
# comment lines like these keep line numbers where they were, which Python 3.13
# bakes into class definitions.)


# Columns read_rows() requires in a results file; the SnpEff and
# allele-fraction columns are optional there.
REQUIRED_LOAD_COLUMNS = [
    "chrom", "pos", "ref", "alt", "gene",
    "DS_AG", "DP_AG", "DS_AL", "DP_AL", "DS_DG", "DP_DG", "DS_DL", "DP_DL",
    "max_score", "source",
]

_INT_COLUMNS = {"pos", "DP_AG", "DP_AL", "DP_DG", "DP_DL"}
_FLOAT_COLUMNS = {"DS_AG", "DS_AL", "DS_DG", "DS_DL", "max_score", "allele_fraction"}


class ResultsFileError(ValueError):
    """Raised when a file passed to read_rows doesn't look like a results file
    this app (CLI or GUI) produced."""


def _fmt(value):
    return "" if value is None else value


def _write_excel_sep_hint(fh, delimiter):
    """Writes an Excel "sep=<char>" directive as the file's very first line,
    but ONLY for comma-delimited output.

    Why this is needed: fields ARE genuinely comma-separated on disk (csv.writer
    handles that correctly) -- but on any system where Windows' regional "list
    separator" isn't a comma (e.g. Turkish/many European locales use ";", since
    "," is their decimal separator instead), Excel auto-detects the WRONG
    delimiter when a .csv is opened by double-click, and dumps every field into
    column A as one string. This is a well-documented, widely-used Excel-only
    override that fixes that regardless of the machine's locale. Tab-delimited
    output is never ambiguous this way (no locale uses tab as a decimal
    separator), so it's left alone. read_rows() skips this line if present.
    """
    if delimiter == ",":
        fh.write("sep=,\r\n")


def write_rows(rows, out_path, delimiter="\t"):
    with open(out_path, "w", newline="") as fh:
        _write_excel_sep_hint(fh, delimiter)
        writer = csv.writer(fh, delimiter=delimiter)
        writer.writerow(FIELDNAMES)
        for row in rows:
            writer.writerow([
                row.chrom, row.pos, row.ref, row.alt, row.gene,
                _fmt(row.DS_AG), _fmt(row.DP_AG), _fmt(row.DS_AL), _fmt(row.DP_AL),
                _fmt(row.DS_DG), _fmt(row.DP_DG), _fmt(row.DS_DL), _fmt(row.DP_DL),
                _fmt(row.max_score), _fmt(row.allele_fraction), row.source,
                _fmt(row.snpeff_transcript), _fmt(row.snpeff_hgvs_c),
            ])


def delimiter_for_path(path):
    return "," if path.lower().endswith(".csv") else "\t"


def write_rows_ordered(rows, out_path, column_specs, delimiter="\t"):
    """Like write_rows, but with an arbitrary caller-supplied column set/order instead
    of the fixed FIELDNAMES schema -- used by the GUI's CSV export so the file matches
    whatever column arrangement the user currently has in the results table (including
    variant_id, which FIELDNAMES/write_rows deliberately doesn't have -- that's the
    CLI's own fixed schema and is left untouched).

    column_specs: list of (header_label, ScoreRow attribute name) pairs, in write order.
    """
    with open(out_path, "w", newline="") as fh:
        _write_excel_sep_hint(fh, delimiter)
        writer = csv.writer(fh, delimiter=delimiter)
        writer.writerow([label for label, _ in column_specs])
        for row in rows:
            writer.writerow([_fmt(getattr(row, attr)) for _, attr in column_specs])


def _parse_cell(raw, name):
    if name in _INT_COLUMNS:
        return None if raw == "" else int(raw)
    if name in _FLOAT_COLUMNS:
        return None if raw == "" else float(raw)
    return raw


def read_rows(path):
    """Reads a results file previously written by this app (via write_rows or
    write_rows_ordered -- CLI output or a GUI export/download) back into a list of
    ScoreRow. Column order in the file doesn't matter (matched by header label, not
    position) and an extra variant_id column, if present, is ignored (it's always
    re-derived from chrom/pos/ref/alt). Raises ResultsFileError, with a message
    fit to show directly to a user, if the file is empty, missing required
    columns, or has a malformed row -- never silently loads partial/garbage data.
    """
    delimiter = delimiter_for_path(path)
    with open(path, newline="", encoding="utf-8") as fh:
        # Skip the Excel "sep=" hint line (see _write_excel_sep_hint) if this
        # file starts with one; otherwise rewind so the header is read normally.
        pos_before_first_line = fh.tell()
        first_line = fh.readline()
        if not first_line.rstrip("\r\n").startswith("sep="):
            fh.seek(pos_before_first_line)

        reader = csv.reader(fh, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration:
            raise ResultsFileError("File is empty.")

        missing = [c for c in REQUIRED_LOAD_COLUMNS if c not in header]
        if missing:
            raise ResultsFileError(
                "This doesn't look like a results file produced by this app -- missing required column(s): "
                f"{', '.join(missing)}."
            )

        col_index = {name: header.index(name) for name in header}
        rows = []
        for line_num, fields in enumerate(reader, start=2):
            if not fields or fields == [""]:
                continue
            if len(fields) != len(header):
                raise ResultsFileError(
                    f"Line {line_num}: expected {len(header)} column(s), found {len(fields)}."
                )
            try:
                values = {name: _parse_cell(fields[idx], name) for name, idx in col_index.items()}
                rows.append(ScoreRow(
                    chrom=values["chrom"], pos=int(values["pos"]), ref=values["ref"], alt=values["alt"],
                    gene=values["gene"],
                    DS_AG=values["DS_AG"], DP_AG=values["DP_AG"],
                    DS_AL=values["DS_AL"], DP_AL=values["DP_AL"],
                    DS_DG=values["DS_DG"], DP_DG=values["DP_DG"],
                    DS_DL=values["DS_DL"], DP_DL=values["DP_DL"],
                    max_score=values["max_score"], source=values["source"],
                    snpeff_transcript=values.get("snpeff_transcript") or None,
                    snpeff_hgvs_c=values.get("snpeff_hgvs_c") or None,
                    allele_fraction=values.get("allele_fraction"),
                ))
            except (ValueError, KeyError) as exc:
                raise ResultsFileError(f"Line {line_num}: could not parse row ({exc}).") from exc

        if not rows:
            raise ResultsFileError("File has a valid header but no data rows.")
        return rows
