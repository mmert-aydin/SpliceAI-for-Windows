"""Dialog for bulk-pasting a gene list, validating it against known gene
symbols, and confirming before it's applied as a results-table filter.

Validated lists can also be saved under a name for reuse across sessions
(and deleted again later) -- see gene_list_storage.py.
"""
from PySide6.QtWidgets import (
    QComboBox, QDialog, QHBoxLayout, QInputDialog, QLabel, QMessageBox, QPlainTextEdit,
    QPushButton, QVBoxLayout,
)

from . import gene_list_storage
from .gene_reference import parse_gene_tokens, validate_genes
from .i18n import t


class GeneImportDialog(QDialog):
    def __init__(self, parent=None, title=None, apply_text=None, genes=None):
        super().__init__(parent)
        self.setWindowTitle(title or t("Import gene list"))
        self.resize(560, 560)
        self._result = None

        layout = QVBoxLayout(self)

        saved_row = QHBoxLayout()
        saved_row.addWidget(QLabel(t("Saved lists:")))
        self.saved_combo = QComboBox()
        saved_row.addWidget(self.saved_combo, stretch=1)
        self.load_saved_button = QPushButton(t("Load"))
        self.load_saved_button.clicked.connect(self._on_load_saved)
        saved_row.addWidget(self.load_saved_button)
        self.save_button = QPushButton(t("Save..."))
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self._on_save_clicked)
        saved_row.addWidget(self.save_button)
        self.delete_saved_button = QPushButton(t("Delete"))
        self.delete_saved_button.clicked.connect(self._on_delete_saved)
        saved_row.addWidget(self.delete_saved_button)
        layout.addLayout(saved_row)

        intro = QLabel(t("Paste gene symbols below, one per line or separated by commas or spaces (for example, copied "
                         "from an HPO phenotype gene list):"))
        intro.setWordWrap(True)         # wraps on a small screen instead of widening the dialog
        layout.addWidget(intro)

        self.paste_box = QPlainTextEdit()
        self.paste_box.setPlaceholderText(t("BRCA1, BRCA2\nANKRD26\nTP53 ..."))
        self.paste_box.textChanged.connect(self._on_text_changed)
        self.paste_box.setMinimumHeight(60)
        layout.addWidget(self.paste_box, stretch=1)

        check_row = QHBoxLayout()
        self.check_button = QPushButton(t("Check"))
        self.check_button.clicked.connect(self._on_check_clicked)
        check_row.addWidget(self.check_button)
        check_row.addStretch(1)
        layout.addLayout(check_row)

        self.summary_label = QLabel("")
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)

        layout.addWidget(QLabel(t("Not recognized (check for typos or non-gene text):")))
        self.invalid_box = QPlainTextEdit()
        self.invalid_box.setReadOnly(True)
        self.invalid_box.setStyleSheet("color: #b00020;")
        self.invalid_box.setMaximumHeight(110)
        self.invalid_box.setMinimumHeight(40)     # two lines: the dialog fits a short screen
        layout.addWidget(self.invalid_box)

        layout.addWidget(QLabel(t("Recognized:")))
        self.valid_box = QPlainTextEdit()
        self.valid_box.setReadOnly(True)
        self.valid_box.setStyleSheet("color: #1a7f37;")
        self.valid_box.setMaximumHeight(110)
        self.valid_box.setMinimumHeight(40)     # two lines: the dialog fits a short screen
        layout.addWidget(self.valid_box)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self.cancel_button = QPushButton(t("Cancel"))
        self.cancel_button.clicked.connect(self.reject)
        button_row.addWidget(self.cancel_button)
        self.apply_button = QPushButton(apply_text or t("Apply Filter"))
        self.apply_button.setEnabled(False)
        self.apply_button.clicked.connect(self.accept)
        button_row.addWidget(self.apply_button)
        layout.addLayout(button_row)

        self._refresh_saved_combo()
        if genes:
            self.paste_box.setPlainText("\n".join(genes))
            self._on_check_clicked()

    def _on_text_changed(self):
        self._result = None
        self.apply_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.summary_label.setText("")
        self.invalid_box.clear()
        self.valid_box.clear()

    def _on_check_clicked(self):
        tokens = parse_gene_tokens(self.paste_box.toPlainText())
        if not tokens:
            self.summary_label.setText(t("No gene symbols detected in the pasted text."))
            self.apply_button.setEnabled(False)
            self.save_button.setEnabled(False)
            return

        try:
            result = validate_genes(tokens)
        except ModuleNotFoundError as exc:
            if exc.name == "spliceai" or (exc.name or "").startswith("spliceai."):
                QMessageBox.critical(
                    self, t("SpliceAI not installed"),
                    t("Checking gene symbols requires SpliceAI to be installed (it provides the reference gene list used for the "
                      "check).\n"
                      "\n"
                      "To install it, click the SpliceAI status button in the main window's top bar.")
                )
            else:
                QMessageBox.critical(self, t("Check failed"), f"{type(exc).__name__}: {exc}")
            return

        self._result = result

        self.summary_label.setText(
            t("{n} entries detected -- {valid} recognized, {invalid} not recognized.", n=len(tokens),
              valid=len(result.valid), invalid=len(result.invalid))
        )
        self.invalid_box.setPlainText("\n".join(result.invalid) if result.invalid else t("(none)"))
        self.valid_box.setPlainText(", ".join(result.valid) if result.valid else t("(none)"))
        self.apply_button.setEnabled(len(result.valid) > 0)
        self.save_button.setEnabled(len(result.valid) > 0)

    def validated_genes(self):
        return list(self._result.valid) if self._result else []

    def _refresh_saved_combo(self, select_name=None):
        self.saved_combo.clear()
        names = sorted(gene_list_storage.load_all().keys())
        self.saved_combo.addItems(names)
        has_any = bool(names)
        self.load_saved_button.setEnabled(has_any)
        self.delete_saved_button.setEnabled(has_any)
        if select_name is not None:
            index = self.saved_combo.findText(select_name)
            if index >= 0:
                self.saved_combo.setCurrentIndex(index)

    def _on_load_saved(self):
        name = self.saved_combo.currentText()
        if not name:
            return
        genes = gene_list_storage.load_all().get(name, [])
        self.paste_box.setPlainText("\n".join(genes))
        self._on_check_clicked()

    def _on_save_clicked(self):
        if not self._result or not self._result.valid:
            return
        name, ok = QInputDialog.getText(self, t("Save gene list"), t("Name for this list:"))
        name = name.strip()
        if not ok or not name:
            return
        if name in gene_list_storage.load_all():
            reply = QMessageBox.question(
                self, t("Overwrite saved list?"),
                t('A saved list named "{name}" already exists. Overwrite it?', name=name)
            )
            if reply != QMessageBox.Yes:
                return
        gene_list_storage.save_list(name, self._result.valid)
        self._refresh_saved_combo(select_name=name)
        QMessageBox.information(self, t("Saved"), t('Saved "{name}" ({n} gene(s)).', name=name, n=len(self._result.valid)))

    def _on_delete_saved(self):
        name = self.saved_combo.currentText()
        if not name:
            return
        reply = QMessageBox.question(self, t("Delete saved list?"), t('Delete saved list "{name}"?', name=name))
        if reply != QMessageBox.Yes:
            return
        gene_list_storage.delete_list(name)
        self._refresh_saved_combo()
