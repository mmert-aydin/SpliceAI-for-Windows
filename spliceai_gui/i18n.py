"""English / Türkçe for the window's texts.

The English text in the code is the key: t("Run") is "Run", or its entry in
tr.TR when the language setting (config "language") is "tr". Texts with
values in them are templates: t("{n} gene(s)", n=3). Only what the window
shows is translated -- written results (CSV/TSV, bundles, reports) stay
English, so files read the same whatever the language.

Texts the pipeline makes (notices, the check's summary, errors) come in
English, made with spliceai_pipeline.texts.tx; t() of such a text finds the
template and values it was made from and fills the Turkish template instead.

The language is read once per run: changing it takes effect at the next
start (module-level help texts are built at import).
"""
from spliceai_pipeline.texts import MADE

from . import config
from .tr import TR

LANGUAGES = {"en": "English", "tr": "Türkçe"}
_language = None


def language():
    global _language
    if _language is None:
        _language = config.load().get("language")
        if _language not in LANGUAGES:
            _language = "en"
    return _language


def t(text, **fields):
    """`text` in the window's language; fields fill its {placeholders}."""
    if language() != "tr":
        return text.format(**fields) if fields else text
    if not fields and text not in TR and text in MADE:
        template, made_with = MADE[text]
        # Values that are themselves made texts (a status inside a notice) are translated too.
        return TR.get(template, template).format(**{k: t(v) if isinstance(v, str) and v != text else v
                                                     for k, v in made_with.items()})
    out = TR.get(text, text)
    # A value the pipeline made (a genotype line, a status) is translated in its place too.
    return out.format(**{k: t(v) if isinstance(v, str) and v in MADE else v for k, v in fields.items()}) \
        if fields else out


def N_(text):
    """Marks a text kept in a constant and translated where it is shown (t(CONSTANT)), so tests/test_i18n.py
    finds it. Returns it unchanged."""
    return text


_qt_translator = None


def install_qt_translation(app):
    """Qt's own texts (Yes/No, the wizard's Next/Back/Cancel, file dialogs) in the window's language, from the
    qtbase_<lang>.qm that comes with PySide6. Call again after the language changes."""
    from PySide6.QtCore import QLibraryInfo, QTranslator

    global _qt_translator
    if _qt_translator is not None:
        app.removeTranslator(_qt_translator)
        _qt_translator = None
    if language() == "en":
        return
    translator = QTranslator(app)
    if translator.load(f"qtbase_{language()}", QLibraryInfo.path(QLibraryInfo.TranslationsPath)):
        app.installTranslator(translator)
        _qt_translator = translator
