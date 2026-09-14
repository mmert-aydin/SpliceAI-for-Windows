"""Per-user, per-machine settings persistence.

Stored under the user's home directory (~/.spliceai_gui/config.json), never
inside the project folder -- so this file is never accidentally shared,
committed, or copied along with the app to another machine/user. On first
launch (or on any machine where this file doesn't exist yet), every field
comes back empty/unset -- there are no bundled defaults for FASTA or
precomputed-dir paths.
"""

import json
import os
import re
from pathlib import Path

CONFIG_PATH = Path.home() / ".spliceai_gui" / "config.json"

DEFAULTS = {
    # Older versions had a single FASTA field; its value is moved to the
    # matching per-build key on load (see _migrate) and then left empty.
    "fasta_path": "",
    "fasta_path_hg19": "",
    "fasta_path_hg38": "",
    "precomputed_mode": None,
    "precomputed_dir": "",
    "build": "hg19",
    "mode": "masked",
    "skip_precomputed": False,
    "last_reference_download_dir": "",
    "column_order": [],
    # Set once _migrate has dropped a column order saved before the columns
    # were reordered; see _migrate.
    "column_order_v2_applied": False,
    # Heights of the three draggable panes, as the user last left them.
    "pane_sizes": [],
    "snpeff_dir": "",
    "mane_dir": "",
    # SnpEff comes with Setup, so annotating is the normal case rather than
    # something to opt into. The window turns it off again if SnpEff can't be
    # found, so this default can never block a run.
    "use_snpeff": True,
    # Set once _migrate has turned the above on for a settings file written
    # before SnpEff shipped with the program; see _migrate.
    "snpeff_default_applied": False,
    "hide_reference_note": False,
}


def fasta_key(build):
    """Settings key of the reference FASTA path for one build ("hg19"/"hg38")."""
    return f"fasta_path_{build}"


_NAME_BUILDS = [
    (re.compile(r"hg19|grch37|hs37|b37", re.IGNORECASE), "hg19"),
    (re.compile(r"hg38|grch38|hs38|b38", re.IGNORECASE), "hg38"),
]


def guess_build_from_name(path):
    """"hg19"/"hg38" if a FASTA's file name names exactly one build (e.g.
    hg38.fa, GRCh37_genomic.fna), else None."""
    name = os.path.basename(path or "")
    builds = {build for pattern, build in _NAME_BUILDS if pattern.search(name)}
    return builds.pop() if len(builds) == 1 else None


def _migrate(settings):
    """Brings a settings file written by an older version up to date:

    - an older single "fasta_path" moves into the per-build key it belongs to,
      by its file name, else the build that was selected with it;
    - "use_snpeff" is turned on once. It used to default to off, because SnpEff
      was something the user installed separately; it now comes with Setup, so
      a saved "false" is almost always just that old default rather than a
      decision. Done once and remembered, so turning it off afterwards sticks.
    """
    legacy = settings.get("fasta_path")
    if legacy and not settings[fasta_key("hg19")] and not settings[fasta_key("hg38")]:
        build = guess_build_from_name(legacy) or (settings.get("build") if settings.get("build") in ("hg19", "hg38") else "hg19")
        settings[fasta_key(build)] = legacy
    if not settings.get("snpeff_default_applied"):
        settings["use_snpeff"] = True
        settings["snpeff_default_applied"] = True
    # A column order saved before the columns were reordered pins every column
    # back where it used to be, so the new order would never be seen. Dropped
    # once; any order set after that is kept.
    if not settings.get("column_order_v2_applied"):
        settings["column_order"] = []
        settings["column_order_v2_applied"] = True
    return settings


def load():
    if not CONFIG_PATH.exists():
        return dict(DEFAULTS)
    try:
        # utf-8-sig: a file saved with a byte-order mark (some Windows editors,
        # Windows PowerShell) would otherwise fail to parse, and every setting
        # would silently fall back to its default.
        data = json.loads(CONFIG_PATH.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return dict(DEFAULTS)
    merged = dict(DEFAULTS)
    merged.update({k: v for k, v in data.items() if k in DEFAULTS})
    return _migrate(merged)


def save(settings):
    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = {k: settings.get(k, DEFAULTS[k]) for k in DEFAULTS}
    CONFIG_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
