"""Persistence for named, reusable gene lists (save/load/delete).

Stored per-user under ~/.spliceai_gui/gene_lists.json -- same location
convention as config.py, never inside the project folder, so it's never
accidentally shared/committed/copied along with the app.
"""

import json
from pathlib import Path

STORAGE_PATH = Path.home() / ".spliceai_gui" / "gene_lists.json"


def load_all():
    """Returns {list_name: [gene_symbol, ...]}."""
    if not STORAGE_PATH.exists():
        return {}
    try:
        # utf-8-sig: also reads a file saved with a byte-order mark (see config.py).
        data = json.loads(STORAGE_PATH.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    return {str(name): list(genes) for name, genes in data.items() if isinstance(genes, list)}


def save_list(name, genes):
    all_lists = load_all()
    all_lists[name] = list(genes)
    _write(all_lists)


def delete_list(name):
    all_lists = load_all()
    all_lists.pop(name, None)
    _write(all_lists)


def _write(all_lists):
    STORAGE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STORAGE_PATH.write_text(json.dumps(all_lists, indent=2), encoding="utf-8")
