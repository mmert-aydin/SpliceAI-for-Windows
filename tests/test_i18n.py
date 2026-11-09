"""The Turkish texts (spliceai_gui/i18n.py, tr.py): every t("...") in the
window's code has a Turkish entry with the same {placeholders}, none is an
f-string (it could never match), the owner's spelling rules hold, and the
language setting switches the texts.

    .venv\\Scripts\\python tests\\test_i18n.py
"""
import ast
import glob
import os
import re
import string

import helpers
from helpers import check

helpers.isolate_home()
from spliceai_gui import config, i18n  # noqa: E402
from spliceai_gui.tr import TR  # noqa: E402

keys, bad = {}, []
gui = glob.glob(os.path.join(str(helpers.REPO), "spliceai_gui", "*.py"))
pipeline = [p for p in glob.glob(os.path.join(str(helpers.REPO), "spliceai_pipeline", "*.py"))
            if not p.endswith("texts.py")]
for path in gui + pipeline:
    tree = ast.parse(open(path, encoding="utf-8").read())
    names = ("tx",) if path in pipeline else ("t", "N_")
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) in names and node.args:
            arg = node.args[0]
            where = f"{os.path.basename(path)}:{node.lineno}"
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                keys.setdefault(arg.value, where)
            elif isinstance(arg, (ast.JoinedStr, ast.BinOp)) or node.func.id in ("N_", "tx"):
                bad.append(where)       # t(name) is fine: the name holds an N_("...") text



def shadowed(path):
    """Scopes (functions, lambdas, comprehensions) that call t(...) while `t` is a local of theirs or of a
    function around them: the call would get that variable, not i18n.t."""
    scopes = (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ListComp, ast.SetComp, ast.DictComp,
              ast.GeneratorExp)
    out = []

    def visit(scope, outer_binds):
        args = getattr(scope, "args", None)
        binds = args is not None and any(a.arg == "t" for a in args.args + args.kwonlyargs + args.posonlyargs)
        calls, nested, stack = [], [], list(ast.iter_child_nodes(scope))
        while stack:                    # this scope only; nested ones are visited with what they see
            node = stack.pop()
            if isinstance(node, scopes):
                nested.append(node)
                continue
            binds |= isinstance(node, ast.Name) and node.id == "t" and isinstance(node.ctx, (ast.Store, ast.Del))
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "t":
                calls.append(node.lineno)
            stack.extend(ast.iter_child_nodes(node))
        if calls and (binds or outer_binds):
            out.append(f"{os.path.basename(path)}:{min(calls)}")
        for n in nested:
            visit(n, binds or outer_binds)

    tree = ast.parse(open(path, encoding="utf-8").read())
    top = [n for n in ast.walk(tree) if isinstance(n, scopes)]
    inner = {id(c) for n in top for c in ast.walk(n) if c is not n and isinstance(c, scopes)}
    for node in top:
        if id(node) not in inner:
            visit(node, False)
    return out


clashes = [c for p in glob.glob(os.path.join(str(helpers.REPO), "spliceai_gui", "*.py")) for c in shadowed(p)]
check("no function both names something `t` and calls t(...)", not clashes, str(clashes))
from spliceai_pipeline.texts import WORDS  # noqa: E402

for w in WORDS:
    keys.setdefault(w, "texts.py")
fields = lambda s: sorted(f for _, f, _, _ in string.Formatter().parse(s) if f)
missing = [f"{w}: {k[:60]}" for k, w in keys.items() if k not in TR]
check("no t(...) of an f-string or a sum of strings (it could never match)", not bad, str(bad[:10]))
check("every text shown has its Turkish entry", not missing, f"{len(missing)} missing: {missing[:10]}")
check("each Turkish entry has the same {placeholders} as its English text",
      all(fields(k) == fields(v) for k, v in TR.items()),
      [k[:60] for k, v in TR.items() if fields(k) != fields(v)][:10])
check("no unused Turkish entries", not set(TR) - set(keys), [k[:60] for k in set(TR) - set(keys)][:10])
words = lambda pattern: [v[:80] for v in TR.values() if re.search(pattern, v)]
check("'ekzon', never 'ekson'", not words(r"(?i)ekson"), words(r"(?i)ekson")[:5])
# "ham veri" (raw data in general, the owner's own wording in the opening notice) isn't the raw mode.
check("masked and raw stay as they are (not 'maskeli', 'ham')", not words(r"(?i)maskeli|\bham\b(?! veri)"),
      words(r"(?i)maskeli|\bham\b(?! veri)")[:5])

TR.setdefault("Run", "Çalıştır")
check("English by default", i18n.language() == "en" and i18n.t("Run") == "Run")
settings = config.load()
settings["language"] = "tr"
config.save(settings)
i18n._language = None
check("Türkçe once chosen", i18n.language() == "tr" and i18n.t("Run") == "Çalıştır")
check("a text without an entry falls back to English", i18n.t("no such text") == "no such text")

# Texts the pipeline makes stay English (files), and show in Türkçe in the window.
from spliceai_pipeline import align, cdna, vcf_bridge  # noqa: E402
from spliceai_pipeline.texts import tx, tx_join  # noqa: E402

try:
    align.parse_region("chrQ:1-2")
    message = ""
except ValueError as exc:
    message = str(exc)
check("a pipeline error is English as made, Türkçe as shown", message.startswith("'chrQ:1-2' isn't a region")
      and i18n.t(message).startswith("'chrQ:1-2', chr7:117,480,025-117,668,665 gibi bir bölge değil"),
      i18n.t(message))
rows, _ = cdna.concordance([(1, 10), (20, 30)], [(1, 12)], "+")
inner = tx_join("; ", [tx("exon {k} ({status})", k=r["ref"], status=r["status"]) for r in rows])
check("made texts inside made texts are translated too (an exon's status in a list)",
      inner == "exon 1 (end shifted +2); exon 2 (missing)"
      and i18n.t(inner) == "ekzon 1 (bitiş +2 kaymış); ekzon 2 (yok)", i18n.t(inner))
g = vcf_bridge.Genotype("0/1", "S1", 1, tx("heterozygous"), tx("The other allele is expected to give the WT transcript."))
check("a genotype line in Türkçe", i18n.t(g.text()) == "0/1 (heterozigot; S1). Diğer alelin WT transkripti vermesi "
      "beklenir.", i18n.t(g.text()))
shown = i18n.t("Genotype:") + " " + i18n.t("{n} gene(s): {genes}", n=1, genes=g.text())
check("...also as a value inside a window text", "heterozigot" in shown and "heterozygous" not in shown, shown)
i18n._language = None
helpers.finish()
