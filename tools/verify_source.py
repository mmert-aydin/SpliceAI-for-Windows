"""Prove a reconstructed .py is the same program as the shipped .pyc.

Compiles each candidate source with this Python 3.13 and compares the code
objects recursively against the originals extracted from the exe:

  EXACT   -- bytecode, constants (incl. docstrings), names, variables, flags,
             exception tables all identical. Only line numbers may differ.
  NOPS    -- identical except for NOP instructions (Python inserts those to
             keep per-line tracing; they depend on how statements are split
             across lines, not on behaviour). Jumps are compared by target
             instruction, so the control flow is still proven identical.
  DIFF    -- anything else; the first differing location is printed.

usage: python verify_source.py <orig_pyz_dir> <src_dir> [module ...]
       (module like spliceai_gui/main; default: every .pyc under orig dir)
"""
import dis
import marshal
import os
import sys
import types

ATTRS = ("co_argcount", "co_posonlyargcount", "co_kwonlyargcount", "co_flags",
         "co_names", "co_varnames", "co_freevars", "co_cellvars", "co_name", "co_qualname")


def load_pyc(path):
    with open(path, "rb") as fh:
        return marshal.loads(fh.read()[16:])


def norm_instructions(co):
    """(opname, arg) list without NOPs; jump args become target indices."""
    ins = [i for i in dis.get_instructions(co) if i.opname not in ("NOP", "CACHE")]
    idx = {i.offset: n for n, i in enumerate(ins)}
    all_ins = list(dis.get_instructions(co))

    def target_index(off):
        # a jump may land on a removed NOP; move forward to the next kept one
        for i in all_ins:
            if i.offset >= off and i.offset in idx:
                return idx[i.offset]
        return len(ins)

    out = []
    for i in ins:
        if i.opcode in dis.hasjrel or i.opcode in dis.hasjabs:
            out.append((i.opname, "->", target_index(i.argval)))
        elif isinstance(i.argval, types.CodeType):
            out.append((i.opname, "<code>"))
        else:
            out.append((i.opname, repr(i.argval)))
    return out


def compare(a, b, path, strict):
    """Return None if equal, else a description of the first difference."""
    for attr in ATTRS:
        if getattr(a, attr) != getattr(b, attr):
            return f"{path}: {attr} {getattr(a, attr)!r} != {getattr(b, attr)!r}"
    if len(a.co_consts) != len(b.co_consts):
        return f"{path}: {len(a.co_consts)} consts != {len(b.co_consts)}"
    for n, (x, y) in enumerate(zip(a.co_consts, b.co_consts)):
        if isinstance(x, types.CodeType) and isinstance(y, types.CodeType):
            r = compare(x, y, f"{path}.{x.co_name}", strict)
            if r:
                return r
        elif type(x) is not type(y) or x != y:
            return f"{path}: const[{n}] {x!r:.120} != {y!r:.120}"
    if strict:
        if a.co_code != b.co_code:
            return f"{path}: co_code differs"
        if a.co_exceptiontable != b.co_exceptiontable:
            return f"{path}: exception table differs"
    else:
        na, nb = norm_instructions(a), norm_instructions(b)
        if na != nb:
            for k, (p, q) in enumerate(zip(na, nb)):
                if p != q:
                    return f"{path}: instr #{k} {p} != {q}"
            return f"{path}: {len(na)} instrs != {len(nb)}"
    return None


def verify(orig_pyc, src_py):
    orig = load_pyc(orig_pyc)
    with open(src_py, encoding="utf-8") as fh:
        source = fh.read()
    try:
        cand = compile(source, orig.co_filename, "exec", dont_inherit=True)
    except SyntaxError as exc:
        return "DIFF", f"SyntaxError line {exc.lineno}: {exc.msg}"
    if compare(orig, cand, "<module>", strict=True) is None:
        return "EXACT", ""
    why = compare(orig, cand, "<module>", strict=False)
    return ("NOPS", "") if why is None else ("DIFF", why)


def main():
    orig_dir, src_dir = sys.argv[1], sys.argv[2]
    mods = sys.argv[3:]
    if not mods:
        for root, _, files in os.walk(orig_dir):
            for f in files:
                if f.endswith(".pyc"):
                    rel = os.path.relpath(os.path.join(root, f), orig_dir)
                    mods.append(rel[:-4].replace(os.sep, "/"))
    counts = {}
    for m in sorted(mods):
        src = os.path.join(src_dir, *m.split("/")) + ".py"
        if not os.path.exists(src):
            status, why = "MISSING", ""
        else:
            status, why = verify(os.path.join(orig_dir, *m.split("/")) + ".pyc", src)
        counts[status] = counts.get(status, 0) + 1
        print(f"{status:7} {m}" + (f"   <- {why}" if why else ""), flush=True)
    print("summary:", counts)


if __name__ == "__main__":
    main()
