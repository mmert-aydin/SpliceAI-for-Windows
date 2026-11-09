"""Texts the window shows that the pipeline writes: always English here (files
and bundles keep them English), and remembered with what they were made from,
so the window can show them in its own language (spliceai_gui/i18n.t).

    tx("Only {same} of {n} exons", same=3, n=5)  ->  "Only 3 of 5 exons"

and MADE["Only 3 of 5 exons"] == ("Only {same} of {n} exons", {"same": 3, "n": 5}).
"""
MADE = {}
LIMIT = 50_000      # ponytail: cleared when full; a text made before that shows in English, raise if it matters


def tx(template, **fields):
    text = template.format(**fields) if fields else template
    if len(MADE) >= LIMIT:
        MADE.clear()
    MADE[text] = (template, fields)
    return text


def tx_join(sep, items):
    """sep.join(items), remembered so that each item (a made text itself) is translated in its place."""
    if len(items) < 2:          # the item itself: remembering it again would replace what it was made from
        return sep.join(items)
    fields = {f"i{k}": v for k, v in enumerate(items)}
    return tx(sep.join("{i%d}" % k for k in range(len(items))), **fields)


# Words the pipeline puts into made texts as values (an event's kind, an anchor's name): the window translates
# them where they stand (spliceai_gui/tr.py has an entry for each).
WORDS = ("deletion", "insertion", "duplication", "inversion", "complex", "snv", "indel",
         "transcript start", "transcript end", "ATG", "stop codon", "acceptor", "donor", "exon", "intron")
