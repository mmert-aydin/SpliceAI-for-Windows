"""bedGraph export of per-base tracks (spec 5.5), for IGV or the UCSC browser.

write() writes one track as a bgzip-compressed bedGraph (Biopython's
BgzfWriter) with its tabix index (.tbi), written here to the tabix format
(https://samtools.github.io/hts-specs/tabix.pdf, UCSC-style 0-based
coordinates, as `tabix -p bed` would). Runs of equal values (at the precision
written) become one line and runs of zero are left out, so a track is small
even for a 2 Mb gene. tabix.TabixFile reads it back.
"""
import struct

import numpy as np
from Bio import bgzf

TBX_UCSC = 0x10000      # tabix format flag: 0-based, half-open coordinates (BED)
LINEAR_SHIFT = 14       # the linear index has one entry per 16 kb


def reg2bin(beg, end):
    """The smallest tabix/BAM bin holding the 0-based half-open [beg, end)."""
    end -= 1
    for shift, first in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return first + (beg >> shift)
    return 0


def runs(values, start, decimals=4):
    """[(begin, end, value)]: 0-based half-open runs of equal values (rounded
    to `decimals`) of a track whose first value is at 1-based position
    `start`; zero runs are left out."""
    v = np.round(np.asarray(values, np.float64), decimals)
    if not len(v):
        return []
    edges = np.nonzero(np.diff(v))[0] + 1
    begins = np.r_[0, edges]
    ends = np.r_[edges, len(v)]
    return [(start - 1 + int(b), start - 1 + int(e), float(v[b])) for b, e in zip(begins, ends) if v[b] != 0]


def write(path, chrom, start, values, decimals=4, description=""):
    """Writes the track (values at 1-based positions start, start + 1, ...) to
    `path` (.bedgraph.gz) and its index to path + ".tbi". Returns the number
    of lines written."""
    lines = runs(values, start, decimals)
    bins, linear = {}, {}
    with bgzf.BgzfWriter(path, "wb") as out:
        if description:
            out.write(f"# {description}\n".encode("ascii", "replace"))
        last_bin = None
        for beg, end, value in lines:
            offset = out.tell()
            out.write(f"{chrom}\t{beg}\t{end}\t{value:.{decimals}f}\n".encode("ascii"))
            after = out.tell()
            b = reg2bin(beg, end)
            if b == last_bin:
                bins[b][-1][1] = after          # the same bin as the line before: one chunk
            else:
                bins.setdefault(b, []).append([offset, after])
            last_bin = b
            for w in range(beg >> LINEAR_SHIFT, ((end - 1) >> LINEAR_SHIFT) + 1):
                linear.setdefault(w, offset)
    _write_index(path + ".tbi", chrom, bins, linear)
    return len(lines)


def _write_index(path, chrom, bins, linear):
    names = chrom.encode("ascii") + b"\x00"
    data = [b"TBI\x01", struct.pack("<8i", 1, TBX_UCSC, 1, 2, 3, ord("#"), 0, len(names)), names]
    data.append(struct.pack("<i", len(bins)))
    for b in sorted(bins):
        data.append(struct.pack("<Ii", b, len(bins[b])))
        for beg, end in bins[b]:
            data.append(struct.pack("<QQ", beg, end))
    # One entry per 16 kb window up to the last one used; an empty window takes the previous offset.
    n = max(linear) + 1 if linear else 0
    entries, previous = [], 0
    for w in range(n):
        previous = linear.get(w, previous)
        entries.append(previous)
    data.append(struct.pack("<i", n))
    data.append(struct.pack(f"<{n}Q", *entries))
    with bgzf.BgzfWriter(path, "wb") as out:
        out.write(b"".join(data))
