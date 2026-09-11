"""Pure-Python tabix (.tbi) index reader with BGZF random access.

cyvcf2/htslib publish no Windows wheels and cannot be built on Windows without a
pre-existing external htslib (verified: cyvcf2's own CMakeLists.txt hard-errors
on CYVCF2_HTSLIB_MODE=BUILTIN when WIN32), so indexed lookups against the
precomputed SpliceAI score files (28-69 GB each) are implemented directly
against the documented tabix index format
(https://samtools.github.io/hts-specs/tabix.pdf) using Bio.bgzf for BGZF
block-level random access -- both cyvcf2 and htslib itself use this exact
virtual-offset (coffset<<16 | uoffset) scheme, which Bio.bgzf.BgzfReader
implements natively via seek()/tell().
"""
import gzip
import struct
from collections import OrderedDict
from dataclasses import dataclass, field

from Bio import bgzf

MAGIC = b"TBI\x01"


def reg2bins(beg, end):
    """Tabix bin ids that can overlap a 0-based half-open region [beg, end)."""
    end -= 1
    bins = [0]
    for start, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        bins.extend(range(start + (beg >> shift), start + (end >> shift) + 1))
    return bins


@dataclass
class _RefIndex:
    bins: dict = field(default_factory=dict)
    intervals: list = field(default_factory=list)


class TabixIndex:
    """Parsed representation of a .tbi index file (itself BGZF-compressed)."""

    def __init__(self, path):
        with gzip.open(path, "rb") as fh:
            data = fh.read()
        if data[:4] != MAGIC:
            raise ValueError(f"{path} does not look like a tabix .tbi index (bad magic)")

        n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", data[4:36])
        offset = 36
        names = data[offset:offset + l_nm].decode("ascii").rstrip("\x00").split("\x00")
        offset += l_nm

        self.format = fmt
        self.col_seq, self.col_beg, self.col_end = col_seq, col_beg, col_end
        self.name_to_id = {name: i for i, name in enumerate(names)}
        self.refs = []

        for _ in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", data, offset)
            offset += 4
            ref = _RefIndex()
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack_from("<Ii", data, offset)
                offset += 8
                chunks = []
                for _ in range(n_chunk):
                    chunk_beg, chunk_end = struct.unpack_from("<QQ", data, offset)
                    offset += 16
                    chunks.append((chunk_beg, chunk_end))
                ref.bins[bin_id] = chunks
            (n_intv,) = struct.unpack_from("<i", data, offset)
            offset += 4
            if n_intv:
                ref.intervals = list(struct.unpack_from(f"<{n_intv}Q", data, offset))
                offset += 8 * n_intv
            self.refs.append(ref)

    def chunks_for_region(self, chrom, beg, end):
        """beg/end are 0-based half-open. Returns sorted (chunk_beg, chunk_end) virtual offsets."""
        ref_id = self.name_to_id.get(chrom)
        if ref_id is None:
            return []
        ref = self.refs[ref_id]

        min_offset = 0
        if ref.intervals:
            iv_idx = beg >> 14
            if iv_idx < len(ref.intervals):
                min_offset = ref.intervals[iv_idx]
            else:
                min_offset = ref.intervals[-1]

        chunks = [
            (chunk_beg, chunk_end)
            for bin_id in reg2bins(beg, end)
            for chunk_beg, chunk_end in ref.bins.get(bin_id, ())
            if chunk_end > min_offset
        ]
        chunks.sort()
        return chunks

    def has_chrom(self, chrom):
        return chrom in self.name_to_id


class TabixFile:
    """Indexed random access to a bgzip-compressed, tabix-indexed tab-delimited file.

    Point queries (single-variant lookups) against a densely-populated file
    still require decompressing and scanning an entire tabix bin/chunk (the
    same is true of real htslib tabix -- there is no finer index
    granularity). Nearby variants in a coordinate-sorted VCF very often land
    in the same chunk, so completed chunks are cached (LRU, bounded) to avoid
    repeatedly re-scanning the same bin for every variant in a dense region.
    """

    def __init__(self, data_path, index_path=None, chunk_cache_size=256):
        self.index = TabixIndex(index_path or data_path + ".tbi")
        self._reader = bgzf.BgzfReader(data_path, "r")
        self._chunk_cache_size = chunk_cache_size
        self._chunk_cache = OrderedDict()

    def has_chrom(self, chrom):
        return self.index.has_chrom(chrom)

    def _read_chunk(self, chunk_beg, chunk_end):
        key = (chunk_beg, chunk_end)
        cached = self._chunk_cache.get(key)
        if cached is not None:
            self._chunk_cache.move_to_end(key)
            return cached

        lines = []
        self._reader.seek(chunk_beg)
        while self._reader.tell() < chunk_end:
            line = self._reader.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))

        self._chunk_cache[key] = lines
        if len(self._chunk_cache) > self._chunk_cache_size:
            self._chunk_cache.popitem(last=False)
        return lines

    def fetch(self, chrom, start, end):
        """1-based inclusive [start, end]. Yields candidate raw lines overlapping the
        query at bin granularity -- callers must still confirm the exact
        position/allele match themselves."""
        beg0, end0 = start - 1, end
        for chunk_beg, chunk_end in self.index.chunks_for_region(chrom, beg0, end0):
            yield from self._read_chunk(chunk_beg, chunk_end)

    def close(self):
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
