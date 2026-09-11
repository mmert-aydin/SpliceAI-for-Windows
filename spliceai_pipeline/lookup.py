"""Indexed lookup of precomputed SpliceAI scores.

The precomputed score files ship as four bgzip+tabix pairs per build
(spliceai_scores.{raw,masked}.{snv,indel}.hg{19,38}.vcf.gz[.tbi]). Each is
tens of gigabytes, so this never reads them sequentially -- it always goes
through the indexed TabixFile.fetch() from tabix.py.
"""
import os

from .tabix import TabixFile

FILENAME_TEMPLATE = "spliceai_scores.{mode}.{kind}.{build}.vcf.gz"


class PrecomputedScores:
    def __init__(self, precomputed_dir, build, mode):
        self.build = build
        self.mode = mode
        snv_path = os.path.join(precomputed_dir, FILENAME_TEMPLATE.format(mode=mode, kind="snv", build=build))
        indel_path = os.path.join(precomputed_dir, FILENAME_TEMPLATE.format(mode=mode, kind="indel", build=build))
        for p in (snv_path, indel_path):
            if not os.path.exists(p):
                raise FileNotFoundError(f"Precomputed score file not found: {p}")
        self._snv = TabixFile(snv_path)
        self._indel = TabixFile(indel_path)
        self._chrom_cache = {}

    def _tabix_for(self, ref, alt):
        return self._snv if len(ref) == 1 and len(alt) == 1 else self._indel

    def _chrom_for(self, tf, chrom):
        key = (id(tf), chrom)
        if key in self._chrom_cache:
            return self._chrom_cache[key]
        candidates = [chrom]
        bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
        candidates.append(bare)
        candidates.append("chr" + bare)
        resolved = next((c for c in candidates if tf.has_chrom(c)), None)
        self._chrom_cache[key] = resolved
        return resolved

    def lookup(self, chrom, pos, ref, alt):
        """Returns a list of 'ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL'
        annotation strings (one per overlapping gene), or None if this exact
        (pos, ref, alt) is not present in the precomputed file."""
        tf = self._tabix_for(ref, alt)
        tchrom = self._chrom_for(tf, chrom)
        if tchrom is None:
            return None

        for line in tf.fetch(tchrom, pos, pos):
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            lpos, lref, lalt, info = int(fields[1]), fields[3], fields[4], fields[7]
            if lpos != pos or lref != ref or lalt != alt:
                continue
            for kv in info.split(";"):
                if kv.startswith("SpliceAI="):
                    return kv[len("SpliceAI="):].split(",")
        return None

    def close(self):
        self._snv.close()
        self._indel.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
