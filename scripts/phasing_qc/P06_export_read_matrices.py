#!/usr/bin/env python3
"""
Export a COMPACT per-region read x CpG matrix from haplotagged ONT data, once,
so that every downstream read-level analysis stops re-streaming multi-GB files.

WHY THIS EXISTS
---------------
`P05_design_effect_real_data.py` streams a 4-5 GB gzipped `modkit extract` file
per sample per window, and its own docstring notes this will not scale
genome-wide. Every read-level analysis in `github/ont_asm_caller` wants the same
thing from that file -- a (read x CpG) 0/1 matrix per region per haplotype --
and there is no reason to pay the streaming cost more than once.

At 10^5 regions x ~40 reads x ~8 CpGs this is ~3 x 10^7 calls: tens of MB as a
packed npz, against gigabytes as TSV. It is small enough to move off the
cluster, which is the actual point -- the methods work needs to iterate on it.

WHAT IT PRODUCES
----------------
`tables/read_matrices/{sample}_{chrom}.npz` holding, for each region that has
at least `MIN_READS` reads on both haplotypes:

    region_id     (R,) int32     index into the region table
    starts/ends   (R,) int32
    positions     concatenated int32 CpG positions, with `pos_offsets`
    m1/m2         concatenated int8 call matrices (0/1, -1 = no call),
                  with `m1_shape`/`m2_shape` per region

Load with `load_read_matrices()` below, which yields objects carrying
`.positions`, `.m1`, `.m2` -- the same duck type `ont_asm_caller.decay`,
`ont_asm_caller.pattern` and `ont_asm_caller.null` already consume, so nothing
downstream needs a real-data code path of its own.

Usage: python3 P06_export_read_matrices.py <SAMPLE> <CHROM> [max_regions]
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

PROJ = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"
SAMTOOLS = "/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools"
MIN_READS = 4
NO_CALL = -1

# --- call confidence -------------------------------------------------------
# modkit extract emits ONE ROW PER (read, position, mod_code), so a CpG appears
# twice: mod_code "m" (5mC) and "h" (5hmC), with mod_qual the probability of
# THAT modification. P(canonical C) = 1 - p_m - p_h. Example from real output:
#
#   ...400581... 0.6855 h        <- 5hmC probable
#   ...400581... 0.3145 m        <- and only 31% 5mC
#
# The earlier `mod_qual >= 0.5` rule on the "m" row alone turns that into a
# confident UNMETHYLATED call, even though the base is almost certainly
# modified. That matters more than it looks: a wrong binary call attenuates
# every correlation measure by (1-2*eps)^2, and the design effect IS a
# correlation measure. Measured on simulation at true DE 6.3:
#   eps=0.02 -> 5.9   eps=0.05 -> 5.3   eps=0.10 -> 4.4   eps=0.20 -> 2.9
# i.e. call noise can move DE across the whole range that decides whether
# region pooling is legal. (decay_bp itself is unaffected -- only amplitude.)
#
# So: take the argmax over {C, m, h}, require it to clear CONF_THRESH, and emit
# a NO-CALL otherwise. A no-call costs data; a wrong call costs calibration.
CONF_THRESH = 0.80
# What to do when 5hmC is the most probable state. "nocall" is the default
# because it is the only option that does not bias the correlation: calling it
# unmethylated is a systematic error concentrated in specific genomic contexts,
# and calling it methylated conflates two marks. Set "methylated" to reproduce
# modkit's --combine-mods behaviour, or "unmethylated" for the old rule.
H_POLICY = "nocall"          # "nocall" | "methylated" | "unmethylated"


def hp_tags(sample, chrom):
    """modkit extract is not haplotype-split; HP comes from the haplotagged BAM."""
    bam = f"{PROJ}/bam/haplotagged/by_chrom/{sample}_{chrom}_haplotagged.bam"
    proc = subprocess.Popen([SAMTOOLS, "view", bam], stdout=subprocess.PIPE, text=True)
    hp = {}
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        for tag in f[11:]:
            if tag.startswith("HP:i:"):
                hp[f[0]] = int(tag[5:])
                break
    proc.wait()
    return hp


def stream_calls(sample, chrom, stats):
    """(read_id, cpg_position, call) for CpG methylation, streamed once.

    `call` is 1 / 0 / None (no-call). `cpg_position` is normalised to the PLUS
    strand -- see below. Column indices verified against the real header:
    read_id(1) ref_position(3) ref_mod_strand(7) mod_qual(11) mod_code(12).

    STRAND NORMALISATION. modkit extract is per-strand and does not combine.
    A CpG occupies (p, p+1) = (C, G) on the plus strand, so its minus-strand C
    sits at p+1 and modkit reports the SAME physical CpG at position p from
    plus-strand reads and p+1 from minus-strand reads. Keying columns on the raw
    ref_position therefore splits every CpG into two columns 1 bp apart that
    share ZERO reads. design_effect survives that (it reduces each read to a
    fraction over whatever columns that read saw), but the decay curve loses
    exactly its short-distance bins, fit_chain/NME/PDM get a chain over twice
    the sites with strand-structured missingness, and n_cpgs -- a null
    stratification variable -- doubles.
    """
    path = f"{PROJ}/modkit/{sample}_{chrom}_modkit_extract.tsv.gz"
    # keep BOTH mod codes: P(canonical) = 1 - p_m - p_h needs both
    awk = ('BEGIN{FS="\\t"} NR>1 && ($12=="m" || $12=="h") && $3>=0 '
           '{print $1"\\t"$3"\\t"$7"\\t"$11"\\t"$12}')
    proc = subprocess.Popen(f"zcat {path} | awk '{awk}'", shell=True,
                            stdout=subprocess.PIPE, text=True)
    pending = {}
    for line in proc.stdout:
        rid, pos, strand, qual, code = line.rstrip("\n").split("\t")
        pos = int(pos) - 1 if strand == "-" else int(pos)
        stats["strand_minus" if strand == "-" else "strand_plus"] += 1
        key = (rid, pos)
        d = pending.setdefault(key, {})
        d[code] = float(qual)
        if len(d) < 2:
            continue                      # wait for the partner row
        del pending[key]
        p_m, p_h = d.get("m", 0.0), d.get("h", 0.0)
        p_c = 1.0 - p_m - p_h
        best = max((p_c, "c"), (p_m, "m"), (p_h, "h"))
        stats["total"] += 1
        if best[0] < CONF_THRESH:
            stats["ambiguous"] += 1
            yield rid, pos, None
        elif best[1] == "m":
            yield rid, pos, True
        elif best[1] == "c":
            yield rid, pos, False
        else:
            stats["hmc"] += 1
            yield rid, pos, (None if H_POLICY == "nocall"
                             else H_POLICY == "methylated")
    # rows whose partner never arrived (shouldn't happen, but don't lose them)
    for (rid, pos), d in pending.items():
        stats["unpaired"] += 1
        yield rid, pos, None
    proc.wait()


def main():
    sample, chrom = sys.argv[1], sys.argv[2]
    max_regions = int(sys.argv[3]) if len(sys.argv) > 3 else None

    regions = pd.read_csv(f"{PROJ}/tables/dasm_betabinom/{sample}/regions_{chrom}.tsv",
                          sep="\t").sort_values("start").reset_index(drop=True)
    if max_regions:
        regions = regions.head(max_regions)
    print(f"{sample} {chrom}: {len(regions)} regions", flush=True)

    hp = hp_tags(sample, chrom)
    print(f"  {len(hp):,} HP-tagged reads", flush=True)

    starts = regions["start"].to_numpy()
    ends = regions["end"].to_numpy()
    # one pass over the calls, bucketed into regions by binary search on start
    buckets = {i: {} for i in range(len(regions))}
    stats = {k: 0 for k in ("total", "ambiguous", "hmc", "unpaired",
                            "strand_plus", "strand_minus")}
    n_rows = n_kept = 0
    for rid, pos, called in stream_calls(sample, chrom, stats):
        n_rows += 1
        h = hp.get(rid)
        if h is None:
            continue
        i = int(np.searchsorted(starts, pos, side="right")) - 1
        if i < 0 or pos > ends[i]:
            continue
        buckets[i].setdefault((h, rid), {})[pos] = called
        n_kept += 1
        if n_rows % 20_000_000 == 0:
            print(f"    {n_rows:,} rows scanned, {n_kept:,} kept", flush=True)
    print(f"  {n_rows:,} call rows scanned, {n_kept:,} assigned to a region", flush=True)
    t = max(stats["total"], 1)
    print(f"  QC: {stats['ambiguous']:,} ({stats['ambiguous']/t:.1%}) below "
          f"CONF_THRESH={CONF_THRESH} -> no-call")
    print(f"      {stats['hmc']:,} ({stats['hmc']/t:.1%}) 5hmC-dominant "
          f"-> H_POLICY={H_POLICY}")
    print(f"      strand: {stats['strand_plus']:,} plus / "
          f"{stats['strand_minus']:,} minus"
          f"{'  <-- both present, normalisation was required' if stats['strand_plus'] and stats['strand_minus'] else ''}")
    if stats["unpaired"]:
        print(f"      WARNING {stats['unpaired']:,} rows had no m/h partner")

    out = {k: [] for k in ("region_id", "start", "end", "positions",
                           "pos_offset", "m1", "m2", "m1_shape", "m2_shape")}
    pos_off = 0
    for i in range(len(regions)):
        cells = buckets[i]
        if not cells:
            continue
        allpos = sorted({p for d in cells.values() for p in d})
        if len(allpos) < 2:
            continue
        pidx = {p: j for j, p in enumerate(allpos)}
        mats = {}
        for h in (1, 2):
            reads = [d for (hh, _), d in cells.items() if hh == h]
            if len(reads) < MIN_READS:
                mats = None
                break
            m = np.full((len(reads), len(allpos)), NO_CALL, np.int8)
            for r, d in enumerate(reads):
                for p, c in d.items():
                    if c is not None:
                        m[r, pidx[p]] = 1 if c else 0
            mats[h] = m
        if mats is None:
            continue
        out["region_id"].append(i)
        out["start"].append(int(regions["start"][i]))
        out["end"].append(int(regions["end"][i]))
        out["positions"].append(np.asarray(allpos, np.int32))
        out["pos_offset"].append(pos_off); pos_off += len(allpos)
        out["m1"].append(mats[1]); out["m2"].append(mats[2])
        out["m1_shape"].append(mats[1].shape); out["m2_shape"].append(mats[2].shape)

    os.makedirs(f"{PROJ}/tables/read_matrices", exist_ok=True)
    path = f"{PROJ}/tables/read_matrices/{sample}_{chrom}.npz"
    np.savez_compressed(
        path,
        region_id=np.asarray(out["region_id"], np.int32),
        start=np.asarray(out["start"], np.int32),
        end=np.asarray(out["end"], np.int32),
        positions=np.concatenate(out["positions"]) if out["positions"] else np.zeros(0, np.int32),
        pos_offset=np.asarray(out["pos_offset"], np.int64),
        m1=np.concatenate([m.ravel() for m in out["m1"]]) if out["m1"] else np.zeros(0, np.int8),
        m2=np.concatenate([m.ravel() for m in out["m2"]]) if out["m2"] else np.zeros(0, np.int8),
        m1_shape=np.asarray(out["m1_shape"], np.int32).reshape(-1, 2),
        m2_shape=np.asarray(out["m2_shape"], np.int32).reshape(-1, 2),
    )
    mb = os.path.getsize(path) / 1e6
    print(f"wrote {path}  ({len(out['region_id'])} regions, {mb:.1f} MB)")


class LoadedRegion:
    """Duck-typed to match simulate_patterns.PatternRegion so every analysis in
    ont_asm_caller consumes real and simulated data through the same path."""
    __slots__ = ("chrom", "start", "end", "positions", "m1", "m2")

    def __init__(self, chrom, start, end, positions, m1, m2):
        self.chrom, self.start, self.end = chrom, int(start), int(end)
        self.positions, self.m1, self.m2 = positions, m1, m2

    @property
    def n_cpgs(self):
        return len(self.positions)


def load_read_matrices(path, chrom="chr1"):
    """Yield LoadedRegion objects. No-calls come back as np.nan, which is the
    convention every consumer in ont_asm_caller already uses."""
    z = np.load(path)
    pos, off = z["positions"], z["pos_offset"]
    m1flat, m2flat = z["m1"], z["m2"]
    s1, s2 = z["m1_shape"], z["m2_shape"]
    o1 = np.concatenate([[0], np.cumsum(s1[:, 0] * s1[:, 1])])
    o2 = np.concatenate([[0], np.cumsum(s2[:, 0] * s2[:, 1])])
    for i in range(len(z["start"])):
        n = s1[i, 1]
        p = pos[off[i]:off[i] + n].astype(float)
        a = m1flat[o1[i]:o1[i + 1]].reshape(s1[i]).astype(float)
        b = m2flat[o2[i]:o2[i + 1]].reshape(s2[i]).astype(float)
        a[a == NO_CALL] = np.nan
        b[b == NO_CALL] = np.nan
        yield LoadedRegion(chrom, z["start"][i], z["end"][i], p, a, b)


if __name__ == "__main__":
    main()
