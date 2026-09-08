#!/usr/bin/env python3
"""
Empirical design-effect measurement on real modkit-extract data, per the
"the number that decides everything" section of
github/ont_asm_caller/docs/2026-09-05_calibration_critique.md: does
region-pooling (region.cluster_cpgs treating every CpG's reads as
independent) actually hold on this project's real data, or are the same
molecules being counted multiple times because ONT reads span every CpG in
a 500bp window?

Scope (pilot, not exhaustive): restricted to a single early window on chr1
(default chr1:1-5,000,000) for tractability -- modkit extract output is not
position-indexed, and this project's chr1 extract files are 4-5GB gzipped
per sample. Uses the SAME region boundaries the production beta-binomial
caller already computed (tables/dasm_betabinom/{sample}/regions_chr1.tsv),
subsampled to the first N_REGIONS within the window, so design_effect is
measured on the actual regions real ASM calls came from, not independently
redefined windows.

Modkit extract itself is NOT haplotype-split (unlike the pileup step) -- HP
tags are recovered separately from the haplotagged BAM via samtools view.

Usage: python3 P05_design_effect_real_data.py <SAMPLE> [n_regions] [window_end]
Output: tables/phasing_qc/design_effect_{SAMPLE}_chr1_pilot.tsv + summary printed
"""
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/u/project/cluo/terencew/claude/project_ideas/asm_lr/github/ont_asm_caller")
from ont_asm_caller.readlevel import design_effect

PROJ = "/u/project/cluo/terencew/claude/project_ideas/asm_lr"
CHROM = "chr1"
SAMTOOLS = "/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools"
MOD_QUAL_THRESH = 0.5


def get_hp_tags(sample, chrom, window_end):
    bam = f"{PROJ}/bam/haplotagged/by_chrom/{sample}_{chrom}_haplotagged.bam"
    cmd = [SAMTOOLS, "view", bam, f"{chrom}:1-{window_end}"]
    hp = {}
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        qname = fields[0]
        tag = None
        for f in fields[11:]:
            if f.startswith("HP:i:"):
                tag = int(f[5:])
                break
        if tag is not None:
            hp[qname] = tag
    proc.wait()
    return hp


def load_extract_window(sample, chrom, window_end):
    """Stream-filter modkit extract for mod_code=='m', ref_position in
    [0, window_end + buffer], with an early exit once we're clearly past
    the window (accepting minor edge loss from very long reads straddling
    the cutoff -- acceptable for a pilot estimate)."""
    path = f"{PROJ}/modkit/{sample}_{chrom}_modkit_extract.tsv.gz"
    buffer = 200_000
    hard_stop = window_end + buffer
    # column indices verified against the header: read_id(1)
    # forward_read_position(2) ref_position(3) chrom(4) mod_strand(5)
    # ref_strand(6) ref_mod_strand(7) fw_soft_clipped_start(8)
    # fw_soft_clipped_end(9) read_length(10) mod_qual(11) mod_code(12) ...
    awk_prog = (
        f'BEGIN{{FS="\\t"}} NR>1 && $12=="m" && $3>=0 {{'
        f'  if ($3 <= {hard_stop}) print $1"\\t"$3"\\t"$11; '
        f'  else {{ n_over++; if (n_over > 200000) exit }} '
        f"}}"
    )
    cmd = f"zcat {path} | awk '{awk_prog}'"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    rows = []
    for line in proc.stdout:
        rid, pos, qual = line.rstrip("\n").split("\t")
        rows.append((rid, int(pos), float(qual)))
    proc.wait()
    df = pd.DataFrame(rows, columns=["read_id", "pos", "mod_qual"])
    df = df[df.pos <= window_end]
    return df


def main():
    sample = sys.argv[1]
    n_regions = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    window_end = int(sys.argv[3]) if len(sys.argv) > 3 else 5_000_000

    print(f"{sample}: loading HP tags for {CHROM}:1-{window_end}...", flush=True)
    hp = get_hp_tags(sample, CHROM, window_end)
    print(f"  {len(hp):,} reads with an HP tag", flush=True)

    print(f"{sample}: streaming modkit extract (mod_code=m) for the same window...", flush=True)
    calls = load_extract_window(sample, CHROM, window_end)
    print(f"  {len(calls):,} (read, CpG) call rows loaded", flush=True)
    calls["called"] = calls["mod_qual"] >= MOD_QUAL_THRESH
    calls["hp"] = calls["read_id"].map(hp)
    calls = calls.dropna(subset=["hp"])
    calls["hp"] = calls["hp"].astype(int)
    print(f"  {len(calls):,} rows after joining HP tags", flush=True)

    regions = pd.read_csv(f"{PROJ}/tables/dasm_betabinom/{sample}/regions_chr1.tsv", sep="\t")
    regions = regions[regions.end < window_end].sort_values("start").head(n_regions)
    print(f"{sample}: measuring design effect on {len(regions)} real regions...", flush=True)

    results = []
    calls_sorted = calls.sort_values("pos")
    pos_arr = calls_sorted["pos"].to_numpy()

    for r in regions.itertuples():
        lo_idx = np.searchsorted(pos_arr, r.start, side="left")
        hi_idx = np.searchsorted(pos_arr, r.end, side="right")
        sub = calls_sorted.iloc[lo_idx:hi_idx]
        if sub.empty:
            continue
        m_list = {1: None, 2: None}
        for hp_val in (1, 2):
            hp_sub = sub[sub.hp == hp_val]
            if hp_sub.empty:
                continue
            piv = hp_sub.pivot_table(index="read_id", columns="pos", values="called", aggfunc="first")
            m_list[hp_val] = piv.to_numpy(dtype=float)
        if m_list[1] is None or m_list[2] is None:
            continue
        if m_list[1].shape[0] < 3 or m_list[2].shape[0] < 3:
            continue
        de = design_effect(m_list[1], m_list[2])
        results.append({
            "chrom": r.chrom, "start": r.start, "end": r.end, "n_cpgs": r.n_cpgs,
            "n_reads1": m_list[1].shape[0], "n_reads2": m_list[2].shape[0],
            "design_effect": de,
        })

    out = pd.DataFrame(results)
    outpath = f"{PROJ}/tables/phasing_qc/design_effect_{sample}_chr1_pilot.tsv"
    out.to_csv(outpath, sep="\t", index=False)

    print(f"\n{sample}: {len(out)} regions with usable design-effect estimates")
    if len(out):
        print(out["design_effect"].describe())
        print(f"  fraction with DE > 2: {(out.design_effect > 2).mean():.3f}")
        print(f"  fraction with DE > 5: {(out.design_effect > 5).mean():.3f}")
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
