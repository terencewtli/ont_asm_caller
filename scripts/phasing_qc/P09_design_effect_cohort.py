#!/usr/bin/env python3
"""
Design effect across donors on a FIXED set of genomic windows, plus the
per-window quantities needed to explain why the two DE estimators disagree.

WHY THIS EXISTS SEPARATELY FROM P08
----------------------------------
P08 reads the first `max_rows` rows of an extract. That was fine for one donor,
but it cannot compare donors: the genomic span covered by N rows depends on
depth, so a 56x donor and an 87x donor would be measured on different loci.
P09 streams the WHOLE chromosome and keeps every K-th window of `window_bp`
(an awk prefilter does the discarding, so Python only sees ~1/K of the rows).
Every donor is then measured on the same windows.

Call logic (strand normalisation, argmax over {C, m, h}, CONF_THRESH, 5hmC
no-call, 1 bp strand-pair merge) is imported from P08 unchanged, so the two
scripts cannot drift apart.

THE ESTIMATOR QUESTION (docs/PRIORITIES.md item 3)
--------------------------------------------------
`decay.implied_design_effect` and `readlevel.design_effect` disagreed ~2x at
500 bp on HG00146 (2.77 vs 1.29). Four candidate explanations differ in what
they predict, so each window records enough to separate them:

  de_split  readlevel.design_effect on two arbitrary read halves (what P08 and
            PRIORITIES.md report)
  de_one    readlevel.design_effect on all reads as one group
  icc_w     variance-weighted pairwise within-read correlation, each CpG
            centred by its own mean -- the same definition decay._accumulate
            pools across windows
  cbar      mean observed calls per read
  de_icc    1 + (cbar - 1) * icc_w          (the implied-DE formula, per window)
  de_het    Var(read fraction) / mean_r( sum_{j in r} mu_j(1-mu_j) / C_r^2 )
            i.e. the design effect against a HETEROGENEOUS-mean baseline

Usage: P09_design_effect_cohort.py <extract.tsv.gz> <out.tsv>
                                   [window_bp=500] [stride=10] [max_rows=0]
"""
import os
import shlex
import subprocess
import sys
from collections import defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))
from P08_decay_from_extract import CONF_THRESH, H_POLICY, to_matrices  # noqa: E402
from ont_asm_caller.readlevel import design_effect  # noqa: E402

# columns: $1 read_id  $3 ref_position  $7 ref_mod_strand  $11 mod_qual  $12 mod_code
AWK = ('NR>1 && $3>=0 && ($12=="m"||$12=="h") '
       '{p=$3; if($7=="-")p=p-1; if(int(p/W)%K==0) print $1"\\t"p"\\t"$11"\\t"$12}')


def stream(path, window_bp, stride, max_rows, stats):
    head = f" | head -n {max_rows + 1}" if max_rows else ""
    cmd = (f"gzip -dc {shlex.quote(path)}{head} | "
           f"awk -F'\\t' -v W={window_bp} -v K={stride} {shlex.quote(AWK)}")
    # pipefail: without it the pipeline reports awk's status, so a TRUNCATED
    # .gz (gzip: unexpected end of file) exits 0 and silently yields a partial
    # chromosome. Found on NA18508's local extract, which returned 2,384 windows
    # against ~16k for complete files.
    proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd],
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    pend = {}
    for line in proc.stdout:
        rid, pos, q, code = line.rstrip("\n").split("\t")
        key = (rid, int(pos))
        d = pend.setdefault(key, {})
        d[code] = float(q)
        if len(d) < 2:
            continue
        del pend[key]
        p_m, p_h = d.get("m", 0.0), d.get("h", 0.0)
        conf, which = max((1.0 - p_m - p_h, "c"), (p_m, "m"), (p_h, "h"))
        stats["calls"] += 1
        if conf < CONF_THRESH:
            stats["ambiguous"] += 1
            continue
        if which == "h":
            stats["hmc"] += 1
            if H_POLICY == "nocall":
                continue
        stats["kept"] += 1
        stats["err_sum"] += 1.0 - conf
        yield rid, key[1], (1 if which in ("m", "h") else 0)
    rc = proc.wait()
    # 141 = SIGPIPE, expected only when head cut the stream short on purpose
    if rc != 0 and not (max_rows and rc == 141):
        raise RuntimeError(f"stream failed (exit {rc}; truncated input?): {cmd}")


def window_metrics(pos, m):
    n_reads, n_cpg = m.shape
    h = n_reads // 2
    obs = np.isfinite(m)
    mu = np.nanmean(m, axis=0)
    c = np.where(obs, m - mu, 0.0)
    # variance-weighted pooled correlation over all CpG pairs (decay._accumulate)
    sxy = sxx = syy = 0.0
    for i in range(n_cpg - 1):
        both = obs[:, i:i + 1] & obs[:, i + 1:]
        sxy += float((c[:, i:i + 1] * c[:, i + 1:] * both).sum())
        sxx += float(((c[:, i:i + 1] ** 2) * both).sum())
        syy += float(((c[:, i + 1:] ** 2) * both).sum())
    icc_w = sxy / np.sqrt(sxx * syy) if sxx > 0 and syy > 0 else np.nan

    ok = obs.sum(axis=1) >= 2
    sub, subobs = m[ok], obs[ok]
    C_r = subobs.sum(axis=1)
    cbar = float(C_r.mean()) if ok.any() else np.nan
    frac = np.nanmean(sub, axis=1) if ok.sum() >= 3 else np.array([])
    v = mu * (1 - mu)
    if len(frac) >= 3:
        exp_het = float(np.mean((subobs * v).sum(axis=1) / C_r ** 2))
        de_het = float(np.var(frac, ddof=1)) / exp_het if exp_het > 0 else np.nan
    else:
        de_het = np.nan
    return dict(
        start=int(pos[0]), n_reads=n_reads, n_cpgs=n_cpg,
        span=int(pos[-1] - pos[0]), meth=float(np.nanmean(m)),
        de_split=design_effect(m[:h], m[h:]),
        de_one=design_effect(m, m[:0]),
        icc_w=icc_w, cbar=cbar,
        de_icc=1.0 + (cbar - 1.0) * icc_w if np.isfinite(icc_w) else np.nan,
        de_het=de_het,
    )


def main():
    path, out = sys.argv[1], sys.argv[2]
    window_bp = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    max_rows = int(float(sys.argv[5])) if len(sys.argv) > 5 else 0
    name = os.path.basename(path).split("_modkit")[0]

    stats = defaultdict(float)
    win = defaultdict(lambda: defaultdict(dict))
    for rid, pos, call in stream(path, window_bp, stride, max_rows, stats):
        win[pos // window_bp][rid][pos] = call
    mats = to_matrices(win)
    del win

    rows = [window_metrics(p, m) for p, m in mats]
    cols = list(rows[0].keys()) if rows else []
    with open(out, "w") as fh:
        fh.write(f"# sample={name} window_bp={window_bp} stride={stride} "
                 f"max_rows={max_rows} conf_thresh={CONF_THRESH} h_policy={H_POLICY}\n")
        fh.write(f"# calls={int(stats['calls'])} ambiguous={int(stats['ambiguous'])} "
                 f"hmc={int(stats['hmc'])} kept={int(stats['kept'])} "
                 f"eps={stats['err_sum'] / max(stats['kept'], 1):.5f}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])
                               for k in cols) + "\n")
    de = np.array([r["de_split"] for r in rows])
    print(f"{name}\twindows={len(rows)}\tmedian_DE={np.median(de):.3f}\t"
          f"eps={stats['err_sum'] / max(stats['kept'], 1):.4f}", flush=True)


if __name__ == "__main__":
    main()
