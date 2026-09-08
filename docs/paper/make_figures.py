"""Figures for docs/paper/method.tex, from the real HG00146 chr15 extract.

Palette: slots 1-3 of the reference categorical order (blue/orange/aqua), used
in fixed order. Print target, so light mode only and no interaction layer.
"""
import os, pickle, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts", "phasing_qc"))
import P08_decay_from_extract as P08
from ont_asm_caller.decay import comethylation_decay

EXTRACT = os.path.join(HERE, "..", "..", "..", "..", "modkit",
                       "HG00146_chr15_modkit_extract.tsv.gz")
CACHE = os.path.join(HERE, "figures", "_cache.pkl")
ROWS = 12_000_000

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#b8b7b2"

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#e8e7e3", "grid.linewidth": 0.5,
})


def de_of(m):
    ok = np.isfinite(m)
    p = np.nanmean(m)
    if ok.sum() < 24 or not (0.02 < p < 0.98) or m.shape[0] < 8:
        return None
    C = ok.sum(axis=1).mean()
    return float(np.var(np.nanmean(m, axis=1), ddof=1) / (p * (1 - p) / C))


def gather():
    if os.path.exists(CACHE):
        return pickle.load(open(CACHE, "rb"))
    out = {}
    P08.MIN_CPGS = 3
    for wbp in (500, 1000, 2000):
        win, confs, _ = P08.build_windows(EXTRACT, ROWS, wbp)
        mats = P08.to_matrices(win)
        des = [d for _, m in mats if (d := de_of(m)) is not None]
        out[wbp] = {"de": np.array(des),
                    "eps": float(np.mean(1 - confs)),
                    "ncpg": float(np.median([len(p) for p, _ in mats]))}
        if wbp == 2000:
            edges = np.array([0, 2, 4, 8, 15, 25, 40, 60, 90, 140,
                              200, 300, 450, 700, 1000, 1500, 2100], float)
            c = comethylation_decay(mats, edges=edges, min_pairs=200)
            out["curve"] = (c.centres, c.corr, c.decay_bp, c.amplitude, c.r2)
        print(f"  {wbp}bp: {len(des)} windows", flush=True)
    pickle.dump(out, open(CACHE, "wb"))
    return out


def fig_curve(d):
    centres, corr, decay, amp, r2 = d["curve"]
    eps = d[2000]["eps"]
    f = (1 - 2 * eps) ** 2
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.grid(True, which="major", axis="y", zorder=0)
    ax.plot(centres, np.minimum(corr / f, 1), color=BLUE, lw=2, zorder=3,
            label="measured", solid_capstyle="round")
    grid = np.geomspace(2, 2100, 200)
    ax.plot(grid, amp / f * np.exp(-grid / decay), color=ORANGE, lw=2, ls="--",
            zorder=2, label=f"single exponential ($r^2$={r2:.2f})")
    ax.set_xscale("log")
    ax.set_xlabel("CpG–CpG distance (bp)")
    ax.set_ylabel("within-read correlation")
    ax.set_ylim(0, 0.85)
    ax.annotate("nucleosome\nrepeat ~167 bp", xy=(167, 0.191 / f),
                xytext=(230, 0.42), fontsize=6.5, color=INK2,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(HERE, "figures", "decay_curve.pdf"))
    fig.savefig(os.path.join(HERE, "figures", "_check_decay.png"), dpi=200)
    plt.close(fig)


def fig_de(d):
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.grid(True, which="major", axis="y", zorder=0)
    for (wbp, col) in ((500, BLUE), (1000, ORANGE), (2000, AQUA)):
        de = np.sort(d[wbp]["de"])
        ax.plot(de, np.arange(1, len(de) + 1) / len(de), color=col, lw=2,
                zorder=3, solid_capstyle="round",
                label=f"{wbp} bp (~{d[wbp]['ncpg']:.0f} CpGs)")
    # reference lines labelled along the TOP: the legend occupies the lower
    # right, which is the only region the ECDFs leave empty
    for x, lab in ((2.40, "pooled FDR 0.07"), (5.49, "pooled FDR 0.70")):
        ax.axvline(x, color=MUTED, lw=0.8, ls=":", zorder=1)
        ax.text(x * 0.93, 0.97, lab, fontsize=6.5, color=INK2, rotation=90,
                ha="right", va="top")
    ax.set_xscale("log")
    ax.set_xlim(0.6, 20)
    ax.set_xlabel("design effect (per region)")
    ax.set_ylabel("cumulative fraction of regions")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, loc="lower right", title="window")
    ax.get_legend().get_title().set_fontsize(7)
    fig.tight_layout(pad=0.3)
    fig.savefig(os.path.join(HERE, "figures", "design_effect.pdf"))
    fig.savefig(os.path.join(HERE, "figures", "_check_de.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    d = gather()
    fig_curve(d)
    fig_de(d)
    print("wrote figures/decay_curve.pdf, figures/design_effect.pdf")
    for w in (500, 1000, 2000):
        de = d[w]["de"]
        print(f"  {w}bp  n={len(de)}  median={np.median(de):.2f}  "
              f"90th={np.percentile(de,90):.2f}  %>5={(de>5).mean():.1%}")
