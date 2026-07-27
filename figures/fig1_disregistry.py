"""
Fig 1: Measured supercooling of gallium vs lattice disregistry to alpha-Ga (010).

Canonical data source: Table_S2_nucleant_screen.csv (this repo) for the carbide/nitride
set, and Table 1 of the manuscript for the Zhang oxide/metal set. The carbide/nitride
undercoolings are READ from the CSV so the figure can never drift from the table again.

  Chakravarty et al. JAP 130:125107 (2021) - carbides/nitrides (filled circles)
  Zhang et al. IJHMT 148:119055 (2020)     - oxides/metals (open squares, Table 1)

alpha-Ga (010) near-square edge: a = 4.523 Ang. Disregistry = |x - 4.523| / 4.523.
Upper-bound undercoolings in the CSV ("<10", "<20") are plotted at the bound.
"""
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, "..", "Table_S2_nucleant_screen.csv")


def parse_uc(s):
    s = s.strip()
    if not s:
        return None
    return float(s.lstrip("<")), s.startswith("<")


# Carbides/nitrides from Table S2 (filled circles); HfN/ScN carry no undercooling.
measured_cn, predictions = [], []
with open(CSV_PATH, encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        if row["class"] not in ("nitride", "carbide"):
            continue
        d = float(row["disregistry_010_percent"])
        uc = parse_uc(row["measured_undercooling_K"])
        name = row["compound"]
        if uc is None:
            if name in ("HfN", "ScN"):
                predictions.append((d, name))
        else:
            measured_cn.append((d, uc[0], uc[1], name))

# Zhang et al. (2020) oxides/metals, manuscript Table 1 (open squares). Not in Table S2.
zhang = [
    (6.3, 38.2, "TeO$_2$"),
    (6.4, 44.6, "CaO"),
    (6.9, 53.7, "MgO"),
    (20.1, 52.5, "Cu"),
    (36.6, 53.8, "Fe"),
]

fig, ax = plt.subplots(figsize=(4.6, 3.1))

ax.axhline(y=67.8, color="k", linestyle="--", linewidth=0.8, zorder=1)
ax.text(39.6, 65.6, "Homogeneous (67.8 K)", fontsize=6.5, va="top", ha="right")

ax.axvspan(5.0, 40.0, alpha=0.05, color="k", zorder=0)
ax.text(30.0, 14, "Saturation region", fontsize=6, ha="center", color="0.55")

for d, uc, lbl in zhang:
    ax.scatter(d, uc, marker="s", s=24, facecolors="none", edgecolors="k",
               linewidths=0.8, zorder=3)
for d, uc, bound, name in measured_cn:
    ax.scatter(d, uc, marker="o", s=24, facecolors="k", edgecolors="k",
               linewidths=0.8, zorder=3)
    if bound:
        ax.annotate("", xy=(d, uc - 4.0), xytext=(d, uc - 0.5),
                    arrowprops=dict(arrowstyle="->", color="0.5", lw=0.7), zorder=2)

# Predictions HfN, ScN: short down-arrows near zero, labels fanned up with leaders.
pred_xy = {"HfN": (0.7, 44.0), "ScN": (3.0, 44.0)}
for d, name in predictions:
    ax.annotate("", xy=(d, 2), xytext=(d, 11),
                arrowprops=dict(arrowstyle="->", color="k", lw=0.8), zorder=3)
    lx, ly = pred_xy[name]
    ax.annotate(name, xy=(d, 11.5), xytext=(lx, ly), fontsize=6.5, ha="center",
                va="center", arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5), zorder=3)

# Point labels: explicit (x, y, ha) so nothing overlaps. beta-Si3N4's CSV name has
# unicode subscripts; render it through mathtext.
si_name = next(n for _, _, _, n in measured_cn if "Si" in n)
cn_labels = {
    "ZrN": (2.3, 10.0, "left"),
    "HfC": (3.4, 20.0, "left"),
    "ZrC": (4.7, 30.0, "left"),
    "NbN": (1.9, 37.0, "right"),
    "TiC": (3.4, 59.5, "right"),
    "TiN": (7.1, 63.0, "left"),
    si_name: (24.23, 73.0, "center"),
}
display = {si_name: r"$\beta$-Si$_3$N$_4$"}
for d, uc, bound, name in measured_cn:
    lx, ly, ha = cn_labels[name]
    ax.text(lx, ly, display.get(name, name), fontsize=6.5, ha=ha, va="center")

zhang_labels = {
    "TeO$_2$": (7.1, 38.2, "left"),
    "CaO": (7.1, 44.6, "left"),
    "MgO": (7.7, 53.7, "left"),
    "Cu": (20.1, 56.5, "center"),
    "Fe": (35.2, 55.8, "right"),
}
for d, uc, lbl in zhang:
    lx, ly, ha = zhang_labels[lbl]
    ax.text(lx, ly, lbl, fontsize=6.5, ha=ha, va="center")

h1 = mpatches.Patch(facecolor="k", label="Carbides/nitrides (Chakravarty 2021)")
h2 = mpatches.Patch(facecolor="none", edgecolor="k", label="Oxides/metals (Zhang 2020)")
ax.legend(handles=[h1, h2], fontsize=6.5, frameon=False, loc="upper left",
          handlelength=1, handletextpad=0.4, borderaxespad=0.4)

ax.set_xlabel("Lattice disregistry to $\\alpha$-Ga (010) (%)", fontsize=8.5)
ax.set_ylabel("Supercooling (K)", fontsize=8.5)
ax.set_xlim(-1, 40)
ax.set_ylim(-2, 80)
ax.tick_params(labelsize=7.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout(pad=0.4)
out = os.path.join(HERE, "fig1_disregistry.png")
plt.savefig(out, dpi=600, bbox_inches="tight")
plt.savefig(out.replace(".png", ".tif"), dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("saved")
