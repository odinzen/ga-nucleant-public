"""Graphical abstract: stability-filtered candidates by disregistry to alpha-Ga (010).

Clean rebuild of the TOC graphic. The earlier version had the legend sitting on top of
the TeO2 bar and a finding-line baked under the plot. This puts the legend in the empty
upper-right (the HfN/ScN bars are tiny, so nothing is there), carries real subscripts,
and bakes no caption or takeaway onto the image; that line belongs in the graphical-
abstract text field, not on the figure.

Disregistry values match Table 2 / Fig. 3.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# (label, disregistry %, measured-on-Ga, untested-prediction)
STABLE = [
    (r"HfN",         0.05, False, True),
    (r"ScN",         0.49, False, True),
    (r"VO$_2$",      1.00, False, False),
    (r"NbC",         1.17, False, False),
    (r"ZrN",         1.19, True,  False),
    (r"TaC",         1.40, False, False),
    (r"TiO$_2$",     1.57, False, False),
    (r"HfC",         2.54, True,  False),
]
REDUCED = [(r"TeO$_2$", 6.30)]

# top-to-bottom = best match first, reduced oxide last
rows = STABLE + [(l, d, None, None) for l, d in REDUCED]
labels = [r[0] for r in rows]
vals = [r[1] for r in rows]
y = list(range(len(rows)))[::-1]          # HfN at top

fig, ax = plt.subplots(figsize=(4.2, 3.0))

for yi, (lab, d, meas, pred) in zip(y, rows):
    reduced = meas is None
    ax.barh(yi, d, height=0.66, zorder=3,
            color="white" if reduced else "0.25",
            edgecolor="k", linewidth=0.7,
            hatch="///" if reduced else None)
    if meas:
        ax.text(d + 0.08, yi, "*", fontsize=11, va="center", ha="left", fontweight="bold")
    if pred:
        ax.text(d + 0.12, yi, "untested", fontsize=6.5, va="center", ha="left",
                color="0.35", style="italic")

# 3% screen threshold, labelled at the top of the line where the short HfN/ScN rows
# leave the space clear (bottom would collide with the x-axis ticks)
ax.axvline(3.0, color="k", linestyle=":", linewidth=0.9, zorder=1)
ax.text(3.15, len(rows) - 0.75, "3% screen threshold", fontsize=6.5,
        ha="left", va="center", color="0.4")

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel(r"Lattice disregistry to $\alpha$-Ga (010)  (%)", fontsize=9)
ax.set_xlim(0, 7.2)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.tick_params(axis="x", labelsize=7)
ax.spines[["top", "right"]].set_visible(False)

# legend in the empty upper-right (HfN/ScN bars are ~0.05-0.5, so the top rows are clear)
handles = [
    mpatches.Patch(facecolor="0.25", edgecolor="k", linewidth=0.7, label="Stable in liquid Ga"),
    mpatches.Patch(facecolor="white", edgecolor="k", linewidth=0.7, hatch="///",
                   label="Reduced by Ga"),
    plt.Line2D([0], [0], marker="$*$", color="k", linestyle="none", markersize=9,
               label="Measured on Ga"),
]
ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right",
          handlelength=1.2, handletextpad=0.5, borderaxespad=0.6)

plt.tight_layout(pad=0.5)
out_png = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ga_nucleant_TOC_graphic.png")
plt.savefig(out_png, dpi=600, bbox_inches="tight")
plt.savefig(out_png.replace(".png", ".tif"), dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("saved", os.path.normpath(out_png))
