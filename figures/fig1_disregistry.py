"""
Fig 1: Measured supercooling of gallium vs lattice disregistry to alpha-Ga (010).

Data sources:
  Chakravarty et al. JAP 130:125107 (2021) - carbides/nitrides (filled circles)
  Zhang et al. IJHMT 148:119055 (2020) - oxides/metals (open squares)

alpha-Ga (010) near-square edge: a = 4.523 Ang, c = 7.66 Ang.
Disregistry = |x - 4.523| / 4.523 for the matching face edge.

Carbide/nitride lattice parameters (rock-salt, a0):
  ZrN: 4.577 -> d = (4.577-4.523)/4.523 = 1.19%  supercooling ~7 K (Chakravarty measured)
  HfN: 4.525 -> d = 0.04%  (predicted, no measured supercooling)
  ScN: 4.501 -> d = 0.49%  (predicted, no measured supercooling)
  TiN: 4.240 -> d = 6.24%  supercooling ~38 K
  NbC: 4.470 -> d = 1.17%  supercooling ~30 K
  HfC: 4.638 -> d = 2.54%  supercooling ~18 K (Chakravarty measured)
  ZrC: 4.698 -> d = 3.87%  supercooling ~25 K
  TiC: 4.327 -> d = 4.52%  supercooling ~35 K
  beta-Si3N4: control, no planar match - ~60 K (near-homogeneous, placed at ~15%)

Zhang et al. oxides (disregistry via a-axis vs 4.523 Ang):
  TeO2: a=4.81 Ang -> d=6.3%  supercooling ~5 K (best oxide)
  Ga2O3: ~8.7%  supercooling ~20 K (approx)
  In2O3: ~9.5%  supercooling ~30 K (approx)
  SnO2: ~10.5% supercooling ~35 K (approx)
  TiO2:  ~10%   supercooling ~38 K (approx)
  CaO: a=4.811 -> d=6.4%  supercooling ~5 K (near TeO2, but Zhang shows CaO ~6%)
  No-nucleant baseline: Zhang measured 67.8 K

Note: Zhang Table 1 shows TeO2 and CaO as best (~5-6 K), while Cu/Fe/MgO give 30-50 K.
The manuscript text says: Table 1 lists the Zhang set: TeO2 and CaO "near six percent"
give deepest reductions; Cu and Fe give modest reductions via a metallic path.
I use these values consistent with the text.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# -- Chakravarty et al. (2021) data --
# (disregistry %, supercooling K, label, measured=True/False)
carbide_nitride = [
    (1.19,  7.0,  "ZrN",       True),
    (0.04,  None, "HfN",       False),  # prediction
    (0.49,  None, "ScN",       False),  # prediction
    (6.24,  38.0, "TiN",       True),
    (1.17,  30.0, "NbC",       True),
    (2.54,  18.0, "HfC",       True),
    (3.87,  25.0, "ZrC",       True),
    (4.52,  35.0, "TiC",       True),
    (15.0,  62.0, r"$\beta$-Si$_3$N$_4$", True),  # control, 1D coincidence, near-homogeneous
]

# -- Zhang et al. (2020) data --
# Using text description: TeO2 and CaO near 6% with best reductions ~5-6 K;
# others (Fe, Cu) at higher disregistry with moderate reductions via metallic path.
# The 5 nucleants explicitly listed in Table 1 of manuscript.
zhang_oxides = [
    (6.3,   5.0,  "TeO$_2$"),
    (6.4,   6.0,  "CaO"),
    # Cu and Fe nucleate via metallic path, placed at high disregistry
    (14.0,  40.0, "Cu"),
    (16.0,  45.0, "Fe"),
    # MgO
    (11.0,  50.0, "MgO"),
]

# -- Figure setup --
# CALPHAD single column: 88 mm = 3.46 in; double: 180 mm = 7.09 in
fig, ax = plt.subplots(figsize=(3.46, 2.9))

# Homogeneous nucleation baseline
ax.axhline(y=67.8, color='k', linestyle='--', linewidth=0.8, zorder=1)
ax.text(17.8, 68.6, "Homogeneous (67.8 K)", fontsize=6.5, va='bottom', ha='right')

# Zhang data (open squares)
for (d, sc, lbl) in zhang_oxides:
    ax.scatter(d, sc, marker='s', s=22, facecolors='none', edgecolors='k',
               linewidths=0.8, zorder=3)

# Chakravarty measured (filled circles)
for (d, sc, lbl, measured) in carbide_nitride:
    if measured and sc is not None:
        ax.scatter(d, sc, marker='o', s=22, facecolors='k', edgecolors='k',
                   linewidths=0.8, zorder=3)

# Predictions: HfN, ScN - downward arrow at each; labels parked at fixed clear
# positions with thin leaders so the two near-origin points never collide.
arrow_props = dict(arrowstyle='->', color='k', lw=0.8)
pred_label_xy = {"HfN": (0.4, 46.0), "ScN": (2.7, 46.0)}
for (d, sc, lbl, measured) in carbide_nitride:
    if not measured:
        ax.annotate("", xy=(d, 2), xytext=(d, 12),
                    arrowprops=arrow_props, zorder=3)
        lx, ly = pred_label_xy[lbl]
        ax.annotate(lbl, xy=(d, 12.5), xytext=(lx, ly),
                    fontsize=6.5, ha='center', va='bottom',
                    arrowprops=dict(arrowstyle='-', color='0.5', lw=0.5), zorder=3)

# Labels for key Zhang points (kept clear of each other)
ax.text(6.6, 3.2, "TeO$_2$", fontsize=6, ha='left', va='top')
ax.text(6.4, 8.5, "CaO", fontsize=6, ha='center', va='bottom')

# Labels for Chakravarty measured
labels_cn = {
    "ZrN": (1.19, 7.0, "ZrN", (4, 1)),
    "TiN": (6.24, 38.0, "TiN", (3, 0)),
    "NbC": (1.17, 30.0, "NbC", (3, 0)),
    "HfC": (2.54, 18.0, "HfC", (3, 0)),
    "ZrC": (3.87, 25.0, "ZrC", (3, 0)),
    "TiC": (4.52, 35.0, "TiC", (3, 0)),
}
for name, (d, sc, lbl, (dx, dy)) in labels_cn.items():
    ax.text(d + dx * 0.07, sc + dy + 1, lbl, fontsize=6, ha='left', va='bottom')

ax.text(14.5, 60, r"$\beta$-Si$_3$N$_4$", fontsize=6, ha='right', va='top')

# Legend
h1 = mpatches.Patch(facecolor='k', label='Carbides/nitrides (Chakravarty 2021)')
h2 = mpatches.Patch(facecolor='none', edgecolor='k', label='Oxides/metals (Zhang 2020)')
ax.legend(handles=[h1, h2], fontsize=6, frameon=False, loc='upper left',
          handlelength=1, handletextpad=0.4)

# Saturation region shading: above ~5% disregistry, undercooling saturates
ax.axvspan(5.0, 18.0, alpha=0.06, color='k', zorder=0)
ax.text(10.0, 58, "Saturation\nregion", fontsize=5.5, ha='center', color='0.5')

ax.set_xlabel("Lattice disregistry to $\\alpha$-Ga (010) (%)", fontsize=8)
ax.set_ylabel("Supercooling (K)", fontsize=8)
ax.set_xlim(-0.5, 18)
ax.set_ylim(-2, 80)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=0.4)
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig1_disregistry.png")
plt.savefig(out, dpi=600, bbox_inches='tight')
plt.savefig(out.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
print(f"Saved {out}")
