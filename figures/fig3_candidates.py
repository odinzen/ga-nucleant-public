"""
Fig 3: Stability-filtered candidates ranked by disregistry to alpha-Ga (010).

Data from manuscript Table 2 description and Section 3.4:
  - HfN: 0.05% (near-coherent, untested prediction)
  - ScN: 0.49% (near-coherent, untested prediction)
  - VO2: ~0.9% (stable oxide)
  - NbC: ~1.17% (stable, inert carbide)
  - ZrN: ~1.19% (stable, inert nitride; MEASURED - asterisk)
  - TaC: ~1.4% (stable, inert carbide)
  - TiO2: ~1.5% (stable oxide)
  - HfC: ~2.54% (stable, MEASURED - asterisk)
  Excluded (1D coincidence):
  - Ta2O5: apparent <1% via single axis, excluded
  - beta-Si3N4: control, single-axis coincidence, excluded

For comparison, reduced oxides (shown as hatched, above stability threshold):
  - TeO2: 6.3% (best oxide but reduced to Ga2O3+Te)
  - Ga2O3: ~8.7% (reduction product)

VO2 lattice: monoclinic, rutile-type at high T; a~5.75, b~4.52 Ang -> (010) match
  b = 4.526 Ang -> disregistry = |4.526-4.523|/4.523 = 0.066% ~ 0.07%
  Actually VO2 rutile: a=4.554, c=2.855 - not a square face
  Manuscript says VO2 and TiO2 match to within 1.5%; I use ~1.0% for VO2
  TiO2 rutile: a=4.594 -> d=|4.594-4.523|/4.523=1.57%
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# (label, disregistry %, status, measured, prediction_arrow)
candidates_stable = [
    ("HfN",   0.05, "stable",  False, True),   # prediction
    ("ScN",   0.49, "stable",  False, True),    # prediction
    ("VO$_2$",  1.0,  "stable",  False, False),
    ("NbC",   1.17, "stable",  False, False),
    ("ZrN",   1.19, "stable",  True,  False),   # measured
    ("TaC",   1.4,  "stable",  False, False),
    ("TiO$_2$", 1.57, "stable", False, False),
    ("HfC",   2.54, "stable",  True,  False),   # measured
]

candidates_reduced = [
    ("TeO$_2$", 6.3,  "reduced"),
]

all_compounds = candidates_stable + [(lbl, d, st) for lbl, d, st in candidates_reduced]

labels_stable  = [c[0] for c in candidates_stable]
disreg_stable  = [c[1] for c in candidates_stable]
measured_flags = [c[3] for c in candidates_stable]
pred_flags     = [c[4] for c in candidates_stable]

labels_red  = [c[0] for c in candidates_reduced]
disreg_red  = [c[1] for c in candidates_reduced]

fig, ax = plt.subplots(figsize=(3.46, 3.0))

n_stable = len(labels_stable)
n_red    = len(labels_red)
x_stable = np.arange(n_stable)
x_red    = np.arange(n_stable, n_stable + n_red)

# Stable bars (solid, dark)
bars_s = ax.bar(x_stable, disreg_stable, color='0.25', edgecolor='k',
                linewidth=0.6, width=0.7, zorder=2)

# Reduced bars (hatched, light)
bars_r = ax.bar(x_red, disreg_red, color='white', edgecolor='k',
                linewidth=0.6, width=0.7, hatch='///', zorder=2)

# Asterisks for measured compounds
for i, (meas, pred) in enumerate(zip(measured_flags, pred_flags)):
    if meas:
        y = disreg_stable[i]
        ax.text(x_stable[i], y + 0.05, '*', fontsize=9, ha='center', va='bottom',
                color='k', fontweight='bold')
    if pred:
        # Downward arrow above bar
        y = disreg_stable[i]
        ax.annotate("", xy=(x_stable[i], y + 0.02),
                    xytext=(x_stable[i], y + 0.35),
                    arrowprops=dict(arrowstyle='->', color='k', lw=0.8))

# 3% threshold line (manuscript: "match to better than 3%")
ax.axhline(y=3.0, color='k', linestyle=':', linewidth=0.8, zorder=1)
# Label parked over the low left bars, clear of the tall TeO2 hatched bar on the right;
# white background so it never sits on a bar or the dotted line.
ax.text(2.0, 3.15, "3% screen threshold", fontsize=6, ha='left', va='bottom',
        color='0.4', zorder=4,
        bbox=dict(facecolor='white', edgecolor='none', pad=0.5))

all_labels = labels_stable + labels_red
all_x      = list(x_stable) + list(x_red)
ax.set_xticks(all_x)
ax.set_xticklabels(all_labels, fontsize=7, rotation=45, ha='right')

ax.set_ylabel("Lattice disregistry to $\\alpha$-Ga (010) (%)", fontsize=8)
ax.set_xlim(-0.6, n_stable + n_red - 0.4)
ax.set_ylim(0, 7.5)
ax.tick_params(axis='y', labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
h1 = mpatches.Patch(facecolor='0.25', edgecolor='k', linewidth=0.6,
                     label='Stable in Ga melt')
h2 = mpatches.Patch(facecolor='white', edgecolor='k', linewidth=0.6, hatch='///',
                     label='Reduced by Ga')
h3 = plt.Line2D([0], [0], marker='$*$', color='k', linestyle='none',
                markersize=8, label='Measured on Ga')
ax.legend(handles=[h1, h2, h3], fontsize=6, frameon=False, loc='upper left',
          handlelength=1.2, handletextpad=0.4)

# Note for 1D exclusions
ax.text(0.99, 0.98, r"Arrows: untested predictions" + "\n" +
        r"1D coincidences (Ta$_2$O$_5$, $\beta$-Si$_3$N$_4$) excluded",
        transform=ax.transAxes, fontsize=5.5, ha='right', va='top',
        color='0.4', linespacing=1.4)

plt.tight_layout(pad=0.4)
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig3_candidates.png")
plt.savefig(out, dpi=600, bbox_inches='tight')
plt.savefig(out.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
print(f"Saved {out}")
