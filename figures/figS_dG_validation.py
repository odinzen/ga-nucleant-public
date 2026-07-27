"""SI figure: exact CALPHAD driving force for solidification of pure gallium vs
the linear Turnbull form, with the percent deviation on a second axis. Data from
the Ga-In-Sn TDB (G_LIQUID - G_ORTHORHOMBIC_GA for pure Ga).

This is the one script here that needs the proprietary Ga-In-Sn CALPHAD database
(not distributed), because its sole purpose is to check the approximation used
everywhere else. The main results use the linear Turnbull driving force, which needs
no database and is reproducible from the paper's equations; this figure shows that
linear form stays within a few percent of the exact CALPHAD driving force over the
relevant undercooling range. The reported results therefore do not depend on the
proprietary database."""
import os
import numpy as np
import matplotlib.pyplot as plt
from pycalphad import Database, calculate

db = Database(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "..", "gainsn-calphad", "databases", "GaInSn.tdb"))
T = np.linspace(200.0, 315.0, 461)
P = 101325.0
Gliq = np.asarray(calculate(db, ["GA", "IN", "SN", "VA"], "LIQUID", T=T, P=P,
                            points=np.array([[1.0, 0.0, 0.0]])).GM).squeeze()
Gsol = np.asarray(calculate(db, ["GA", "VA"], "ORTHORHOMBIC_GA", T=T, P=P,
                            points=np.array([[1.0]])).GM).squeeze()
dG = Gliq - Gsol
Tm = float(np.interp(0.0, dG[::-1], T[::-1]))
dSf = -float(np.interp(Tm, T, np.gradient(dG, T)))

dTs = np.linspace(0.5, 74.0, 300)
exact = np.interp(Tm - dTs, T, dG)
linear = dSf * dTs
dev = 100.0 * (linear - exact) / exact

fig, ax = plt.subplots(figsize=(3.46, 2.7))
ax.plot(dTs, exact, "k-", lw=1.3, label="Exact CALPHAD")
ax.plot(dTs, linear, color="0.55", ls="--", lw=1.1, label="Linear Turnbull")
ax.axvline(67.8, color="k", lw=0.6, ls=":")
ax.text(67.0, 120, "homogeneous\nundercooling", fontsize=5.5, ha="right",
        va="bottom", color="0.35", linespacing=1.3)
ax.set_xlabel("Undercooling $\\Delta T$ (K)", fontsize=8)
ax.set_ylabel("Driving force $\\Delta G$ (J mol$^{-1}$)", fontsize=8)
ax.set_xlim(0, 75)
ax.set_ylim(0, 1400)
ax.tick_params(labelsize=7)
ax.spines["top"].set_visible(False)
ax.legend(fontsize=6.5, frameon=False, loc="upper left", handlelength=1.6)

ax2 = ax.twinx()
ax2.plot(dTs, dev, color="0.55", ls=":", lw=1.0)
ax2.set_ylabel("Linear over exact, deviation (%)", fontsize=7, color="0.4")
ax2.set_ylim(0, 5)
ax2.tick_params(labelsize=7, colors="0.4")
ax2.spines["top"].set_visible(False)
ax2.text(74, dev[-1] + 0.15, f"{dev[-1]:.1f}%", fontsize=6, ha="right",
         va="bottom", color="0.4")

plt.tight_layout(pad=0.4)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figS_dG_validation.png")
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"Saved {out}  (Tm={Tm:.2f} K, dSf={dSf:.2f} J/mol-K, max dev={dev[-1]:.2f}%)")
