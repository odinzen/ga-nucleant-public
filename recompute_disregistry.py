"""Two-axis (010) planar disregistry from EXPERIMENTAL lattice constants.

Reviewer point 3: the screen computed a single best edge (1D) while the text claimed a
planar (two-axis) match. This applies the honest two-axis rule and changes ONLY that; the
lattice source stays experimental, as in the original screen. Disregistry for a sub-1%
lattice match must not be computed from PBE-relaxed (Materials Project) cells, whose ~1%
error is comparable to the quantity itself.

Cubic and tetragonal candidates present a square (010)/(001) face, so the two-axis value
uses their single experimental lattice constant a; that a is recovered from the original
ICSD-based single-edge disregistry (a = 4.523 * (1 + sign * d_1D/100), sign from whether the
phase is larger or smaller than alpha-Ga), which ties every value to the paper's own verified
data. The three DFT nitrides cross-check against dft/build_qe_inputs.py (HfN 4.525, ScN 4.501,
ZrN 4.577). Low-symmetry phases (Ta2O5, Nb2O5, beta-Si3N4, and the reduction products) have no
square face and are taken from Materials Project exactly as the original broader search did;
all are excluded or are poor templates regardless.

    python recompute_disregistry.py            # print numbers
    python recompute_disregistry.py --write     # also rewrite Table_S2_nucleant_screen.csv
"""
import argparse, csv, itertools, os
import numpy as np
from scipy.stats import spearmanr, pearsonr
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mp_api.client import MPRester

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = os.environ.get("MP_API_KEY")  # a Materials Project API key; never write it into a file
GA_MEAN, A_GA, C_GA = 4.523, 4.5197, 4.5257     # alpha-Ga (010) face: near-square, two edges

# cubic/tetragonal: (original 1D disregistry, sign of (a - 4.523)), a square (010)/(001) face.
# a is recovered as 4.523*(1 + sign*d1D/100); two-axis d = max mismatch to the two Ga edges.
CUBIC_TET = {  # compound: (class, d_1D_orig, sign, density, measured_uc_K, note)
    "HfN":  ("nitride", 0.04, +1, 13.8, None, ""),
    "ScN":  ("nitride", 0.49, -1, 4.29, None, ""),
    "VO2":  ("oxide",   0.69, +1, 4.57, None, "Stable polymorph (rutile)"),
    "NbC":  ("carbide", 1.17, -1, 7.82, None, ""),
    "ZrN":  ("nitride", 1.22, +1, 7.29, 10.0, "Upper bound"),
    "TaC":  ("carbide", 1.48, -1, 14.5, None, ""),
    "TiO2": ("oxide",   1.57, +1, 4.25, None, "Stable polymorph (rutile)"),
    "HfC":  ("carbide", 2.54, +1, 12.7, 20.0, "Upper bound"),
    "NbN":  ("nitride", 2.90, -1, 8.47, 36.5, ""),
    "ZrC":  ("carbide", 3.87, +1, 6.73, 30.0, ""),
    "TiC":  ("carbide", 4.31, -1, 4.93, 59.5, ""),
    "TiN":  ("nitride", 6.21, -1, 5.22, 63.0, ""),
    "TeO2": ("oxide",   6.35, +1, 6.02, 38.0, "Reduced in Ga melt (Zhang 2020)"),
    "VN":   ("nitride", 8.49, -1, 6.13, None, ""),
}
LOWSYM = {  # low-symmetry: full cell from Materials Project (as the original broader search)
    "Ta2O5": ("oxide",   "mp-1539317", 8.20, None, "No planar match, excluded"),
    "Nb2O5": ("oxide",   "mp-680944",  4.60, None, "No planar match, excluded"),
    "Si3N4": ("nitride", "mp-2245",    3.19, 68.9, "Failed control (Chakravarty 2021)"),
}
# Zhang measured set for Table 1 / Fig 1 (all cubic): (d_1D_orig, sign, measured_uc)
TABLE1 = {"TeO2": (6.35, +1, 38.2), "CaO": (6.4, +1, 44.6), "MgO": (6.9, -1, 53.7),
          "Cu": (20.1, -1, 52.5), "Fe": (36.6, -1, 53.8)}
PRODUCTS = {"Te": None, "GaTe": None, "Ga2O3": None}   # Section 2.3, from MP


def two_axis_from_a(a):
    return max(abs(a - A_GA) / A_GA, abs(a - C_GA) / C_GA) * 100.0


def recover_a(d1d, sign):
    return GA_MEAN * (1 + sign * d1d / 100.0)


def two_axis_from_cell(lat, mult=(1,), ang_tol=8.0):
    a, b, c = lat.a, lat.b, lat.c
    best = np.inf
    for (e1, e2), phi in [((a, b), lat.gamma), ((a, c), lat.beta), ((b, c), lat.alpha)]:
        if abs(phi - 90.0) > ang_tol:
            continue
        for m, n in itertools.product(mult, mult):
            for x, y in ((m * e1, n * e2), (n * e2, m * e1)):
                best = min(best, max(abs(x - A_GA) / A_GA, abs(y - C_GA) / C_GA) * 100.0)
    return best


def potency(d_pct):
    if not np.isfinite(d_pct):
        return 1.0
    theta = np.radians(180.0 * min(d_pct / 100.0, 0.15) / 0.15)
    return (2 - 3 * np.cos(theta) + np.cos(theta) ** 3) / 4


def mp_cell(mpr, formula, mpid=None):
    docs = mpr.materials.summary.search(formula=formula,
        fields=["material_id", "structure", "symmetry", "energy_above_hull"])
    if mpid:
        docs = [d for d in docs if str(d.material_id) == mpid]
    d = min(docs, key=lambda z: z.energy_above_hull)
    return SpacegroupAnalyzer(d.structure).get_conventional_standard_structure().lattice


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    screen = []
    for name, (cls, d1d, sign, dens, uc, note) in CUBIC_TET.items():
        a = recover_a(d1d, sign)
        d = two_axis_from_a(a)
        assert abs(abs(a - GA_MEAN) / GA_MEAN * 100 - d1d) < 0.01, f"{name} 1D cross-check failed"
        screen.append(dict(name=name, cls=cls, d=d, f=potency(d), dens=dens, uc=uc,
                           note=note, src="ICSD (experimental)"))
    if not KEY:
        raise SystemExit("Set MP_API_KEY to a Materials Project API key first.")
    with MPRester(KEY) as mpr:
        for name, (cls, mpid, dens, uc, note) in LOWSYM.items():
            d = two_axis_from_cell(mp_cell(mpr, name, mpid))
            screen.append(dict(name=name, cls=cls, d=d, f=potency(d), dens=dens, uc=uc,
                               note=note, src=mpid))
        t1 = {n: (two_axis_from_a(recover_a(d1d, s)), uc) for n, (d1d, s, uc) in TABLE1.items()}
        prod = {n: two_axis_from_cell(mp_cell(mpr, n)) for n in PRODUCTS}

    print("=== SCREEN (Table 2 / Table S2) ===")
    for r in sorted(screen, key=lambda z: z['d']):
        ds = f"{r['d']:6.2f}" if np.isfinite(r['d']) else "  excl"
        print(f"  {r['name']:6} d={ds}  f={r['f']:.3f}  {r['cls']:<7} {r['src']}  {r['note']}")
    print("\n=== Table 1 (Zhang) ===")
    for n, (d, uc) in t1.items():
        print(f"  {n:5} d={d:6.2f}  f={potency(d):.3f}  measured {uc} K")
    print("\n=== Section 2.3 reduction products ===")
    for n, d in prod.items():
        print(f"  {n:6} -> {d:.1f}%" if np.isfinite(d) and d < 1e3 else f"  {n:6} -> no near-square planar match")

    print("\n=== Table S3 correlations ===")
    cn = [(r['d'], r['uc']) for r in screen if r['cls'] in ('carbide', 'nitride') and r['uc'] is not None and np.isfinite(r['d'])]
    ox = [(t1["TeO2"]), (t1["CaO"]), (t1["MgO"])]
    zh = list(t1.values())
    stats = []  # (dataset label as in SI Table S3, n, rho, p_rho, R2, p_R2)
    for label, si_label, data in [
            ("Cubic carbides/nitrides", "Cubic carbides and nitrides (Chakravarty et al. 2021)", cn),
            ("Ceramics+oxides", "Lattice-matched ceramics and oxides (combined)", cn + ox),
            ("All measured", "All measured phases", cn + zh),
            ("Oxides only (Zhang)", "Oxides only (Zhang et al. 2020)", ox)]:
        x = [p[0] for p in data]; y = [p[1] for p in data]
        rho, pr = spearmanr(x, y); r, pp = pearsonr(x, y)
        print(f"  {label:<26} n={len(data):2d}  Spearman {rho:.2f} (p={pr:.4f})  R2 {r**2:.2f} (p={pp:.3f})")
        stats.append((si_label, len(data), rho, pr, r**2, pp))
    xa = np.log10([p[0] for p in cn + zh]); ya = np.array([p[1] for p in cn + zh])
    r, pp = pearsonr(xa, ya)
    print(f"  {'Log-disregistry, all':<26} n={len(cn+zh):2d}  R2 {r**2:.2f} (p={pp:.4f})")
    stats.append(("Log-disregistry fit, all phases", len(cn + zh), None, None, r**2, pp))

    if args.write:
        disp = {"VO2": "VO₂ (rutile)", "TiO2": "TiO₂ (rutile)", "TeO2": "TeO₂",
                "Ta2O5": "Ta₂O₅", "Nb2O5": "Nb₂O₅", "Si3N4": "β-Si₃N₄"}
        path = os.path.join(HERE, "Table_S2_nucleant_screen.csv")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["compound", "class", "structure_source", "disregistry_010_percent",
                        "density_g_cm3", "predicted_potency_f", "measured_undercooling_K", "notes"])
            for r in sorted(screen, key=lambda z: z['d']):
                nm = disp.get(r['name'], r['name'])
                dv = f"{r['d']:.2f}" if np.isfinite(r['d']) else "excluded"
                uc = "" if r['uc'] is None else (f"<{int(r['uc'])}" if 'Upper' in r['note'] else f"{r['uc']:g}")
                w.writerow([nm, r['cls'], r['src'], dv, r['dens'], f"{r['f']:.3f}", uc, r['note']])
        print(f"\nwrote {path}")
        # Table S3 at full precision; the SI rounds these. With n = 3 (oxides only) the
        # p-values carry no weight and the SI prints them as "-".
        path = os.path.join(HERE, "Table_S3_validation_statistics.csv")
        sig = lambda v: "" if v is None else f"{v:.2g}"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["dataset", "n", "spearman_rho", "spearman_p", "pearson_r2", "pearson_p"])
            for si_label, n, rho, pr, r2, pp in stats:
                few = n < 4  # p-values blank, as the SI prints them
                w.writerow([si_label, n, "" if rho is None else f"{rho:.2f}", "" if few else sig(pr),
                            f"{r2:.2f}", "" if few else sig(pp)])
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
