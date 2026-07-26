"""
Build Quantum ESPRESSO (pw.x) inputs for the first-principles nucleant-potency
study: the alpha-Ga(010) // nucleant(001) coherent interface, its two isolated
slabs, and the two bulk references. From these you get the work of adhesion
W_ad = (E_slab_Ga + E_slab_MN - E_interface) / A, which replaces the manuscript's
calibrated potency surrogate.

Runs anywhere ASE is installed (build the inputs locally, scp the tree to Sol).
It does NOT run QE. Verify the geometry (dft/_check/*.png) before submitting.

Convergence is YOUR job: ecutwfc/ecutrho, k-points, slab thickness, and vacuum
below are sensible STARTING values, not converged ones. Test them.
"""
import os
import numpy as np
from ase.spacegroup import crystal
from ase.build import bulk, surface
from ase.io.espresso import write_espresso_in

# --- lattice parameters (experimental; QE vc-relax refines them) ---
A_MN = {"HfN": 4.525, "ScN": 4.501, "ZrN": 4.577}   # rocksalt a (Angstrom)
NUCLEANTS = ["HfN", "ScN", "ZrN"]                    # ZrN is the measured anchor

# --- slab / interface geometry (starting values -- converge) ---
NLAY_SUB, NLAY_GA, VAC, GAP = 5, 4, 12.0, 2.35

# --- SSSP-efficiency PBE pseudopotentials: EDIT to match your SSSP install ---
PSEUDO = {
    "Ga": "Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "N":  "N.pbe-n-radius_5.UPF",
    "Hf": "Hf-sp.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Sc": "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Zr": "Zr-sp.pbe-spn-kjpaw_psl.1.0.0.UPF",
}
ECUTWFC, ECUTRHO = 60.0, 480.0    # Ry -- CONVERGE THESE

def ga_bulk():
    return crystal("Ga", [(0.0, 0.1549, 0.0810)], spacegroup=64,
                   cellpar=[4.5197, 7.6633, 4.5257, 90, 90, 90])

def slabs(mn, a):
    sub = surface(bulk(mn, "rocksalt", a=a, cubic=True), (0, 0, 1), NLAY_SUB, vacuum=None)
    sub = sub[sub.positions[:, 2].argsort()]
    ga = surface(ga_bulk(), (0, 1, 0), NLAY_GA, vacuum=None)
    ga = ga[ga.positions[:, 2].argsort()]
    axy = sub.cell.lengths()[:2]
    sc = ga.cell.copy()
    sc[0] *= axy[0] / ga.cell.lengths()[0]
    sc[1] *= axy[1] / ga.cell.lengths()[1]
    ga.set_cell(sc, scale_atoms=True)          # coherent: strain Ga to the substrate
    inter = sub.copy()
    g2 = ga.copy()
    g2.positions[:, 2] += sub.positions[:, 2].max() + GAP - ga.positions[:, 2].min()
    inter += g2
    for s in (sub, ga, inter):
        c = s.cell.copy()
        c[2, 2] = s.positions[:, 2].max() + VAC
        s.set_cell(c); s.center(axis=2); s.pbc = True
    return sub, ga, inter, float(axy[0] * axy[1])

def base_input(calc, atoms):
    return {
        "control": {"calculation": calc, "restart_mode": "from_scratch",
                    "tprnfor": True, "tstress": True, "pseudo_dir": "./pseudo",
                    "outdir": "./out", "prefix": "pw"},
        "system": {"ibrav": 0, "nat": len(atoms),
                   "ntyp": len(set(atoms.get_chemical_symbols())),
                   "ecutwfc": ECUTWFC, "ecutrho": ECUTRHO,
                   "occupations": "smearing", "smearing": "cold", "degauss": 0.01},
        "electrons": {"conv_thr": 1e-7, "mixing_beta": 0.3},
        "ions": {"ion_dynamics": "bfgs"},
        "cell": {"cell_dofree": "all"},
    }

def kmesh(atoms, planar):
    L = atoms.cell.lengths()
    k = [max(1, int(round(30.0 / L[i]))) for i in range(3)]
    if planar:
        k[2] = 1
    return k

def write_job(path, atoms, calc, planar):
    os.makedirs(path, exist_ok=True)
    data = base_input(calc, atoms)
    if calc != "vc-relax":
        data.pop("cell")
    with open(os.path.join(path, "pw.in"), "w", encoding="utf-8") as fd:
        write_espresso_in(fd, atoms, input_data=data, pseudopotentials=PSEUDO,
                          kpts=kmesh(atoms, planar))

def main():
    root = os.path.dirname(os.path.abspath(__file__))
    # shared bulk Ga reference
    write_job(os.path.join(root, "bulk_Ga"), ga_bulk(), "vc-relax", planar=False)
    areas = {}
    for mn in NUCLEANTS:
        a = A_MN[mn]
        write_job(os.path.join(root, mn, "bulk_MN"),
                  bulk(mn, "rocksalt", a=a, cubic=True), "vc-relax", planar=False)
        sub, ga, inter, area = slabs(mn, a)
        write_job(os.path.join(root, mn, "slab_MN"), sub, "relax", planar=True)
        write_job(os.path.join(root, mn, "slab_Ga"), ga, "relax", planar=True)
        write_job(os.path.join(root, mn, "interface"), inter, "relax", planar=True)
        areas[mn] = area
        print(f"{mn}: interface {len(inter)} atoms, area {area:.3f} A^2 "
              f"-> {mn}/{{bulk_MN,slab_MN,slab_Ga,interface}}/pw.in")
    with open(os.path.join(root, "areas.txt"), "w") as fd:
        for mn, ar in areas.items():
            fd.write(f"{mn} {ar:.6f}\n")
    print("wrote bulk_Ga/ and per-nucleant trees; areas -> areas.txt")

if __name__ == "__main__":
    main()
