# First-principles nucleant potency (Quantum ESPRESSO)

Quantum ESPRESSO inputs and the analysis script behind the work-of-adhesion numbers in
the paper (Section 2.4). Computed on the Sol supercomputer at Arizona State University.
The multi-GB raw `pw.out`/wavefunction files are not archived here; the inputs plus
`analyze.py` reproduce the values from a fresh run.

Goal: replace the manuscript's *calibrated* potency surrogate with a first-principles
**work of adhesion** at the alpha-Ga(010) // nucleant(001) coherent interface,

    W_ad = ( E_slab_MN + E_slab_Ga - E_interface ) / A            [J/m^2]

A higher W_ad means a more potent nucleant (lower solid-substrate interfacial energy,
smaller contact angle, smaller barrier). Feed W_ad into the same CNT model the paper
already runs; nothing is fit.

## Results (`results_W_ad.csv`)

| nucleant | W_ad (J m^-2) | (010) disregistry | role |
|---|---|---|---|
| ScN | 0.884 | 0.48% | prediction |
| HfN | 0.740 | 0.07% | prediction |
| ZrN | 0.609 | 1.19% | measured anchor (Chakravarty et al.) |

Both untested candidates bind alpha-Ga more strongly than the measured-good ZrN anchor,
which is the out-of-sample support for the screen. Lattice match and adhesion disagree on
the top candidate (HfN matches best, ScN adheres strongest), direct evidence for the
paper's thesis that lattice match alone does not set potency.

**Scope, stated plainly.** These are a proof-of-concept ranking at one coherent interface
registry, with plane-wave cutoffs (60/480 Ry), k-points, slab thickness, and vacuum set
for a consistent ranking across the three nitrides, not for a converged absolute value.
All three interfaces are BFGS-relaxed to 0.005 Ry/bohr. Treat the numbers as a relative
ranking, not a publication-grade absolute work of adhesion. See the convergence notes below.

**Validation logic (why this is not circular).** Compute W_ad for the *measured* anchor
**ZrN** (Chakravarty et al.) alongside the untested predictions **HfN** and **ScN**. If the
DFT ordering reproduces the measured supercooling ranking for the anchor(s), the HfN/ScN
numbers are a genuine out-of-sample prediction. Extend to HfC, TiN, NbC, ... for a full
ranking (add them to `NUCLEANTS`/`A_MN` in `build_qe_inputs.py`).

## Layout
```
build_qe_inputs.py   generator (ASE); writes the pw.in tree, verified geometry
areas.txt            interface area A per nucleant (A^2), used by analyze.py
job_pw.slurm         Sol SLURM template (edit account/partition/modules)
analyze.py           parse pw.out energies -> W_ad table
bulk_Ga/pw.in                         vc-relax, alpha-Ga reference
<MN>/bulk_MN/pw.in                    vc-relax, rocksalt nitride
<MN>/slab_MN/pw.in    <MN>/slab_Ga/   fixed-cell relax, isolated slabs (common in-plane cell)
<MN>/interface/pw.in                  fixed-cell relax, the joined interface (72 atoms)
```

## Run order (on Sol)
1. Put SSSP-efficiency PBE pseudopotentials in `pseudo/`; edit the `PSEUDO` filenames in
   `build_qe_inputs.py` (and `pseudo_dir` if not `./pseudo`). If your QE build rejects the
   empty `&FCP`/`&RISM` namelists ASE writes, delete those two lines from each `pw.in`.
2. `sbatch job_pw.slurm` in each leaf directory (or loop, see below).
3. `python analyze.py` to print the W_ad table.

Copy-paste loop (adjust to your scheduler):
```bash
for d in bulk_Ga */bulk_MN */slab_MN */slab_Ga */interface; do
  (cd "$d" && sbatch --job-name "$(echo $d|tr / _)" ../../job_pw.slurm)   # depth-adjust the path
done
```

## Convergence is YOUR job (the numbers below are starting guesses, not converged)
- **ecutwfc / ecutrho** (60 / 480 Ry): raise until total energy / W_ad is stable to your tol.
- **k-points**: generated as ~30/L; densify the in-plane mesh until W_ad converges.
- **slab thickness** `NLAY_SUB/NLAY_GA` (5/4) and **vacuum** (12 A): thicken until each
  isolated-slab surface energy converges and the two surfaces are decoupled.
- **interface registry**: the generator places one stacking. Scan lateral registries (and
  the **gap**, 2.35 A) for the minimum-energy interface; that minimum is the physical W_ad.
- **strain reference**: the Ga slab is strained to the substrate (coherent). For near-coherent
  HfN/ScN this is negligible; for larger mismatches account for the coherency strain energy.

## Honest caveats
- W_ad here is a **solid** alpha-Ga/substrate adhesion, a proxy for the **liquid**-Ga
  nucleation potency (the Turnbull-Bramfitt premise that crystallographic coherency sets
  potency). It gives a defensible **relative ranking + a physical interfacial energy**, not
  an absolute measured supercooling. The fully rigorous version is AIMD of liquid Ga wetting
  the substrate, which is far heavier and out of scope here.
- alpha-Ga is metallic and orthorhombic (covalent Ga2 dimers); use cold/MV smearing and
  check k-convergence carefully.

## Order-of-magnitude cost
~13 jobs for the 3-system proof of concept; interface cells are 72 atoms (metallic). Plan
a few thousand core-hours once converged; the near-coherent HfN/ScN are the cheapest cells.
