# Supercooling nucleants in gallium: data and code

Data and reproducibility code for:

> **Stability-filtered lattice matching for the discovery of supercooling nucleants in
> gallium and its low-melting alloys.** M. E. Bustamante, G. Bustamante, K. Lilova.
> *Computational Materials Science* (submitted).

Gallium supercools by tens of kelvin before it freezes, which makes it unreliable as a
phase-change medium. The usual way to pick a heterogeneous nucleant is to match its
lattice to the freezing solid. The paper shows that a lattice match alone is not enough
in gallium, and that a nucleant has to satisfy four requirements together: a close match
to the alpha-Ga (010) plane, stability against reduction by the melt, wetting, and a
density near that of liquid gallium so it stays dispersed. No phase in the screen
satisfies all four.

## Contents

| File | What it is |
|---|---|
| `Table_S2_nucleant_screen.csv` | The full candidate screen, seventeen compounds: class, structure source, (010) disregistry, density, predicted potency factor, and measured undercooling where one exists |
| `Table_S3_validation_statistics.csv` | Rank and linear correlations of measured undercooling against disregistry, per dataset |
| `pcm_cell_model.py` | The lumped phase-change cell model behind Figure 4. Run it to reproduce the junction-temperature history for the three cases |
| `figures/` | The scripts that generate Figures 1-4, Figure S1, and the graphical abstract |
| `dft/` | Quantum ESPRESSO inputs and analysis for the first-principles work of adhesion (Section 2.4), computed on the Sol supercomputer at Arizona State University. Inputs and `analyze.py` reproduce the W_ad values in `dft/results_W_ad.csv`; raw outputs are not archived |

## Reproducing

```
pip install numpy matplotlib
python pcm_cell_model.py            # Figure 4
python figures/fig1_disregistry.py  # and the rest
```

## Notes on the data

Disregistry is the planar match to the alpha-Ga (010) face; single-axis coincidences are
not admitted. Densities are experimental crystallographic values, set against liquid
gallium at 6.05 g cm^-3; TeO2 is paratellurite. Measured undercoolings marked `<` are
upper bounds. Sources for each structure are given in the `structure_source` column.

The first-principles interface calculations were run on the Sol supercomputer at Arizona
State University; their Quantum ESPRESSO inputs and analysis are archived here under `dft/`.
The CALPHAD database and the nucleation engine used in the paper are proprietary, are not
included here, and are available from the corresponding author on reasonable request.

## Licence

The code (`pcm_cell_model.py` and everything under `figures/`) is under the MIT License,
see `LICENSE`. The two CSV data tables are under CC BY 4.0, see `LICENSE-DATA`. Cite the
paper if you use either.

## Contact

Michael E. Bustamante, Odinzen LLC, michaelbusta@odinzen.io
