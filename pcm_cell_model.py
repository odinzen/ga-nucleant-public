"""Lumped phase-change cell model behind Figure 4.

A gallium phase-change buffer of mass m is bonded to a chip. The junction sits at the
buffer temperature, which the high conductivity of gallium justifies for a small cell.
Under load the buffer warms, then melts at a fixed temperature and holds the junction
near the melt point; at rest a cold plate pulls the heat back out. On the cooling stroke
the buffer releases its latent heat only if it nucleates above the cold plate, so the
recovery turns on a single inequality,

    dT_residual  <  Tm - Ts

which is the design rule the paper reports. Three cases are run: a near-coherent
nucleant (2 K residual, refreezes every cycle), the best measured oxide (38 K), and bulk
gallium (58 K). The last two cannot nucleate above the plate, so after the first cycle
they carry sensible heat only and the junction runs past its limit.

Reproduces Figure 4 of "Stability-filtered lattice matching for the discovery of
supercooling nucleants in gallium and its low-melting alloys".

    python pcm_cell_model.py          # prints the protection times, writes figure_4.png
"""
import os
import numpy as np
import matplotlib.pyplot as plt

# --- cell parameters (Table S5) -------------------------------------------------
M = 10.0        # buffer mass, g
L = 80.0        # latent heat, J/g
CP = 0.37       # specific heat, J/(g K)
Q_LOAD = 5.0    # chip load, W
TS = 20.0       # cold plate, degC
TM = 29.76      # gallium melt point, degC
T_LIMIT = 60.0  # junction cap, degC
R_COOL = 1.0 / 0.30   # buffer-to-plate resistance = 1/UA, K/W (UA = 0.30 W/K, Table S5)
T_ON, T_OFF = 180.0, 430.0  # load-on / load-off, s (Table S5)
DT = 0.5        # timestep, s
N_CYCLES = 6
LCAP = M * L    # total latent capacity, J


def simulate(residual_undercooling, n_cycles=N_CYCLES):
    """Junction history for a given residual supercooling (K below the melt point).

    A single bounded variable, `stored`, tracks the latent heat held in the solid
    (0 = fully liquid, LCAP = fully solid). The buffer melts at Tm when heated and,
    if it can nucleate above the cold plate, freezes at the nucleation temperature when
    cooled. If the residual supercooling puts the nucleation temperature below the cold
    plate, the buffer never re-solidifies and carries sensible heat only after the
    first melt.
    """
    t_nuc = TM - residual_undercooling
    can_freeze = t_nuc > TS

    t, temp, stored = 0.0, TS, LCAP        # start fully solid at the plate temperature
    times, temps = [t], [temp]
    cycle, load_on, phase_start = 0, True, 0.0

    while cycle < n_cycles:
        phase_end = phase_start + (T_ON if load_on else T_OFF)
        while t < phase_end:
            dQ = (Q_LOAD if load_on else -max(0.0, (temp - TS) / R_COOL)) * DT

            if dQ >= 0:                     # heating
                if stored > 0:              # solid present: warm to Tm, then melt
                    if temp < TM:
                        temp += dQ / (M * CP)
                        if temp > TM:
                            stored -= (temp - TM) * (M * CP)
                            temp = TM
                    else:
                        stored -= dQ
                    if stored < 0:          # fully melted, spill into sensible heat
                        temp = TM + (-stored) / (M * CP)
                        stored = 0.0
                else:                       # fully liquid
                    temp += dQ / (M * CP)
            else:                           # cooling
                if can_freeze and stored < LCAP:   # liquid present: cool to t_nuc, freeze
                    if temp > t_nuc:
                        temp += dQ / (M * CP)
                        if temp < t_nuc:
                            stored += (t_nuc - temp) * (M * CP)
                            temp = t_nuc
                    else:
                        stored += -dQ
                    if stored > LCAP:       # fully frozen, spill into sensible cooling
                        temp = t_nuc - (stored - LCAP) / (M * CP)
                        stored = LCAP
                else:                       # cannot nucleate, or already solid
                    temp += dQ / (M * CP)

            temp = max(TS, temp)
            t += DT
            times.append(t)
            temps.append(temp)

        load_on = not load_on
        if load_on:
            cycle += 1
        phase_start = t

    return np.array(times), np.array(temps)


def main():
    cases = [(2.0, "Near-coherent (HfN/ScN)", "k-"),
             (38.0, "Best oxide (TeO$_2$, 38 K)", "k--"),
             (58.0, "Bulk Ga (58 K)", "k:")]

    fig, ax = plt.subplots(figsize=(3.46, 2.8))
    for residual, label, style in cases:
        times, temps = simulate(residual)
        ax.plot(times / 60, temps, style, linewidth=1.0, label=label)
        refreezes = (TM - residual) > TS
        print(f"{label:28} residual {residual:4.0f} K   peak {temps.max():6.1f} degC   "
              f"refreezes every cycle: {'yes' if refreezes else 'no'}")

    for y in (TM, T_LIMIT, TS):
        ax.axhline(y, color="k", linewidth=0.5, alpha=0.4)
    ax.set_xlabel("Time (min)", fontsize=8)
    ax.set_ylabel("Junction temperature (°C)", fontsize=8)
    ax.set_ylim(15, 80)
    ax.legend(fontsize=6, frameon=True, framealpha=1.0, edgecolor="0.75", loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout(pad=0.4)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure_4.png")
    plt.savefig(out, dpi=600, bbox_inches="tight")
    plt.savefig(out.replace(".png", ".tif"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    print("\nwrote figure_4.png and figure_4.tif")


if __name__ == "__main__":
    main()
