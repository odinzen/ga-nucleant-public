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
R_COOL = 0.73   # buffer-to-plate resistance, K/W
T_ON, T_OFF = 180.0, 60.0   # load-on / load-off, s
DT = 0.5        # timestep, s
N_CYCLES = 6


def simulate(residual_undercooling, n_cycles=N_CYCLES):
    """Junction history for a given residual supercooling (K below the melt point).

    The buffer can only refreeze if its nucleation temperature sits above the cold
    plate; otherwise it stays liquid for the rest of the run and the latent capacity
    is lost after the first melt.
    """
    t_nucleate = TM - residual_undercooling
    can_freeze = t_nucleate > TS

    t, temp, phase, latent = 0.0, TS, "solid", 0.0
    times, temps = [t], [temp]
    cycle, load_on, phase_start = 0, True, 0.0

    while cycle < n_cycles:
        phase_end = phase_start + (T_ON if load_on else T_OFF)
        while t < phase_end:
            q = Q_LOAD if load_on else -max(0.0, (temp - TS) / R_COOL)

            if phase == "solid":
                temp += q * DT / (M * CP)
                if temp >= TM and load_on:
                    temp, phase, latent = TM, "melting", M * L
            elif phase == "melting":
                temp = TM
                latent -= q * DT
                if latent <= 0:
                    latent, phase = 0.0, "liquid"
            elif phase == "liquid":
                temp += q * DT / (M * CP)
                if not load_on and can_freeze and temp <= t_nucleate:
                    temp, phase, latent = t_nucleate, "freezing", M * L
            elif phase == "freezing":
                temp = t_nucleate
                latent += q * DT          # q is negative while cooling
                if latent <= 0:
                    latent, phase = 0.0, "solid"

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
    plt.savefig("figure_4.png", dpi=600, bbox_inches="tight")
    print("\nwrote figure_4.png")


if __name__ == "__main__":
    main()
