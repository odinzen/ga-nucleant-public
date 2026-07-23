"""
Fig 4: Junction temperature of a Ga PCM cell over 6 load cycles.

Lumped thermal model (from manuscript Section 3.7):
  - 10 g Ga buffer, latent heat 80 J/g -> L = 800 J total
  - 5 W continuous load
  - Cold plate at Ts = 20 C
  - Ga melt point: Tm = 29.76 C
  - Junction limit: 60 C

Three cases:
  1. Near-coherent nucleant (HfN/ScN): dT_residual ~2 K -> refreezes at ~27.8 C, inside
     the 9.8 K gap (Tm - Ts = 9.76 K). Buffer available every cycle.
  2. TeO2 (best oxide): manuscript says 38 K supercooling, but Zhang measured ~5 K.
     The manuscript Section 3.7 says "TeO2 at 38 K" does not clear the 9.8 K gap.
     Wait - re-read: "best oxide TeO2 at 38 K and bulk gallium at 58 K do not".
     So the manuscript assigns 38 K to TeO2 for the device model context.
     This is the supercooling FROM the melt point, so it freezes at 29.76 - 38 = -8.2 C,
     which is below the cold plate (20 C), so it NEVER refreezes in the device window.
     Actually: if cold plate is 20 C and TeO2 needs 38 K undercooling, it would try to
     freeze at 29.76 - 38 = -8.24 C, but the cold plate only goes to 20 C. So it cannot
     nucleate and the buffer stays liquid = no latent heat recovery.
  3. Bulk Ga: ~58 K supercooling (manuscript: "TeO2 at 38 K and bulk gallium at 58 K").
     Also cannot freeze above cold plate.

Timing:
  Manuscript: "protection time of 180 s" for near-coherent case; "30 s" for others.
  10 g * 80 J/g = 800 J; at 5 W that's 160 s of latent heat phase, plus sensible heat.
  Sensible heat from 20 C to Tm (29.76 C): c_p(Ga) ~ 0.37 J/(g K) * 10 g * 9.76 K = 36 J
  -> 7 s to reach melt. Then latent: 800 J / 5 W = 160 s at melt point.
  After melt absorbed, temperature rises at 5 W / (10 g * 0.37 J/gK) = 1.35 C/s.
  From 29.76 to 60 C = 30.2 C / 1.35 C/s = 22 s. Total cycle ~ 189 s ~ 180 s. Good.

  For the near-coherent case (refreezes each cycle):
  - Cooling (rest) phase: with cold plate at 20 C and thermal resistance:
    R = (Tm - Ts) / Q_cool. We need the buffer to refreeze and release 800 J + sensible.
    Assume symmetric cycle: on-time = 180 s, off-time = ~60 s (matches 6 cycles in plot).
    During cooling: buffer freezes at Tm - 2 K = ~27.76 C, releases latent heat to cold plate.
    Cooling power = (Tm - Ts) / R. If cycle time ~240 s total, off-time = 60 s.
    Q out = 800 J in 60 s -> 13.3 W cooling during off time.
    R = (Tm - Ts) / Q_cool = 9.76 K / 13.3 W = 0.73 K/W. Reasonable for a thermal spreader.

  For TeO2 / bulk case:
  - First cycle: sensible heat to Tm (7 s) + latent (160 s) + sensible to 60 C (22 s) = 189 s
  - Buffer stays liquid after (no nucleation). On-time 2+: only sensible heat from 20 C to 60 C.
    Q_sensible_total = 10 g * 0.37 J/(gK) * 40 K = 148 J -> 148/5 = 30 s. Matches "30 s" in ms.

Simulation approach: explicit time-stepping with state machine.
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
m   = 10.0      # g
L   = 80.0      # J/g  -> total latent = 800 J
cp  = 0.37      # J/(g K)
Q_load = 5.0    # W
Ts  = 20.0      # C (cold plate)
Tm  = 29.76     # C (Ga melting point)
T_lim = 60.0    # C (junction limit)

# Thermal resistance (cold plate to buffer) - tuned to give ~180 s protection
R_cool = 0.73   # K/W  -> Q_cool_max = (Tm - Ts) / R = 13.4 W

# Cycle structure: alternate load (on) and rest (off)
t_on  = 180.0   # s per on-phase (approximately - simulation will find actual)
t_off = 60.0    # s per off-phase

n_cycles = 6
dt = 0.5        # time step in s

def simulate(dT_nucleation, label, n_cycles=6):
    """
    dT_nucleation: undercooling required to nucleate (K below Tm).
    If Tm - dT_nucleation < Ts, the buffer can never nucleate (stays liquid).
    """
    T_nucleate = Tm - dT_nucleation  # freezing onset temp
    can_freeze = T_nucleate > Ts

    T_j = Ts        # start at cold plate temp
    t = 0.0
    latent_stored = m * L   # J available in solid state initially (solid buffer)
    phase = 'solid'  # 'solid', 'melting', 'liquid', 'freezing'
    # latent_remaining tracks how much latent heat is still in the buffer
    # positive = energy that must be added to melt; after full melt = 0
    latent_remaining = 0.0  # buffer starts solid, latent_remaining = 0 means solid

    times = [t]
    temps = [T_j]

    cycle = 0
    in_on_phase = True
    phase_start_t = 0.0

    # We'll integrate with on/off alternation
    while cycle < n_cycles:
        # Determine current phase duration
        t_phase_end = phase_start_t + (t_on if in_on_phase else t_off)

        while t < t_phase_end:
            if in_on_phase:
                Q_net = Q_load
            else:
                # cooling: Q_out = (T_j - Ts) / R_cool, limited by available driving force
                Q_out = max(0, (T_j - Ts) / R_cool)
                Q_net = -Q_out

            # State machine
            if phase == 'solid':
                # sensible heat of solid buffer
                dT = Q_net * dt / (m * cp)
                T_j += dT
                if T_j >= Tm and in_on_phase:
                    T_j = Tm
                    phase = 'melting'
                    latent_remaining = m * L

            elif phase == 'melting':
                # clamp at Tm, consume latent
                T_j = Tm
                latent_remaining -= Q_net * dt
                if latent_remaining <= 0:
                    latent_remaining = 0
                    phase = 'liquid'

            elif phase == 'liquid':
                dT = Q_net * dt / (m * cp)
                T_j += dT
                # Check nucleation on cooling stroke
                if not in_on_phase and can_freeze and T_j <= T_nucleate:
                    T_j = T_nucleate
                    phase = 'freezing'
                    latent_remaining = m * L

            elif phase == 'freezing':
                # clamp near nucleation temp, release latent
                T_j = T_nucleate
                latent_remaining += Q_net * dt  # Q_net is negative (cooling)
                if latent_remaining <= 0:
                    latent_remaining = 0
                    phase = 'solid'
                    T_j = T_nucleate

            # Clamp junction temperature
            T_j = max(Ts, T_j)

            t += dt
            times.append(t)
            temps.append(T_j)

        # Switch on/off
        if in_on_phase:
            in_on_phase = False
        else:
            in_on_phase = True
            cycle += 1
        phase_start_t = t

    return np.array(times), np.array(temps)

t1, T1 = simulate(dT_nucleation=2.0,  label='Near-coherent (HfN/ScN)')
t2, T2 = simulate(dT_nucleation=38.0, label='TeO2 (38 K residual)')
t3, T3 = simulate(dT_nucleation=58.0, label='Bulk Ga (58 K residual)')

# Normalize time to seconds
fig, ax = plt.subplots(figsize=(3.46, 2.8))

ax.plot(t1 / 60, T1, 'k-',   linewidth=1.0, label='Near-coherent (HfN/ScN)')
ax.plot(t2 / 60, T2, 'k--',  linewidth=1.0, label='Best oxide (TeO$_2$, 38 K)')
ax.plot(t3 / 60, T3, 'k:',   linewidth=1.2, label='Bulk Ga (58 K)')

# Reference lines
ax.axhline(y=Tm,    color='k', linestyle='-',  linewidth=0.5, alpha=0.4)
ax.axhline(y=T_lim, color='k', linestyle='-',  linewidth=0.5, alpha=0.4)
ax.axhline(y=Ts,    color='k', linestyle='-',  linewidth=0.5, alpha=0.4)

ax.set_xlabel("Time (min)", fontsize=8)
ax.set_ylabel("Junction temperature (°C)", fontsize=8)
xmax = max(t1[-1], t2[-1], t3[-1]) / 60
ax.set_xlim(0, xmax)
# The failing cases run away to ~260 degC in this un-throttled lumped model. No axis
# shows that without crushing the 20-60 degC operating band, so the band is resolved and
# the failing curves exit the top as runaway spikes, labelled below.
ax.set_ylim(15, 80)
ax.tick_params(labelsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Reference-line labels on the left, each on a white background so the failing-curve
# spikes that cross the left of the plot cannot run through the text.
_lbl = dict(fontsize=6, ha='left', color='0.4', zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax.text(0.2, Tm + 0.6, f"$T_m$ = {Tm} °C", va='bottom', **_lbl)
ax.text(0.2, T_lim + 0.5, f"Junction limit ({T_lim:.0f} °C)", va='bottom', **_lbl)
ax.text(0.2, Ts - 0.6, f"Cold plate ({Ts:.0f} °C)", va='top', **_lbl)

# Mark where the failing curves actually go, so the top clip reads as runaway not cut-off.
# Text parked in the clear first-cycle gap (x<3, above the near-coherent 47 degC peak).
peak = max(T2.max(), T3.max())
ax.annotate(f"peak ≈ {peak:.0f} °C",
            xy=(t2[np.argmax(T2 > 79)] / 60, 79.3), xytext=(0.4, 71),
            fontsize=6, color='0.35', ha='left', va='center',
            arrowprops=dict(arrowstyle='->', color='0.5', lw=0.7))

# Legend inside, upper-right, white background over the sparse spike region
ax.legend(fontsize=6, frameon=True, framealpha=1.0, edgecolor='0.75',
          loc='upper right', ncol=1, handlelength=1.8, handletextpad=0.5,
          borderpad=0.5)

plt.tight_layout(pad=0.4)
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig4_thermal.png")
plt.savefig(out, dpi=600, bbox_inches='tight')
plt.savefig(out.replace('.png', '.tif'), dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'})
print(f"Saved {out}")
