"""
Work of adhesion from the relaxed pw.out files:
    W_ad = (E_slab_MN + E_slab_Ga - E_interface) / A     [J/m^2]
Run after the QE jobs finish. Reads areas.txt written by build_qe_inputs.py.
Higher W_ad = more potent nucleant. Compare the measured anchor (ZrN) against the
predictions (HfN, ScN); a DFT ordering that matches the measured supercooling ranking
validates the predictions.
"""
import os, re

RY_TO_EV = 13.605693
EV_A2_TO_J_M2 = 16.021766          # eV/Angstrom^2 -> J/m^2

def final_energy_ry(pwout):
    """Final '! total energy' from a pw.x run, but only if the run really finished.

    Taking the last '!' line unconditionally is what produced two false harvests:
    pw.x writes one every ionic step, and an SCF that exhausted electron_maxstep
    still leaves a plausible-looking number behind. Require JOB DONE, and require
    that nothing after the energy we are about to use reports a failed SCF.
    Returns (energy_or_None, reason).
    """
    if not os.path.exists(pwout):
        return None, "no pw.out"
    with open(pwout, encoding="utf-8", errors="replace") as fd:
        text = fd.read()
    if "JOB DONE." not in text:
        return None, "no JOB DONE (running, killed, or out of walltime)"
    energy, end = None, -1
    for m in re.finditer(r"^!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", text, re.M):
        energy, end = float(m.group(1)), m.end()
    if energy is None:
        return None, "no total energy line"
    if "convergence NOT achieved" in text[end:]:
        return None, "SCF did not converge"
    return energy, "ok"

def main():
    root = os.path.dirname(os.path.abspath(__file__))
    areas = {}
    ap = os.path.join(root, "areas.txt")
    if os.path.exists(ap):
        for ln in open(ap):
            mn, a = ln.split()
            areas[mn] = float(a)
    nucleants = sorted(d for d in os.listdir(root)
                       if os.path.isdir(os.path.join(root, d))
                       and os.path.exists(os.path.join(root, d, "interface", "pw.in")))
    print(f"{'nucleant':>8} | {'W_ad (J/m^2)':>13} | status")
    print("-" * 64)
    rows = []
    for mn in nucleants:
        terms, bad = {}, []
        for part in ("slab_MN", "slab_Ga", "interface"):
            e, why = final_energy_ry(os.path.join(root, mn, part, "pw.out"))
            terms[part] = e
            if e is None:
                bad.append(f"{part}: {why}")
        A = areas.get(mn)
        if A is None:
            bad.append("no interface area in areas.txt")
        if bad:
            print(f"{mn:>8} | {'--':>13} | {'; '.join(bad)}")
            continue
        w_ry = (terms["slab_MN"] + terms["slab_Ga"] - terms["interface"]) / A
        w = w_ry * RY_TO_EV * EV_A2_TO_J_M2             # J/m^2
        rows.append((mn, w))
        print(f"{mn:>8} | {w:13.3f} | ok")
    if len(rows) > 1:
        rows.sort(key=lambda r: -r[1])
        order = " > ".join(f"{m} ({w:.2f})" for m, w in rows)
        print("\nadhesion ranking (strongest first):", order)
        print("check: does this order track the measured supercooling ranking for the "
              "anchor(s)? If yes, the HfN/ScN predictions are validated.")

if __name__ == "__main__":
    main()
