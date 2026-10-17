"""Generic closure test: can ANY twin-sector equation of state rescue a no-dark-matter universe?"""
import sys
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from likelihoods import SN, chi2_cc, chi2_bao, ZG

OMB = 0.049
sn = SN(sys.argv[1] if len(sys.argv) > 1 else "pantheon")
N_DATA = 15 + 13 + len(sn.z)
REF_LCDM, REF_JANUS = 1422.50, 1739.05  # combined-fit totals (Pantheon+) from combined_fit.py


def make_E(ok, w0, wa):
    ot = 1 - OMB - ok
    def E(z):
        a = 1 / (1 + z)
        f = (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
        return np.sqrt(np.clip(OMB * (1 + z) ** 3 + ok * (1 + z) ** 2 + ot * f, 1e-30, None))
    e2 = OMB * (1 + ZG) ** 3 + ok * (1 + ZG) ** 2 + ot * ((1 + ZG) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * ZG / (1 + ZG)))
    return E, ot, bool(np.all(e2 > 0))


def total(p, sign):
    H0, ok, w0, wa, rd = p
    E, ot, valid = make_E(ok, w0, wa)
    if not valid or (sign > 0 and ot <= 0) or (sign < 0 and ot >= 0):
        return 1e9
    return chi2_cc(E, H0) + chi2_bao(E, ok, H0, rd) + sn.chi2(E, ok)


BOUNDS = [(50, 90), (-1.0, 1.5), (-3.0, 1.0), (-6.0, 6.0), (100, 200)]
lines = [f"Closure test: baryons only (Omega_b = {OMB}) + curvature + twin sector with w(a) = w0 + wa(1-a); r_d free",
         f"data: 15 CC + 13 DESI DR2 BAO + {sn.name}", "",
         f"reference: flat LCDM total chi2 = {REF_LCDM}, Janus = {REF_JANUS}", "",
         f"{'twin density':16s} {'chi2':>9s} {'dchi2 vs LCDM':>14s} {'H0':>6s} {'Om_k':>7s} {'Om_twin':>8s} {'w0':>7s} {'wa':>7s} {'rd':>6s}   effective w(z=0), w(z=1)"]
best = {}
for label, sign in [("negative (Janus)", -1), ("positive", +1), ("free sign", 0)]:
    de = differential_evolution(total, BOUNDS, args=(sign,), seed=1, popsize=20, maxiter=300, tol=1e-8, polish=False)
    res = minimize(total, de.x, args=(sign,), method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-7))
    p = res.x if res.fun < de.fun else de.x
    c = min(res.fun, de.fun)
    H0, ok, w0, wa, rd = p
    ot = 1 - OMB - ok
    best[label] = (c, p)
    lines.append(f"{label:16s} {c:9.2f} {c-REF_LCDM:+14.1f} {H0:6.1f} {ok:7.3f} {ot:8.3f} {w0:7.3f} {wa:7.3f} {rd:6.1f}   {w0:+.2f}, {w0 + wa*0.5:+.2f}")
lines += ["", "Interpretation key: a twin sector with Om_twin ~ +0.95 and w ~ -0.7 is 'dark matter + dark energy in one fluid';",
          "Om_twin < 0 is the Janus-type negative-density sector. k = 5 parameters vs 3 for flat LCDM;",
          "a dchi2 above about +10 means the extra freedom could not buy back the fit."]
text = "\n".join(lines)
print("\n" + text)
open("results/twin_eos_fit.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(8, 5))
zz = np.linspace(0, 2.4, 300)
E_l = lambda z: np.sqrt(0.304 * (1 + z) ** 3 + 0.696)
for (label, (c, p)), col in zip(best.items(), ["C1", "C2", "C3"]):
    E, ot, _ = make_E(p[1], p[2], p[3])
    ax.plot(zz, E(zz) / E_l(zz) * (68.1 / p[0]) ** 0 , color=col, label=f"twin {label}: χ² = {c:.0f}")
ax.axhline(1, color="C0", label="flat LCDM (χ² = 1422)")
Ej = lambda z: (1 + z) * np.sqrt(1 + 0.17 * z)
ax.plot(zz, Ej(zz) / E_l(zz), "C7--", label="Janus best (χ² = 1739)")
ax.set_xlabel("z"); ax.set_ylabel("E(z) / E_LCDM(z)"); ax.legend(fontsize=8)
ax.set_title("Best no-dark-matter backgrounds with a free twin-sector equation of state")
fig.savefig("results/twin_eos_fit.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_twin_eos_fit.png", dpi=120, bbox_inches="tight")
print("\nplot saved to results/twin_eos_fit.png")
