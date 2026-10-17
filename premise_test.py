import sys
import numpy as np
from scipy.optimize import minimize
from likelihoods import SN, chi2_cc, chi2_bao

sn = SN(sys.argv[1] if len(sys.argv) > 1 else "pantheon")


def total(H0, q0, rd):
    E = lambda z: (1 + z) * np.sqrt(np.clip(1 + 2 * q0 * z, 1e-15, None))
    ok = 1 - 2 * q0
    return chi2_cc(E, H0) + chi2_bao(E, ok, H0, rd) + sn.chi2(E, ok)


def profile(q0):
    res = minimize(lambda p: total(p[0], q0, p[1]), [60, 154], method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-6, maxiter=5000))
    return res.fun, res.x


best = minimize(lambda p: total(*p), [59, 0.08, 154], method="Nelder-Mead",
                options=dict(xatol=1e-7, fatol=1e-7, maxiter=10000))
lines = [f"Janus model, CC + DESI DR2 BAO + {sn.name}", "",
         f"best fit: H0 = {best.x[0]:.2f}, q0 = {best.x[1]:+.4f}, rd = {best.x[2]:.1f}, chi2 = {best.fun:.2f}", ""]
for label, q in [("q0 = 0 (coasting, zero net mass)", 0.0),
                 ("q0 = -0.087 (Petit & d'Agostini 2018)", -0.087)]:
    c, p = profile(q)
    d = c - best.fun
    lines.append(f"{label:40s} chi2 = {c:9.2f}  delta chi2 = {d:7.1f}  (~{np.sqrt(d):.1f} sigma)  H0 = {p[0]:.1f}, rd = {p[1]:.1f}")
lines += ["", "Within the Janus background, Omega_m = 2 q0: q0 > 0 means positive net mass,",
          "the opposite of the model's premise that negative net mass drives acceleration."]
text = "\n".join(lines)
print(text)
open(f"results/premise_test_{sys.argv[1] if len(sys.argv) > 1 else 'pantheon'}.txt", "w").write(text + "\n")
