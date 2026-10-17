import numpy as np
from scipy.optimize import minimize
from likelihoods import comoving, transverse, ZG, C_KMS, chi2_cc, chi2_bao

tab = np.genfromtxt("data/Pantheon+SH0ES.dat", names=True, dtype=None, encoding=None)
with open("data/Pantheon+SH0ES_STAT+SYS.cov") as f:
    n = int(f.readline())
    cov = np.loadtxt(f).reshape(n, n)
calib = tab["IS_CALIBRATOR"] == 1
use = calib | (tab["zHD"] > 0.01)
z, zhel, mb, ceph, cal = tab["zHD"][use], tab["zHEL"][use], tab["m_b_corr"][use], tab["CEPH_DIST"][use], calib[use]
icov = np.linalg.inv(cov[np.ix_(use, use)])
print(f"{use.sum()} supernovae, of which {cal.sum()} in Cepheid-calibrated hosts")


def chi2_sn_abs(E, ok, H0, M):
    chi = comoving(E)
    if chi is None:
        return np.inf
    with np.errstate(all="ignore"):
        mu = 5 * np.log10((1 + zhel) * np.interp(z, ZG, transverse(chi, ok)) * C_KMS / H0) + 25
    mu = np.where(cal, ceph, mu)
    r = mb - M - mu
    return float(r @ icov @ r) if np.all(np.isfinite(r)) else np.inf


bg = {"Janus": lambda s: ((lambda x: (1 + x) * np.sqrt(np.clip(1 + 2 * s * x, 1e-15, None))), 1 - 2 * s),
      "Flat LCDM": lambda s: ((lambda x: np.sqrt(s * (1 + x) ** 3 + 1 - s)), 0.0)}
start = {"Janus": 0.0, "Flat LCDM": 0.3}

lines = ["Pantheon+ with SH0ES Cepheid calibration (absolute H0)", ""]
for name in bg:
    def sn_only(p):
        E, ok = bg[name](p[2])
        return chi2_sn_abs(E, ok, p[0], p[1])

    def combined(p):
        E, ok = bg[name](p[2])
        return chi2_sn_abs(E, ok, p[0], p[1]) + chi2_cc(E, p[0]) + chi2_bao(E, ok, p[0], p[3])

    a = minimize(sn_only, [73, -19.25, start[name]], method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-6, fatol=1e-6))
    b = minimize(combined, [70, -19.3, start[name], 148], method="Nelder-Mead", options=dict(maxiter=40000, xatol=1e-6, fatol=1e-6))
    shape = "q0" if name == "Janus" else "Om"
    lines.append(f"{name}:")
    lines.append(f"  calibrated SNe alone:      H0 = {a.x[0]:.2f}, M = {a.x[1]:.3f}, {shape} = {a.x[2]:+.3f}, chi2 = {a.fun:.1f}")
    lines.append(f"  + BAO + CC (all combined):  H0 = {b.x[0]:.2f}, M = {b.x[1]:.3f}, {shape} = {b.x[2]:+.3f}, rd = {b.x[3]:.1f}, chi2 = {b.fun:.1f}")
    lines.append(f"  chi2 cost of combining = {b.fun - a.fun:.1f}")
    lines.append("")
text = "\n".join(lines)
print(text)
open("results/h0_test.txt", "w").write(text + "\n")
