import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Sagredo, Nesseris & Sapone 2018 (arXiv:1806.10822), Table I: z, fs8, sigma, Om_ref
D = np.array([
    [0.02, 0.428, 0.0465, 0.3], [0.02, 0.398, 0.065, 0.3], [0.02, 0.314, 0.048, 0.266], [0.10, 0.370, 0.130, 0.3],
    [0.15, 0.490, 0.145, 0.31], [0.17, 0.510, 0.060, 0.3], [0.18, 0.360, 0.090, 0.27], [0.38, 0.440, 0.060, 0.27],
    [0.25, 0.3512, 0.0583, 0.25], [0.37, 0.4602, 0.0378, 0.25], [0.32, 0.384, 0.095, 0.274], [0.59, 0.488, 0.060, 0.307115],
    [0.44, 0.413, 0.080, 0.27], [0.60, 0.390, 0.063, 0.27], [0.73, 0.437, 0.072, 0.27], [0.60, 0.550, 0.120, 0.3],
    [0.86, 0.400, 0.110, 0.3], [1.40, 0.482, 0.116, 0.27], [0.978, 0.379, 0.176, 0.31], [1.23, 0.385, 0.099, 0.31],
    [1.526, 0.342, 0.070, 0.31], [1.944, 0.364, 0.106, 0.31]])
zd, fs8, sig, omref = D.T
C = np.diag(sig**2)
C[12:15, 12:15] = 1e-3 * np.array([[6.400, 2.570, 0.000], [2.570, 3.969, 2.540], [0.000, 2.540, 5.184]])
C[18:22, 18:22] = 1e-2 * np.array([[3.098, 0.892, 0.329, -0.021], [0.892, 0.980, 0.436, 0.076],
                                   [0.329, 0.436, 0.490, 0.350], [-0.021, 0.076, 0.350, 1.124]])
iC = np.linalg.inv(C)
xi_mm = fs8 @ iC @ fs8
OMB = 0.049  # ordinary matter = baryons in Janus (no dark matter)


def hda(E, ok, z):
    """H(z) * d_A(z) in units independent of H0."""
    zg = np.linspace(0, z.max() * 1.001, 4001)
    chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
    dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return E(z) * np.interp(z, zg, dm) / (1 + z)


E_ref = lambda om: (lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om))
HDA_REF = np.array([hda(E_ref(o), 0.0, np.array([z]))[0] for z, o in zip(zd, omref)])


def chi2_marg(nu, E, ok):
    """chi2 with sigma8,0 marginalised analytically (sigma8 > 0 enforced)."""
    mu = nu / (hda(E, ok, zd) / HDA_REF)
    xmn, xnn = fs8 @ iC @ mu, mu @ iC @ mu
    if xmn <= 0:
        return xi_mm, 0.0
    return xi_mm - xmn**2 / xnn, xmn / xnn


# ---- flat LCDM: single fluid ----
def lcdm(om):
    E = E_ref(om)
    dlnE = lambda x: -1.5 * om * np.exp(-3 * x) / (om * np.exp(-3 * x) + 1 - om)
    def rhs(x, y):
        return [y[1], -(2 + dlnE(x)) * y[1] + 1.5 * om * np.exp(-3 * x) / E(np.exp(-x) - 1) ** 2 * y[0]]
    ai = 1 / 51
    s = solve_ivp(rhs, [np.log(ai), 0], [ai, ai], dense_output=True, rtol=1e-8, atol=1e-12)
    d1 = s.sol(0)[0]
    nu = np.array([s.sol(np.log(1 / (1 + z)))[1] for z in zd]) / d1
    return chi2_marg(nu, E, 0.0)


# ---- Janus: two coupled fluids ----
IC = {"in phase (d+ = d-)": (1, 1, 0, 0), "in phase, growing": (1, 1, 1, 1), "anti-phase (d+ = -d-)": (1, -1, 0, 0),
      "twin matter only": (0, 1, 0, 0), "ordinary matter only": (1, 0, 0, 0)}


def janus(q0, ic, return_curve=False):
    om_p, om_m = OMB, OMB - 2 * q0
    if om_m < 0:
        return np.inf, 0.0, None
    E = lambda z: (1 + z) * np.sqrt(np.clip(1 + 2 * q0 * z, 1e-15, None))
    zmax = 1 / (2 * abs(q0)) if q0 < 0 else 1e9
    zi = min(0.9 * zmax, 30.0)

    def rhs(x, y):
        a = np.exp(x); z = 1 / a - 1
        dlnE = -1 - q0 * (1 + z) / (1 + 2 * q0 * z)
        src = 1.5 / (E(z) ** 2 * a**3)
        dp, vp, dm, vm = y
        return [vp, -(2 + dlnE) * vp + src * (om_p * dp - om_m * dm),
                vm, -(2 + dlnE) * vm + src * (om_m * dm - om_p * dp)]

    dp0, dm0, vp0, vm0 = ic
    s = solve_ivp(rhs, [np.log(1 / (1 + zi)), 0], [dp0, vp0, dm0, vm0], dense_output=True, rtol=1e-8, atol=1e-12)
    d1 = s.sol(0)[0]
    if abs(d1) < 1e-12:
        return np.inf, 0.0, None
    nu = np.array([s.sol(np.log(1 / (1 + z)))[1] for z in zd]) / d1
    c2, s8 = chi2_marg(nu, E, 1 - 2 * q0)
    curve = None
    if return_curve:
        zz = np.linspace(0, min(2.0, 0.89 * zmax), 200)
        curve = (zz, s8 * np.array([s.sol(np.log(1 / (1 + z)))[1] for z in zz]) / d1)
    return c2, s8, curve


lines = [f"Growth of structure: {len(zd)} fs8 measurements (Gold compilation, Sagredo et al. 2018), sigma8 marginalised", ""]
rl = minimize_scalar(lambda o: lcdm(o)[0], bounds=(0.05, 0.9), method="bounded")
c_l, s8_l = lcdm(rl.x)
lines.append(f"Flat LCDM (Om free): chi2 = {c_l:.2f} for {len(zd)-2} dof, Om = {rl.x:.3f}, sigma8 = {s8_l:.3f}")
lines += ["", f"Janus, ordinary matter = baryons (Om+ = {OMB}), twin matter = Om+ - 2 q0:"]
best = (np.inf, None, None)
for name, ic in IC.items():
    c, s8, _ = janus(-0.087, ic)
    lines.append(f"  q0 = -0.087, start '{name:24s}': chi2 = {c:8.2f}, sigma8 = {s8:.3f}")
    qs = np.linspace(-0.2, 0.024, 57)
    cs = [janus(q, ic)[0] for q in qs]
    i = int(np.argmin(cs))
    if cs[i] < best[0]:
        best = (cs[i], qs[i], name)
lines += ["", f"Best Janus over all q0 < 0.0245 and starting conditions: chi2 = {best[0]:.2f} (q0 = {best[1]:.3f}, start '{best[2]}')",
          f"chi2(best Janus) - chi2(flat LCDM) = {best[0] - c_l:.2f}"]
text = "\n".join(lines)
print(text)
open("results/growth_fit.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(8, 5))
ax.errorbar(zd, fs8, sig, fmt="o", color="k", label="fσ8 data (Gold 2018)")
zz = np.linspace(0, 2, 200)
om = rl.x
for name, col in [("in phase (d+ = d-)", "C1"), ("anti-phase (d+ = -d-)", "C3"), (best[2], "C2")]:
    q = -0.087 if col != "C2" else best[1]
    _, _, curve = janus(q, IC[name], return_curve=True)
    if curve is not None:
        ax.plot(*curve, color=col, label=f"Janus q0 = {q:.3f}, {name}")
# LCDM curve
E = E_ref(om)
dlnE = lambda x: -1.5 * om * np.exp(-3 * x) / (om * np.exp(-3 * x) + 1 - om)
s = solve_ivp(lambda x, y: [y[1], -(2 + dlnE(x)) * y[1] + 1.5 * om * np.exp(-3 * x) / E(np.exp(-x) - 1) ** 2 * y[0]],
              [np.log(1 / 51), 0], [1 / 51, 1 / 51], dense_output=True, rtol=1e-8)
ax.plot(zz, s8_l * np.array([s.sol(np.log(1 / (1 + z)))[1] for z in zz]) / s.sol(0)[0], color="C0", label=f"flat LCDM, Om = {om:.2f}")
ax.set_xlabel("z"); ax.set_ylabel("fσ8"); ax.legend(fontsize=8)
fig.savefig("results/growth_fit.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_growth_fit.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/growth_fit.txt and results/growth_fit.png (copy on your desktop)")
