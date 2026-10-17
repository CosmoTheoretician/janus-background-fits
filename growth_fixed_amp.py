import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S8_PRIOR = (0.78, 0.05)  # weak-lensing-style prior on sigma8 today

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
OMB = 0.049


def hda(E, ok, z):
    zg = np.linspace(0, z.max() * 1.001, 4001)
    chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
    dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return E(z) * np.interp(z, zg, dm) / (1 + z)


E_ref = lambda om: (lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om))
HDA_REF = np.array([hda(E_ref(o), 0.0, np.array([z]))[0] for z, o in zip(zd, omref)])


def chi2_with_prior(nu, E, ok):
    """chi2 minimised over sigma8 with a Gaussian prior (closed form: quadratic in sigma8)."""
    mu = nu / (hda(E, ok, zd) / HDA_REF)
    m0, s0 = S8_PRIOR
    A = mu @ iC @ mu + 1 / s0**2
    B = fs8 @ iC @ mu + m0 / s0**2
    s8 = B / A
    r = fs8 - s8 * mu
    return float(r @ iC @ r + ((s8 - m0) / s0) ** 2), s8


def lcdm(om):
    E = E_ref(om)
    dlnE = lambda x: -1.5 * om * np.exp(-3 * x) / (om * np.exp(-3 * x) + 1 - om)
    rhs = lambda x, y: [y[1], -(2 + dlnE(x)) * y[1] + 1.5 * om * np.exp(-3 * x) / E(np.exp(-x) - 1) ** 2 * y[0]]
    s = solve_ivp(rhs, [np.log(1 / 51), 0], [1 / 51, 1 / 51], dense_output=True, rtol=1e-8, atol=1e-12)
    d1 = s.sol(0)[0]
    nu = lambda z: s.sol(np.log(1 / (1 + z)))[1] / d1
    c2, s8 = chi2_with_prior(np.array([nu(z) for z in zd]), E, 0.0)
    return c2, s8, nu


IC = {"in phase (d+ = d-)": (1, 1, 0, 0), "in phase, growing": (1, 1, 1, 1), "anti-phase (d+ = -d-)": (1, -1, 0, 0),
      "twin matter only": (0, 1, 0, 0), "ordinary matter only": (1, 0, 0, 0)}


def janus(q0, ic):
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
    nu = lambda z: s.sol(np.log(1 / (1 + z)))[1] / d1
    c2, s8 = chi2_with_prior(np.array([nu(z) for z in zd]), E, 1 - 2 * q0)
    return c2, s8, nu


lines = [f"Growth of structure, 22 fs8 points, with prior sigma8 = {S8_PRIOR[0]} +/- {S8_PRIOR[1]}", ""]
rl = minimize_scalar(lambda o: lcdm(o)[0], bounds=(0.05, 0.9), method="bounded")
c_l, s8_l, nu_l = lcdm(rl.x)
lines.append(f"Flat LCDM (Om free): chi2 = {c_l:.2f}, Om = {rl.x:.3f}, sigma8 = {s8_l:.3f}")
lines += ["", "Janus (ordinary matter = baryons, twin = Om+ - 2 q0):",
          f"  {'start':26s} {'chi2 @ q0=-0.087':>17s} {'s8':>6s}   {'best chi2':>9s} {'q0':>7s} {'s8':>6s}"]
qs = np.linspace(-0.2, 0.024, 57)
curves = []
overall = (np.inf, None, None, None)
for name, ic in IC.items():
    c_p, s_p, _ = janus(-0.087, ic)
    cs = [janus(q, ic)[0] for q in qs]
    i = int(np.argmin(cs))
    c_b, s_b, nu_b = janus(qs[i], ic)
    lines.append(f"  {name:26s} {c_p:17.2f} {s_p:6.3f}   {c_b:9.2f} {qs[i]:7.3f} {s_b:6.3f}")
    curves.append((name, qs[i], s_b, nu_b))
    if c_b < overall[0]:
        overall = (c_b, qs[i], name, s_b)
lines += ["", f"Best Janus: chi2 = {overall[0]:.2f} (q0 = {overall[1]:.3f}, start '{overall[2]}', sigma8 = {overall[3]:.3f})",
          f"chi2(best Janus) - chi2(flat LCDM) = {overall[0] - c_l:.2f}",
          f"chi2(Janus q0=-0.087, best start) - chi2(flat LCDM) = {min(janus(-0.087, ic)[0] for ic in IC.values()) - c_l:.2f}"]
text = "\n".join(lines)
print(text)
open("results/growth_fixed_amp.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(8, 5))
ax.errorbar(zd, fs8, sig, fmt="o", color="k", label="fσ8 data (Gold 2018)")
zz = np.linspace(0, 2, 200)
ax.plot(zz, [s8_l * nu_l(z) for z in zz], color="C0", label=f"flat LCDM, Om = {rl.x:.2f}, σ8 = {s8_l:.2f}")
for (name, q, s8, nu), col in zip(curves, ["C1", "C2", "C3", "C4", "C5"]):
    zm = 1 / (2 * abs(q)) if q < 0 else 2.5
    zq = zz[zz < 0.89 * zm]
    ax.plot(zq, [s8 * nu(z) for z in zq], color=col, lw=1, label=f"Janus q0 = {q:.3f}, {name}, σ8 = {s8:.2f}")
ax.set_xlabel("z"); ax.set_ylabel("fσ8"); ax.set_title(f"σ8 prior {S8_PRIOR[0]} ± {S8_PRIOR[1]}"); ax.legend(fontsize=7)
fig.savefig("results/growth_fixed_amp.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_growth_fixed_amp.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/growth_fixed_amp.txt and results/growth_fixed_amp.png")

prof = np.array([min(janus(q, ic)[0] for ic in IC.values()) for q in qs])
ok = qs[prof <= prof.min() + 1]
open("results/growth_q0_range.txt", "w").write(f"{ok.min():.4f} {ok.max():.4f}\n")
print(f"growth (sigma8 prior) q0 68% window: {ok.min():.3f} to {ok.max():.3f}")
