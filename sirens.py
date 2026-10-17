import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = 299792.458
GW_H0 = (68.0, 8.0, 6.0)  # LVK GWTC-3, GW170817 + dark sirens: value, +err, -err
GW170817 = (70.0, 12.0, 8.0)

E_janus = lambda q: (lambda z: (1 + z) * np.sqrt(np.clip(1 + 2 * q * z, 1e-15, None)))
E_lcdm = lambda om: (lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om))
MODELS = [  # name, H0 from the fits, E(z), Omega_k, z limit
    ("Flat LCDM (combined fit)", 68.1, E_lcdm(0.304), 0.0, None),
    ("Janus q0 = -0.087 (published), H0 from fit", 63.7, E_janus(-0.087), 1.174, 5.75),
    ("Janus q0 = +0.085 (combined best)", 59.0, E_janus(0.085), 0.83, None),
]


def dl(H0, E, ok, z):
    zg = np.linspace(0, max(z) * 1.001, 20001)
    chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
    dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return (1 + z) * np.interp(z, zg, dm) * C / H0


lines = ["Standard sirens", "", "1. Siren H0 (no supernova calibration) versus the H0 each model needs to fit CC + BAO + SN:", ""]
for name, H0, *_ in MODELS:
    for lab, (h, p, m) in [("GWTC-3", GW_H0), ("GW170817 alone", GW170817)]:
        err = m if H0 < h else p
        lines.append(f"  {name:44s} H0 = {H0:5.1f}   vs {lab:15s} {h:.0f} +{p:.0f}/-{m:.0f}:  {(H0 - h) / err:+.1f} sigma")
lines += ["", "2. Luminosity distance D_L(z) [Mpc] and ratio to flat LCDM (future sirens with redshifts measure this directly):", "",
          f"  {'z':>5s}" + "".join(f" {n.split('(')[0].strip()[:22]:>24s}" for n, *_ in MODELS)]
ZS = np.array([0.01, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0])
D = [dl(H0, E, ok, ZS) for _, H0, E, ok, _ in MODELS]
for i, z in enumerate(ZS):
    row = f"  {z:5.2f}"
    for j, (name, H0, E, ok, zl) in enumerate(MODELS):
        if zl is not None and z >= zl:
            row += f" {'n/a':>24s}"
        else:
            row += f" {D[j][i]:9.0f} ({D[j][i] / D[0][i]:.3f})" + " " * 7
    lines.append(row)
lines += ["", "Same-H0 comparison (shape only, H0 = 68 for all): ratio of D_L to flat LCDM", ""]
for name, _, E, ok, zl in MODELS[1:]:
    d = dl(68.0, E, ok, ZS)
    lines.append(f"  {name:44s} " + " ".join(f"z={z:.1f}: {'n/a' if zl and z >= zl else f'{d[i] / dl(68.0, MODELS[0][2], 0.0, ZS)[i]:.3f}'}" for i, z in enumerate(ZS) if z >= 0.3))
lines += ["", "Einstein Telescope / Cosmic Explorer forecasts: binary neutron stars with counterparts to z ~ 2-3,",
          "D_L to a few percent per event, percent-level H(z) after ~100 events. A 5-10% shape difference",
          "at z ~ 1-2 is therefore a decisive test independent of all supernova and BAO systematics."]
text = "\n".join(lines)
print(text)
open("results/sirens.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(8, 5))
zz = np.linspace(0.01, 3, 300)
ref = dl(68.0, MODELS[0][2], 0.0, zz)
for name, _, E, ok, zl in MODELS[1:]:
    zp = zz if zl is None else zz[zz < zl]
    ax.plot(zp, dl(68.0, E, ok, zp) / ref[: len(zp)], label=name.split(", H0")[0])
ax.axhline(1, color="C0", label="flat LCDM")
ax.axvspan(0, 0.05, color="gold", alpha=0.3, label="current sirens (GW170817 at z = 0.01)")
ax.axvspan(0.3, 3, color="gray", alpha=0.12, label="ET / CE reach with counterparts")
ax.set_xlabel("z"); ax.set_ylabel("D_L / D_L(flat LCDM), same H0"); ax.legend(fontsize=8)
fig.savefig("results/sirens.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_sirens.png", dpi=120, bbox_inches="tight")
