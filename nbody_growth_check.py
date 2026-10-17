"""Compare the growth of the lowest box modes in a 3-D twin-matter N-body run with linear two-fluid theory."""
import sys, json, glob, os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def mixed_expansion_runs():
    out = []
    for a in sorted(glob.glob(os.path.expanduser("~/scraper/sim_gpu2/out/*/assumptions.txt"))):
        r = os.path.dirname(a)
        if "switched on" not in open(a).read():
            continue
        p = sorted(glob.glob(r + "/pts_*.npz"))
        if p and np.load(p[0])["neg"].shape[0] > 0:
            out.append((r, len(p)))
    return out

if len(sys.argv) > 1:
    RUN = sys.argv[1]
else:
    cands = mixed_expansion_runs()
    print("mixed expansion runs (snapshots):", [(os.path.basename(r), n) for r, n in cands])
    RUN = cands[-1][0]
print("using run:", RUN)
cfg = json.load(open(RUN + "/config.json"))
cfg = {**cfg.get("derived", {}), **cfg.get("config", cfg)}
print("config keys:", sorted(cfg.keys()))
box = float(np.load(sorted(glob.glob(RUN + "/pts_*.npz"))[0])["box"])
N_FULL = int(cfg.get("npos", 4_000_000))
MASS_RATIO = float(cfg.get("mass_ratio", 64))
TJ = 0.2216                       # t_J,neg in code units (assumption A4)
T0 = float(cfg.get("t0_exp", 0.1)) * TJ
G_RHO_NEG = np.pi * 0.2821**2 / (0.0625 * box) ** 2 * (1.0 if box == 1.0 else 1.0)  # from lambda_J = sigma sqrt(pi/(G rho))
G_RHO_POS = G_RHO_NEG / MASS_RATIO
print(f"box = {box}, N = {N_FULL}, t0 = {T0:.4f}, G rho_neg = {G_RHO_NEG:.2f}, G rho_pos = {G_RHO_POS:.3f}")

# ---- read snapshots ----
files = sorted(glob.glob(RUN + "/pts_*.npz"))
d0 = np.load(files[0])
print("npz keys:", {k: d0[k].shape for k in d0.files})
arr3 = [k for k in d0.files if d0[k].ndim == 2 and d0[k].shape[1] == 3]
kpos = next((k for k in arr3 if "pos" in k.lower() or k.lower().startswith("p") or "plus" in k.lower()), None)
kneg = next((k for k in arr3 if k != kpos and ("neg" in k.lower() or k.lower().startswith("n") or "minus" in k.lower())), None)
if kpos is None or kneg is None:
    sys.exit(f"could not identify position arrays among {arr3}; edit kpos/kneg by hand")
tkey = next((k for k in d0.files if k.lower() in ("t", "time")), None)
snap_every = float(cfg.get("snap_every", -1))
print(f"using positions '{kpos}' (ordinary) and '{kneg}' (twin); time key: {tkey or f'index * snap_every={snap_every} t_J'}")

KV = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]], float)


def power(x):
    n = len(x)
    ph = 2 * np.pi * (x / box) @ KV.T
    dk = (np.exp(-1j * ph)).mean(axis=0)
    return np.abs(dk) ** 2 - 1.0 / n  # subsample shot noise removed


ts, Pp, Pn = [], [], []
for i, f in enumerate(files):
    d = np.load(f)
    t = float(d[tkey]) if tkey else i * snap_every * TJ
    ts.append(t)
    Pp.append(power(d[kpos].astype(float)).mean())
    Pn.append(power(d[kneg].astype(float)).mean())
ts, Pp, Pn = map(np.array, (ts, Pp, Pn))

# ---- linear theory on the same background ----
def a_of(t):
    return 1.0 if t < T0 else (t / T0) ** (2 / 3)


def H_of(t):
    return 0.0 if t < T0 else 2 / (3 * t)


def rhs(t, y):
    dp, vp, dn, vn = y
    a, H = a_of(t), H_of(t)
    s = 4 * np.pi / a**3
    return [vp, -2 * H * vp + s * (G_RHO_POS * dp - G_RHO_NEG * dn),
            vn, -2 * H * vn + s * (G_RHO_NEG * dn - G_RHO_POS * dp)]


tt = np.linspace(1e-6, ts.max(), 2000)
resp = {}
for name, y0 in [("from_pos", [1, 0, 0, 0]), ("from_neg", [0, 0, 1, 0])]:
    s = solve_ivp(rhs, [tt[0], tt[-1]], y0, t_eval=tt, rtol=1e-9, atol=1e-12, max_step=T0 / 20)
    resp[name] = s.y
P0 = 1.0 / N_FULL  # Poisson seed power per mode, each population
Pp_th = P0 * (resp["from_pos"][0] ** 2 + resp["from_neg"][0] ** 2)
Pn_th = P0 * (resp["from_pos"][2] ** 2 + resp["from_neg"][2] ** 2)

lines = [f"Run {os.path.basename(RUN)}: growth of the 9 lowest box modes vs linear two-fluid theory (no free parameters)", "",
         f"{'t/tJ':>6s} {'a':>7s} {'twin: sim':>11s} {'theory':>9s} {'ratio':>6s}   {'ordinary: sim':>14s} {'theory':>9s} {'ratio':>6s}"]
for t, pp, pn in zip(ts, Pp, Pn):
    th_n, th_p = np.interp(t, tt, Pn_th), np.interp(t, tt, Pp_th)
    sn, sp = np.sqrt(max(pn, 0)), np.sqrt(max(pp, 0))
    lines.append(f"{t/TJ:6.2f} {a_of(t):7.2f} {sn:11.2e} {np.sqrt(th_n):9.2e} {sn/np.sqrt(th_n):6.2f}   {sp:14.2e} {np.sqrt(th_p):9.2e} {sp/np.sqrt(th_p):6.2f}")
lines += ["", "delta_rms is sqrt of the noise-subtracted power averaged over the 9 modes; the subsample's own shot noise",
          f"(1/{len(d0[kpos])} per mode) is subtracted, so early values below ~1e-3 are noise-dominated and unreliable.",
          "Linear theory is pressureless (box = 16 twin Jeans lengths, pressure correction 1/256) and valid until",
          "nonlinear clumps couple into the box modes (clump counts start rising at ~8 t_J)."]
text = "\n".join(lines)
print("\n" + text)
open("results/nbody_growth_check.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(tt / TJ, np.sqrt(Pn_th), "C1", label="twin matter, linear theory")
ax.plot(tt / TJ, np.sqrt(Pp_th), "C0", label="ordinary matter, linear theory")
ax.plot(ts / TJ, np.sqrt(np.clip(Pn, 0, None)), "C1o", ms=4, label="twin matter, N-body")
ax.plot(ts / TJ, np.sqrt(np.clip(Pp, 0, None)), "C0s", ms=4, label="ordinary matter, N-body")
ax.axhline(1 / np.sqrt(len(d0[kpos])), color="gray", ls=":", label="subsample noise floor")
ax.set_yscale("log"); ax.set_xlabel("t / t_J,neg"); ax.set_ylabel("δ_rms at box scale (|k| = 1, √2)")
ax.legend(fontsize=8); ax.set_title("3-D expansion run vs linear two-fluid growth, no free parameters")
fig.savefig("results/nbody_growth_check.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_nbody_growth_check.png", dpi=120, bbox_inches="tight")
print("\nplot saved to results/nbody_growth_check.png")
