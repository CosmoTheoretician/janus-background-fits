import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g_lo, g_hi = map(float, open("results/growth_q0_range.txt").read().split())
JWST = 1 / (2 * 14.4)  # |q0| must be below this for the dust era to start before MoM-z14 (z = 14.4)

ROWS = [  # label, lo, hi, colour
    ("Cosmic chronometers", -0.075, 0.100, "C7"),
    ("Pantheon+ SNe", -0.034, -0.008, "C0"),
    ("DES-SN5YR SNe", -0.040, -0.015, "C0"),
    ("DESI DR2 BAO (r_d free)", 0.098, 0.104, "C3"),
    ("Combined CC+BAO+SN (MCMC)", 0.079, 0.091, "C3"),
    ("Growth fσ8 + σ8 prior", g_lo, g_hi, "C2"),
    ("JWST galaxies inside dust era", -JWST, 0.6, "C8"),
    ("Janus premise (q0 < 0)", -0.25, 0.0, "C1"),
    ("Petit & d'Agostini 2018", -0.102, -0.072, "k"),
]
fig, ax = plt.subplots(figsize=(9, 5))
for i, (lab, lo, hi, c) in enumerate(ROWS):
    y = len(ROWS) - i
    ax.barh(y, min(hi, 0.2) - max(lo, -0.25), left=max(lo, -0.25), height=0.5, color=c, alpha=0.6)
    ax.text(-0.255, y, lab, ha="right", va="center", fontsize=9)
ax.axvline(0, color="gray", lw=0.8, ls=":")
ax.axvline(-0.087, color="k", lw=0.8, ls="--")
ax.set_xlim(-0.25, 0.2); ax.set_yticks([]); ax.set_xlabel("q0 (68% windows; bars clipped to the plotted range)")
ax.set_title("Janus deceleration parameter q0: what each test allows")
fig.savefig("results/q0_map.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_q0_map.png", dpi=120, bbox_inches="tight")

lines = ["Allowed q0 windows (68%):"] + [f"  {lab:32s} {lo:+.3f} to {hi:+.3f}" for lab, lo, hi, _ in ROWS]
obs = [r for r in ROWS if "premise" not in r[0] and "Petit" not in r[0]]
lo, hi = max(r[1] for r in obs), min(r[2] for r in obs)
lines += ["", f"Intersection of all observational windows: {'EMPTY' if lo > hi else f'{lo:+.3f} to {hi:+.3f}'}",
          f"Intersection of SN + BAO alone: {'EMPTY' if max(-0.034, 0.098) > min(-0.008, 0.104) else 'non-empty'}",
          f"Intersection of growth + JWST + premise: {'EMPTY' if max(g_lo, -JWST) > min(g_hi, 0.0) else f'{max(g_lo, -JWST):+.3f} to {min(g_hi, 0.0):+.3f}'}"]
text = "\n".join(lines)
print(text)
open("results/q0_map.txt", "w").write(text + "\n")
