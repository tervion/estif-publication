"""
test_desi_wz_consistency.py  —  ESTIF v6.0  —  March 2026

TEST 1: DESI w(z) Consistency Check
ESTIF pre-existing prediction w_eff ≈ -1.08 vs DESI DR1 + DR2
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
import urllib.request

import estif_ec_gr_constants as const
import estif_ec_gr_model as estif

# ── constants ────────────────────────────────────────────────────────────────
H0  = const.H_0
c   = const.c
x0  = (c / H0) / 4.4e26
RD  = 147.09          # Planck 2018 sound horizon [Mpc]
MPC = 3.085677581e22  # metres per Mpc

# DESI DR2 best-fit w0wa  (arXiv:2503.14738 Table 7, DESI+CMB+Union3)
W0_DR2, WA_DR2, W0_ERR, WA_ERR = -0.73, -0.66, 0.10, 0.25
# DESI DR1 best-fit (arXiv:2404.03002, DESI+CMB+Union3)
W0_DR1, WA_DR1 = -0.65, -1.27

print("=" * 70)
print("TEST 1 — DESI w(z) CONSISTENCY CHECK")
print("ESTIF pre-existing prediction vs DESI DR1 + DR2")
print("=" * 70)
print(f"\n  ESTIF prediction (made before DR1): w_eff ≈ -1.08")
print(f"  Sound horizon: rd = {RD} Mpc (Planck 2018, unchanged by ESTIF)")
print(f"  DESI DR2 w0 (DESI+CMB+Union3): {W0_DR2} ± {W0_ERR}\n")

# ── load data ────────────────────────────────────────────────────────────────
BASE1 = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master"
BASE2 = f"{BASE1}/desi_bao_dr2"

URLS = {
    'dr1_mean': f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
    'dr1_cov':  f"{BASE1}/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
    'dr2_mean': f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    'dr2_cov':  f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt",
}
CACHE = {
    'dr1_mean': '/home/claude/desi_dr1_mean.txt',
    'dr1_cov':  '/home/claude/desi_dr1_cov.txt',
    'dr2_mean': '/home/claude/desi_dr2_mean.txt',
    'dr2_cov':  '/home/claude/desi_dr2_cov.txt',
}


def load_mean(key):
    p = CACHE[key]
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[key], p)
    rows = []
    for line in open(p):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[1]), parts[2]))
    return rows


def load_cov(key, n):
    p = CACHE[key]
    if not (os.path.exists(p) and os.path.getsize(p) > 100):
        urllib.request.urlretrieve(URLS[key], p)
    return np.loadtxt(p).reshape(n, n)


dr1 = load_mean('dr1_mean')
dr2 = load_mean('dr2_mean')
cov1 = load_cov('dr1_cov', len(dr1))
cov2 = load_cov('dr2_cov', len(dr2))
obs1, obs2 = np.array([r[1] for r in dr1]), np.array([r[1] for r in dr2])
err1, err2 = np.sqrt(np.diag(cov1)), np.sqrt(np.diag(cov2))
print(f"  Loaded: DR1={len(dr1)} bins, DR2={len(dr2)} bins\n")

# ── model predictions ────────────────────────────────────────────────────────

def H_lcdm(z):
    return H0 * np.sqrt(estif.OMEGA_M*(1+z)**3 + estif.OMEGA_LAMBDA)


def DM(z, Hf):
    if z <= 0:
        return 0.0
    val, _ = quad(lambda zp: c / (Hf(zp) * MPC), 0, z, limit=200)
    return val


def DH(z, Hf):
    return c / (Hf(z) * MPC)


def DV(z, Hf):
    return (z * DH(z, Hf) * DM(z, Hf)**2) ** (1/3)


def predict(rows, Hf):
    out = []
    for z, _, qty in rows:
        if qty == 'DV_over_rs':
            out.append(DV(z, Hf) / RD)
        elif qty == 'DM_over_rs':
            out.append(DM(z, Hf) / RD)
        elif qty == 'DH_over_rs':
            out.append(DH(z, Hf) / RD)
    return np.array(out)


print("  Computing predictions (takes ~30s)...", flush=True)
p_lcdm1  = predict(dr1, H_lcdm)
p_estif1 = predict(dr1, estif.H_estif)
p_lcdm2  = predict(dr2, H_lcdm)
p_estif2 = predict(dr2, estif.H_estif)
print("  Done.\n")

# ── chi-squared ──────────────────────────────────────────────────────────────

def chi2(pred, obs, cov):
    d = obs - pred
    try:
        return float(d @ np.linalg.inv(cov) @ d)
    except np.linalg.LinAlgError:
        return float(np.sum((d / np.sqrt(np.diag(cov)))**2))


c2_lcdm1  = chi2(p_lcdm1,  obs1, cov1)
c2_estif1 = chi2(p_estif1, obs1, cov1)
c2_lcdm2  = chi2(p_lcdm2,  obs2, cov2)
c2_estif2 = chi2(p_estif2, obs2, cov2)
n1, n2 = len(dr1), len(dr2)

print("=" * 70)
print("CHI-SQUARED RESULTS")
print("=" * 70)
print(f"  {'Dataset':<10} {'N':>4}  {'LCDM chi2/N':>13}  {'ESTIF chi2/N':>13}  {'Delta chi2':>12}")
print(f"  {'-'*58}")
print(f"  {'DESI DR1':<10} {n1:>4}  {c2_lcdm1/n1:>13.3f}  {c2_estif1/n1:>13.3f}  {c2_lcdm1-c2_estif1:>+12.3f}")
print(f"  {'DESI DR2':<10} {n2:>4}  {c2_lcdm2/n2:>13.3f}  {c2_estif2/n2:>13.3f}  {c2_lcdm2-c2_estif2:>+12.3f}")
print(f"\n  Delta chi2 > 0: ESTIF fits DESI better than LCDM")
print(f"  Delta chi2 < 0: LCDM fits DESI better than ESTIF")

# ── per-bin pulls ─────────────────────────────────────────────────────────────
pulls_l2 = [(p_lcdm2[i]  - obs2[i]) / err2[i] for i in range(n2)]
pulls_e2 = [(p_estif2[i] - obs2[i]) / err2[i] for i in range(n2)]

print("\n" + "=" * 70)
print("PER-BIN — DESI DR2")
print("=" * 70)
print(f"  {'z':<7} {'qty':<14} {'DESI DR2':>10} {'err':>7} {'LCDM':>10} {'ESTIF':>10} {'pull_L':>8} {'pull_E':>8}")
print("  " + "-" * 76)
for i, (z, obs_v, qty) in enumerate(dr2):
    fl = "ok" if abs(pulls_e2[i]) < 1 else ("!" if abs(pulls_e2[i]) < 2 else "!!")
    print(f"  {z:<7.3f} {qty:<14} {obs_v:>10.4f} {err2[i]:>7.4f} "
          f"{p_lcdm2[i]:>10.4f} {p_estif2[i]:>10.4f} "
          f"{pulls_l2[i]:>+8.3f} {pulls_e2[i]:>+8.3f}  {fl}")

n_1s = sum(abs(p) < 1 for p in pulls_e2)
n_2s = sum(abs(p) < 2 for p in pulls_e2)
print(f"\n  Within 1sigma: ESTIF {n_1s}/{n2}   LCDM {sum(abs(p)<1 for p in pulls_l2)}/{n2}")
print(f"  Within 2sigma: ESTIF {n_2s}/{n2}   LCDM {sum(abs(p)<2 for p in pulls_l2)}/{n2}")

# ── effective w(z) ────────────────────────────────────────────────────────────

def w_eff(z, dz=0.01):
    z = max(z, dz)
    hi = estif.omega_tilt(z + dz)
    lo = estif.omega_tilt(max(z - dz, 1e-4))
    dlnde_dz = (np.log(hi + 1e-30) - np.log(lo + 1e-30)) / (2*dz)
    return -1.0 + (1.0 + z) / 3.0 * dlnde_dz


w_z0   = w_eff(0.05)
w_z05  = w_eff(0.50)
w_z1   = w_eff(1.00)

print("\n" + "=" * 70)
print("EFFECTIVE w(z) FROM ESTIF")
print("=" * 70)
print(f"  w_eff(z=0.05) = {w_z0:.4f}   ESTIF prediction ≈ -1.08")
print(f"  w_eff(z=0.50) = {w_z05:.4f}")
print(f"  w_eff(z=1.00) = {w_z1:.4f}")
print(f"\n  DESI DR2 w0 = {W0_DR2} +/- {W0_ERR}")
print(f"  Pull of ESTIF w_z0 from DR2 w0: {(w_z0-W0_DR2)/W0_ERR:+.2f}sigma")
print(f"  Pull of -1.08 from DR2 w0:      {(-1.08-W0_DR2)/W0_ERR:+.2f}sigma")

# ── visualisation ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 12))
fig.suptitle(
    'ESTIF vs DESI DR1 + DR2 — Test 1: w(z) Consistency\n'
    f'Pre-existing prediction w_eff ≈ -1.08  |  '
    f'DESI DR2 w0 = {W0_DR2}±{W0_ERR} (DESI+CMB+Union3)',
    fontsize=13, fontweight='bold'
)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

z_s = np.linspace(0.1, 2.5, 150)
DM_l = np.array([DM(z, H_lcdm)        / RD for z in z_s])
DM_e = np.array([DM(z, estif.H_estif) / RD for z in z_s])
DH_l = np.array([DH(z, H_lcdm)        / RD for z in z_s])
DH_e = np.array([DH(z, estif.H_estif) / RD for z in z_s])

# panel 1: DM/rd
ax1 = fig.add_subplot(gs[0, 0:2])
ax1.plot(z_s, DM_l, 'royalblue', lw=2.0, label='LCDM')
ax1.plot(z_s, DM_e, 'tomato',    lw=2.5, ls='--', label='ESTIF')
for rows, errs, obs, col, mk, lbl in [
    (dr1, err1, obs1, 'gray',  'o', 'DESI DR1'),
    (dr2, err2, obs2, 'black', 's', 'DESI DR2'),
]:
    pts = [(rows[i][0], obs[i], errs[i]) for i in range(len(rows)) if rows[i][2]=='DM_over_rs']
    if pts:
        zd, vd, ed = zip(*pts)
        ax1.errorbar(zd, vd, yerr=ed, fmt=mk, color=col, ms=7, capsize=4,
                     label=lbl, zorder=5, alpha=0.85 if col=='gray' else 1.0)
ax1.set_xlabel('Redshift z', fontsize=11); ax1.set_ylabel('DM / rd', fontsize=11)
ax1.set_title('Transverse Comoving Distance DM/rd', fontsize=11, fontweight='bold')
ax1.legend(fontsize=9); ax1.grid(alpha=0.3)

# panel 2: DH/rd
ax2 = fig.add_subplot(gs[0, 2])
ax2.plot(z_s, DH_l, 'royalblue', lw=2.0, label='LCDM')
ax2.plot(z_s, DH_e, 'tomato',    lw=2.5, ls='--', label='ESTIF')
for rows, errs, obs, col, mk, lbl in [
    (dr1, err1, obs1, 'gray',  'o', 'DR1'),
    (dr2, err2, obs2, 'black', 's', 'DR2'),
]:
    pts = [(rows[i][0], obs[i], errs[i]) for i in range(len(rows)) if rows[i][2]=='DH_over_rs']
    if pts:
        zd, vd, ed = zip(*pts)
        ax2.errorbar(zd, vd, yerr=ed, fmt=mk, color=col, ms=6, capsize=4,
                     label=lbl, alpha=0.85 if col=='gray' else 1.0)
ax2.set_xlabel('Redshift z', fontsize=11); ax2.set_ylabel('DH / rd', fontsize=11)
ax2.set_title('Hubble Distance DH/rd', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9); ax2.grid(alpha=0.3)

# panel 3: w(z)
ax3 = fig.add_subplot(gs[1, 0])
z_wp = np.linspace(0.05, 2.0, 120)
w_e  = [w_eff(z) for z in z_wp]
w_d2 = [W0_DR2 + WA_DR2*z/(1+z) for z in z_wp]
w_d1 = [W0_DR1 + WA_DR1*z/(1+z) for z in z_wp]
w_hi = [(W0_DR2+W0_ERR) + (WA_DR2+WA_ERR)*z/(1+z) for z in z_wp]
w_lo = [(W0_DR2-W0_ERR) + (WA_DR2-WA_ERR)*z/(1+z) for z in z_wp]
ax3.fill_between(z_wp, w_lo, w_hi, alpha=0.25, color='green', label=f'DESI DR2 ±1s')
ax3.plot(z_wp, w_d2, 'green',  lw=2,   label=f'DR2 CPL (w0={W0_DR2})')
ax3.plot(z_wp, w_d1, 'purple', lw=1.5, ls=':', label=f'DR1 CPL')
ax3.plot(z_wp, w_e,  'tomato', lw=2.5, label='ESTIF w_eff(z)')
ax3.axhline(-1,    color='black', lw=1.5, ls='--', label='w=-1 (LCDM)')
ax3.axhline(-1.08, color='tomato', lw=1, ls=':', alpha=0.6, label='prediction -1.08')
ax3.set_xlabel('z', fontsize=10); ax3.set_ylabel('w(z)', fontsize=10)
ax3.set_title('Dark Energy EoS w(z)', fontsize=10, fontweight='bold')
ax3.legend(fontsize=7.5); ax3.grid(alpha=0.3); ax3.set_ylim(-2.2, 0.3)

# panel 4: pull bars (DR2)
ax4 = fig.add_subplot(gs[1, 1])
x_pos = np.arange(n2)
cols  = ['green' if abs(p)<1 else ('orange' if abs(p)<2 else 'red') for p in pulls_e2]
ax4.bar(x_pos, pulls_e2, color=cols, alpha=0.85, edgecolor='black', lw=0.7)
ax4.axhline(0, color='black', lw=1.5)
for lv, cl, lb in [(1,'gray','±1s'),(-1,'gray',''),  (2,'orange','±2s'),(-2,'orange','')]:
    ax4.axhline(lv, color=cl, lw=1, ls='--' if abs(lv)==1 else ':', alpha=0.7,
                label=lb if lb else None)
ax4.set_xticks(x_pos)
ax4.set_xticklabels([f"z={r[0]:.3f}\n{r[2]}" for r in dr2], fontsize=6.5,
                     rotation=35, ha='right')
ax4.set_ylabel('(model - data) / sigma', fontsize=10)
ax4.set_title(f'ESTIF Pull per Bin (DESI DR2)\n{n_1s}/{n2} within 1s, {n_2s}/{n2} within 2s',
              fontsize=10, fontweight='bold')
ax4.legend(fontsize=8); ax4.grid(axis='y', alpha=0.3)

# panel 5: chi2 bars
ax5 = fig.add_subplot(gs[1, 2])
lbls = ['DR1 LCDM','DR1 ESTIF','DR2 LCDM','DR2 ESTIF']
vals = [c2_lcdm1/n1, c2_estif1/n1, c2_lcdm2/n2, c2_estif2/n2]
brs  = ax5.bar(lbls, vals,
               color=['royalblue','tomato','royalblue','tomato'],
               alpha=0.85, edgecolor='black', lw=0.9)
ax5.axhline(1.0, color='green',  lw=2, ls='--', label='ideal chi2/N=1')
ax5.axhline(2.0, color='orange', lw=1.5, ls=':', label='tension chi2/N=2')
for b, v in zip(brs, vals):
    ax5.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}',
             ha='center', va='bottom', fontsize=9, fontweight='bold')
ax5.set_ylabel('chi2 / N', fontsize=10)
ax5.set_title('Goodness of Fit\nESTIF vs LCDM', fontsize=10, fontweight='bold')
ax5.legend(fontsize=9); ax5.grid(axis='y', alpha=0.3)
ax5.set_ylim(0, max(vals)*1.35)

plt.savefig('desi_wz_consistency.png', dpi=150, bbox_inches='tight')
plt.close()
print("\nPlot saved: desi_wz_consistency.png")

# ── final assessment ──────────────────────────────────────────────────────────
chi2_ok   = c2_estif2/n2 < 2.0
chi2_good = c2_estif2/n2 < 1.5
bins_2s   = n_2s == n2
bins_1s   = n_1s >= int(n2 * 0.7)
w_ok      = abs(w_z0 - W0_DR2) < 2*W0_ERR
better    = c2_estif2 <= c2_lcdm2

all_good = chi2_ok and bins_2s and w_ok
mixed    = chi2_ok and not (bins_1s and w_ok)

result = "GOOD" if all_good else ("MIXED" if mixed else "BAD")

print("\n" + "=" * 70)
print("FINAL ASSESSMENT")
print("=" * 70)
print(f"""
  Results: {result}

  Criteria:
    chi2/N < 2.0  (no gross tension):       {'PASS' if chi2_ok   else 'FAIL'}
    chi2/N < 1.5  (good fit):               {'PASS' if chi2_good else 'FAIL'}
    All bins within 2sigma (DR2):           {'PASS' if bins_2s   else 'FAIL'}
    70%+ bins within 1sigma (DR2):          {'PASS' if bins_1s   else 'FAIL'}
    ESTIF w_eff(z~0) within 2sigma of DR2:  {'PASS' if w_ok      else 'FAIL'}
    ESTIF chi2 <= LCDM chi2 (DR2):          {'PASS' if better    else 'FAIL'}

  Pre-existing prediction w_eff ≈ -1.08:
    DESI DR2 w0 = {W0_DR2} ± {W0_ERR}
    Pull of -1.08 from DR2 w0: {(-1.08-W0_DR2)/W0_ERR:+.2f}sigma
    Status: {'CONSISTENT (within 2sigma)' if abs(-1.08-W0_DR2) < 2*W0_ERR else 'TENSION'}

  Project status:
""")

if all_good:
    print("    The ESTIF cosmological sector passes the DESI DR2 test.")
    print("    The pre-existing w_eff prediction is consistent with new data.")
    print("    The model is worthy of continued development.")
elif mixed:
    print("    ESTIF is broadly consistent with DESI DR2 but shows tension in")
    print("    specific bins. The cosmological sector needs attention but the")
    print("    project should not be abandoned — refine the high-z behaviour.")
else:
    print("    ESTIF shows significant tension with DESI DR2.")
    print("    The cosmological sector requires revision before proceeding.")

print(f"""
  Recommended action: {'CONTINUE to Test 2' if all_good or mixed else 'RE-EXAMINE cosmological sector'}
""")
print("=" * 70)
