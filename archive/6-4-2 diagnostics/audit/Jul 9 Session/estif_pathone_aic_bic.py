"""
PATH ONE - Information-criterion parity check (AIC / BIC)
========================================================
"Frozen eddy ties LCDM at chi2/N = 1.92" is a raw goodness-of-fit statement.
This script makes it RIGOROUS with AIC and BIC, which penalise free parameters,
and answers three questions on real DESI DR2 data:

  Q1  Is ESTIF-Core (frozen eddy, Omega_m = x0) statistically distinguishable
      from LCDM?  -> if |Delta AIC| < 2 they are indistinguishable.
  Q2  Is the RETIRED self-consistent tilt worse by the criteria (justifying
      its retirement)?
  Q3  The fitted evolving-w (CPL) fits better -- but it costs 2 free parameters.
      By how much must a DERIVED (zero-parameter) Path Two thawing beat the
      frozen eddy to actually win on AIC?  This quantifies the Path Two target.

Parameter counting (free parameters fit to THIS dataset, dark-energy sector):
  LCDM                : k = 0   (Omega_m from Planck; w = -1 fixed)
  ESTIF-Core frozen   : k = 0   (Omega_m = x0 predicted; w = -1 derived)
  Retired tilt (SC)   : k = 0 for DESI (N_MAX, B calibrated to EHT/Lambda, not
                        DESI); k = 2 shown too, since they are not free globally
  CPL (fitted)        : k = 2   (w0, wa fit to DESI here)

AIC = chi2 + 2k ;  BIC = chi2 + k*ln(N).  Lower is better. Delta < 2 ~ tie.

Run:  python3 estif_pathone_aic_bic.py
Deps: numpy, scipy, internet (DESI fetch, cached beside script).
"""
import os, urllib.request
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

CACHE_DIR = os.path.dirname(os.path.abspath(__file__))
MPC = 3.085677581e22
c = 2.99792458e8
H0 = 67.66 * 1000.0 / MPC
OMEGA_M_PLANCK = 0.3111
RD = 147.09
R_UNIV = 4.4e26
x0 = (c / H0) / R_UNIV
N_MAX, B = 33.265, 15.429

def observable(x):
    x = np.asarray(x, float); n = N_MAX * np.exp(-B * x)
    val = np.where(x > 0, x ** (2.0 * n), 0.0)
    beta = np.where(val >= 1.0, 0.0, np.sqrt(np.maximum(0.0, 1.0 - val)))
    return np.sqrt(beta)
OBS_NOW = float(observable(x0))

def H_lcdm(z): return H0 * np.sqrt(OMEGA_M_PLANCK*(1+z)**3 + (1-OMEGA_M_PLANCK))
def H_core(z): return H0 * np.sqrt(x0*(1+z)**3 + (1-x0))          # frozen eddy, Om=x0
def H_tilt(z):                                                    # retired self-consistent tilt
    matter = OMEGA_M_PLANCK*(1+z)**3; h = np.sqrt(matter + (1-OMEGA_M_PLANCK))
    for _ in range(400):
        x = x0*(1+z)/h; oz = float(observable(x))
        om = (1-OMEGA_M_PLANCK)*(OBS_NOW/oz)**2 if oz>0 else (1-OMEGA_M_PLANCK)
        hn = np.sqrt(matter+om)
        if abs(hn-h)<1e-13: h=hn; break
        h = 0.5*h+0.5*hn
    return H0*h
def H_cpl(z, w0, wa):
    fde = (1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z))
    return H0*np.sqrt(OMEGA_M_PLANCK*(1+z)**3 + (1-OMEGA_M_PLANCK)*fde)

def DH(z,Hf): return c/(Hf(z)*MPC)
def DM(z,Hf):
    if z<=0: return 0.0
    v,_=quad(lambda zp: c/(Hf(zp)*MPC),0,z,limit=200); return v
def DV(z,Hf): return (z*DH(z,Hf)*DM(z,Hf)**2)**(1/3)
def predict(rows,Hf):
    return np.array([(DV if q=='DV_over_rs' else DM if q=='DM_over_rs' else DH)(z,Hf)/RD
                     for z,_,q in rows])
def chi2(pred,obs,cov):
    d=obs-pred
    try: return float(d@np.linalg.inv(cov)@d)
    except np.linalg.LinAlgError: return float(np.sum((d/np.sqrt(np.diag(cov)))**2))

BASE2="https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2"
URLS={'m':f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_mean.txt",
      'c':f"{BASE2}/desi_gaussian_bao_ALL_GCcomb_cov.txt"}
def cache(k):
    p=os.path.join(CACHE_DIR,f"aicbic_{k}.txt")
    if not (os.path.exists(p) and os.path.getsize(p)>100):
        urllib.request.urlretrieve(URLS[k],p)
    return p

def main():
    rows=[]
    for line in open(cache('m')):
        line=line.strip()
        if line and not line.startswith('#'):
            p=line.split(); rows.append((float(p[0]),float(p[1]),p[2]))
    N=len(rows); cov=np.loadtxt(cache('c')).reshape(N,N); obs=np.array([r[1] for r in rows])
    lnN=np.log(N)

    # fit CPL (2 free params) to DESI
    def negfit(p): return chi2(predict(rows, lambda z: H_cpl(z,p[0],p[1])), obs, cov)
    fit = minimize(negfit, x0=[-0.85,-0.45], method="Nelder-Mead",
                   options={'xatol':1e-3,'fatol':1e-3,'maxiter':400})
    w0b, wab = fit.x; chi2_cpl = fit.fun

    chi2_lcdm = chi2(predict(rows,H_lcdm),obs,cov)
    chi2_core = chi2(predict(rows,H_core),obs,cov)
    chi2_tilt = chi2(predict(rows,H_tilt),obs,cov)

    print("="*74)
    print("PATH ONE - AIC / BIC PARITY CHECK vs REAL DESI DR2")
    print("="*74)
    print(f"  N data points = {N}   ln N = {lnN:.3f}")
    print()
    rows_out = [
        ("LCDM (Om=Planck, w=-1)",            chi2_lcdm, 0),
        ("ESTIF-Core frozen eddy (Om=x0)",    chi2_core, 0),
        ("Retired tilt SC  [k=0, DESI]",      chi2_tilt, 0),
        ("Retired tilt SC  [k=2, if counted]",chi2_tilt, 2),
        (f"CPL fitted (w0={w0b:+.2f},wa={wab:+.2f})", chi2_cpl, 2),
    ]
    print(f"  {'model':<38}{'chi2':>8}{'k':>4}{'chi2/N':>9}{'AIC':>9}{'BIC':>9}")
    print("  "+"-"*70)
    recs={}
    for name,ch,k in rows_out:
        aic=ch+2*k; bic=ch+k*lnN
        recs[name]=(ch,k,aic,bic)
        print(f"  {name:<38}{ch:>8.2f}{k:>4}{ch/N:>9.3f}{aic:>9.2f}{bic:>9.2f}")
    print()

    aic_lcdm=recs["LCDM (Om=Planck, w=-1)"][2]
    aic_core=recs["ESTIF-Core frozen eddy (Om=x0)"][2]
    aic_tilt0=recs["Retired tilt SC  [k=0, DESI]"][2]
    aic_cpl=recs[f"CPL fitted (w0={w0b:+.2f},wa={wab:+.2f})"][2]

    print("="*74); print("READING (Delta AIC: <2 tie, 4-7 considerable, >10 decisive)")
    print("="*74)
    print(f"  Q1  ESTIF-Core vs LCDM:   Delta AIC = {aic_core-aic_lcdm:+.2f}")
    if abs(aic_core-aic_lcdm)<2:
        print("      -> INDISTINGUISHABLE. Frozen eddy = LCDM statistically. The")
        print("         geometric Omega_m = x0 and the eddy reinterpretation cost")
        print("         nothing in fit quality. Path One parity is rigorous.")
    print()
    print(f"  Q2  Retired tilt vs frozen eddy: Delta AIC = {aic_tilt0-aic_core:+.2f} (k=0 count)")
    print(f"      (k=2 count: Delta AIC = {recs['Retired tilt SC  [k=2, if counted]'][2]-aic_core:+.2f})")
    if aic_tilt0-aic_core>10:
        print("      -> DECISIVELY WORSE. The information criteria confirm retiring the")
        print("         tilt cosmology: it is worse than the plain cosmological constant")
        print("         under it, by any parameter count.")
    print()
    dchi2 = chi2_lcdm - chi2_cpl
    print(f"  Q3  Fitted CPL vs LCDM (FIXED params): Delta chi2 = {dchi2:.2f} -> "
          f"~{np.sqrt(max(dchi2,0)):.1f} sigma")
    print("      *** CAUTION (Erratum v6.3.1): this holds H0, Om, rd FIXED to Planck.")
    print("      BAO measures D/rd (degenerate with rd) and depends on Om, so fixing")
    print("      them INFLATES the significance. The marginalized comparison follows.")
    print()

    # --- Marginalized comparison (Correction 4): Om and rd-scale free ---
    print("="*74); print("MARGINALIZED COMPARISON (Om and rd-scale free for every model)")
    print("="*74)
    def chi2_scaled(Hf, s):
        d = obs - predict(rows, Hf)/s
        try: return float(d@np.linalg.inv(cov)@d)
        except np.linalg.LinAlgError: return float(np.sum((d/np.sqrt(np.diag(cov)))**2))
    def H_lcdm_om(z, Om): return H0*np.sqrt(Om*(1+z)**3 + (1-Om))
    def H_cpl_om(z, Om, w0, wa):
        fde=(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z))
        return H0*np.sqrt(Om*(1+z)**3 + (1-Om)*fde)
    fL = minimize(lambda p: chi2_scaled(lambda z:H_lcdm_om(z,p[0]), p[1]), [0.31,1.0],
                  method="Nelder-Mead", options={'xatol':1e-4,'fatol':1e-4,'maxiter':900})
    fC = minimize(lambda p: chi2_scaled(lambda z:H_cpl_om(z,p[0],p[1],p[2]), p[3]),
                  [0.31,-0.85,-0.45,1.0], method="Nelder-Mead",
                  options={'xatol':1e-4,'fatol':1e-4,'maxiter':2500})
    chi2_L_m, chi2_C_m = fL.fun, fC.fun
    dchi2_m = chi2_L_m - chi2_C_m
    print(f"  LCDM  (Om,rd free):        chi2/N = {chi2_L_m/N:.3f}   (Om={fL.x[0]:.3f})")
    print(f"  CPL   (Om,w0,wa,rd free):  chi2/N = {chi2_C_m/N:.3f}   "
          f"(Om={fC.x[0]:.3f}, w0={fC.x[1]:+.2f}, wa={fC.x[2]:+.2f})")
    print(f"  Delta chi2 (LCDM-CPL) = {dchi2_m:.2f} for 2 params -> "
          f"~{np.sqrt(max(dchi2_m,0)):.1f} sigma  (vs ~{np.sqrt(max(dchi2,0)):.1f} sigma fixed)")
    print("  NOTE: the marginalized CPL best-fit often runs to a non-physical corner;")
    print("  BAO-alone does NOT robustly prefer a sensible thawing. The strong, sensible")
    print("  evolving-DE signal comes from DESI+CMB+SNe, which ESTIF must eventually face.")
    print()
    print("="*74); print("THE PATH TWO TARGET (honest, marginalized)")
    print("="*74)
    print(f"  The real bar is NOT 'reach chi2/N ~ 0.66' (a fixed-parameter artifact).")
    print(f"  Against marginalized LCDM (chi2/N = {chi2_L_m/N:.3f}), a DERIVED (k=0)")
    print(f"  thawing must clear the marginalized ~{np.sqrt(max(dchi2_m,0)):.1f} sigma BAO gap --")
    print(f"  an EASIER bar than advertised, but a SMALLER prize.")
    print("  A real Path Two claim must ALSO survive DESI + CMB + SNe with full")
    print("  marginalization (CLASS/Cobaya-level) -- beyond BAO-only, fixed-param chi2.")
    print("  If the derived thawing only matches the frozen eddy, Path One already has")
    print("  that (and ties LCDM) -- for free.")
    print("="*74)

if __name__ == "__main__":
    main()
