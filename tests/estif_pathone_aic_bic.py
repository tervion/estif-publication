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
    print(f"  Q3  Fitted CPL vs LCDM:  Delta chi2 = {dchi2:.2f} for 2 extra params")
    print(f"      -> evolving DE preferred at ~{np.sqrt(max(dchi2,0)):.1f} sigma-equivalent")
    print(f"         (consistent with DESI's own evolving-DE preference).")
    print(f"      Fitted CPL AIC = {aic_cpl:.2f}  vs frozen-eddy AIC = {aic_core:.2f}")
    print(f"      -> fitted thawing beats frozen by Delta AIC = {aic_core-aic_cpl:+.2f}.")
    print()
    print("="*74); print("THE PATH TWO TARGET (what a DERIVED thawing must achieve)")
    print("="*74)
    print("  A DERIVED thawing correction has k=0 (no free parameters), so its AIC")
    print("  is just its chi2. To beat the frozen eddy it must reach:")
    print(f"      chi2  <  {aic_core:.2f}   (i.e. chi2/N < {aic_core/N:.3f})")
    print("  To beat the FITTED CPL outright (and justify the whole Path Two effort)")
    print("  a derived k=0 model would need:")
    print(f"      chi2  <  {aic_cpl:.2f}   (i.e. chi2/N < {aic_cpl/N:.3f})")
    print("  For reference the fitted CPL reaches chi2/N = {:.3f} WITH 2 free params.".format(chi2_cpl/N))
    print("  So Path Two is worth it only if the vorticity stress tensor DERIVES a")
    print("  thawing that lands near chi2/N ~ 0.7 with zero free parameters. If it")
    print("  can only match the frozen eddy, Path One already has that -- for free.")
    print("="*74)

if __name__ == "__main__":
    main()
