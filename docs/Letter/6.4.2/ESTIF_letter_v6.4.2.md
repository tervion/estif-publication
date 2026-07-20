# Geometric Construction of the MOND Critical Acceleration from 4D Hypersurface Tilt Geometry: Horizon-Scale Origin, Prediction, and SPARC Validation

**Peter Angelov**
Independent Researcher · tervion@gmail.com
March 2026 (v6.2) · updated July 2026 (v6.3) · updated July 2026 (v6.4) · updated July 2026 (v6.4.2)

---

## Abstract

The MOND critical acceleration $a_0 \approx 1.2 \times 10^{-10}\,\mathrm{m/s^2}$ has been empirically successful at predicting galactic rotation curves for 40 years, yet its physical origin has never been derived from first principles. We present a motivated geometric construction within the Emergent Spacetime from Inward Flow (ESTIF) framework that reproduces the MOND normalization without internal tuning, using only the Planck 2018 values of $H_0$ and $\Omega_m$ as independently measured inputs. The four-step construction gives $a_0 = H_0 c x_0 / \sqrt{3} = 1.179 \times 10^{-10}\,\mathrm{m/s^2}$, agreeing with the MOND empirical value to 1.72% with zero free parameters in the construction itself. Under the ESTIF horizon reading (RHAC-010), the derived content is the horizon-scale form $a_0 \propto cH_0$ — the de Sitter surface-gravity scale — from which mass-independence, asymptotically constant rotation curves, and the baryonic Tully-Fisher scaling already follow; the construction supplies a candidate $O(1)$ coefficient, $x_0/\sqrt{3}$, whose first-principles derivation remains open. We prove that $a_0$ is exactly constant across cosmic time: in the local non-expanding frame of bound systems, which is the physically correct frame for galaxy dynamics, the Hubble parameter $H(z)$ cancels algebraically, giving $a_0 = c^2 / (r_\mathrm{universe}\sqrt{3}) = \text{constant}$. This constancy emerges from the frame definition, not from an assumption. It is confirmed consistent with high-redshift Tully-Fisher observations at $z \approx 0.75\text{--}2.2$. Validation against 87 quality-1 SPARC galaxies yields RMS = 15.6%, within the observed scatter of the baryonic Tully-Fisher relation. The formula is robust across all published cosmological datasets: every combination of $H_0 \in [65, 75]\,\mathrm{km/s/Mpc}$ and $\Omega_m \in [0.27, 0.33]$ remains within SPARC observational scatter.

**Code:** https://github.com/tervion/estif-publication

---

## 1. Introduction

Modified Newtonian Dynamics (MOND; Milgrom 1983) has been empirically successful at predicting galactic rotation curves for over 40 years. Its central claim is that below a critical acceleration $a_0 \approx 1.2 \times 10^{-10}\,\mathrm{m/s^2}$, galaxy dynamics transition to a regime that naturally produces flat rotation curves. Despite this empirical success, the physical origin of $a_0$ has remained unexplained: it is fitted, not derived. Its numerical proximity to $cH_0$ has been noted repeatedly (e.g. Milgrom 1999; Famaey & McGaugh 2012), but no first-principles derivation of the exact numerical coefficient has been established.

Several approaches have attempted to give $a_0$ a physical footing. Entropic gravity proposals (Verlinde 2016) reproduce MOND-like behaviour but require additional assumptions about dark energy entanglement. Covariant formulations (Bekenstein 2004; Skordis & Zlosnik 2021) take $a_0$ as an input parameter. No approach has derived the specific value from independently measured cosmological quantities without fitting.

Here we present a geometric construction within the ESTIF framework (Angelov 2026) that reproduces $a_0$ without internal tuning. We present the construction, prove redshift constancy, demonstrate cosmological robustness, and validate against 87 quality-1 SPARC galaxies. We state explicitly what has been established and what remains as future theoretical work.

---

## 2. The ESTIF Framework

### 2.1 Background

ESTIF proposes that 3D space is a hypersurface moving through 4D space, with a background eddy whose cosmological curvature ratio $x_0 = R_H / r_\mathrm{universe}$ identifies with the matter density $\Omega_m$ to 0.12% (a consistency relation: $r_\mathrm{universe}$ here is the $\Lambda$CDM particle horizon, itself dependent on $\Omega_m$). The underlying tilt formula has been calibrated separately to strong-field observational data; in this paper we use only the background eddy scale at cosmic curvature $x = x_0$.

> The tilt formula $n(x) = N_\mathrm{MAX} \times \exp(-B \times x)$ with $N_\mathrm{MAX} = 33.265$ and $B = 15.429$ simultaneously satisfies EHT M87\* shadow (0.00σ tension), Planck $\Lambda$ (ratio 1.0000), and the predicted LISA GW delay (491 µs, S/N = 49σ) with no free parameters after calibration. **The LISA deviation-from-GR figure is conditional on the eddy-stress sector remaining as currently specified (C1); the ESTIF vacuum solution is exactly Schwarzschild, so any deviation from GR enters only through that sector.** The strong-field results are independent of the MOND construction; we mention them only to establish that the framework is internally consistent at multiple scales.

### 2.2 The gravity–eddies identity

The eddy spin rate is $\omega(x) = H_0 \times x^{n(x)}$. Gravitational acceleration equals the gradient of eddy spin energy:

$$a = -\frac{c^2}{2} \nabla(\omega/H_0)^2 = \frac{GM}{r^2} \quad \text{(exact at the GR crossover } n = 1/2\text{)}$$

This force law is verified numerically to 10 decimal places and is the starting point for the MOND construction.

---

## 3. Construction of $a_0$

We construct $a_0$ in four steps, stating at each step what is derived, what is motivated, and what is assumed.

### Step 1 — Force law (exact)

The gradient of the eddy spin energy recovers Newton exactly, with no approximation:

$$a = -\frac{c^2}{2} \nabla(\omega/H_0)^2 = \frac{GM}{r^2}$$

This is the established ESTIF force law, not an assumption introduced for the MOND construction.

### Step 2 — Cosmological velocity scale

The eddy background velocity at the Hubble radius, weighted by the matter fraction $x_0 = R_H / r_\mathrm{universe} = \Omega_m$, gives:

$$v_\mathrm{flow} = c \times x_0 = c \times \Omega_m$$

Two independent routes confirm this: (i) $x_0$ is the matter fraction, so the matter-weighted flow speed is $cx_0$; (ii) the tilt formula evaluated at cosmic curvature $x = x_0$ gives the same result. $H_0$ and $x_0$ are taken from Planck 2018 as independently measured inputs; they are not adjusted to improve agreement with MOND. Under the horizon reading (Section 6.1), $cH_0$ is the load-bearing scale; $x_0$ enters as part of the $O(1)$ coefficient.

### Step 3 — 3D isotropic projection

The eddy background velocity $v_\mathrm{flow}$ is isotropic in 3D space. Projecting onto one spatial dimension:

$$v_\mathrm{3D} = v_\mathrm{flow} / \sqrt{3}$$

This is motivated by 3D spatial isotropy and consistent with the equipartition theorem applied to 3 spatial dimensions: $\langle v^2 \rangle = \langle v_x^2 \rangle + \langle v_y^2 \rangle + \langle v_z^2 \rangle$, giving $v_\mathrm{1D} = v_\mathrm{rms}/\sqrt{3}$. The same factor appears in kinetic theory ($c_s = v_\mathrm{rms}/\sqrt{3}$), the Jeans criterion, and the ESTIF parameter $B = L/3$. A complete kinetic theory of the eddy background — formally specifying its distribution function and equation of state — is identified as future theoretical work; the factor is physically motivated, not formally derived.

The uniqueness of $1/\sqrt{3}$ was confirmed against 12 candidate geometric factors ($1$, $1/\sqrt{2}$, $1/2$, $1/\pi$, $1/\sqrt{2\pi}$, $x_0$, $\sqrt{x_0}$, $2/3$, $1/\sqrt{4\pi/3}$, $1/e$, $1/\sqrt{4}$). It is the only factor giving < 5% agreement with the MOND empirical value and the only one with an independent physical derivation. This is a confirmatory test, not a mathematical uniqueness proof from first principles; the latter requires the complete kinetic theory noted above.

### Step 4 — MOND acceleration threshold

The natural acceleration scale associated with $v_\mathrm{3D}$ and $H_0$ is:

$$a_0 = v_\mathrm{3D} \times H_0 = H_0 \times c \times x_0 / \sqrt{3}$$

### Numerical result

Using Planck 2018 values ($H_0 = 67.66\,\mathrm{km/s/Mpc}$, $x_0 = 0.3107 = \Omega_m$ to 0.12%):

$$a_0 = 1.179 \times 10^{-10}\,\mathrm{m/s^2} \quad (-1.72\% \text{ from MOND empirical } 1.200 \times 10^{-10}\,\mathrm{m/s^2})$$

**Free parameters in the construction itself: zero.** $H_0$ and $x_0 = \Omega_m$ are independently measured cosmological quantities, not adjusted to improve agreement. (The underlying tilt formula has separately calibrated parameters; those are not used here.)

---

## 4. Redshift Constancy of $a_0$

### 4.1 Frame definition and the origin of constancy

The MOND acceleration scale $a_0$ is inferred from galaxy rotation curves and therefore corresponds to a quantity defined in the local rest frame of gravitationally bound systems. Such systems are decoupled from the Hubble expansion; their internal dynamics are governed by physical (non-expanding) coordinates rather than comoving coordinates that scale with the cosmic expansion factor.

A naive redshift dependence $a_0 \propto H(z)$ implicitly assumes that the relevant length scales evolve with the background expansion. This assumption is not applicable to galaxies, whose sizes and internal kinematics remain approximately constant in physical coordinates across cosmic time.

In the ESTIF framework, the dimensionless geometric factor $x(z)$ is defined in the same local physical frame in which galaxy dynamics are measured. When expressed consistently in this frame, the redshift dependence of the Hubble parameter is compensated by the corresponding scaling of the geometric ratio, leading to an exact cancellation. As a result, the predicted acceleration scale remains invariant:

$$a_0(z) = \text{constant}$$

This invariance is not imposed — it follows from evaluating both cosmological and dynamical quantities in the physically relevant non-expanding frame of bound systems. Any alternative that applies cosmological expansion directly to bound galactic systems would predict strong redshift evolution of rotation curves, which is ruled out by observations.

### 4.2 Algebraic proof

In the comoving frame, $r_\mathrm{universe}$ is constant. The ESTIF generalisation of $x_0$ to redshift $z$ is:

$$x(z) = R_H(z) / r_\mathrm{universe,comoving} = c / [H(z) \times r_\mathrm{universe,comoving}]$$

Substituting:

$$a_0(z) = H(z) \times c \times x(z) / \sqrt{3} = c^2 / (r_\mathrm{universe,comoving} \times \sqrt{3})$$

The $H(z)$ cancels exactly. Numerical verification across $z = 0$ to $z = 10$ with both $H_{\Lambda\mathrm{CDM}}(z)$ and $H_\mathrm{ESTIF}(z)$ shows the deviation from constancy is $2.22 \times 10^{-16}$, floating-point machine epsilon. This is an algebraic identity. The non-comoving $x(z) = x_0 \times (1+z) \times H_0 / H(z)$, appropriate for expansion-history calculations, would give $a_0(z) = (1+z) \times a_0(0)$, rising 3–4× by $z = 2\text{--}3$ and clearly ruled out by observations. Using the comoving frame for galaxy dynamics is the physically correct choice.

### 4.3 Observational confirmation

Di Teodoro et al. (2021), Übler et al. (2017), and Tiley et al. (2019) find the baryonic Tully-Fisher normalisation consistent with constant $a_0$ at $z = 0.75\text{--}2.2$, with deviations ≤ 2σ. ESTIF predicts exactly zero evolution.

---

## 5. Validation Against SPARC

### 5.1 Data and method

We test $v_\mathrm{flat} = (G \times M_\mathrm{bar} \times a_0)^{1/4}$ against the SPARC catalog (Lelli et al. 2016, AJ 152, 157), 175 disk galaxies with Spitzer 3.6 µm photometry and high-quality rotation curves. Baryonic masses are $M_\mathrm{bar} = 0.50 \times L_{3.6} + 1.33 \times M_\mathrm{HI}$ with $\Upsilon_* = 0.50\,M_\odot/L_\odot$ from McGaugh & Schombert (2014), independent of MOND. No parameters are adjusted: $a_0$ is the value constructed in Section 3.

### 5.2 Results

The quality-1 sample (87 galaxies with the most reliable rotation curves) gives RMS = 15.6%, consistent with the observed scatter of the baryonic Tully-Fisher relation (15–20% in $v_\mathrm{flat}$; Lelli et al. 2016). 82% of galaxies fall within 20%, and 97% within 30%.

| Sample      | N   | RMS error | Within 20%    | Pass? |
|-------------|-----|-----------|---------------|-------|
| Quality-1   | 87  | 15.6%     | 82% (71/87)   | ✓     |
| Quality-1+2 | 129 | 18.4%     | 84% (108/129) | ✓     |

*Table 1. SPARC BTFR results using $a_0 = H_0 c x_0/\sqrt{3}$ (zero free parameters in the construction).*

### 5.3 Bias analysis

A systematic mean bias of −7.6% is present. Investigation shows this correlates with morphological type, gas fraction, and surface brightness before stellar mass correction, but all correlations become non-significant ($p > 0.16$) after correcting $\Upsilon_*$ to $0.85\,M_\odot/L_\odot$. The bias originates in the stellar mass calibration, not in the force law structure. The literature range $\Upsilon_* = 0.60\text{--}0.70$ for 3.6 µm reduces the bias to 4–5%. We adopt $\Upsilon_* = 0.50$ as the standard benchmark and note the offset explicitly.

### 5.4 Parameter robustness

Testing $a_0$ across 3,600 combinations of $H_0 \in [65, 75]\,\mathrm{km/s/Mpc}$ and $\Omega_m \in [0.27, 0.33]$ shows every combination within ±20% SPARC scatter. Eight published datasets (Planck, WMAP, SH0ES, DES, KiDS, SPT, ACT, H0LiCOW) all pass. The Planck–SH0ES Hubble tension shifts $a_0$ by only 4.1%. The formula is not critically dependent on any particular cosmological measurement.

---

## 6. Discussion

### 6.1 What has been established

The quantity $a_0 = H_0 c x_0 / \sqrt{3}$ reproduces the MOND critical acceleration to 1.72% without internal tuning. Under the horizon reading (RHAC-010), the derived content is the scale: $a_0$ is a horizon-scale acceleration, $a_0 \approx cH$ — the de Sitter surface gravity — and mass-independence, asymptotically constant rotation curves, the baryonic Tully-Fisher scaling $v^4 = GMa_0$, and the correct order of magnitude follow from this alone. The specific $O(1)$ coefficient $x_0/\sqrt{3} \approx 0.18$ is this construction's proposal: the $1/\sqrt{3}$ factor is motivated by 3D spatial isotropy and $x_0$ by the matter-fraction weighting, but the precise number is not claimed as derived and remains open. The redshift constancy is proved algebraically and emerges from evaluating the formula in the locally non-expanding frame natural to galaxy dynamics. The formula is robust across cosmological parameter uncertainty and consistent with SPARC. These are the specific claims.

### 6.2 Comparison to prior work

The proximity $a_0 \approx cH_0$ has been noted since Milgrom (1999). The specific form $H_0 c x_0/\sqrt{3} = H_0 c \Omega_m / \sqrt{3}$ differs from prior attempts by providing a geometric construction of the numerical coefficient through the isotropy argument, embedded within a framework in which $\Omega_m = x_0$ holds to 0.12% as a consistency relation. Verlinde (2016) obtains MOND-like behaviour from entropic gravity but does not predict the specific value of $a_0$ without additional assumptions. No prior work derives the exact coefficient from an independently motivated geometric argument.

### 6.3 Honest limitations

The following are explicitly open:

(i) **The coefficient is proposed, not derived.** *(v6.4.2 update, per the July 2026 horizon doctrine, RHAC-010.)* $a_0$ is a horizon-scale acceleration, $a_0 \approx cH$ (de Sitter surface gravity): the scale is derived, and mass-independence, asymptotically constant rotation curves, the BTFR, and the correct magnitude follow from it alone. The precise $O(1)$ prefactor — here $x_0/\sqrt{3}$ — is horizon-set but open: the $1/\sqrt{3}$ projection is motivated by 3D spatial isotropy, and a complete kinetic theory of the flow background would be required to derive it.

(ii) **The tilt sector governs strong-field and cosmological regimes; galactic dynamics enter through the separately constructed $a_0$. The interpolation function $\mu(a/a_0)$ remains open — the local-$n$ reading of the tilt force law collapses to zero force at galactic accelerations (verified numerically, `mu_extraction.py`), and no per-system evaluation rule for $n(x)$ is yet derived. The RAR shape is therefore inherited from the deep-MOND normalization, not derived from the tilt sector.** *(v6.3 update)*

(iii) The identity $\Omega_m = x_0$ is confirmed numerically to 0.12% but not yet derived from the 4D stress-energy tensor $T_{\mu\nu}$ projection. The derivability question is now closed in the negative for principle P ($\Omega_m = R_H/r_p$): the equality holds only near $a \approx 1$ (the two quantities cross once, today), so it cannot follow from the time-symmetric axioms; $\Omega_m = x_0$ stands as a consistency relation and P as a predictive postulate. *(v6.4.2 update, RHAC-009.)*

(iv) Galactic halo structure requires N-body simulation with the ESTIF force law.

(v) The ESTIF cosmological dark energy sector has been fully re-examined since v6.3. Path Two (a dynamical, thawing dark-energy sector) was pre-registered and closed: the derivation returns $w = -1$ exactly, an honorable null -- no dynamical dark energy is found, and none is claimed. Path One (frozen-eddy, constant-$\Lambda$ limit) ties $\Lambda$CDM on DESI DR2, $\chi^2/N = 1.965$ (13-bin pipeline). No cosmological claims beyond this are made here. *(v6.4 update, supersedes the v6.3 note.)*

(vi) ESTIF reproduces the normalization of $a_0$ but inherits MOND's known cluster-scale missing-mass problem; this is not addressed here.

(vii) **Weak-lensing BTFR consistency and residual $a_0$-normalization gap.** *(v6.3 addition)* Inverting the baryonic Tully-Fisher relation for the $a_0$ it prefers yields +1.17σ tension from the kinematic BTFR (Lelli+2019 binned) and +1.52σ from the weak-lensing BTFR (Mistele et al. 2024, ApJL 969 L3, primary sample $R < 1000\,\mathrm{kpc}$), both measured against the bootstrap value $a_0 = 1.192\times10^{-10}\,\mathrm{m/s^2}$ used internally by these two scripts -- not the paper's headline constructed value $a_0 = 1.179\times10^{-10}\,\mathrm{m/s^2}$ (Section 3); relative to the headline value the tensions are marginally larger (kinematic ≈ +1.21σ). Both figures are fully accounted for by the paper-adopted 0.1 dex $M_*/L$ systematic treated as fully correlated across mass bins. A prior report of +2.67σ / +2.12σ (`btfr_lensing.py`) incorrectly propagated the correlated systematic per-bin; the corrected inversion (`a0_tension_corrected.py`) supersedes it. The residual gap between RAR-calibrated $a_0$ (0.66% agreement) and BTFR-normalized $a_0$ is inherited from the deep-MOND construction, is partially $H_0$-convention-dependent for the lensing sample and $h$-independent for the kinematic sample, and remains unresolved.

(viii) **Bullet Cluster / cluster-scale lensing-baryon offset.** *(v6.4 addition, flagged open -- not resolved)* The Bullet Cluster is the sharpest instance of the cluster-scale problem in (vi). No ESTIF-specific analysis of this system has been carried out; this is flagged here as an explicit open item rather than left unaddressed by omission.

---

## 7. Conclusions

We have presented a four-step geometric construction reproducing the MOND critical acceleration $a_0$ from the Planck 2018 values of $H_0$ and $\Omega_m$, with zero free parameters in the construction itself. The result $a_0 = H_0 c x_0 / \sqrt{3} = 1.179 \times 10^{-10}\,\mathrm{m/s^2}$ agrees with the MOND empirical value to 1.72%. Under the horizon reading, the derived content is the $cH_0$ scale — the de Sitter surface gravity of the horizon — and the precise $O(1)$ coefficient remains open.

We have proved that $a_0$ is exactly constant across cosmic time by an algebraic identity: in the local non-expanding frame of bound systems — the physically correct frame for galaxy dynamics — the Hubble parameter $H(z)$ cancels exactly. This constancy is not assumed; it emerges from the frame definition. Any alternative that applies cosmological expansion directly to bound galactic systems would predict strong redshift evolution of rotation curves, which is ruled out by observations.

Validation against 87 quality-1 SPARC galaxies gives RMS = 15.6%, within the observed baryonic Tully-Fisher scatter. The formula is robust across all published cosmological datasets; the full Hubble tension shifts $a_0$ by only 4.1%.

The specific claims are: (1) $a_0$ is reproduced without fitting; (2) the $1/\sqrt{3}$ factor is unique among simple geometric candidates and motivated by isotropy — under the horizon doctrine the derived content is the $cH_0$ scale, and the precise $O(1)$ coefficient is open; (3) redshift constancy follows algebraically from the correct frame; (4) cosmological robustness is confirmed across 3,600 parameter combinations; (5) the result is consistent with SPARC. Open questions are explicitly stated.

Code, tests, and the full validation suite are available at https://github.com/tervion/estif-publication.

---

## Acknowledgements

The author thanks the SPARC collaboration (Lelli, McGaugh, Schombert) for making the catalog publicly available. This research made use of the VizieR catalog service (Ochsenbein et al. 2000).

---

## References

- Begeman K. G., Broeils A. H., Sanders R. H., 1991, MNRAS, 249, 523
- Bekenstein J. D., 2004, Phys. Rev. D, 70, 083509
- Di Teodoro E. M. et al., 2021, A&A, 655, A82
- Event Horizon Telescope Collaboration, 2019, ApJL, 875, L1
- Famaey B., McGaugh S. S., 2012, Living Rev. Rel., 15, 10
- Lelli F., McGaugh S. S., Schombert J. M., 2016, AJ, 152, 157
- McGaugh S. S., Schombert J. M., 2014, AJ, 802, 18
- Meidt S. E. et al., 2014, ApJ, 788, 144
- Milgrom M., 1983, ApJ, 270, 365
- Milgrom M., 1999, Phys. Lett. A, 253, 273
- Mistele T., McGaugh S., Lelli F., Schombert J., Li P., 2024, ApJL, 969, L3 *(v6.3 addition)*
- Planck Collaboration, 2020, A&A, 641, A6
- Riess A. G. et al., 2022, ApJL, 934, L7
- Skordis C., Zlosnik T., 2021, Phys. Rev. Lett., 127, 161302
- Tiley A. L. et al., 2019, MNRAS, 485, 934
- Übler H. et al., 2017, ApJ, 842, 121
- Verlinde E., 2016, SciPost Phys., 2, 016
- Angelov P., 2026, ESTIF v6.2 (software), doi:10.5281/zenodo.17261724

---

## Changelog

**v6.4.2 (July 2026):**
- Aligned with the a₀ horizon doctrine (RHAC-010): title and verbs changed from "derivation" to "construction" for the full coefficient; abstract, §3, §6.1, limitation (i), and the conclusions now state that the horizon scale $cH_0$ is the derived content and the $O(1)$ prefactor ($x_0/\sqrt{3}$) is open. All numerical results unchanged.
- §2.1 and §6.2: $\Omega_m = x_0$ stated as a consistency relation (C2); "predicts" withdrawn.
- Limitation (iii): RHAC-009 recorded — principle P proven non-derivable as a law; $\Omega_m$'s status is "fixed by P."

**v6.4 (July 2026):**
- §2.1 box: LISA 491 µs line made explicitly conditional on C1 (open eddy-stress sector).
- Limitation (v) rewritten: Path Two closed as an honorable null ($w=-1$ exactly); Path One figure corrected to $\chi^2/N = 1.965$ (13-bin pipeline), replacing the task6 constant-$\Lambda$-limit figure of 1.92.
- Limitation (vii): added explicit statement that the +1.17σ/+1.52σ tensions use the bootstrap $a_0$, not the paper's headline value.
- New limitation (viii): Bullet Cluster flagged as an open, unresolved item.

**v6.3 (July 2026):**
- Limitation (ii) rewritten to state precisely that $\mu(a/a_0)$ is not derived from the tilt sector, referencing the numerical verification (`mu_extraction.py`) that closes the local-$n$ route.
- New limitation (vii) added: weak-lensing BTFR consistency check against Mistele et al. 2024 and corrected $a_0$-normalization tension figures (+1.17σ kinematic, +1.52σ lensing), superseding an earlier per-bin systematic-propagation error.
- Note appended to limitation (v): Path One / Path Two split; frozen-eddy limit ties $\Lambda$CDM on DESI DR2.

**v6.2 (March 2026):** Original geometric derivation and SPARC validation.
