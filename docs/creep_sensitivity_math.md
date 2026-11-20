# Creep Model and Sensitivity Comparison: Mathematical Description

## Rest-shape (bend–twist) creep model

The cable is simulated with orientation-based elastic elements (Cosserat-rod-like segment orientations). For a bend–twist constraint between consecutive segment orientations \(q_i, q_{i+1} \in \mathbb{S}^3\) (unit quaternions), the **current relative orientation** is
\[
q^{\mathrm{cur}}_{i,i+1} = q_i^{-1}\, q_{i+1}, \qquad \|q^{\mathrm{cur}}_{i,i+1}\|=1,
\]
with quaternion sign/hemisphere consistency enforced so that \(q^{\mathrm{rest}}\cdot q^{\mathrm{cur}} \ge 0\).

Each such constraint stores a rest Darboux/relative orientation \(q^{\mathrm{rest}}_{i,i+1}\). With creep disabled, \(q^{\mathrm{rest}}\) remains fixed (up to initialization). With creep enabled, the simulator performs two stages controlled by counters, thresholds, and rates:

1. **Activation / rest-state capture.** After `activationDelaySteps` simulation substeps (initial settling under the configured forces, including gravity-only settling in this setup), the current relative orientation is captured as the new rest state:
   \[
   q^{\mathrm{rest}}_{i,i+1} \leftarrow q^{\mathrm{cur}}_{i,i+1}.
   \]
   This gives the deformed, settled configuration a reference state without changing stretch rest lengths.

2. **Post-activation evaluation delay.** An additional `evaluationDelaySteps` (CLI: `--creep-settle-steps`) may elapse before unattended evaluation/auditing begins, allowing the captured/rested state to evolve under the chosen creep law.

3. **Creep update.** After activation, at every `stepInterval` substeps, if the mismatch exceeds `activationThreshold`, the rest orientation is relaxed toward the current orientation:
   \[
   q^{\mathrm{rest}}_{i,i+1} \leftarrow \operatorname{normalize}\!\left((1-\gamma)\,q^{\mathrm{rest}}_{i,i+1} + \gamma\, q^{\mathrm{cur}}_{i,i+1}\right),
   \]
   where \(\gamma \in [0,1]\) is the creep rate (`--creep`, `RestShapeCreepConfig::creep`). The mismatch gate is a squared coefficient difference
   \[
   m = \|q^{\mathrm{rest}}_{i,i+1} - q^{\mathrm{cur}}_{i,i+1}\|^2 > \tau,
   \]
   with threshold \(\tau =\) `activationThreshold`.

Special cases: \(\gamma=0\) means no post-capture creep relaxation (fixed captured rest state); larger \(\gamma\) means faster relaxation of bend–twist rest shape toward the instantaneous deformed shape (reduced long-term elastic memory / resistance to sustained deformation, within this model).

All variants in sensitivity runs use the same lowered activation threshold (e.g. \(10^{-4}\)) so that capture/creep can activate under gravity-only settling; with a too-large threshold the mismatch may never exceed \(\tau\) and nominal/high-\(\gamma\) variants can become indistinguishable from the no-creep baseline.

## Sensitivity comparison methodology

The goal is to quantify how conclusions depend on \(\gamma\) (and, separately, on tangential contact response), not to claim a uniquely calibrated material parameter.

### Creep sweep

Compare at least the nominal and bounding/stronger settings on the same seed sets/catalog slices and with identical solver/audit settings except \(\gamma\):
- \(\gamma = 0\) (no creep / captured-fixed baseline in effect),
- \(\gamma = 0.02\) (nominal),
- \(\gamma = 0.05\) (moderate/stronger, physically reasonable within the same phenomenological model),
- optionally \(\gamma = 0.2\) (high, included as an upper-stress sensitivity bound where computationally feasible).

For each condition, run the relevant pipeline/audit mode (topology+residual-penetration audit on a selected seed subset and/or full catalog phases as available) and record:

- topology outcomes: total/audited seeds and phases, steps, violation counts/types, seeds with violations, centerline/surface-contact counts, first-violation metadata;
- contact/solver outcomes: candidate and penetrating contact counts, pre-correction \(\max/\)mean penetration, residual candidate/penetrating counts and \(\max/\)mean residual penetration, max ground penetration, capped-correction and max-applied-correction statistics, normal-flip counts;
- task/geometry outcomes when available from the dataset pipeline (e.g. resolved rates, cuts-to-zero/heurstic-oracle comparisons, configuration-level consistency) to assess stability of reported conclusions.

### Comparison logic

Sensitivity is reported as condition-wise differences/ratios relative to nominal (\(\gamma=0.02\)) and to baseline (\(\gamma=0\)): changes in violation presence/severity, shifts in residual penetration distributions/extrema, and any qualitative change in strategy/outcome conclusions. Conclusions stated as robust should remain directionally and materially stable across \(\gamma \in \{0,0.02,0.05\}\) (and, when reported, under the \(\gamma=0.2\) stress bound); any conclusion that reverses or depends on a particular \(\gamma\) must be scoped accordingly.

### Tangential contact response (damping vs Coulomb friction)

The current tangential interaction is modeled as **velocity-proportional damping** opposing relative tangential motion at contacts, not Coulomb (rate-independent, stick–slip, normal-load-proportional) friction. Let \(\mathbf{v}_t\) be tangential relative velocity and \(b_t \ge 0\) the tangential damping coefficient (`--tangential-damping`); the damping response scales with \(\mathbf{v}_t\) (and the implemented contact weighting), so it regularizes sliding/jitter but cannot represent static friction, sticking thresholds, or load-dependent friction cones.

Sensitivity to this modeling choice is assessed by sweeping \(b_t\) over a small set of values (including reduced/nominal/increased damping) while holding other parameters fixed, and comparing the same topology/penetration/outcome metrics. Physical interpretations that would rely on true Coulomb friction, stick–slip transitions, or friction-limited rearrangements are explicitly out of scope for the present results except insofar as they remain stable under these damping variations.
