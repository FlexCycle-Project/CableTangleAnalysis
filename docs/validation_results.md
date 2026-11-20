# Numerical/topology validation and sensitivity results

This document records the executed audits and sensitivity comparisons produced with the diagnostics-enabled `MyCableDemo` binary (extended topology-audit outputs with pre-correction penetration, post-correction/residual penetration, capped-correction, normal-flip, and ground-penetration diagnostics; CLI controls `--dt`, `--collision-iterations`, `--collision-compliance`, `--collision-max-correction`, `--tangential-damping`, `--normal-damping`, `--hash-cell-size`).

**Scope.** All topology/penetration audits are **discrete-substep observations**. The simulation uses finite $\Delta t$, finite self-collision iterations per substep, compliance-softened and capped/averaged (Jacobi-style) corrections, and finite contact passes. Reported absence of violations or low residual penetration is evidence at audited states only, not a continuous-time or continuous-collision no-tunneling/no-penetration guarantee. See `docs/topology_audit_math.md`.

**Creep model.** Bend–twist rest orientations are captured after settling and optionally relaxed as $q^{\mathrm{rest}}\leftarrow \operatorname{normalize}((1-\gamma)q^{\mathrm{rest}}+\gamma q^{\mathrm{cur}})$ at activated intervals when mismatch exceeds threshold; `--creep-settle-steps` controls post-capture evaluation delay. See `docs/creep_sensitivity_math.md`.

**Tangential contact.** The implemented tangential response is velocity-proportional damping (coefficient $b_t$, `--tangential-damping`), not Coulomb friction (no stick–slip, no normal-load-dependent friction cone). Physical interpretations requiring Coulomb/stick–slip behavior are out of scope except where stable under damping sweeps.

## 1. Baseline (existing, legacy schema)

From `results/full_dataset_20261106_112203/topology_audit_summary.txt` (current/nominal creep, 1000 seeds, 300 audited substeps/phase, 3 heuristic cuts, contact $\epsilon=10^{-4}$; no penetration-diagnostic columns):

| seeds | phases | steps | violations | seeds w/ violations | centerline contacts | surface overlaps |
|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 2977 | 893100 | 2 | 2 | 0 | 9369048 |

Existing creep-variant pipeline outcomes (`results/full_dataset_20261106_112203/{zero,current,high}/summary.txt`, 1000 seeds each, settle 300, threshold $10^{-4}$):

| Variant | $\gamma$ | heuristic resolved | heuristic cuts mean | $\Delta T$ heuristic mean |
|---|---:|---:|---:|---:|
| zero | 0.00 | 1.000 | 2.307 | 12.606 |
| current | 0.02 | 1.000 | 2.306 | 12.575 |
| high | 0.20 | 1.000 | 2.306 | 12.624 |

These aggregate task outcomes are nearly invariant across $\gamma\in\{0,0.02,0.2\}$, indicating low sensitivity of the reported heuristic/cuts-to-zero performance summaries to this creep range in that dataset (modest evidence; not a full distributional equivalence proof).

## 2. Diagnostics-enabled nominal and spot-check audits

Common settings unless noted: catalog `results/full_dataset_20261106_112203/current/catalog.tsv`, 300 steps/phase, 3 heuristic cuts, $\epsilon_c=10^{-4}$, settle 300, threshold $10^{-4}$, $\Delta t=0.005$, 8 self-collision iterations, compliance $10^{-10}$, max correction $0.08$, normal damping $0.15$, tangential damping $b_t=0.05$, hash cell $0.2$, $\gamma=0.02$. “Residual” penetration is measured after the final self-collision iteration of each substep, before the ground clamp.

| Run | Seeds | Phases | Steps | Viol. | Seeds w/ viol. | Max pen. | Resid. max | Resid. mean* | Max ground | Capped | Max corr. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Nominal (diagnostics, 100-seed) | 100 | 287 | 86100 | 0 | 0 | 0.1534 | 0.3351 | 0.01770 | 0.0662 | 1228 | 0.0800 |
| $\Delta t=0.0025$ spot | 20 | 57 | 17100 | 0 | 0 | 0.0893 | 0.0424 | 0.002664 | 0.0291 | 0 | 0.0447 |
| Iterations $=16$ spot | 20 | 57 | 17100 | 0 | 0 | 0.0769 | 0.00694 | 0.000234 | 0.0485 | 0 | 0.0468 |

\*`residual_mean_penetration` is the run-aggregate mean reported by the auditor over residual penetrating contacts.

**Nominal status:** `results/validation_20261118/nominal` is now fully merged (24/24 shards: 100 seeds, 287 phase rows, 86100 audited substeps, 0 event rows) in `topology_audit_{summary,seeds,events}.tsv`. The 6 previously stalled shards (`part04,11,12,13,16,19`) were rerun and incorporated; earlier partial artifacts remain as `topology_audit_{summary,seeds,events}.partial.tsv/.txt`.

**Observations:** No topology/crossing violations or centerline-contact events were observed in the completed diagnostic runs above. However, nonzero **residual** self-collision penetrations (and, in nominal/creep sweeps, nonzero capped correction counts and corrections reaching the cap) were observed at audited states, indicating that the finite, compliance/cap-limited solver does not eliminate all interpenetration after the final iteration of a substep under this finite, compliance/cap-limited solver configuration. The $16$-iteration spot check substantially reduced residual penetration relative to the 8-iteration configurations (residual max $0.00694$, residual mean $2.34\times 10^{-4}$, no caps), consistent with improved per-substep convergence. The smaller-$\Delta t$ spot check showed lower pre/post residual maxima than the 20-seed $\gamma=0.02$ reference in this subset but a larger residual max than the 16-iteration run; these are subset-specific spot checks, not exhaustive trends. The full nominal 100-seed run exhibits a larger residual maximum ($0.3351$, residual mean $0.01770$) than the 20-seed $\gamma=0.02$ subset, underscoring that penetration residuals are configuration- and case-dependent and that subset results should not be extrapolated as global bounds.

## 3. Creep sensitivity audits (20-seed subsets)

Same settings as nominal diagnostics except $\gamma$; 20 seeds, 58 phases, 17400 steps each. No topology violations observed.

| $\gamma$ | Max pen. | Resid. max | Resid. mean | Max ground | Capped | Max corr. |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.1440 | 0.1136 | 0.006655 | 0.0695 | 394 | 0.0800 |
| 0.02 | 0.1534 | 0.1157 | 0.006332 | 0.0662 | 714 | 0.0800 |
| 0.05 | 0.1536 | 0.0943 | 0.005442 | 0.0708 | 282 | 0.0800 |
| 0.20 | 0.1230 | 0.0482 | 0.004137 | 0.0619 | 6 | 0.0800 |

**Interpretation:** Over this 20-seed audited subset, residual penetration and capped-correction burden vary with $\gamma$ (lowest residual max/mean and fewest caps at $\gamma=0.20$ here), while crossing-topology violations remain absent in all tested conditions. Together with the near-invariant 1000-seed pipeline task metrics for $\gamma\in\{0,0.02,0.2\}$ (Section 1), this suggests the headline task outcomes are not strongly driven by creep in the tested range, but contact-tightness/residual-penetration diagnostics are **not** invariant and can be materially affected by $\gamma$ and by solver limits (caps). Claims depending on precise contact geometry or penetration-free states must be qualified accordingly.

## 4. Tangential-damping sensitivity audits

Subset audits (20 seeds, nominal $\Delta t$/iterations, $\gamma=0.02$) compare $b_t\in\{0.0,0.05,0.2\}$ as **velocity-proportional damping** sensitivity only.

| $b_t$ | Seeds | Phases | Viol. | Max pen. | Resid. max | Resid. mean | Max ground | Capped | Max corr. |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 20 | 56 | 0 | 0.1182 | 0.1006 | 0.007685 | 0.0728 | 15 | 0.0800 |
| 0.05 (nominal, = $\gamma=0.02$ creep run) | 20 | 58 | 0 | 0.1534 | 0.1157 | 0.006332 | 0.0662 | 714 | 0.0800 |
| 0.20 | 20 | 58 | 0 | 0.1058 | 0.0635 | 0.005212 | 0.0619 | 0 | 0.0681 |

$b_t=0.00$ completed with 56 audited phase rows (16800 steps) versus 58 for the other conditions in this subset; values are reported as observed.

**Scope note:** $b_t$ controls viscous tangential damping, not Coulomb friction. Therefore these runs bound sensitivity to tangential damping magnitude but do **not** validate behavior under static friction, stick–slip, or load-dependent friction; conclusions requiring those mechanisms remain out of scope. No topology violations were observed across the tested damping values, while residual/capped-correction diagnostics varied (e.g. lower residual max/mean and no caps at $b_t=0.20$ in this subset), so contact-quality conclusions remain conditioned on the damping model and solver settings.

## 5. TeX tables

`tables/topology_audit_tables.tex` contains two summary tables corresponding to Sections 2–4 (audit/spot checks; creep and tangential-damping sensitivity) with the currently observed values filled in, including the fully merged 100-seed nominal diagnostics run.

## 6. Reproducibility artifacts

- Binary: `bin/MyCableDemo` (rebuilt with diagnostics/topology-audit hooks)
- Runner: `scripts/run_topology_audit.sh` (supports `--dt`, `--collision-iterations`, `--collision-compliance`, `--collision-max-correction`, `--tangential-damping`, `--normal-damping`, `--hash-cell-size`, `--max-seeds`, correct max/sum/config summary merge)
- Outputs: `results/validation_20261118/{smoke,smoke2,nominal,dt_spot,iter_spot,creep_000,creep_005,creep_020,creep_200,damp_000,damp_020}/...` (nominal 100-seed audit fully merged; prior nominal partial artifacts retained as `*.partial.*` where present)
- Math docs: `docs/topology_audit_math.md`, `docs/creep_sensitivity_math.md`
