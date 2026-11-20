# Topology and Penetration Audit: Mathematical Description

## Scope and interpretation

This audit is an **observed discrete-substep** check of the simulated cable state. The solver advances with finite time steps, a finite number of contact/collision correction passes, compliance-softened responses, and capped/averaged (Jacobi-style) corrections. Therefore the audit provides evidence that no audited crossing violation or residual penetration above the reported thresholds was observed at the audited states; it is **not** a formal continuous-time or continuous-collision no-tunneling/no-penetration guarantee.

The audit is run after the initial settle/creep-activation phase and, when cuts are requested, after each heuristic cut. A cut is an intentional topology change, so the auditor re-baselines immediately on the rebuilt geometry after the cut and only judges subsequent motion.

## 2.5-D cable model and XZ crossing diagram

Gravity acts along \(-Y\), and the cable rests on/above the ground plane. For a pair of cable centerline points with the same planar coordinates \(x,z\), the vertical coordinate \(y\) determines which strand is above the other. This makes the over/under structure observable in the planar projection onto \(XZ\), with \(Y\) as the height/overpass coordinate.

Let each cable piece be a polyline \(P^{(p)} = (p^{(p)}_0,\ldots,p^{(p)}_{n_p})\), \(p^{(p)}_i=(x_i,y_i,z_i)\in\mathbb{R}^3\). For a non-adjacent edge pair
\[
e_A = [p^{(a)}_i,p^{(a)}_{i+1}],\quad e_B=[p^{(b)}_j,p^{(b)}_{j+1}],
\]
project to \(XZ\): \(\pi(p)=(x,z)\in\mathbb{R}^2\). The projected segments \([\pi(p^{(a)}_i),\pi(p^{(a)}_{i+1})]\) and \([\pi(p^{(b)}_j),\pi(p^{(b)}_{j+1})]\) cross if there exist parameters \(u,v\in(0,1)\) (excluding endpoints to avoid shared-particle double counting and degenerate adjacency) such that
\[
\pi(p^{(a)}_i)+u\bigl(\pi(p^{(a)}_{i+1})-\pi(p^{(a)}_i)\bigr)=
\pi(p^{(b)}_j)+v\bigl(\pi(p^{(b)}_{j+1})-\pi(p^{(b)}_j)\bigr).
\]

At that projected crossing point, define the interpolated heights
\[
Y_A(u)=y^{(a)}_i + u\,(y^{(a)}_{i+1}-y^{(a)}_i),\qquad
Y_B(v)=y^{(b)}_j + v\,(y^{(b)}_{j+1}-y^{(b)}_j),
\]
and the signed height separation (over/under indicator)
\[
\Delta h = Y_A(u)-Y_B(v).
\]
Because the points coincide in \(x,z\), the 3D centerline distance at the crossing is
\[
d_c = \lVert (x,z,Y_A)-(x,z,Y_B)\rVert = |\Delta h|.
\]

A crossing is recorded with its edge identity, projected location, and \(\Delta h\). Adjacent edges of the same piece are excluded.

## Temporal matching and topology violations

Let \(C_k\) be the set of crossings at audited substep \(k\), keyed by the non-adjacent edge pair (and piece indices). Crossings present in both \(C_{k-1}\) and \(C_k\) are **persistent**; crossings only in \(C_k\) are **born**; crossings only in \(C_{k-1}\) are **dead**.

With centerline contact threshold \(\varepsilon_c\) (`--topology-audit-contact-eps`), the auditor classifies:

- **`CROSSING_FLIP`** (topology violation): a persistent crossing changes over/under sign, \(\operatorname{sign}(\Delta h_{k-1}) \neq \operatorname{sign}(\Delta h_k)\), with the usual convention that a sign change through/near zero indicates strands passing through each other in the audited transition.
- **`BIRTH_AT_CONTACT`**: a crossing appears while \(|\Delta h_k| \le \varepsilon_c\), i.e. born through centerline contact rather than a separated planar rearrangement.
- **`DEATH_AT_CONTACT`**: a crossing disappears while \(|\Delta h_{k-1}| \le \varepsilon_c\).
- **`CENTERLINE_CONTACT`**: a crossing remains within \(\varepsilon_c\), \(|\Delta h_k| \le \varepsilon_c\), indicating touching or centerline-interpenetrating strands at that crossing.

Crossings that appear/disappear with \(|\Delta h|\) well above \(\varepsilon_c\) are ordinary projection rearrangements (plan-view Reidemeister-type moves) and are not reported as violations. For a ground-resting tangle of circular cross-section with physical rope radius \(r\), two resting strands in contact have centerline separation \(\approx 2r\); \(\varepsilon_c\) is chosen far below \(2r\) so that normal resting contact does not trigger a centerline-contact/topology event.

## Self-collision penetration and residual diagnostics

Self-collision is handled on non-adjacent cable segments by closest-point queries and corrective responses over a finite number of Jacobi-style iterations per substep, with compliance, damping, and a maximum correction cap.

For a candidate segment pair, let \(d\) be the closest-point distance between segment centerlines and \(r\) the collision radius (configured rope/collision radius). Define nominal penetration
\[
p = 2r - d.
\]
A candidate is **penetrating** when \(p > 0\). The solver records, per iteration and accumulated per substep:

- candidate pair count and penetrating contact count,
- pre-correction penetration statistics \(\max(p),\, \sum p\) over penetrating candidates,
- counts of capped/applied endpoint corrections and \(\max\) applied correction magnitude,
- normal-flip/side-switch indicators when the contact normal relationship is inconsistent with a stable separating correction.

After the final collision iteration of a substep (and before the ground clamp), a **residual** closest-point pass measures post-correction state: residual candidate count, residual penetrating contact count, and residual penetration statistics \(\max(p_{\mathrm{res}}),\, \sum p_{\mathrm{res}}\). Ground interaction reports pre-clamp maximum penetration below the ground plane \(y=y_g\) as \(\max(0, y_g-y)\) over constrained particles.

These diagnostics distinguish *detected and corrected* penetrations from *remaining/residual* penetrations at audited states. A zero (or negligible, relative to scale/tolerances) residual is evidence against observed pass-through at audited states, but finite \(\Delta t\), finite iterations, compliance, and correction caps mean it does not exclude tunneling between audited states.

## Outputs

The audit produces per-phase seed rows (`*_seeds.tsv`), per-violation event rows (`*_events.tsv`), and an aggregate summary (`*_summary.txt`). Key fields include step/phase counts, crossing/violation/contact totals, first-violation metadata, candidate/penetrating/residual contact counts, max/mean penetration and residual penetration, max applied correction and cap/normal-flip counts, max ground penetration, and the solver/configuration parameters (\(\Delta t\), collision iterations, compliance, max correction, normal/tangential damping, hash cell size, creep) used for the run.
