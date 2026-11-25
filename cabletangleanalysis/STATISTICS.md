# Statistical methods in CableTangleAnalysis

This note documents the statistics used by the Julia analysis. It is kept
outside version control so the numerical assumptions and test choices remain
easy to inspect while the analysis is being developed.

## Unit of analysis and pairing

The unit is one simulated tangle, identified by `seed`. Every strategy is
evaluated on the same seed, so strategy comparisons are **paired by seed**.
For random strategies, the per-seed value is the mean over the random trials
(`--eval-random-trials`, `--ctz-random-trials`).

For a strategy value `x_i` and the heuristic value `h_i` on seed `i`, the
paired difference is

```text
d_i = x_i - h_i
```

The direction that counts as better depends on the metric: for cuts needed,
negative `d_i` is better; for threading reduction and complexity reduction,
positive `d_i` is better.

Because many seeds can have exactly the same value under two strategies, the
paired differences often contain a large mass at zero. The **median of all
paired differences** can therefore be zero even when the two marginal
distributions differ visibly (for example, heuristic median 3 cuts vs. oracle
median 2 cuts). The reports therefore also include:

- `median_nonzero_diff`: median of the paired differences after excluding ties
- `n_improved`, `n_tied`, `n_worse`: outcome counts relative to the heuristic
- `improved_rate`, `tied_rate`, `worse_rate`: the corresponding fractions

The `better` label is derived from the direction of the effect — the sign of
the matched-pairs rank-biserial effect size, with the sign-test majority and
the mean difference as fallbacks — rather than from the sign of `median_diff`.
This avoids reporting “tie” when a strategy improves a large minority of seeds
while just over half are unchanged.

## Threading outcomes

- **Cuts to zero**: number of cuts needed to remove all threading.
- **Single-cut threading reduction**: `T_before - T_after` after one cut.

### Oracle scope

The complete cuts-to-zero oracle (`oracle_cuts`) reports the shortest sequence
that reaches zero threading using only **hoop-midpoint cuts** — the same action
set as the random-hoop, biggest-hoop, and smallest-hoop strategies. It may use
more cuts than a heuristic that is free to cut any edge, and it is not a lower
bound over all possible cable cuts. Optimality within the hoop-midpoint action
set is certified by `oracle_exhaustive = 1`, which means every shorter depth was
exhausted before the winning sequence was found. Intermediate threading is
allowed to increase; only the final state is constrained to be zero.

Oracle stage rows contain the replayed winning path only, so the resolved
violin uses one trajectory per seed rather than one row per expanded search
node. When comparing the oracle with the heuristic, negative differences still
mean the oracle used fewer cuts, but a positive difference is possible and
reflects the narrower action set rather than a search failure.

Threading reduction and cuts needed are also shown as violin distributions in
two variants each:

- **first cut**: threading removed by stage 1 relative to stage 0.
- **resolved**: the final cuts-to-zero stage, restricted to trials where
  threading actually reached zero.

For cuts needed, one variant includes every trial (unresolved trials count the
cuts that were attempted) and one is restricted to resolved trials. Tangles
whose initial state has no threading are excluded from the plots and paired
tests, since they cannot need a cut. For cuts needed, negative differences
mean the strategy used fewer cuts. Significance markers are described in the
next section.

## Significance markers in the cut-evaluation plots

The four cut violins — `threading_reduction_first_cut.pdf`,
`threading_reduction_resolved.pdf`, `cuts_needed_all.pdf`, and
`cuts_needed_resolved.pdf` — compare each strategy with the heuristic using
the paired Wilcoxon signed-rank test. Each has a split-violin companion
(`*_paired.pdf`) that shows the paired differences themselves: the left half
contains every `strategy − heuristic` difference, including exact zeros, and
the right half contains only the non-zero differences that enter the
conditional test. The tie and discordant counts printed above each strategy
are descriptive; the asterisks still come from the conditional Wilcoxon test
described below. The procedure is:

1. **Paired differences.** For each seed present in both series, compute
   `d = strategy − heuristic`. Seeds where either value is non-finite are
   dropped; a row is tested only when at least two such seeds remain.
2. **Zero handling.** Differences with `|d| ≤ 1e-12` are treated as exact
   ties and removed before the Wilcoxon test. They are still counted in the
   improved/tied/worse rates and in the reported medians. This is deliberate:
   the signed-rank and sign tests test symmetry about zero, or equivalently
   equal probability of a positive and a negative discordant difference. Under
   that null a zero pair is unchanged by a sign flip, so it carries no
   directional information. The p-value is therefore conditional on the
   observed number of discordant pairs; adding or removing ties does not change
   the evidence about direction.
3. **Wilcoxon signed-rank test.** Two-sided test of zero median on the
   non-zero differences. `HypothesisTests.SignedRankTest` uses the exact null
   distribution when there are at most 15 non-zero differences, or at most 50
   with no tied absolute ranks; otherwise it uses the normal approximation
   with tie correction and a 0.5 continuity correction. When every difference
   is zero the test is skipped and `p = 1`.
4. **Multiple testing.** Benjamini-Hochberg adjustment is applied within each
   metric, across all strategies compared with the heuristic for that metric
   (one family per violin/table; the heuristic itself and rows with fewer than
   two paired seeds or a non-finite p-value are excluded from the family). The
   conditional and unconditional p-values are adjusted in separate families.
5. **Marker threshold.** A strategy gets a marker only if
   `p_adj < 0.05`; otherwise it is drawn without a marker but remains in the
   paired-test TSV. The boundary values are strict, so `p_adj = 0.05` is not
   marked.

The asterisk count is:

| adjusted p-value | marker |
|---|---|
| `p_adj < 0.001` | `***` |
| `0.001 ≤ p_adj < 0.01` | `**` |
| `0.01 ≤ p_adj < 0.05` | `*` |
| `p_adj ≥ 0.05` or non-finite | none |

Every marked violin also prints its effect size above the asterisk as
`r = ...`, where `r` is the matched-pairs rank-biserial
`(W+ − W−)/(W+ + W−)` computed on the average ranks of the absolute non-zero
differences. The sign follows `strategy − heuristic`, so positive `r` means
the strategy is better for metrics where higher is better, and worse for
metrics where lower is better (cuts needed). `r` is conditional on the
discordant pairs, as described under
[Conditional effect versus overall effect](#conditional-effect-versus-overall-effect).

Marker colour encodes direction only:

- **red** = the strategy is better than the heuristic
- **green** = the heuristic is better than the strategy
- **gray** = no direction (tie); this case is not expected in the current data

Direction comes from the `better` label, which is evaluated in this order:

1. the sign of the matched-pairs rank-biserial effect size, computed as
   `(W+ − W−)/(W+ + W−)` on the average ranks of the absolute non-zero
   differences;
2. if the effect size is non-finite or its absolute value is `≤ 1e-12`, the
   sign-test majority: more positive than negative differences (also with a
   `1e-12` tolerance);
3. if still tied, the sign of the mean paired difference;
4. if still tied, the sign of the median of the non-zero paired differences;
5. otherwise the result is a tie.

This is why the median of all paired differences is not used for direction:
when many seeds are exact ties, that median can be zero even though the
effect size and every non-tied seed favour one strategy. The asterisks
therefore report significance and the colour reports direction.

### Conditional effect versus overall effect

The conditional test deliberately ignores ties, so the rank-biserial effect
size is conditional as well: `effect_size = −1` means **every discordant pair**
favours the strategy, not that the strategy wins on every tangle. Every paired
table therefore also reports unconditional summaries:

- `n_discordant = n_improved + n_worse`
- `probability_of_superiority = (n_improved + 0.5 × n_tied) / n`, where ties
  count as half a win (`0.5` = no effect, `1` = the strategy is always better)
- `net_improvement_rate = (n_improved − n_worse) / n`

As a sensitivity check, each table also contains `unconditional_p_value` and
`unconditional_p_adjusted`: a one-sided exact binomial test of
`H0: P(strategy better) ≤ 0.5` over all `n` pairs, with ties counted as
failures. This answers the stricter question “does the strategy improve more
than half of all cases?” and is Benjamini-Hochberg adjusted separately within
the same family. It is **not** used for the plot markers. The markers keep the
conditional signed-rank result because that is the test for a directional
difference among pairs where the strategies actually differ. When the two
disagree, report both and state which question is being answered.

For example, in the cuts-needed comparison Oracle has
`322 improved / 421 tied / 0 worse` out of `n = 743`:

- conditional Wilcoxon: `p_adj = 1.7e-60`, effect size `−1.0`
- `probability_of_superiority = 0.717`, `net_improvement_rate = 0.433`
- unconditional majority test: `p_adj = 1.0`

The first statement is “among the 322 tangles where Oracle and the heuristic
differ, Oracle always uses fewer cuts”; the second is “Oracle is better in
43.3% of all tangles and worse in none”. Both are true, and the table now
contains the information needed to make either statement.

For completeness, the TSV also reports a 95% percentile bootstrap confidence
interval for the median paired difference (10,000 resamples, base seed
`20240505` plus the strategy index, 2.5th/97.5th percentiles) and a
randomization symmetry diagnostic (1,000 sign flips, same seed scheme,
two-sample KS distance), which is Benjamini-Hochberg adjusted separately and
summarized as `symmetry_ok` when its adjusted p-value is `≥ 0.05`. Neither
affects the markers.

The complexity-reduction violins use the same test, adjustment, thresholds,
and colour convention.

## Complexity outcomes

The complexity metrics are computed at every cutting stage. The complexity
score is a PCA projection of the five standardized metrics:

```text
segmentCount, selfCollisionCount, closePairCount,
totalCurvature, outsideKnotCurvature
```

Two scores are reported:

- **positive PC1**: absolute PC1 loadings, so every metric increases complexity.
- **signed PC1**: the actual signed first-principal-component projection.

Complexity reduction is `complexity_before - complexity_after`, so positive
means less complexity. The complexity-reduction figures show the signed PC1
score only. Strategies whose paired Wilcoxon comparison with the heuristic is
significant after Benjamini-Hochberg correction are marked with red asterisks
when the strategy is better and green asterisks when the heuristic is better.
Each marked violin also prints its matched-pairs rank-biserial effect size as
`r = ...` (strategy − Heuristic, conditional on the discordant pairs).

Two cut outcomes are compared:

- **first cut**: stage 1 relative to stage 0.
- **resolved**: the final cuts-to-zero stage, restricted to trials where
  threading actually reached zero.

## Paired Wilcoxon signed-rank test

The main strategy comparison uses the two-sided Wilcoxon signed-rank test on
the paired differences. Zero differences are discarded.

The test assumes independent pairs and, for the usual median-shift
interpretation, that the paired differences are symmetrically distributed
about their median. It does **not** assume normality. The sample sizes here
are large (`n` is typically several hundred), so the asymptotic approximation
is reasonable, but symmetry is still checked.

## Symmetry diagnostic

For every paired comparison, a randomization test checks whether the paired
differences are symmetric about their median:

1. Center the differences: `c_i = d_i - median(d)`.
2. Use the two-sample Kolmogorov-Smirnov distance between `c` and `-c` as the
   test statistic.
3. Build the null distribution by randomly flipping the sign of each `c_i`
   (1000 permutations).
4. Report `p = (1 + #{statistic* >= statistic_obs}) / (1 + B)`.

The p-values are Benjamini-Hochberg adjusted within each metric family (all
strategy-vs-heuristic comparisons for one outcome). The column `symmetry_ok`
is false when the adjusted p-value is below 0.05. This is a diagnostic, not a
prerequisite for reporting the Wilcoxon result.

## Sign test

A distribution-free sign test is reported as a robustness check. It uses the
number of positive and negative paired differences (zeros are discarded) and
tests the null that the median difference is zero. It does not require
symmetry but is generally less powerful than the Wilcoxon test.

## Confidence interval for the median difference

The 95% confidence interval is a percentile bootstrap interval:

- resample the paired differences with replacement,
- compute the median of each resample (10,000 resamples),
- report the 2.5% and 97.5% quantiles.

This interval does not assume normality or symmetry.

## Effect size

The matched-pairs rank-biserial correlation is

```text
(W+ - W-) / (W+ + W-)
```

where `W+` and `W-` are the summed ranks of positive and negative paired
differences. It ranges from -1 to 1; the sign follows the median difference.

## Multiple testing: why Benjamini–Hochberg?

Each cut-evaluation figure compares several strategies against the same
reference strategy (the heuristic). A figure with `m` non-reference strategies
therefore performs `m` two-sided tests on the same outcome. Without any
adjustment, the probability of at least one false positive under a global null
is approximately `1 - (1 - 0.05)^m` (about 23% for `m = 5` and 26% for
`m = 6`, assuming independence). Reporting unadjusted stars would then
overstate how many strategies genuinely differ from the heuristic, especially
for borderline comparisons.

We control the **false discovery rate (FDR)** rather than the family-wise
error rate (FWER) with the Benjamini–Hochberg (BH) procedure:

1. Sort the `m` p-values: `p_(1) <= ... <= p_(m)`.
2. Find the largest `k` with `p_(k) <= (k / m) * 0.05`.
3. Call the hypotheses with the `k` smallest p-values significant.
   Equivalently, the adjusted p-value is
   `p_adj_(k) = min_{j >= k} (m / j * p_(j))`.

BH is used here because:

- It controls the expected proportion of false discoveries among the
  comparisons that are called significant, which matches how these figures are
  read: "which strategies look different?" rather than "is *any* strategy
  different?".
- It is less conservative than Bonferroni when several effects are real, which
  matters because a single figure can contain two to four genuine differences.
- The tests within a figure share the reference strategy and the same seeds,
  so their test statistics are positively correlated. BH keeps its FDR
  guarantee under positive regression dependence (Benjamini & Yekutieli,
  2001), which is the plausible dependence structure here.

Caveats:

- BH does not guarantee FDR control under arbitrary dependence. For a
  distribution-free guarantee the Benjamini–Yekutieli (BY) procedure would be
  required. BY is more conservative and would only strengthen the current
  conclusions, whose adjusted p-values are already many orders of magnitude
  below 0.05.
- Adjustment is performed **within one outcome metric**, i.e. within one
  figure. The six figures (signed PC1 first cut, signed PC1 resolved,
  threading first cut, threading resolved, cuts needed for all trials, cuts
  needed for resolved trials) are separate families because they answer
  different questions. They are not pooled into one global family. A global BH
  correction across all strategy comparisons would be more conservative; it
  would not change the current stars because every starred comparison has
  `p_adjusted << 0.05`.
- Adjustment is applied to the Wilcoxon p-values before the star levels are
  assigned: `*` means `p_adjusted < 0.05`, `**` means `< 0.01`, and `***`
  means `< 0.001`.
- The symmetry diagnostic and the sign test are also BH-adjusted within their
  own families. Adjusting the symmetry p-values makes the precondition check
  more conservative: more comparisons can be classified as
  `symmetry_ok = false` and therefore receive a black rather than a red star.
- Multiplicity correction cannot compensate for dependence between different
  outcome metrics, small samples, or a comparison chosen after looking at the
  data. Those issues must be handled by design (pre-specified strategies and
  outcomes).
- If the comparisons were strictly confirmatory and pre-specified one at a
  time, reporting unadjusted p-values together with effect sizes could be
  defensible. These figures scan all available strategies, so adjustment is
  the safer default.

The adjusted values are written as `p_adjusted`, `symmetry_p_adjusted`, and
`sign_p_adjusted`. In the current dataset the correction does not change any
conclusion: all starred comparisons have adjusted p-values far below 0.05.

### Why not Bonferroni?

Bonferroni controls the FWER (the probability of *any* false positive) at
`0.05 / m`. It is valid under arbitrary dependence but very conservative when
tests are correlated and several effects are real. With `m = 6`, Bonferroni
requires `p < 0.0083` for every discovery, whereas BH uses a data-driven
threshold that approaches `0.05` when many hypotheses are truly non-null.
Because the analysis screens several strategies instead of making a single
confirmatory claim, FDR control is the more appropriate and more powerful
target.

## PCA

The PCA uses the five metrics above. They are standardized with the population
standard deviation before the SVD. The loadings, eigenvalues, and explained
variance are written to `pca_results.tsv`; per-seed PC scores are written to
`pca_scores.tsv`.

The component signs are arbitrary. The complexity score fixes the PC1 sign so
the summed loadings are positive, then either keeps the signed projection or
uses the absolute loadings (positive orientation). The score is normalized to
[0, 1] over the catalog rows and also reported as a percentile.

A separate model is fitted per creep variant. When multiple variants are
analyzed together, a pooled model is also fitted.

## Descriptive statistics

- The cuts histogram and the violin plots print per-strategy means and
  medians; the paired-test tables carry the complete descriptive summary.
- The knot-type distribution includes every catalog row. Labels such as
  `none`, `unknown`, `??`, and empty values are grouped as `Undefined`.
- Tangles with an undefined knot type are excluded from the metric,
  correlation, PCA, complexity, and cut analyses because their outside-knot
  geometry is not defined. The knot-type distribution is the only output that
  keeps them.

## Plot filtering

Some cut plots show only the heuristic and strategies with a statistically
significant Wilcoxon difference from the heuristic after Benjamini-Hochberg
correction. All strategies remain in the TSV reports.

## Reproducibility

- Bootstrap seed: `20240505 + index`.
- Symmetry permutation seed: `20240505 + index`.
- Bootstrap repetitions: 10,000.
- Symmetry permutations: 1,000.

## Software

The analysis uses Julia with:

```text
HypothesisTests, MultipleTesting, MultivariateStats,
Distributions, StatsBase, Polynomials, CairoMakie
```

## LaTeX tables

Generate the paper tables with:

```bash
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/make_latex_tables.jl \
  results/dataset_final/analysis_julia/single \
  results/dataset_final/latex
```

The script writes:

```text
table_pca.tex
table_complexity_first_cut.tex
table_complexity_resolved.tex
tables.tex
```
