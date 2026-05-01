# biomAid

[![R-CMD-check](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml)
[![Codecov](https://app.codecov.io/gh/DrJ001/biomAid/graph/badge.svg)](https://app.codecov.io/gh/DrJ001/biomAid)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)

Biometric tools to support the analysis of comparative experiment data fitted
with **ASReml-R V4** mixed models. Watch this space, there are a lot more functions coming.

---

## Vignettes

| Vignette | Description |
|----------|-------------|
| [Wald Tests on Fixed-Effect Contrasts](https://DrJ001.github.io/biomAid/waldTest.html) | Mathematical framework and worked examples for `waldTest()` and `plot_waldTest()` |
| [Multi-treatment Random Regression](https://DrJ001.github.io/biomAid/randomRegress.html) | Conditioning schemes, efficiency/responsiveness decomposition, and all plot types for `randomRegress()` and `plot_randomRegress()` |
| [Multi-treatment Fixed-Effects Regression](https://DrJ001.github.io/biomAid/fixedRegress.html) | OLS conditioning schemes, efficiency/response index decomposition, and plot types for `fixedRegress()` and `plot_fixedRegress()` |
| [Extracting and Padding Field Trial Layouts](https://DrJ001.github.io/biomAid/padTrial.html) | Step-by-step guide to guard-row removal, missing-plot padding, and Before/After visualisation with `padTrial()` and `plot_padTrial()` |
| [Multiple Comparison Criteria](https://DrJ001.github.io/biomAid/compare.html) | HSD, LSD, and Bonferroni criteria, by-group comparisons, and all three plot types for `compare()` and `plot_compare()` |
| [BLUP Accuracy in Multi-Environment Trials](https://DrJ001.github.io/biomAid/accuracy.html) | Mrode accuracy and Cullis H², supported random structures, and all six plot types for `accuracy()` and `plot_accuracy()` |

## Function reference

Full documentation for every exported function is available on the package website:

👉 **[https://DrJ001.github.io/biomAid/reference/](https://DrJ001.github.io/biomAid/reference/)**

---

## Installation

```r
# From GitHub (requires remotes or pak)
pak::pkg_install("DrJ001/biomAid")

# or
remotes::install_github("DrJ001/biomAid")
```

> **Note:** The core modelling functions require an [ASReml-R V4](https://vsni.co.uk/software/asreml-r/)
> licence. `simTrialData()` and all `plot_*()` functions are fully standalone.

---

## Functions

### `compare()` — Pairwise comparison criteria

Computes **HSD**, **LSD**, or **Bonferroni**-corrected LSD for predicted values
from an ASReml-R V4 model, optionally within subgroups. Two predicted values in
the same group are significantly different when their absolute difference exceeds
the returned criterion value.

```r
compare(model, term, by = NULL,
        type  = c("HSD", "LSD", "Bonferroni"),
        pev   = TRUE,
        alpha = 0.05,
        ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | Classify string passed to `predict.asreml()` |
| `by` | Column(s) to split comparisons by. `NULL` = one group |
| `type` | `"HSD"` (default), `"LSD"`, or `"Bonferroni"` |
| `pev` | `TRUE` (default) uses prediction error variance; `FALSE` uses posterior variance |
| `alpha` | Significance level. Default `0.05` |

---

### `waldTest()` — Wald / F-tests on contrasts

Tests linear contrasts of predicted values using prediction error information
from `predict.asreml()`. Supports pairwise, custom contrast matrix, and joint
zero tests with optional p-value adjustment.

```r
waldTest(pred, cc, by = NULL,
         test     = c("Wald", "F"),
         df_error = NULL,
         adjust   = c("none", "bonferroni", "holm", "fdr", "BH", "BY"))
```

| Argument | Description |
|----------|-------------|
| `pred` | List returned by `predict(model, vcov = TRUE)` |
| `cc` | Named list of test specifications (`coef`, `type`, `comp`, `group`) |
| `by` | Column(s) to run tests within. `NULL` = single group |
| `test` | `"Wald"` (default, χ² statistic) or `"F"` (requires `df_error`) |
| `df_error` | Denominator degrees of freedom for F-tests (e.g. `model$nedf`) |
| `adjust` | P-value adjustment method. Default `"none"` |

---

### `plot_waldTest()` — Forest plot for Wald test contrasts

Produces a publication-ready forest plot from the output of `waldTest()`.
Each contrast is shown as a filled circle with horizontal confidence interval
bars. Points are coloured by **−log₁₀(p)**: non-significant results appear in
grey with a warm gradient from gold through to dark red as evidence strengthens.
The raw p-value is printed beside each row.

```r
plot_waldTest(res,
              facet       = TRUE,
              ci_level    = 0.95,
              alpha       = 0.05,
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `waldTest()` |
| `facet` | `TRUE` (default): one panel per `by`-group; `FALSE`: single panel |
| `ci_level` | Confidence level for CI arms. Default `0.95` |
| `alpha` | Significance threshold for colour-scale break. Default `0.05` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `randomRegress()` — Random regression (BLUP-based)

Decomposes Genotype × Environment × (Treatment/Trait) variety BLUPs from an
ASReml-R V4 model into **efficiency and responsiveness indices** using a natural
genetic regression derived from Gaussian conditional distribution theory.
Supports four conditioning schemes.

```r
randomRegress(model, Env = "TSite:Variety", levs = NULL,
              type = "baseline", cond = NULL,
              sep = "-", pev = TRUE, ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `Env` | Classify string identifying the G-matrix term. Default `"TSite:Variety"` |
| `levs` | Character vector of treatment levels to decompose |
| `type` | `"baseline"`, `"sequential"`, `"partial"`, or `"custom"` |
| `cond` | User-supplied conditioning list when `type = "custom"` |
| `sep` | Separator between site and treatment in composite labels. Default `"-"` |
| `pev` | `TRUE` (default) uses PEV; `FALSE` uses posterior variance |

---

### `plot_randomRegress()` — Visualise random regression results

Generates ggplot2 visualisations from `randomRegress()` output. Three plot
types are available, each returned as a ggplot object that can be further
customised with `+`.

```r
plot_randomRegress(res,
                   type        = c("regress", "quadrant", "gmat"),
                   treatments  = NULL,
                   highlight   = "default",
                   centre      = FALSE,
                   theme       = ggplot2::theme_bw(),
                   return_data = FALSE,
                   ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `randomRegress()` |
| `type` | `"regress"`, `"quadrant"`, or `"gmat"` |
| `treatments` | Character vector to restrict conditioning pairs plotted. `NULL` = all |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` = no highlighting |
| `centre` | `TRUE` adds back within-site means (useful for demo data). Default `FALSE` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `fixedRegress()` — Fixed regression (BLUE-based)

The fixed-effects analogue of `randomRegress()`. Regresses treatment BLUEs via
OLS within each group and returns **efficiency and response indices** for every
genotype. The same four conditioning schemes are available.

```r
fixedRegress(model, term = "Treatment:Genotype",
             by = NULL, levs = NULL,
             type = "baseline", cond = NULL,
             min_obs = NULL, ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | Classify string passed to `predict.asreml()`. Default `"Treatment:Genotype"` |
| `by` | Column(s) defining groups for separate regressions |
| `levs` | Character vector of treatment levels to decompose |
| `type` | `"baseline"`, `"sequential"`, `"partial"`, or `"custom"` |
| `cond` | User-supplied conditioning list when `type = "custom"` |
| `min_obs` | Minimum common genotypes required to fit a regression. `NULL` = auto |

---

### `plot_fixedRegress()` — Visualise fixed regression results

Generates ggplot2 visualisations from `fixedRegress()` output. Two plot types
are available, each returned as a ggplot object that can be further customised
with `+`.

```r
plot_fixedRegress(res,
                  type        = c("regress", "quadrant"),
                  treatments  = NULL,
                  highlight   = "default",
                  centre      = TRUE,
                  theme       = ggplot2::theme_bw(),
                  return_data = FALSE,
                  ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `fixedRegress()` |
| `type` | `"regress"` or `"quadrant"` |
| `treatments` | Character vector to restrict conditioning pairs plotted. `NULL` = all |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` = no highlighting |
| `centre` | `TRUE` (default) subtracts within-group means; `FALSE` = raw BLUEs |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `padTrial()` — Extract and pad a field trial layout

Extracts a rectangular sub-trial from a field layout data frame based on a
plot-type classification, trims guard rows outside the bounding box of the
target plot types, and pads any missing interior grid positions with blank rows.
Useful for preparing irregular trial layouts for spatial analysis.

```r
padTrial(data,
         pattern    = "Row:Column",
         match      = "DH",
         split      = "Block",
         pad        = TRUE,
         keep       = split,
         fill_value = "Blank",
         type_col   = "Type",
         verbose    = FALSE)
```

| Argument | Description |
|----------|-------------|
| `data` | Data frame containing the trial layout |
| `pattern` | Colon-separated names of the spatial coordinate columns. Default `"Row:Column"` |
| `match` | Plot type(s) defining the target sub-trial (e.g. `"DH"`, `c("DH","Check")`) |
| `split` | Column(s) to process independently (e.g. `"Block"`); `NULL` = whole dataset |
| `pad` | `TRUE` (default) inserts blank rows for missing grid cells |
| `keep` | Column(s) whose values are carried into padded rows. Defaults to `split` |
| `fill_value` | String written into character/factor columns of padded rows. Default `"Blank"` |
| `type_col` | Name of the plot-type column. Default `"Type"` |
| `verbose` | `TRUE` prints a per-group summary message |

---

### `plot_padTrial()` — Before/after field layout tile map

Visualises the effect of `padTrial()` as a pair of tile maps — **Before** on
top, **After** below. Tiles are coloured by plot type; inserted (padded) cells
appear in light grey. When `data` is supplied the Before panel shows the
complete original layout including guard rows.

```r
plot_padTrial(result,
              data        = NULL,
              type_col    = "Type",
              pattern     = "Row:Column",
              split       = NULL,
              label       = NULL,
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```

| Argument | Description |
|----------|-------------|
| `result` | Data frame returned by `padTrial()` |
| `data` | Original data passed to `padTrial()`. `NULL` (default) reconstructs Before from result |
| `type_col` | Name of the plot-type column used to colour tiles. Default `"Type"` |
| `pattern` | Colon-separated names of the spatial coordinate columns. Default `"Row:Column"` |
| `split` | Grouping column(s) matching the `split` used in `padTrial()` |
| `label` | Column whose values are printed inside each tile. `NULL` = no labels |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `fast()` — Factor Analytic Selection Tools

Implements the **FAST** (Smith & Cullis 2018) and **iClass** (Smith et al. 2021)
approaches for summarising variety performance from an FA mixed model.

```r
fast(model, term = "fa(Site, 4):Genotype",
     type   = c("all", "FAST", "iClass"),
     ic.num = 2L,
     ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | FA model term string. Default `"fa(Site, 4):Genotype"` |
| `type` | `"all"` (default), `"FAST"`, or `"iClass"` |
| `ic.num` | Number of interaction class factors (1 to k). Default `2` |

---

### `simTrialData()` — Simulate trial data

Generates a balanced MET or split-plot dataset with a realistic
**factor-analytic genetic covariance structure** across environments.
Variety BLUPs are drawn from `G = ΛΛᵀ + Ψ` where `Λ` is an `n_fa`-factor
loadings matrix — so fitted FA models have genuine signal to recover.

Set `treatments = NULL` for a pure MET (simple RCB, no Treatment column).
Supply `treatments = c("T0", "T1", ...)` for a split-plot design; the FA
structure then operates over all Treatment × Site (`TSite`) combinations.

```r
simTrialData(nvar        = 20L,
             nsite       = 10L,
             treatments  = NULL,
             nrep        = 2L,
             n_fa        = 2L,
             seed        = NULL,
             verbose     = TRUE,
             sim.options = list())
```

| Argument | Description |
|----------|-------------|
| `nvar` | Number of varieties. Default `20` |
| `nsite` | Number of sites. Default `10` |
| `treatments` | Character vector of treatment labels, or `NULL` for MET-only. Default `NULL` |
| `nrep` | Number of replicates per site. Default `2` |
| `n_fa` | Number of FA factors in the data-generating model. Default `2` |
| `seed` | Random seed. `NULL` = no fixed seed |
| `verbose` | Print design summary and suggested model. Default `TRUE` |
| `sim.options` | Named list of optional controls (see below) |

Key `sim.options` elements (all have built-in defaults):

| Element | Description | Default |
|---------|-------------|---------|
| `site_mean` / `site_sd` | Grand mean and SD of site means | `4500` / `600` |
| `sigma_genetic` | Target mean genetic SD per group | `250` |
| `loading_min` / `loading_max` | Range of FA loadings before scaling | `0.5` / `2.5` |
| `specific_pct` | Specific variance as fraction of mean `diag(ΛΛᵀ)` | `0.20` |
| `treat_effects` | Fixed treatment effects vector (multi-treatment only) | auto-spaced |
| `error_sd` / `rep_sd` / `row_sd` / `col_sd` | Error SD components | `350`/`150`/`80`/`60` |
| `sep` | Separator for `TSite` labels | `"-"` |
| `variety_prefix` / `site_prefix` | Label prefixes | `"Var"` / `"Env"` |
| `outfile` | Optional CSV output path | `NULL` |

```r
# MET-only: 30 varieties, 8 sites, FA(2) genetic structure
out <- simTrialData(nvar = 30, nsite = 8, n_fa = 2, seed = 42)
head(out$data)
round(cov2cor(out$params$G), 2)   # true genetic correlations

# Multi-treatment with custom error SD
out2 <- simTrialData(nvar = 20, nsite = 6,
                     treatments  = c("T0", "T1", "T2"),
                     n_fa        = 2,
                     seed        = 1,
                     sim.options = list(treat_effects = c(0, 150, 350),
                                        error_sd      = 200))
```

Returns a list with `$data` (the field layout data frame) and `$params`
(the true `G`, `Lambda`, `Psi`, `site_means`, and for multi-treatment:
`treat_effects`, `g_arr`).

---

### `accuracy()` — Model-based BLUP accuracy

Computes per-environment mean BLUP accuracy from a fitted ASReml-R V4
mixed model. Two complementary metrics are available:

- **Mrode accuracy** — `r = mean[sqrt(1 - PEV / G_jj)]` per variety,
  averaged over environments.
- **Cullis H²** — `1 - mean(SED²) / (2 G_jj)`, a variety-comparison
  reliability analogue.

Supports `fa()`, `diag()`, `corgh()`, `corh()`, `us()`, and single-environment
`id()` random structures. The `only =` argument to `predict.asreml()` is used
internally to avoid fixed-effect inflation of the PEV.

```r
accuracy(model,
         classify   = NULL,
         metric     = c("accuracy", "cullis"),
         pworkspace = "2gb",
         by_variety = FALSE)
```

| Argument | Description |
|----------|-------------|
| `model` | Fitted ASReml-R V4 model object |
| `classify` | Override the auto-detected `predict()` classify string |
| `metric` | One or both of `"accuracy"` and `"cullis"`. Default = both |
| `pworkspace` | Passed to `predict.asreml()`. Default `"2gb"` |
| `by_variety` | `TRUE` returns one row per variety × environment |

```r
# Group-level summaries (one row per environment)
acc <- accuracy(model_fa, metric = c("accuracy", "cullis"))

# Per-variety accuracies
acc_bv <- accuracy(model_fa, by_variety = TRUE)

# Cullis H² only
accuracy(model_fa, metric = "cullis")
```

---

### `plot_accuracy()` — Visualise BLUP accuracy

Six plot types for visualising and comparing accuracy results.
Pass a second accuracy object via `res2` for head-to-head model comparisons.

```r
plot_accuracy(res,
              res2        = NULL,
              type        = c("lollipop", "violin", "heatmap",
                              "dumbbell", "scatter", "diff"),
              metric      = c("accuracy", "cullis"),
              label1      = "Model 1",
              label2      = "Model 2",
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```

| Type | Input | Description |
|------|-------|-------------|
| `"lollipop"` | Group-level | Mean accuracy per environment with ±SD error bars |
| `"violin"` | `by_variety` | Distribution of per-variety accuracies per environment |
| `"heatmap"` | `by_variety` | Variety × environment accuracy grid (viridis fill) |
| `"dumbbell"` | Group-level + `res2` | Paired dots per environment; segment colour = improvement direction |
| `"scatter"` | `by_variety` + `res2` | Per-variety model1 vs model2 accuracy scatter |
| `"diff"` | Group-level + `res2` | Signed accuracy gain (model2 − model1) per environment |

When `metric = c("accuracy", "cullis")` the plot is automatically faceted
into two panels. All types support `return_data = TRUE`.

```r
acc_diag <- accuracy(m_diag, metric = c("accuracy", "cullis"))
acc_fa2  <- accuracy(m_fa2,  metric = c("accuracy", "cullis"))

# Single model
plot_accuracy(acc_fa2, type = "lollipop")
plot_accuracy(acc_fa2_bv, type = "violin")

# Model comparison
plot_accuracy(acc_diag, acc_fa2, type = "dumbbell",
              label1 = "Diag", label2 = "FA(2)")
plot_accuracy(acc_diag, acc_fa2, type = "diff",
              label1 = "Diag", label2 = "FA(2)")
```

---

## References

Smith, A.B. & Cullis, B.R. (2018). Plant breeding selection tools built on
factor analytic mixed models for multi-environment trial data.
*Euphytica*, **214**, 143.

Smith, A.B., Norman, A., Kuchel, H. & Cullis, B.R. (2021). Plant variety
selection using interaction classes derived from factor analytic linear mixed
models. *Frontiers in Plant Science*, **12**, 737462.

---

## License

MIT © J Taylor
