# PAM-Gower Pipeline

**A reusable, dataset-independent R pipeline for k-medoids cluster analysis of
mixed-type data, with the diagnostics needed to defend the solution.**

Version 1.0.0 · Code MIT · Documentation CC BY 4.0
Last updated 09-13-26

---

## What this is

Most applied cluster analyses in the social and health sciences face the same
two problems at once. The variables are of mixed type, some continuous, some
ordinal, some categorical, some binary, and the number of groups is not known
in advance. This pipeline handles both: it builds a Gower dissimilarity matrix
that treats each variable in its own terms, partitions the cases with
partitioning around medoids (PAM), and then assembles the evidence that a
reviewer will ask for and that most published cluster analyses do not report.

It is written to be adapted. Everything specific to a dataset lives in one
configuration block; the rest of the script runs unchanged. It is heavily
annotated, and the annotation explains not only what each section does but why
it is there and what decision it exposes.

**Scope.** The pipeline derives and defends a partition. It does not analyse
what you do with the partition afterwards and exports the cluster assignments in
a form any downstream script can read.

**Before you use it**, read [`docs/decision_guide.md`](docs/decision_guide.md).
It exists to help you decide whether this method suits your data and your
question at all, and it will sometimes tell you it does not.

---

## Quick start

The pipeline ships with a synthetic dataset, so it runs before you adapt it.

1. Download the repository and keep the folder structure intact.
2. Open `pam-gower-pipeline.Rproj` by double-clicking it. Do **not** open the
   `.R` file directly (see [Why the .Rproj file](#why-the-rproj-file-matters)).
3. Install the packages (see [Requirements](#requirements)).
4. Open `pam_gower_pipeline.R` from the Files pane and run it from the top.

It will generate a synthetic community-survey dataset with five known groups,
run the full pipeline on it, and write about fifty tables and nine figures to
`outputs/`. Start reading at `outputs/reporting_checklist.md`.

Then open section 2 of the script and start replacing the configuration with
your own.

---

## What is in this repository

Every file, and what it is for.

| File | What it does |
|---|---|
| `README.md` | This file. |
| `pam_gower_pipeline.R` | **The pipeline.** The only file you edit, and only section 2 of it. Sections 0–1 set up; section 2 is the configuration; sections 3–5 import, screen and build the dissimilarity; section 6 chooses the number of clusters; sections 7–9 describe, validate and report the retained solution. |
| `make_example_data.R` | Generates the synthetic demonstration dataset with five known groups, mixed variable types, and three kinds of missingness. Called automatically when `DATA_PATH` is `"example"`. Run it on its own to save the data to disk. Delete it once you are working with your own data — but keep it if you intend to redistribute this pipeline, since it is what makes the deposit runnable by others. |
| `pam-gower-pipeline.Rproj` | RStudio project file. Opening it sets the working directory so that every relative path in the script resolves. |
| `docs/decision_guide.md` | **Read this before running anything.** A decision tree and comparison table for choosing between k-medoids with Gower, k-means, hierarchical clustering, latent class and latent profile analysis, and k-prototypes — plus the prior question of whether a clustering approach fits your aim at all. |
| `docs/decision_guide.html` | The same guide as a self-contained web page, for reading and sharing. Identical content. |
| `docs/reporting_checklist_template.md` | A blank version of the checklist the pipeline generates, for planning an analysis or reviewing someone else's. The filled-in version appears in `outputs/` after a run. |
| `CITATION.cff` | Machine-readable citation metadata. GitHub and Zenodo read it automatically and offer a formatted citation. |
| `LICENSE` | MIT License, covering the code. |
| `LICENSE-docs` | Creative Commons Attribution 4.0, covering the documentation. |
| `CHANGELOG.md` | What changed between versions. |
| `outputs/` | Created at run time. Emptied at the start of every run unless you set `CLEAR_OUTPUTS <- FALSE`. |

---

## Requirements

**R 4.1 or later.** The pipeline uses no feature specific to a recent version.

**Core packages**, checked at the top of the run, which stops with a single
message listing whatever is missing. Nothing is installed automatically.

```r
install.packages(c("cluster", "fpc", "dplyr", "tidyr", "tibble",
                   "purrr", "readr", "ggplot2"))
```

**Optional packages**, needed only for particular modules:

| Package | Needed for |
|---|---|
| `naniar` | Little's MCAR test (section 4.2). Skipped with a warning if absent. |
| `mice` | Multiple imputation (`MISSING_METHOD <- "imputation"`). |
| `haven` | Reading SPSS `.sav` files. |
| `readxl` | Reading Excel files. |
| `MASS` | Non-metric MDS, if you extend section 8.1. |

```r
install.packages(c("naniar", "mice", "haven", "readxl"))
```

### Why the `.Rproj` file matters

Every path in the script is relative to the project root. Opening the `.Rproj`
file sets R's working directory to that root. Opening the `.R` file on its own
leaves the working directory wherever R happened to start, and the run stops in
section 3 with a message telling you so. If you would rather not use the
project file, set the working directory by hand first:

```r
setwd("path/to/the/downloaded/folder")
```

---

## Adapting it to your own data

Section 2 of the script is the whole interface. Work through it in order.

### 2.1 Point it at your data

```r
DATA_PATH   <- "data/my_survey.rds"
DATA_OBJECT <- NULL          # only for .RData files holding several objects
ID_VAR      <- "participant_id"
```

Recognised formats: `.RData`, `.rda`, `.rds`, `.csv`, `.tsv`, `.sav` (SPSS),
`.xlsx`. You need a column that uniquely identifies each case; every exclusion
log and every export is keyed to it.

### 2.3 Declare the clustering variables and their measurement types

This is the decision that defines the typology. Each variable is declared with
`cvar()`:

```r
CLUSTER_VARS <- list(
  cvar("depression_score", type = "interval",
       note = "PHQ-9 sum, 0-27"),

  cvar("education", type = "ordinal",
       levels = 1:4,
       labels = c("None", "Secondary", "Undergraduate", "Postgraduate")),

  cvar("region", type = "nominal",
       levels = 1:5,
       labels = c("North", "East", "South", "West", "Central")),

  cvar("employed", type = "symm_binary", levels = c(0, 1)),

  cvar("hospitalised", type = "asymm_binary", levels = c(0, 1))
)
```

The five types, and when each applies:

| Type | Use it for | How Gower treats it |
|---|---|---|
| `interval` | Continuous measures; sums or means of several items | Absolute difference, divided by the variable's observed range |
| `ordinal` | Ordered categories with no assumption of equal spacing: single Likert items, education level, severity stage | Rank-based extension of Gower (Podani, 1999) |
| `nominal` | Unordered categories | 0 if the two cases share a category, 1 otherwise |
| `symm_binary` | Two categories, both substantive: agree/disagree, urban/rural | Two cases who both answer "no" count as similar |
| `asymm_binary` | Presence of a rare attribute or event | Shared *presence* counts; shared *absence* does not |

Two rules that matter more than they look:

- **Never enter one construct twice.** Five one-hot indicators derived from one
  single-choice question are *one nominal variable*. Entering all five makes
  that question count five times in every distance while everything else counts
  once. Section 5.2 tries to detect this, but it cannot read your codebook.
- **Keep your validation variables out.** A cluster difference on a variable
  that helped build the clusters is guaranteed and means nothing. Put those in
  `EXTERNAL_VARS` instead.

Declaring the type is not a technical formality. `daisy()` infers measurement
type from the *class* of a column, not from its content: a nominal variable
stored as 1–5 is silently treated as interval, so category 1 sits closer to
category 2 than to category 5 — an ordering that does not exist. No warning is
produced. Section 5.5 audits what `daisy()` actually applied and stops the run
if it differs from what you declared.

### 2.4 Declare scale scores, if you have item-level data

```r
SCALE_DEFS <- list(
  cscale(target      = "depression_score",
         items       = paste0("phq_", 1:9),
         score       = "sum",
         max_missing = 0.12,              # at most 1 item of 9
         reverse     = character(0),      # see the warning below
         item_min    = 0, item_max = 3)
)
```

This does two things. It puts the scoring rule in the script, so it is
reproducible and auditable; and it makes **proration** possible, which is what
lets the "layers" missing-data route retain a respondent who skipped one item
out of nine instead of deleting them.

> **The commonest silent error in this pipeline** is reverse-scoring items that
> your source file already stores reversed. Nothing errors, every descriptive
> statistic looks plausible, and the clustering is wrong. Leave `reverse` empty
> unless you are certain the file holds raw responses. Section 4.6
> cross-checks the recomputed score against any stored score of the same name
> and warns if they disagree.

### 2.5 Eligibility

```r
ELIGIBILITY_FILTER <- quote(consent == 1 & age >= 18 & completed == 1)
```

Eligibility (who belongs in the study at all) is applied *before* any
missing-data rule. The order matters: running the MCAR diagnostic after a
missingness rule evaluates the randomness of missingness on a sample from which
the most incomplete cases have already been removed, which is circular.

### 2.6 Choose a missing-data route

See the [comparison below](#the-three-missing-data-routes).

### 2.8–2.9 Choose k, in two passes

The pipeline is designed to be run twice.

**Pass 1** — leave `K_FINAL <- NULL`. Sections 1–6 run, every model-selection
diagnostic is written to `outputs/tables/`, and the script stops with a summary.

**Pass 2** — read the diagnostics, set `K_FINAL`, write your reasons in
`K_FINAL_RATIONALE`, run again. Sections 7–9 now execute.

The stop at the end of pass 1 is an intentional halt, not a failure. It exists
because choosing k is the decision everything else rests on, and no index
chooses it well enough to be trusted with it (Milligan & Cooper, 1985). Making
the run pause turns the choice into something you wrote down.

`K_FINAL_RATIONALE` is reproduced verbatim in the generated reporting
checklist. If you cannot write that paragraph, you have not finished choosing.

### 2.10 Declare external validation variables

```r
EXTERNAL_VARS <- list(
  cvar("service_utilisation", type = "interval"),
  cvar("employment_status",   type = "nominal", levels = 1:3)
)
```

Variables that did **not** enter the clustering, on which the clusters should
differ if the typology means anything. This is where validity comes from
(Rapkin & Luke, 1993). Everything before section 7.7 establishes that the
partition is internally coherent, stable and not an artefact of noise — all
necessary, none of it evidence that the typology corresponds to anything
outside itself.

---

## The three missing-data routes

PAM needs a complete dissimilarity matrix. `daisy()` will compute Gower
distances in the presence of `NA`s by rescaling over the variables each *pair*
happens to share, so it returns something rather than an error — but two cases
compared on four variables and two compared on six are then placed on the same
scale, and a case with many omissions drifts closer to everything. The
partition then partly reflects the missingness pattern.

So the data entering section 5 must be complete. How they got that way is an
analytic decision with consequences for whom your typology describes.

| | `complete_cases` | `layers` | `imputation` |
|---|---|---|---|
| **What it does** | Deletes any case with a missing value on any clustering variable | Ordered exclusion rules of decreasing severity, then prorates scale scores from the items each case answered | Multiple imputation by chained equations, then clustering |
| **Assumes** | MCAR | That answered items fairly represent the scale | MAR, and a correctly specified imputation model |
| **Adds uncertainty?** | No | No | Yes (and a random seed) |
| **Typical loss** | Largest | Substantially smaller | None |
| **Needs item-level data?** | No | Yes | No |
| **Main risk** | Biases the sample toward complete responders — who are often not the ones populating your smaller, more distinct clusters | A tolerance set too high prorates a score from too little information | Partitions cannot be pooled across imputations (see below) |

### Why imputation is awkward here, and what the pipeline does about it

Multiple imputation works by pooling *estimates* across imputations with
Rubin's rules. A cluster partition is not an estimate: cluster 3 in imputation 1
and cluster 3 in imputation 2 are not the same cluster, they are not guaranteed
to correspond at all, and there is no accepted rule for averaging partitions.
Software that offers "clustering with multiple imputation" usually clusters one
completed dataset and says nothing about the rest.

The pipeline offers two honest alternatives, set with `MI_MODE`:

- **`"sensitivity"` (recommended).** The reported analysis runs on one completed
  dataset. The other imputations answer a question that *can* be answered: how
  much does the partition depend on the imputation? Each imputation is clustered
  separately and the partitions are compared pairwise by adjusted Rand index,
  and case by case by how often two cases land together (section 7.9). High
  agreement means the imputation is not driving the typology. Low agreement
  means it partly is, and you report that.
- **`"consensus"`.** All *m* partitions are combined into a co-membership matrix
  and the final partition is derived from it. Uses all the data, at the cost
  that the result is no longer a PAM solution and has no medoids.

Either way, section 7.8 separately reports what the route bought: it re-clusters
the subset of cases that would also have survived listwise deletion and compares
that partition with the retained one. A high agreement means the missing-data
decision was inconsequential and can be reported in a sentence; a low one means
it was consequential and must be reported as a choice with an effect.

---

## What the script produces

```
outputs/
├── reporting_checklist.md     every decision the run made, with the value used
├── output_inventory.csv       a list of every file below, with descriptions
├── gower_dist.rds             the dissimilarity matrix, saved so it need not
│                              be recomputed
├── session_info.txt           R and package versions for this run
├── tables/                    ~50 .csv files
└── figures/                   9 .png files
```

`outputs/reporting_checklist.md` is the place to start. It is generated from the
configuration and the results, so it cannot drift out of step with what was
actually run, and it is ordered so that it maps onto a methods section. It is a
starting point, not a finished methods section — edit it.

`outputs/output_inventory.csv` describes every other file, so you do not need a
key to this repository beyond the run itself.

### The main figures

| Figure | Shows |
|---|---|
| `model_selection_indices.png` | Six validity indices plotted against k. The *shape* of these curves carries information the numbers do not: a flat silhouette curve means the index is not identifying an optimum, and the later diagnostics decide. |
| `permutation_reference.png` | The observed silhouette against what the same procedure produces on data with the same marginal distributions and no multivariate structure. |
| `stability_bootstrap.png` | Bootstrap Jaccard similarity of every cluster in every candidate solution, with the 0.60 / 0.75 / 0.85 conventions marked. |
| `mds_2d.png`, `mds_pairwise_projections.png` | The solution projected into two and three dimensions, labelled with how much of the distance structure the projection actually preserves. |
| `final_silhouette.png` | Every case's silhouette width, sorted within cluster — the shape of the classification, not just its average. |
| `final_profile_interval.png`, `final_profile_categorical.png` | What distinguishes the clusters. |
| `final_distance_to_medoid.png` | Internal compactness: a long upper tail means a cluster held together by the absence of alternatives rather than by resemblance. |

---

## Decisions the pipeline exposes

These are the choices you are obliged to report. The pipeline records each one
in the generated checklist with the value you actually used.

| Decision | CONFIG setting | Section |
|---|---|---|
| Which cases are eligible at all | `ELIGIBILITY_FILTER` | 4.1 |
| Which variables define the clustering | `CLUSTER_VARS` | 2.3 |
| What measurement type each variable has | `CLUSTER_VARS` (`type=`) | 2.3, 5.1, 5.5 |
| Whether variables are weighted equally | `CLUSTER_VARS` (`weight=`) | 5.2, 5.4 |
| How missing data are handled | `MISSING_METHOD` | 4.3 |
| Whether extreme values are trimmed | `WINSORIZE_INTERVAL` | 5.3 |
| Which values of k are examined | `K_RANGE` | 6.1 |
| Which solutions are compared in depth | `K_CANDIDATES` | 6.4 |
| How stability is quantified | `STABILITY_*` | 6.8 |
| Which k is retained, and why | `K_FINAL`, `K_FINAL_RATIONALE` | 7.1 |
| What the clusters are called | `CLUSTER_LABELS` | 7.5 |
| What evidence validates the solution | `EXTERNAL_VARS` | 7.7 |

---

## Built-in checks

The script verifies its own assumptions rather than continuing silently with
results that would be wrong. Some checks **stop** the run; others **warn** and
let it continue, so that the diagnostic outputs are available for inspection.

| Check | Section | Behaviour |
|---|---|---|
| Required packages installed | 1.1 | stops |
| Configuration internally consistent | 2.12 | stops |
| Data file present, working directory correct | 3.1 | stops |
| Configured columns present in the data | 3.2 | stops |
| Identifier unique and complete | 3.2 | stops |
| Declared types consistent with the data | 3.3 | stops |
| Asymmetric binary coded 0/1 | 5.1 | stops |
| `daisy()` applied the declared types | 5.5 | stops |
| Recomputed scores match stored scores | 4.6 | warns |
| No missing values remain before clustering | 4.8 | stops |
| Sample large enough for the top of `K_RANGE` | 4.8 | warns |
| Highly associated variable pairs | 5.2 | warns |
| Binary variables that are secretly one-hot | 5.2 | warns |
| Interval range driven by extreme values | 5.3 | warns |
| Dissimilarities nearly constant | 5.6 | warns |
| Many tied (identical) profiles | 5.6 | warns |
| Clusters below `MIN_CLUSTER_SIZE` | 6.6 | warns |
| Bootstrap precision target reached | 6.8 | warns |
| Any k exceeds the no-structure reference | 6.11 | warns |
| Duplicate cluster labels | 7.5 | stops |
| Any external variable related to the clusters | 7.7 | warns |

A warning is not something to scroll past. If the score cross-check or the
permutation reference fires, something about the analysis needs tracing before
the outputs are used.

---

## Messages you should expect

**"QUICK_RUN is on."** You set `QUICK_RUN <- TRUE`. Resampling parameters are
reduced so that a first pass finishes quickly. Do not report results from a
quick run.

**"'naniar' is not installed, so the MCAR diagnostic was skipped."** Install
`naniar`, or set `RUN_MCAR_TEST <- FALSE`.

**"Pass 1 complete: set K_FINAL in section 2.9, then run again."** The intended
halt described above, not a failure.

**"k = N: reached STABILITY_B_MAX ..."** A candidate solution has a cluster
whose bootstrap distribution is too dispersed to pin down within the runtime
ceiling — usually because it is bimodal. That is evidence about the solution,
not a problem with the run.

**Little's MCAR test failing to converge**, if you run it on ordinal or binary
items. It assumes multivariate normality; ordinal items produce a rank-deficient
covariance matrix. Report it as indicative at best.

Anything beginning with `Error` other than the pass-1 halt means the run
genuinely stopped and something needs looking at.

---

## Runtime

The shipped example — 731 analysed cases, six clustering variables, every
module on — takes **about two and a half minutes**. Runtime grows roughly with
the square of the number of cases, because the dissimilarity matrix does:
expect something like fifteen to twenty-five minutes at 2,000 cases.

| Module | Setting | Share of the run |
|---|---|---|
| Bootstrap stability | `RUN_STABILITY` | roughly half |
| Permutation reference | `RUN_PERMUTATION` | roughly a quarter |
| Split-half replication | `RUN_SPLIT_HALF` | under a tenth |
| Algorithm robustness | `RUN_ALGO_ROBUSTNESS` | under a tenth |
| Everything else | | seconds |

`QUICK_RUN <- TRUE` lowers every resampling parameter for a first pass.

**Above roughly 10,000 cases** the dissimilarity matrix alone is several
gigabytes and PAM becomes impractical. Use `cluster::clara()`, which samples,
or reconsider the method. See the decision guide.

---

## What this pipeline does not do

Said plainly, because choosing a method partly means knowing what it excludes.

- **It is not model-based.** There is no likelihood, no information criterion,
  no probability of membership. If you want those, latent class or latent
  profile analysis is the better tool, and the decision guide says when.
- **It does not scale.** See above.
- **It does not select variables.** Which variables define the typology is a
  theoretical question, and the pipeline will cluster whatever you give it.
- **It does not handle nested or longitudinal structure.** Cases are assumed
  independent and measured once.
- **It cannot make a partition meaningful.** Every diagnostic here can pass on
  data that answer no interesting question.

---

## How to cite

If this pipeline contributed to your analysis, please cite it. The citation
metadata is in `CITATION.cff`; GitHub and Zenodo render it as a formatted
citation automatically.

```
St-James, D.E. (2026). PAM-Gower Pipeline: k-medoids cluster analysis of
mixed-type data in R (Version 1.0.0) [Computer software].
https://doi.org/[Zenodo DOI]
```

Please also cite the methods themselves, rather than only the software that 
applies them: Gower (1971) for the dissimilarity, Kaufman and Rousseeuw (1990) 
for PAM, Hennig (2007) for the stability assessment.
---

## Licence

Two licences, because code and prose are usually reused differently.

- **`pam_gower_pipeline.R` and `make_example_data.R`** are released under the
  [MIT License](LICENSE). You may use, modify and redistribute them, including
  commercially, provided the copyright notice and licence text travel with the
  code. There is no warranty.
- **Everything else** — this README, the decision guide, the checklist template
  — is released under
  [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/),
  recorded in `LICENSE-docs`. You may reuse and adapt it, including
  commercially, provided you give credit and indicate changes.

Both permit reuse. The difference is mostly conventional: software licences
address warranty and liability, which content licences do not, and Creative
Commons itself recommends against using CC licences for software 
(https://creativecommons.org/faq/#can-i-apply-a-creative-commons-license-to-software).

---

## Contributing

Issues and pull requests are welcome, particularly: additional measurement
types, support for larger samples, and corrections to the methodological
claims in the annotation. If you use the pipeline in a published analysis, an
issue saying so is genuinely useful: it helps establish what requires better
documentation.

---

## References

Gower, J. C. (1971). A general coefficient of similarity and some of its
properties. *Biometrics, 27*(4), 857–871. https://doi.org/10.2307/2528823

Hennig, C. (2007). Cluster-wise assessment of cluster stability.
*Computational Statistics & Data Analysis, 52*(1), 258–271.
https://doi.org/10.1016/j.csda.2006.11.025

Hennig, C., & Liao, T. F. (2013). How to find an appropriate clustering for
mixed-type variables with application to socio-economic stratification.
*Journal of the Royal Statistical Society: Series C (Applied Statistics),
62*(3), 309–369. https://doi.org/10.1111/j.1467-9876.2012.01066.x

Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
Classification, 2*(1), 193–218. https://doi.org/10.1007/BF01908075

Kaufman, L., & Rousseeuw, P. J. (1990). *Finding groups in data: An
introduction to cluster analysis*. Wiley.
https://doi.org/10.1002/9780470316801

Luke, D. A. (2005). Getting the big picture in community science: Methods that
capture context. *American Journal of Community Psychology, 35*(3–4), 185–200.
https://doi.org/10.1007/s10464-005-3397-z

Milligan, G. W., & Cooper, M. C. (1985). An examination of procedures for
determining the number of clusters in a data set. *Psychometrika, 50*(2),
159–179. https://doi.org/10.1007/BF02294245

Podani, J. (1999). Extending Gower's general coefficient of similarity to
ordinal characters. *Taxon, 48*(2), 331–340. https://doi.org/10.2307/1224438

Rapkin, B. D., & Luke, D. A. (1993). Cluster analysis in community research:
Epistemology and practice. *American Journal of Community Psychology, 21*(2),
247–277. https://doi.org/10.1007/BF00941623

van Buuren, S. (2018). *Flexible imputation of missing data* (2nd ed.). CRC
Press. https://stefvanbuuren.name/fimd/
