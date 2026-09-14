# PAM-GOWER PIPELINE ======================================================
#  k-medoids cluster analysis of mixed-type data, with Gower dissimilarity
#  A reusable, dataset-independent analytic pipeline for R
#
# 0.0 License -------------------------------------------------------------
#
#  Version   : 1.0.0
#  Code      : MIT License          (see LICENSE)
#  Documents : CC BY 4.0            (see LICENSE-docs)
#  Citation  : see CITATION.cff in the repository root
#  Repository: https://github.com/DaemSt-James/pam-gower-pipeline

#  0.1  WHAT THIS SCRIPT DOES ---------------------------------------------
#  To derive and defend a partition in a given sample this script takes a 
#  dataset whose variables of interest are of MIXED TYPE (continuous, ordinal, 
#  nominal and/or binary) and :
#  1) partitions the CASES into k groups (k = number of clusters in the dataset)
#  2) produce the diagnostic evidence needed to defend the choice of k, the 
#  choice of method, and the claim that the groups are worth interpreting (not
#  due to chance).
#
#  Method: partitioning around medoids (PAM; Kaufman & Rousseeuw, 1990) on a
#  Gower dissimilarity matrix (Gower, 1971).
#
#  Section 7.6 exports the cluster assignments for downstream script.
#
#  0.2  WHAT THIS SCRIPT ASSUMES ABOUT YOU -----------------------------------
#  That you have decided, on substantive grounds, that a person-oriented
#  question is the right question; that the variables entering the clustering
#  were chosen by theory rather than by availability; and that you are willing
#  to report the decisions you made rather than only the solution you kept.
#
#  If you are not yet sure that clustering, or this particular clustering,
#  suits your data and your aims, read docs/decision_guide.md BEFORE running
#  anything here.
#
#  0.3  THE ONE THING YOU EDIT ------------------------------------------------
#  SECTION 2 (CONFIGURATION). Nothing else. Every downstream section reads its
#  settings from there, and section 2.10 refuses to continue if the settings
#  are internally inconsistent.
#
#  To see the pipeline run before you adapt it, leave the configuration exactly
#  as shipped: DATA_PATH is set to "example" and a synthetic mixed-type dataset
#  with known structure is generated on the fly.
#
#  0.4  HOW THE ANNOTATION WORKS-----------------------------------------------
#  Every numbered section carries a four-line header:
#
#    WHAT    what the section computes
#    WHY     why the pipeline needs it, and what goes wrong without it
#    CHOICE  the decision the section requieres, and the CONFIG setting that
#            controls it (omitted when the section requires no choice)
#    STATUS  REQUIRED, or OPTIONAL with the condition under which it matters
#
#  Read the WHY lines first if you are deciding what to run. Read the CHOICE
#  lines when you are writing your method section. These are the decisions you 
#  are obliged to report when publishing or sharing your results. 
#
#  0.5  DECISION MAP (where each decision is made) ---------------------------
#
#   Decision                                     CONFIG setting        Section
#   -------------------------------------------  --------------------  -------
#   Which cases are eligible at all              ELIGIBILITY_FILTER    4.1
#   Which variables define the clustering        CLUSTER_VARS          2.3
#   What measurement type each variable has      CLUSTER_VARS (type=)  2.3 / 5.1
#   Whether variables are weighted equally       CLUSTER_VARS (weight=) 5.4
#   How missing data are handled                 MISSING_METHOD        4.3
#   Whether extreme values are trimmed           WINSORIZE_INTERVAL    5.3
#   Which values of k are examined               K_RANGE               6.1
#   Which solutions are compared in detail       K_CANDIDATES          6.4
#   How stability is quantified                  STABILITY_*           6.8
#   Which k is retained, and on what grounds     K_FINAL               7.1
#   What the clusters are called                 CLUSTER_LABELS        7.5
#   What evidence validates the solution         EXTERNAL_VARS         7.7
#
#  0.6  METHODOLOGICAL COMMITMENTS---------------------------------------------
#  Four commitments are built into the pipeline in accordance with the 
#  methodological literature on cluster analysis in the social and community 
#  sciences, principally Rapkin and Luke (1993) and Luke (2005).
#
#  (1) A CLUSTERING ALGORITHM ALWAYS RETURNS CLUSTERS.
#      PAM will partition random noise into k tidy groups and report a
#      silhouette width for them. The question is never "did the algorithm
#      find groups" but "is this partition distinguishable from what the same
#      algorithm produces on data with no multivariate structure". Section
#      6.11 answers that question with a permutation reference distribution.
#
#  (2) CLUSTERS ARE CONSTRUCTED, NOT DISCOVERED.
#      The solution is a joint product of the variables you entered, the
#      dissimilarity you chose, the algorithm, and k. Change any one and the
#      typology changes. This is not a flaw to be apologised for, but it is why
#      all four should be reported, which section 9.1 does automatically.
#
#  (3) INTERNAL INDICES ARE EVIDENCE, NOT VERDICTS.
#      Silhouette, Dunn and Calinski-Harabasz measure geometric separation. 
#      They cannot tell you whether a cluster is substantively distinct. 
#      Section 6 therefore assembles seven criteria but will not rank the 
#      candidates for you.
#
#  (4) VALIDITY IS EXTERNAL.
#      A partition earns interpretation by relating to variables that did not
#      help produce it. Testing clusters on the variables used to build them
#      is circular and always "significant". Section 7.3 blocks the circular
#      test; section 7.7 runs the external one.
#
#  0.7  REFERENCES CITED IN THIS SCRIPT ---------------------------------------
#  Full list, with DOIs, in README.md. Short forms used in the comments:
#
#    Gower (1971)              general coefficient of similarity, mixed data
#    Kaufman & Rousseeuw (1990) PAM, silhouettes, "Finding Groups in Data"
#    Rapkin & Luke (1993)      cluster analysis in community research
#    Luke (2005)               context-capturing methods in community science
#    Hennig (2007)             bootstrap cluster stability, Jaccard criterion
#    Hennig & Liao (2013)      choosing dissimilarity and weights, mixed data
#    Milligan & Cooper (1985)  comparison of number-of-clusters indices
#    Podani (1999)             extension of Gower to ordinal variables
#    van Buuren (2018)         multiple imputation by chained equations
#
#
#  1.  SETUP ==================================================================
#  1.1  Dependency check ------------------------------------------------------
#  WHAT    Verifies that every package the run needs is installed, and stops
#          with a single message listing whatever is missing.
#  WHY     A run that fails eight minutes in, on a missing package, wastes the
#          eight minutes. Checking first costs a second.
#  CHOICE  None. Note that nothing is installed automatically: the script does
#          not modify your library without your say-so.
#  STATUS  REQUIRED

PKGS_CORE <- c(
  "cluster",   # daisy(), pam(), silhouette(), agnes(), diana(), fanny()
  "fpc",       # cluster.stats(), clusterboot()
  "dplyr", "tidyr", "tibble", "purrr", "readr",  # data handling and export
  "ggplot2"    # figures
)

# Checked only if the module that needs them is switched on (see 2.9).
PKGS_OPTIONAL <- c(
  naniar      = "missingness summaries and Little's MCAR test  (section 4.2)",
  mice        = "multiple imputation                           (section 4.3C)",
  haven       = "reading SPSS .sav files                       (section 3.1)",
  readxl      = "reading Excel files                           (section 3.1)",
  MASS        = "non-metric MDS and its stress value           (section 8.1)"
)

check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "\n  Missing package(s): ", paste(missing, collapse = ", "),
      "\n  Install them with:\n\n    install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), "))\n",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

check_packages(PKGS_CORE)

library(cluster)
library(dplyr)
library(ggplot2)

# Small helper used throughout: fail with a readable message when an optional
# package is needed by a module the user switched on.
need_pkg <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("\n  '", pkg, "' is required for ", what,
         ".\n  Install it with: install.packages(\"", pkg, "\")",
         "\n  Or switch that module off in section 2.", call. = FALSE)
  }
  invisible(TRUE)
}


#  1.2  Variable-specification helper ------------------------------------------
#  WHAT    cvar() builds one entry of the variable specification used in 2.3.
#  WHY     Gower does not infer measurement type from content; it infers it
#          from the R class of the column. A nominal variable stored as 1..5
#          is silently treated as interval, which makes category 1 closer to
#          category 2 than to category 5 -- an ordering that does not exist.
#          Declaring the type explicitly, once, removes the guesswork and
#          gives section 5.5 something to audit the result against.
#  CHOICE  The type of each variable. This is a measurement decision, not a
#          technical one, and it belongs to you rather than to the software.
#  STATUS  REQUIRED
#
#  type =
#    "interval"      numeric, differences meaningful. Gower rescales each
#                    interval variable by its observed range, so the variable
#                    contributes a value in [0, 1] to the distance. See 5.3
#                    for what a single extreme value does to that range.
#
#    "ordinal"       ordered categories with no assumption of equal spacing
#                    (Likert responses, education level, stage). Handled by
#                    the rank-based extension of Gower (Podani, 1999). Use
#                    this rather than "interval" for a single Likert item;
#                    use "interval" for a SUM or MEAN of several items.
#
#    "nominal"       unordered categories. Contributes 0 if two cases share a
#                    category and 1 otherwise (simple matching).
#
#    "symm_binary"   two categories, both of which are substantive positions.
#                    "Agree / disagree", "urban / rural". Two cases who both
#                    answer "no" are counted as similar.
#
#    "asymm_binary"  two categories where only one is informative -- presence
#                    of a rare attribute, occurrence of an event. Two cases
#                    who both lack the attribute are NOT counted as similar,
#                    because "neither has it" says little about resemblance.
#                    Must be coded 0 = absent, 1 = present (enforced in 2.10).
#
#  levels / labels
#    levels  the codes as they appear in your data, in the order you want.
#    labels  optional human-readable names, same length as levels. Used in
#            profile tables and figures. For binary variables labels are
#            ignored: the first level becomes 0 and the second becomes 1,
#            because daisy() coerces declared binary columns to numeric.
#
#  weight
#    Relative contribution to the distance. Default 1 = all variables equal.
#    Read section 5.4 before changing it: the default is already a decision.
#
#  note
#    Free text. Reproduced in the codebook (3.4) and in the reporting
#    checklist (9.1). Use it to record WHY the variable is in the model.

cvar <- function(name,
                 type   = c("interval", "ordinal", "nominal",
                            "symm_binary", "asymm_binary"),
                 levels = NULL,
                 labels = NULL,
                 weight = 1,
                 note   = "") {
  type <- match.arg(type)
  if (!is.null(labels) && !is.null(levels) &&
      length(labels) != length(levels)) {
    stop("cvar('", name, "'): 'labels' and 'levels' differ in length.",
         call. = FALSE)
  }
  list(name = name, type = type, levels = levels,
       labels = labels, weight = weight, note = note)
}

# Accessors, so later sections never index into the list structure by hand.
spec_names   <- function(spec) vapply(spec, function(v) v$name,   character(1))
spec_types   <- function(spec) vapply(spec, function(v) v$type,   character(1))
spec_weights <- function(spec) vapply(spec, function(v) v$weight, numeric(1))
spec_get     <- function(spec, nm) spec[[which(spec_names(spec) == nm)]]


#  1.3  Scale-definition helper------------------------------------------------
#  WHAT    cscale() declares that a clustering variable is a SCALE SCORE
#          computed from a set of items, rather than a variable read directly
#          from the file.
#  WHY     Two reasons, and the second is the important one.
#          First, the scoring rule then lives in this script rather than in
#          whatever software produced the file, so it is reproducible and
#          auditable.
#          Second, it makes PRORATION possible (4.3, route B): a respondent
#          who skipped one item out of twenty-four can still be scored from
#          the twenty-three they answered, instead of being deleted. Without
#          item-level information the pipeline can only delete or impute.
#  CHOICE  Whether your clustering variables are scale scores at all. If they
#          are single measured variables, leave SCALE_DEFS empty.
#  STATUS  OPTIONAL - required only for route B of section 4.3.
#
#  items            character vector of item column names.
#  score            "sum" or "mean". Prorated scores use the mean of answered
#                   items, multiplied by the number of items when "sum".
#  max_missing      largest proportion of items a case may omit and still be
#                   scored. 0 means the scale tolerates no missing item.
#  reverse          items to reverse-score. LEAVE EMPTY if your file already
#                   stores them reversed -- reversing twice restores the
#                   original and silently distorts the scale. Checked in 4.4.
#  item_min/max     the item response range, needed for reversal.
#  target           the name of the resulting clustering variable. It must
#                   appear in CLUSTER_VARS.

cscale <- function(target, items, score = c("sum", "mean"),
                   max_missing = 0, reverse = character(0),
                   item_min = NA_real_, item_max = NA_real_, note = "") {
  score <- match.arg(score)
  if (length(reverse) > 0 && (is.na(item_min) || is.na(item_max))) {
    stop("cscale('", target, "'): reverse-scoring requires item_min and ",
         "item_max.", call. = FALSE)
  }
  list(target = target, items = items, score = score,
       max_missing = max_missing, reverse = reverse,
       item_min = item_min, item_max = item_max, note = note)
}


#  1.4  Console logging --------------------------------------------------------
#  WHAT    Uniform section banners and status lines.
#  WHY     A long run should be legible while it is running, and the console
#          transcript should be usable as a record of what happened.
#  STATUS  REQUIRED

banner <- function(txt) {
  cat("\n", strrep("=", 74), "\n  ", txt, "\n", strrep("=", 74), "\n", sep = "")
}
step <- function(txt) cat("\n--- ", txt, " ", strrep("-", max(0, 62 - nchar(txt))),
                          "\n", sep = "")
ok   <- function(...) cat("   [ok]   ", ..., "\n", sep = "")
note <- function(...) cat("   [note] ", ..., "\n", sep = "")
warn <- function(...) cat("   [WARN] ", ..., "\n", sep = "")

#  1.5  Type coercion and Gower helpers ---------------------------------------
#  WHAT    Turns each declared variable into the R class that makes daisy()
#          treat it as intended, and assembles the arguments daisy() needs.
#  WHY     These live here, before the configuration, because two separate
#          sections need them: the imputation route (4.5) must coerce types
#          BEFORE imputing, so that mice chooses the right elementary model,
#          and section 5.1 coerces them again for the distance computation.
#          Section 5.1 is where the reasoning is explained; this is only the
#          machinery.
#  STATUS  REQUIRED
#
#  The mapping from declared type to R class, and to the type code daisy()
#  reports in attr(d, "Types"):
#
#     declared        R class passed to daisy()      daisy code
#     ________        __________________________     ___________
#     interval        numeric                        I
#     ordinal         ordered factor                 O
#     nominal         unordered factor (>2 levels)   N
#     symm_binary     numeric 0/1, listed in symm    S
#     asymm_binary    numeric 0/1, listed in asymm   A
#
#  Binary variables are passed as numeric rather than as labelled factors:
#  daisy() coerces declared binary columns to numeric anyway, and character
#  labels produce a coercion warning while yielding an identical matrix.

#  The function is IDEMPOTENT: applying it twice is the same as applying it
#  once. This is not a nicety. The imputation route coerces types before
#  imputing (4.5) and section 5.1 coerces again before computing distances, so
#  a second pass happens by design. Without the guards below, re-levelling a
#  column whose values are already labels would map every one of them to NA,
#  and the failure would surface far downstream as an empty contingency table
#  rather than as anything recognisable.

coerce_gower_column <- function(x, v) {
  lv <- v$levels
  if (is.null(lv)) lv <- sort(unique(x[!is.na(x)]))
  lb <- if (is.null(v$labels)) as.character(lv) else as.character(v$labels)

  if (v$type == "interval") return(as.numeric(x))

  if (v$type %in% c("ordinal", "nominal")) {
    want_ordered <- v$type == "ordinal"
    if (is.factor(x) && identical(levels(x), lb) &&
        is.ordered(x) == want_ordered) return(x)          # already coerced
    return(factor(x, levels = lv, labels = lb, ordered = want_ordered))
  }

  # Binary: the FIRST declared level becomes 0, the second becomes 1. For
  # asymmetric variables this ordering is substantive -- 1 must be the
  # category whose shared presence signals resemblance.
  if (is.numeric(x) && all(stats::na.omit(x) %in% c(0, 1)) &&
      !all(as.character(lv) %in% c("0", "1"))) return(as.numeric(x))
  as.numeric(match(x, lv) - 1)
}

prepare_gower_frame <- function(data, spec) {
  out <- data[, spec_names(spec), drop = FALSE]
  for (v in spec) out[[v$name]] <- coerce_gower_column(out[[v$name]], v)
  out
}

# daisy() is told which columns are binary, and how. Everything else it reads
# from the column classes set above.
build_daisy_type <- function(spec) {
  nm  <- spec_names(spec)
  ty  <- spec_types(spec)
  out <- list()
  if (any(ty == "symm_binary"))  out$symm  <- nm[ty == "symm_binary"]
  if (any(ty == "asymm_binary")) out$asymm <- nm[ty == "asymm_binary"]
  out
}

expected_daisy_codes <- function(spec) {
  codes <- c(interval = "I", ordinal = "O", nominal = "N",
             symm_binary = "S", asymm_binary = "A")
  stats::setNames(unname(codes[spec_types(spec)]), spec_names(spec))
}

# Adjusted Rand index between two partitions, corrected for chance agreement
# (Hubert & Arabie, 1985). Implemented here rather than taken from a package
# so that every partition comparison in the pipeline -- nesting, split-half,
# algorithm robustness, imputation sensitivity -- uses the same function.
#   1.0  identical partitions
#   0.0  the agreement expected by chance
#   < 0  less agreement than chance
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n   <- sum(tab)
  if (n < 2) return(NA_real_)
  comb2   <- function(x) sum(choose(x, 2))
  index   <- comb2(as.vector(tab))
  a_sum   <- comb2(rowSums(tab))
  b_sum   <- comb2(colSums(tab))
  expect  <- a_sum * b_sum / choose(n, 2)
  maximum <- (a_sum + b_sum) / 2
  if (isTRUE(all.equal(maximum, expect))) return(NA_real_)
  (index - expect) / (maximum - expect)
}

#  2.  CONFIGURATION  ***THE ONLY SECTION YOU EDIT*** ==========================
#
#  The pipeline is designed to be run TWICE.
#
#    Pass 1   Leave K_FINAL as NULL (2.9). The script runs sections 1 to 6,
#             writes every model-selection diagnostic to outputs/tables/, and
#             stops with a summary of what it found. You then read those
#             diagnostics and decide how many clusters to retain.
#
#    Pass 2   Set K_FINAL to the number you decided on, write your reasons in
#             K_FINAL_RATIONALE, and run again. Sections 7 to 9 now execute.
#
#  This is deliberate. Choosing k is the decision on which everything else
#  rests, and no index chooses it well enough to be trusted with it (Milligan
#  & Cooper, 1985). Forcing the run to pause makes the choice visible, dated
#  and attributable instead of implicit in a line of code.
#
#  The configuration below is filled in for the synthetic example dataset, so
#  that the pipeline runs end to end before you adapt it.
#
#  2.1  Data source -----------------------------------------------------------
#  DATA_PATH     path to your file, relative to the project root.
#                Recognised: .RData/.rda, .rds, .csv, .tsv, .sav (SPSS, needs
#                'haven'), .xlsx/.xls (needs 'readxl').
#                The literal string "example" generates the synthetic
#                demonstration dataset instead of reading a file.
#  DATA_OBJECT   for .RData files containing more than one object, the name of
#                the one to use. NULL when the file holds a single object.
#  ID_VAR        a column that uniquely identifies each case. Required: every
#                exclusion log, every export and every cross-check is keyed to
#                it. If you have no such column, create one before running.

DATA_PATH   <- "example"
DATA_OBJECT <- NULL
ID_VAR      <- "id"

#  2.2  Output location -------------------------------------------------------
#  CLEAR_OUTPUTS = TRUE empties OUTPUT_DIR at the start of every run, so the
#  inventory printed at the end lists exactly what this run produced and
#  nothing left over from a previous version. Move anything you edited by hand
#  out of that folder first.

OUTPUT_DIR    <- "outputs"
CLEAR_OUTPUTS <- TRUE

#  2.3  Variables entering the cluster analysis <- the core decision -----------
#  This list defines the space in which distance is measured, and therefore
#  defines the typology. Nothing later in the pipeline can repair a badly
#  chosen variable set; adding a variable is not a neutral act, and neither is
#  omitting one (Rapkin & Luke, 1993).
#
#  Three rules worth holding to:
#
#  (a) Enter variables because theory says they distinguish types, not because
#      they were collected. A variable that is conceptually irrelevant does
#      not "wash out" -- it adds noise to every pairwise distance.
#
#  (b) Do not enter the same construct twice. Five one-hot indicators derived
#      from one single-choice question are ONE nominal variable, not five
#      binary ones; entering all five makes that construct count five times.
#      Section 5.2 tries to detect this, but it cannot read your codebook.
#
#  (c) Do not enter the variables you intend to validate the solution on.
#      Keep those for EXTERNAL_VARS (2.10). A cluster difference on a variable
#      that helped build the clusters is guaranteed and means nothing.
#
#  See cvar() in section 1.2 for the meaning of each type.

CLUSTER_VARS <- list(

  cvar("efficacy_total", type = "interval",
       note = "Collective efficacy, sum of 8 items (see SCALE_DEFS)."),

  cvar("social_support", type = "interval",
       note = "Perceived social support, 0-100."),

  cvar("econ_strain", type = "ordinal",
       levels = 1:5,
       labels = c("None", "Slight", "Moderate", "Considerable", "Severe"),
       note   = "Ordered severity; spacing between categories not assumed."),

  cvar("neighborhood_type", type = "nominal",
       levels = 1:4,
       labels = c("Urban core", "Inner suburb", "Outer suburb", "Rural"),
       note   = "Unordered. Entered once, not as one-hot indicators."),

  cvar("civic_member", type = "symm_binary",
       levels = c(0, 1),
       note   = "Member of a civic association. Symmetric: non-membership is
                 itself a position, so two non-members count as similar."),

  cvar("service_use", type = "asymm_binary",
       levels = c(0, 1),
       note   = "Used a specialised service in the past year. Asymmetric:
                 the vast majority did not, and shared non-use carries little
                 information about resemblance.")
)

#  2.4  Scale scores computed from items ---------------------------------------
#  Leave as list() if your clustering variables are read directly from the
#  file. See cscale() in section 1.3.
#
#  Any target named here is RECOMPUTED from its items and overwrites whatever
#  column of that name the file contained. Section 4.5 cross-checks the
#  recomputed values against the stored ones on complete cases: if the two
#  disagree, your scoring rule differs from the one used originally and the
#  discrepancy must be resolved before you go further.

SCALE_DEFS <- list(
  cscale(target      = "efficacy_total",
         items       = paste0("eff_", 1:8),
         score       = "sum",
         max_missing = 0.125,   # at most 1 of 8 items may be missing
         reverse     = character(0),
         item_min    = 1, item_max = 5,
         note        = "Prorated from answered items when 1 item is omitted.")
)

#  2.5  Eligibility screening -------------------------------------------------
#  An expression evaluated against the imported data. Cases for which it is
#  not TRUE are removed FIRST, before any missing-data rule is applied.
#
#  Keep eligibility and missingness separate and in this order. Eligibility is
#  about who belongs in the study at all (consent, age, completion); it is a
#  design criterion. Missingness is about how much of the instrument a valid
#  respondent left blank. Mixing the two makes the attrition table
#  uninterpretable, and -- more seriously -- applying an eligibility rule
#  after a missingness rule means the MCAR diagnostic in 4.2 is computed on a
#  sample from which the most incomplete cases have already been removed,
#  which is circular.
#
#  Set to NULL to skip.

ELIGIBILITY_FILTER <- quote(finished == 1 & consent == 1 & age >= 18)

#  2.6  Missing data (see 4.3 for detail) --------------------------------------
#  MISSING_METHOD, one of:
#
#    "complete_cases"  Delete any case with a missing value on any clustering
#                      variable. Simplest and fully transparent. Defensible
#                      when losses are small (say under 5%) and the MCAR
#                      diagnostic is not rejected. Biased otherwise.
#
#    "layers"          Ordered exclusion rules of decreasing severity, then
#                      proration of scale scores from the items each case did
#                      answer. Deletes far fewer cases than listwise deletion,
#                      introduces no model-based uncertainty, and needs no
#                      random seed. Requires SCALE_DEFS.
#
#    "imputation"      Multiple imputation by chained equations (van Buuren,
#                      2018), then clustering. The statistically principled
#                      option under MAR, but the one that interacts worst with
#                      clustering -- see 4.3C, which explains why and what the
#                      pipeline does about it.
#
#  OVERALL_MISSING_MAX  Used by "layers". A case missing more than this
#                       proportion across ALL clustering variables is removed
#                       before the per-scale rules apply. NULL skips the layer.

MISSING_METHOD      <- "layers"
OVERALL_MISSING_MAX <- 0.20

#  Settings for MISSING_METHOD == "imputation" only.
#
#  MI_MODE = "sensitivity"  the main analysis runs on one completed dataset
#                           (MI_REFERENCE), and the remaining imputations are
#                           used to quantify how much the partition depends on
#                           the imputation. Recommended, and the easier of the
#                           two to report.
#          = "consensus"    all m partitions are pooled into a co-membership
#                           matrix and the final partition is derived from it.
#                           Uses all the information, but the resulting object
#                           is no longer a PAM solution and has no medoids.
#
#  MI_AUXILIARY_VARS  variables NOT used for clustering that help predict the
#                     missing values. Including good auxiliaries is the single
#                     most effective way to make the MAR assumption plausible.

MI_N_IMPUTATIONS  <- 5
MI_MODE           <- "sensitivity"
MI_REFERENCE      <- 1
MI_AUXILIARY_VARS <- c("age")

#  2.7  Preparation of interval variables --------------------------------------
#  Gower divides each interval variable by its OBSERVED RANGE. One extreme
#  value therefore stretches the denominator and compresses every other
#  difference on that variable towards zero -- the variable quietly stops
#  contributing. Winsorising caps the extremes before the range is taken.
#
#  This is a real analytic decision, not a cleaning step: it changes the
#  distance matrix. Section 5.3 reports what changed so you can judge it.
#  Leave FALSE unless 5.3 shows a variable whose range is driven by a handful
#  of cases.

WINSORIZE_INTERVAL <- FALSE
WINSORIZE_PROBS    <- c(0.01, 0.99)

#  2.8  Model selection --------------------------------------------------------
#  K_RANGE          values of k to evaluate. Start at 2. The upper end should
#                   be a number you could still interpret substantively.
#  K_CANDIDATES     the solutions examined in depth in 6.5 to 6.10. NULL means
#                   "the three highest average silhouette widths in K_RANGE".
#                   Override when theory points at particular values.
#  MIN_CLUSTER_SIZE the size below which a cluster is flagged as too small to
#                   carry an interpretation or a downstream comparison. A flag,
#                   not a rule: nothing is excluded automatically.

K_RANGE          <- 2:10
K_CANDIDATES     <- NULL
MIN_CLUSTER_SIZE <- 30

#  2.9  The retained solution --------------------------------------------------
#  K_FINAL            NULL on pass 1. On pass 2, the number of clusters you
#                     decided to retain after reading section 6's diagnostics.
#  K_FINAL_RATIONALE  your reasons, in prose. Reproduced verbatim in the
#                     reporting checklist (9.1). Not decoration: if you cannot
#                     write this paragraph, you have not finished choosing.
#  CLUSTER_LABELS     NULL, or a character vector of length K_FINAL giving a
#                     substantive name to each cluster, in cluster order.
#                     Section 7.5 checks the labels against the profiles.

K_FINAL <- 5

K_FINAL_RATIONALE <- "
  Example dataset; five clusters retained.

  The average silhouette is effectively tied between k = 5 (.580) and k = 6
  (.579), so it does not decide between them, and on stability k = 6 is if
  anything the stronger of the two. The decision rests on the profiles in
  6.7. The sixth cluster is obtained by splitting the low-efficacy urban
  cluster into service users and non-users; the two halves are otherwise
  indistinguishable -- collective efficacy, social support, economic strain,
  neighbourhood type and civic membership are all within sampling noise of
  each other. That is a subdivision on one indicator, not a distinct profile,
  so k = 5 is retained. Both solutions clear the permutation reference in
  6.11 and replicate across halves in 6.9.

  (Replace this text with your own reasoning. It is reproduced verbatim in
  the reporting checklist and is part of what you deposit.)
"

CLUSTER_LABELS <- NULL

#  2.10  External validation variables -----------------------------------------
#  Variables that did NOT enter the clustering, on which the clusters should
#  differ if the typology means anything. This is the evidence that separates
#  a typology from an arbitrary partition, and the reason it is worth keeping
#  a substantively interesting variable OUT of the clustering.
#
#  Declared with cvar(), same types. Leave as list() to skip -- but read 7.7
#  before you do.

EXTERNAL_VARS <- list(
  cvar("distress", type = "interval",
       note = "Psychological distress, 0-40. Not used to build the clusters."),
  cvar("employment_status", type = "nominal",
       levels = 1:3,
       labels = c("Employed", "Unemployed", "Not in labour force"))
)

#  2.11  Modules and run parameters --------------------------------------------
#  Each module below can be switched off without affecting anything else. The
#  estimates are for roughly 900 cases and six clustering variables; runtime
#  grows with the square of the number of cases.
#
#    RUN_MCAR_TEST        ~1 s     Little's test and missingness summaries
#    RUN_STABILITY        2-5 min  bootstrap Jaccard per candidate solution
#    RUN_SPLIT_HALF       ~1 min   replication across random halves
#    RUN_PERMUTATION      2-4 min  reference distribution under no structure
#    RUN_ALGO_ROBUSTNESS  ~1 min   hierarchical and fuzzy alternatives
#    RUN_MDS              ~5 s     projection figures
#
#  QUICK_RUN = TRUE lowers every resampling parameter for a first pass. Never
#  report results from a quick run.

RUN_MCAR_TEST       <- TRUE
RUN_STABILITY       <- TRUE
RUN_SPLIT_HALF      <- TRUE
RUN_PERMUTATION     <- TRUE
RUN_ALGO_ROBUSTNESS <- TRUE
RUN_MDS             <- TRUE

QUICK_RUN <- FALSE

SEED <- 2026

# Bootstrap stability (6.8). Replications are added in blocks until the mean
# Jaccard of every cluster is estimated to within TARGET_MCSE. This is a
# precision rule -- it depends on the spread of the bootstrap distribution,
# never on the cluster means -- so it cannot become a stop-when-favourable
# rule. Set STABILITY_ADAPTIVE to FALSE for a fixed B = STABILITY_B_MIN.
STABILITY_ADAPTIVE    <- TRUE
STABILITY_TARGET_MCSE <- 0.010
STABILITY_BLOCK       <- 100
STABILITY_B_MIN       <- 100
STABILITY_B_MAX       <- 1000

PERM_N  <- 50    # permutation replications (6.11)
SPLIT_N <- 25    # random splits (6.9)

if (isTRUE(QUICK_RUN)) {
  STABILITY_B_MIN <- 25; STABILITY_B_MAX <- 50; STABILITY_BLOCK <- 25
  PERM_N <- 10; SPLIT_N <- 5; MI_N_IMPUTATIONS <- 2
}
#
# ======= END OF CONFIGURATION (nothing below this line needs editing) =========
#
#  2.12  Configuration validation ---------------------------------------------
#  WHAT    Checks the settings above against each other and stops on the first
#          inconsistency, naming the setting at fault.
#  WHY     Most configuration errors are silent rather than fatal: a variable
#          declared nominal but absent from the data, a scale whose target is
#          not in CLUSTER_VARS, an asymmetric binary coded 1/2 instead of 0/1.
#          Each produces a run that completes and is wrong. Catching them here
#          costs nothing; catching them in the manuscript costs a correction.
#  STATUS  REQUIRED

validate_config <- function() {
  errs <- character(0)
  add  <- function(...) errs <<- c(errs, paste0(...))

  nm <- spec_names(CLUSTER_VARS)
  if (length(nm) < 2)          add("CLUSTER_VARS: at least two variables are needed.")
  if (anyDuplicated(nm))       add("CLUSTER_VARS: duplicated variable name(s): ",
                                   paste(unique(nm[duplicated(nm)]), collapse = ", "))
  if (any(spec_weights(CLUSTER_VARS) <= 0))
    add("CLUSTER_VARS: weights must be strictly positive.")

  for (s in SCALE_DEFS) {
    if (!s$target %in% nm)
      add("SCALE_DEFS: target '", s$target,
          "' is not declared in CLUSTER_VARS.")
    if (length(s$items) < 2)
      add("SCALE_DEFS('", s$target, "'): needs at least two items.")
    if (!all(s$reverse %in% s$items))
      add("SCALE_DEFS('", s$target, "'): reverse items not in the item list.")
    if (s$max_missing < 0 || s$max_missing >= 1)
      add("SCALE_DEFS('", s$target, "'): max_missing must be in [0, 1).")
  }

  if (!MISSING_METHOD %in% c("complete_cases", "layers", "imputation"))
    add("MISSING_METHOD: must be 'complete_cases', 'layers' or 'imputation'.")
  if (MISSING_METHOD == "layers" && length(SCALE_DEFS) == 0)
    add("MISSING_METHOD = 'layers' requires at least one entry in SCALE_DEFS. ",
        "With no item-level information there is nothing to prorate; use ",
        "'complete_cases' or 'imputation' instead.")
  if (MISSING_METHOD == "imputation") {
    if (!MI_MODE %in% c("sensitivity", "consensus"))
      add("MI_MODE: must be 'sensitivity' or 'consensus'.")
    if (MI_N_IMPUTATIONS < 2)
      add("MI_N_IMPUTATIONS: at least 2 imputations are needed.")
    if (MI_MODE == "sensitivity" &&
        (MI_REFERENCE < 1 || MI_REFERENCE > MI_N_IMPUTATIONS))
      add("MI_REFERENCE: must be between 1 and MI_N_IMPUTATIONS.")
  }

  if (min(K_RANGE) < 2) add("K_RANGE: must start at 2 or higher.")
  if (!is.null(K_CANDIDATES) && !all(K_CANDIDATES %in% K_RANGE))
    add("K_CANDIDATES: every candidate must also appear in K_RANGE.")
  if (!is.null(K_FINAL)) {
    if (!K_FINAL %in% K_RANGE)
      add("K_FINAL: must appear in K_RANGE.")
    if (!is.null(CLUSTER_LABELS) && length(CLUSTER_LABELS) != K_FINAL)
      add("CLUSTER_LABELS: must have exactly K_FINAL entries.")
    if (!nzchar(trimws(K_FINAL_RATIONALE)))
      add("K_FINAL_RATIONALE: must not be empty when K_FINAL is set.")
  }

  ext <- spec_names(EXTERNAL_VARS)
  overlap <- intersect(ext, nm)
  if (length(overlap) > 0)
    add("EXTERNAL_VARS: ", paste(overlap, collapse = ", "),
        " also appear(s) in CLUSTER_VARS. A variable cannot both build the ",
        "clusters and validate them.")

  if (length(errs) > 0) {
    stop("\n  Configuration problems:\n",
         paste0("   - ", errs, collapse = "\n"), "\n", call. = FALSE)
  }
  invisible(TRUE)
}

banner("PAM-GOWER PIPELINE")
validate_config()
ok("Configuration is internally consistent.")
if (isTRUE(QUICK_RUN))
  warn("QUICK_RUN is on. Resampling parameters are reduced; do not report ",
       "these results.")

set.seed(SEED)

# Output folders.
dir_tables  <- file.path(OUTPUT_DIR, "tables")
dir_figures <- file.path(OUTPUT_DIR, "figures")

if (isTRUE(CLEAR_OUTPUTS) && dir.exists(OUTPUT_DIR)) {
  unlink(list.files(OUTPUT_DIR, full.names = TRUE, recursive = TRUE))
}
for (d in c(OUTPUT_DIR, dir_tables, dir_figures)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Every table and figure in the run passes through these two helpers, so the
# inventory in section 9.2 is complete by construction rather than by hand.
OUTPUT_LOG <- new.env(parent = emptyenv())
OUTPUT_LOG$files <- character(0)

save_table <- function(x, filename, description = "") {
  path <- file.path(dir_tables, filename)
  readr::write_csv(as.data.frame(x), path)
  OUTPUT_LOG$files <- c(OUTPUT_LOG$files,
                        paste(file.path("tables", filename), description,
                              sep = "\t"))
  invisible(path)
}

save_figure <- function(plot_obj, filename, description = "",
                        width = 8, height = 5.5, dpi = 300) {
  path <- file.path(dir_figures, filename)
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
  OUTPUT_LOG$files <- c(OUTPUT_LOG$files,
                        paste(file.path("figures", filename), description,
                              sep = "\t"))
  invisible(path)
}

save_base_figure <- function(expr, filename, description = "",
                             width = 2400, height = 1650, res = 240) {
  path <- file.path(dir_figures, filename)
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  OUTPUT_LOG$files <- c(OUTPUT_LOG$files,
                        paste(file.path("figures", filename), description,
                              sep = "\t"))
  invisible(path)
}


#  3.  IMPORT AND STRUCTURAL CHECKS ===========================================

banner("3.  IMPORT AND STRUCTURAL CHECKS")

#  3.1  Import -----------------------------------------------------------------
#  WHAT    Reads the dataset, dispatching on the file extension.
#  WHY     Nothing conceptual; it is here so that changing file format is a
#          one-line change in CONFIG rather than an edit to the pipeline.
#  CHOICE  DATA_PATH, DATA_OBJECT (2.1)
#  STATUS  REQUIRED
#
#  A note on SPSS files: haven returns labelled vectors, which carry value
#  labels that daisy() does not understand. Labels are stripped on import and
#  the measurement type is taken from CLUSTER_VARS instead, which is where it
#  belongs.

import_data <- function(path, object_name = NULL) {

  if (identical(path, "example")) {
    gen <- "make_example_data.R"
    if (!file.exists(gen)) {
      stop("\n  DATA_PATH is \"example\" but make_example_data.R was not found",
           "\n  in the working directory. Either restore that file, or set",
           "\n  DATA_PATH to your own dataset in section 2.1.", call. = FALSE)
    }
    source(gen, local = TRUE)
    note("Using the synthetic example dataset (make_example_data.R).")
    return(make_example_data(seed = SEED))
  }

  if (!file.exists(path)) {
    stop("\n  Data file not found: ", path,
         "\n  Paths are relative to the project root. If you opened the .R ",
         "file\n  directly rather than the .Rproj file, the working directory ",
         "is\n  probably wrong. Current working directory:\n    ", getwd(),
         call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  out <- switch(
    ext,
    "rdata" = ,
    "rda"   = {
      env    <- new.env()
      loaded <- load(path, envir = env)
      if (is.null(object_name)) {
        if (length(loaded) != 1L)
          stop("\n  ", path, " contains ", length(loaded), " objects (",
               paste(loaded, collapse = ", "), ").",
               "\n  Name the one you want in DATA_OBJECT (section 2.1).",
               call. = FALSE)
        object_name <- loaded
      }
      get(object_name, envir = env)
    },
    "rds" = readRDS(path),
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "tsv" = readr::read_tsv(path, show_col_types = FALSE),
    "sav" = { need_pkg("haven", "reading SPSS files")
              haven::zap_labels(haven::read_sav(path)) },
    "xlsx" = ,
    "xls"  = { need_pkg("readxl", "reading Excel files")
               readxl::read_excel(path) },
    stop("\n  Unsupported file extension: '", ext, "'.",
         "\n  Supported: .RData .rda .rds .csv .tsv .sav .xlsx .xls",
         call. = FALSE)
  )
  as.data.frame(out, stringsAsFactors = FALSE)
}

data_raw <- import_data(DATA_PATH, DATA_OBJECT)
ok("Imported ", nrow(data_raw), " cases and ", ncol(data_raw), " columns.")

#  3.2  Required columns and identifier ----------------------------------------
#  WHAT    Confirms that every column the configuration refers to exists, and
#          that the identifier is usable.
#  WHY     A misspelled variable name is the most common configuration error
#          and the least visible: dplyr would drop it, the clustering would run
#          on five variables instead of six, and nothing would look wrong.
#  STATUS  REQUIRED

scale_targets <- if (length(SCALE_DEFS)) {
  vapply(SCALE_DEFS, function(s) s$target, character(1))
} else character(0)
scale_items   <- unique(unlist(lapply(SCALE_DEFS, function(s) s$items)))
direct_vars   <- setdiff(spec_names(CLUSTER_VARS), scale_targets)
external_vars <- spec_names(EXTERNAL_VARS)

# Columns that must be present in the imported file.
required_cols <- unique(c(ID_VAR, direct_vars, scale_items, external_vars,
                          if (MISSING_METHOD == "imputation") MI_AUXILIARY_VARS))

missing_cols <- setdiff(required_cols, names(data_raw))
if (length(missing_cols) > 0) {
  stop("\n  These columns are named in the configuration but absent from the ",
       "data:\n    ", paste(missing_cols, collapse = ", "),
       "\n  Check for typographic differences against names(data_raw).",
       call. = FALSE)
}
ok("All ", length(required_cols), " configured columns are present.")

# The identifier is converted to character: numeric IDs are otherwise liable
# to be reformatted in exports, and IDs are only ever matched, never compared.
data_raw[[ID_VAR]] <- as.character(data_raw[[ID_VAR]])
if (anyDuplicated(data_raw[[ID_VAR]])) {
  stop("\n  ID_VAR ('", ID_VAR, "') is not unique: ",
       sum(duplicated(data_raw[[ID_VAR]])), " duplicate value(s).",
       "\n  Every exclusion log and every export is keyed to it.", call. = FALSE)
}
if (anyNA(data_raw[[ID_VAR]]))
  stop("\n  ID_VAR ('", ID_VAR, "') contains missing values.", call. = FALSE)
ok("Identifier '", ID_VAR, "' is unique and complete.")

#  3.3  Measurement-type audit -------------------------------------------------
#  WHAT    Compares the type declared for each variable with what the column
#          actually contains, and reports anything that cannot be reconciled.
#  WHY     Three failure modes, all silent:
#          - a variable declared nominal whose codes are not in 'levels' loses
#            those cases to NA during coercion;
#          - a binary variable that turns out to have three values is not
#            binary, and Gower will not say so;
#          - a constant variable contributes zero to every distance while
#            still counting in the denominator, diluting every other variable.
#  STATUS  REQUIRED

audit_types <- function(data, spec) {
  purrr::map_dfr(spec, function(v) {
    x  <- data[[v$name]]
    xn <- x[!is.na(x)]
    u  <- unique(xn)
    tibble::tibble(
      variable       = v$name,
      declared_type  = v$type,
      storage        = class(x)[1],
      n_distinct     = length(u),
      n_missing      = sum(is.na(x)),
      pct_missing    = round(100 * mean(is.na(x)), 2),
      observed_range = if (is.numeric(xn) && length(xn))
                         paste0(signif(min(xn), 4), " to ", signif(max(xn), 4))
                       else paste(utils::head(sort(as.character(u)), 6),
                                  collapse = ", "),
      declared_levels = if (is.null(v$levels)) ""
                        else paste(v$levels, collapse = ", "),
      weight          = v$weight,
      note            = gsub("\\s+", " ", v$note)
    )
  })
}

type_audit <- audit_types(data_raw, CLUSTER_VARS)
print(as.data.frame(type_audit[, 1:7]), row.names = FALSE)

problems <- character(0)
for (v in CLUSTER_VARS) {
  if (v$name %in% scale_targets) next   # computed later, in 4.4
  x  <- data_raw[[v$name]]
  xn <- x[!is.na(x)]
  u  <- unique(xn)

  if (length(u) < 2)
    problems <- c(problems, paste0("'", v$name, "' is constant; it cannot ",
                                   "distinguish cases and should be removed."))

  if (v$type %in% c("symm_binary", "asymm_binary") && length(u) > 2)
    problems <- c(problems, paste0("'", v$name, "' is declared binary but has ",
                                   length(u), " distinct values."))

  if (v$type %in% c("nominal", "ordinal", "symm_binary", "asymm_binary") &&
      !is.null(v$levels)) {
    unmatched <- setdiff(as.character(u), as.character(v$levels))
    if (length(unmatched) > 0)
      problems <- c(problems, paste0("'", v$name, "': value(s) ",
        paste(utils::head(unmatched, 5), collapse = ", "),
        " are not in the declared levels and would become NA on coercion."))
  }

  if (v$type == "interval" && !is.numeric(x))
    problems <- c(problems, paste0("'", v$name, "' is declared interval but ",
                                   "stored as ", class(x)[1], "."))
}

if (length(problems) > 0) {
  stop("\n  Measurement-type problems:\n",
       paste0("   - ", problems, collapse = "\n"),
       "\n  Fix the data or the declaration in CLUSTER_VARS (2.3).",
       call. = FALSE)
}
ok("Declared types are consistent with the data.")

#  3.4  Codebook ---------------------------------------------------------------
#  WHAT    Exports the variable specification as a table.
#  WHY     This is the first thing a reader of your deposit needs, and the
#          first thing you will want when you return to the analysis in a
#          year. It is generated rather than written so it cannot drift out of
#          step with the configuration.
#  STATUS  REQUIRED

save_table(type_audit, "codebook_clustering_variables.csv",
           "Variables entering the clustering: declared type, range, notes")

if (length(EXTERNAL_VARS) > 0)
  save_table(audit_types(data_raw, EXTERNAL_VARS),
             "codebook_external_variables.csv",
             "Validation variables, not used to build the clusters")


#  4.  MISSING DATA ===========================================================
#
#  WHY THIS SECTION IS LONG :
#  PAM needs a complete dissimilarity matrix. daisy() will compute Gower
#  distances in the presence of NAs by rescaling over the variables a given
#  PAIR happens to share, so the function returns something rather than an
#  error, but two cases compared on four variables and two cases compared on
#  six are then placed on the same scale, and a case with many omissions ends
#  up closer to everything. The resulting partition partly reflects the
#  missingness pattern rather than the constructs.
#
#  So the data entering section 5 must be complete on the clustering
#  variables. The only question is how they got that way, and that is an
#  analytic decision with consequences for who your sample represents. 

banner("4.  MISSING DATA")

attrition <- tibble::tibble(step = character(), criterion = character(),
                            n_before = integer(), n_excluded = integer(),
                            n_after = integer())

log_step <- function(label, criterion, n_before, n_after) {
  attrition <<- dplyr::bind_rows(attrition, tibble::tibble(
    step = label, criterion = criterion,
    n_before = as.integer(n_before),
    n_excluded = as.integer(n_before - n_after),
    n_after = as.integer(n_after)))
}

excluded_ids <- list()

#  4.1  Eligibility screening -------------------------------------------------
#  WHAT    Applies ELIGIBILITY_FILTER and records who it removed.
#  WHY     See the note at 2.5: eligibility must precede every missingness
#          rule, or the MCAR diagnostic below is computed on a sample already
#          purged of its most incomplete cases.
#  CHOICE  ELIGIBILITY_FILTER (2.5)
#  STATUS  OPTIONAL - skipped when the filter is NULL

step("4.1  Eligibility")

data_eligible <- data_raw
if (!is.null(ELIGIBILITY_FILTER)) {
  keep <- eval(ELIGIBILITY_FILTER, envir = data_raw)
  if (!is.logical(keep))
    stop("ELIGIBILITY_FILTER must evaluate to a logical vector.", call. = FALSE)
  keep[is.na(keep)] <- FALSE           # a missing criterion is not eligibility
  excluded_ids[["Layer 1 (eligibility)"]] <- data_raw[[ID_VAR]][!keep]
  data_eligible <- data_raw[keep, , drop = FALSE]
  log_step("Layer 1", paste("Eligibility:",
           paste(deparse(ELIGIBILITY_FILTER), collapse = " ")),
           nrow(data_raw), nrow(data_eligible))
  ok("Eligible: ", nrow(data_eligible), " of ", nrow(data_raw), " cases.")
} else {
  log_step("Layer 1", "No eligibility filter applied",
           nrow(data_raw), nrow(data_raw))
  note("No eligibility filter (ELIGIBILITY_FILTER is NULL).")
}
if (nrow(data_eligible) < 20)
  stop("Fewer than 20 eligible cases remain. Check ELIGIBILITY_FILTER.",
       call. = FALSE)

# The columns on which missingness is assessed: the items of every declared
# scale, plus the clustering variables that are read directly from the file.
missing_scope <- unique(c(direct_vars, scale_items))

#  4.2  Extent, pattern and mechanism of missingness --------------------------
#  WHAT    Missingness by variable, by case and overall; the pattern of
#          co-occurrence; and Little's test of the MCAR hypothesis.
#  WHY     The choice between deletion, proration and imputation depends on
#          how much is missing and on whether the missingness is plausibly
#          unrelated to the values themselves. Making that choice without
#          looking is guessing.
#
#          On Little's test, two cautions that matter here. It assumes
#          multivariate normality of the variables it examines, which ordinal
#          and binary items do not satisfy; treat a result on such items as
#          indicative only. And it is a test of a null hypothesis: a
#          non-significant result is weak evidence of MCAR, not proof of it,
#          and with a large n even trivial departures will reject.
#  CHOICE  RUN_MCAR_TEST (2.11)
#  STATUS  OPTIONAL, but STRONGLY recommended

step("4.2  Extent and mechanism of missingness")

miss_block <- data_eligible[, missing_scope, drop = FALSE]

missing_by_variable <- tibble::tibble(
  variable    = missing_scope,
  n_missing   = vapply(miss_block, function(x) sum(is.na(x)), integer(1)),
  pct_missing = round(100 * vapply(miss_block, function(x) mean(is.na(x)),
                                   numeric(1)), 2)
) %>% dplyr::arrange(dplyr::desc(pct_missing))

missing_by_case <- tibble::tibble(
  id          = data_eligible[[ID_VAR]],
  n_missing   = as.integer(rowSums(is.na(miss_block))),
  pct_missing = round(100 * rowMeans(is.na(miss_block)), 2)
)

overall_missing <- tibble::tibble(
  n_cases            = nrow(miss_block),
  n_variables        = ncol(miss_block),
  n_cells            = nrow(miss_block) * ncol(miss_block),
  n_missing_cells    = sum(is.na(miss_block)),
  pct_missing_cells  = round(100 * mean(is.na(miss_block)), 3),
  n_complete_cases   = sum(stats::complete.cases(miss_block)),
  pct_complete_cases = round(100 * mean(stats::complete.cases(miss_block)), 2)
)

print(as.data.frame(overall_missing), row.names = FALSE)
cat("\n  Variables with any missing data:\n")
print(as.data.frame(dplyr::filter(missing_by_variable, n_missing > 0)),
      row.names = FALSE)

save_table(missing_by_variable, "missingness_by_variable.csv",
           "Missing values per variable, eligible sample")
save_table(missing_by_case, "missingness_by_case.csv",
           "Missing values per case, eligible sample")
save_table(overall_missing, "missingness_overall.csv",
           "Overall missingness, eligible sample")

# The pattern matters as much as the amount: 5% spread evenly across cases is
# a different problem from 5% concentrated in 5% of cases.
pattern_counts <- miss_block %>%
  is.na() %>% as.data.frame() %>%
  dplyr::group_by(dplyr::across(dplyr::everything())) %>%
  dplyr::tally(name = "n_cases") %>%
  dplyr::arrange(dplyr::desc(n_cases)) %>%
  utils::head(25)
save_table(pattern_counts, "missingness_patterns.csv",
           "25 most frequent patterns of missingness (TRUE = missing)")

mcar_result <- NULL
if (isTRUE(RUN_MCAR_TEST) && !requireNamespace("naniar", quietly = TRUE)) {
  warn("'naniar' is not installed, so the MCAR diagnostic was skipped. ",
       "Install it with install.packages(\"naniar\"), or set ",
       "RUN_MCAR_TEST <- FALSE in section 2.11 to stop asking for it.")
  RUN_MCAR_TEST <- FALSE
}
if (isTRUE(RUN_MCAR_TEST)) {
  numeric_block <- as.data.frame(lapply(miss_block, function(x)
    if (is.numeric(x)) x else as.numeric(as.factor(x))))
  mcar_result <- tryCatch(
    naniar::mcar_test(numeric_block),
    error = function(e) {
      warn("Little's MCAR test could not be computed: ", conditionMessage(e))
      warn("This is common with ordinal or binary items, which produce a ",
           "rank-deficient covariance matrix. It is not a failure of the run.")
      NULL
    })
  if (!is.null(mcar_result)) {
    cat("\n  Little's MCAR test:\n")
    print(as.data.frame(mcar_result), row.names = FALSE)
    save_table(mcar_result, "mcar_test_little.csv",
               "Little's test of the MCAR hypothesis, eligible sample")
    if (mcar_result$p.value < .05)
      note("MCAR is rejected (p < .05). Listwise deletion is likely to bias ",
           "the sample; prefer 'layers' or 'imputation' in 2.6.")
    else
      note("MCAR is not rejected. This is weak evidence, not proof; see the ",
           "cautions in the section header.")
  }
}

#  4.3  The three routes to complete data -------------------------------------
#  WHAT    Explains the options. The code is in 4.5; read this first.
#  CHOICE  MISSING_METHOD (2.6)
#  STATUS  REQUIRED - one route must be chosen
#
#  ROUTE A - "complete_cases"
#  Delete every case with a missing value on any clustering variable.
#
#    For     It is transparent, needs no assumptions beyond MCAR, adds no
#            model-based uncertainty and is trivially reproducible.
#    Against Under anything other than MCAR it biases the sample, and the bias
#            is towards the cases who answer everything. In clustering this is
#            worse than in regression: the cases who omit items are often
#            precisely the ones who would populate the smaller, more distinct
#            clusters, so deletion can remove the structure you are looking
#            for.
#    Use it  when the loss is small (a few percent), when MCAR is not
#            rejected, and when you report the loss.
#
#  ROUTE B - "layers" (exclusion in ordered layers, then proration)
#  Successive rules of decreasing severity, then scale scores computed from
#  the items each surviving case actually answered.
#
#    For     It separates "this case told us almost nothing" from "this case
#            skipped one item of twenty-four", which listwise deletion treats
#            identically. It typically retains substantially more cases than
#            route A while adding no model, no seed and no uncertainty: a
#            prorated score is an arithmetic mean over answered items, and it
#            reduces exactly to the ordinary score when nothing is missing.
#    Against Proration assumes the answered items are a fair sample of the
#            scale's content - reasonable for a homogeneous scale with one
#            item missing, unreasonable for a heterogeneous one with a third
#            missing. The tolerance you set is doing real work, so set it low.
#    Use it  when your clustering variables are scale scores and you have the
#            item-level data. Requires SCALE_DEFS.
#
#    The layers, in order:
#      Layer 2  overall: a case missing more than OVERALL_MISSING_MAX of all
#               the relevant columns is removed. Catches abandonment.
#      Layer 3  zero tolerance: any missing value on a variable read directly
#               from the file, or on a scale whose max_missing is 0.
#      Layer 4  proportional: a scale with max_missing > 0 may lose up to that
#               proportion of its items; beyond it, the case is removed.
#
#    Because the layers are intersecting criteria, the FINAL n does not depend
#    on the order in which they are applied. Only the intermediate counts do,
#    which is why the attrition table reports them as a cascade rather than as
#    independent exclusions.
#
#  ROUTE C - "imputation" (multiple imputation, then clustering)
#  Impute m completed datasets by chained equations, then cluster.
#
#    For     Under MAR it is the statistically principled option, it uses all
#            observed information, and good auxiliary variables make the MAR
#            assumption more plausible than the MCAR assumption route A needs.
#    Against Multiple imputation works by pooling ESTIMATES across imputations 
#            using Rubin's rules. A cluster partition is not an estimate: 
#            cluster 3in imputation 1 and cluster 3 in imputation 2 are not the 
#            same cluster, they are not even guaranteed to correspond, and there 
#            is no accepted rule for averaging partitions. Software that offers
#            "clustering with multiple imputation" usually clusters one
#            completed dataset and says nothing about the rest.
#
#    What this pipeline does instead, and why it is defensible:
#
#      MI_MODE = "sensitivity"  (recommended)
#        The reported analysis runs on ONE completed dataset (MI_REFERENCE).
#        The other imputations are used to answer a question you can actually
#        answer: how much does the partition depend on the imputation? Each
#        imputation is clustered separately at the retained k and the
#        partitions are compared pairwise by adjusted Rand index, and
#        case by case by how often two cases land together. A high agreement
#        means the imputation is not driving the typology; a low one means the
#        typology is partly an artefact of the imputation model and should be
#        reported as such. This is section 7.9.
#
#      MI_MODE = "consensus"
#        All m partitions are combined into a co-membership matrix -- the
#        proportion of imputations in which each pair of cases falls in the
#        same cluster -- and the final partition is obtained by clustering
#        1 minus that matrix. It uses all m datasets, at the cost that the
#        result is no longer a PAM solution: it has no medoids, and the
#        silhouette widths are computed in a different space from the one the
#        model selection used.
#
#    Use it  when missingness is substantial, plausibly MAR, and you have
#            auxiliary variables worth conditioning on. Report the sensitivity
#            analysis either way.
#
#  WHAT THIS PIPELINE WILL NOT DO
#  It will not impute a variable it is about to cluster on using the cluster
#  membership, and it will not carry NAs into daisy(). Both produce numbers.
#  Neither produces evidence.

#  4.4  Scoring rule and proration ---------------------------------------------
#  WHAT    Computes each declared scale score from its items.
#  WHY     Putting the scoring rule in the script makes it reproducible and
#          auditable, and it is what makes proration possible at all. The
#          alternative (trusting a score column computed years ago in other
#          software) is exactly the situation 4.6 exists to detect.
#  CHOICE  SCALE_DEFS (2.4): items, sum or mean, tolerance, reverse-scoring
#  STATUS  OPTIONAL - no-op when SCALE_DEFS is empty
#
#  ON REVERSE-SCORING. The commonest silent error in this whole pipeline is
#  reversing items that the source file already stores reversed. The result is
#  a scale that runs the right way for some items and the wrong way for
#  others; nothing errors, every descriptive statistic looks plausible, and
#  the clustering is wrong. Leave 'reverse' empty unless you are certain the
#  file holds raw responses, and let 4.6 confirm it.

compute_scale_scores <- function(data, defs, prorate = TRUE) {
  for (s in defs) {
    m <- as.matrix(data[, s$items, drop = FALSE])
    storage.mode(m) <- "double"

    if (length(s$reverse) > 0) {
      j <- match(s$reverse, s$items)
      m[, j] <- (s$item_min + s$item_max) - m[, j]
    }

    prop_missing <- rowMeans(is.na(m))
    item_mean    <- rowMeans(m, na.rm = TRUE)
    value        <- if (s$score == "sum") item_mean * length(s$items) else item_mean

    tolerance <- if (isTRUE(prorate)) s$max_missing else 0
    value[prop_missing > tolerance] <- NA_real_
    value[prop_missing == 1]        <- NA_real_

    data[[s$target]] <- value
  }
  data
}

#  4.5  Route execution -------------------------------------------------------
#  WHAT    Runs the route selected in 2.6 and produces data_analytic: one row
#          per retained case, complete on every clustering variable.
#  STATUS  REQUIRED

step(paste0("4.5  Applying missing-data route: ", MISSING_METHOD))

cluster_var_names <- spec_names(CLUSTER_VARS)
carry_cols <- unique(c(ID_VAR, cluster_var_names, external_vars))
imputed_frames <- NULL

if (MISSING_METHOD == "complete_cases") {

  scored <- compute_scale_scores(data_eligible, SCALE_DEFS, prorate = FALSE)
  n0 <- nrow(scored)
  cc <- stats::complete.cases(scored[, cluster_var_names, drop = FALSE])
  excluded_ids[["Layer 2 (listwise deletion)"]] <- scored[[ID_VAR]][!cc]
  data_analytic <- scored[cc, carry_cols, drop = FALSE]
  log_step("Layer 2", "Any missing value on a clustering variable",
           n0, nrow(data_analytic))
  ok("Listwise deletion retained ", nrow(data_analytic), " of ", n0, " cases.")

} else if (MISSING_METHOD == "layers") {

  current <- data_eligible
  layer   <- 2L

  ## --- Layer 2: overall person-level missingness ----------------------------
  if (!is.null(OVERALL_MISSING_MAX)) {
    n0   <- nrow(current)
    pct  <- rowMeans(is.na(current[, missing_scope, drop = FALSE]))
    keep <- pct <= OVERALL_MISSING_MAX     # strict: exactly at the cap is kept
    excluded_ids[[paste0("Layer ", layer, " (overall missingness)")]] <-
      current[[ID_VAR]][!keep]
    current <- current[keep, , drop = FALSE]
    log_step(paste("Layer", layer),
             paste0("More than ", OVERALL_MISSING_MAX * 100,
                    "% missing across all relevant columns"),
             n0, nrow(current))
    layer <- layer + 1L
  }

  ## --- Layer 3: zero tolerance ---------------------------------------------
  strict_cols <- c(
    direct_vars,
    unlist(lapply(SCALE_DEFS, function(s) if (s$max_missing == 0) s$items))
  )
  strict_cols <- unique(strict_cols[!is.na(strict_cols)])

  if (length(strict_cols) > 0) {
    n0   <- nrow(current)
    keep <- rowSums(is.na(current[, strict_cols, drop = FALSE])) == 0
    excluded_ids[[paste0("Layer ", layer, " (zero tolerance)")]] <-
      current[[ID_VAR]][!keep]
    current <- current[keep, , drop = FALSE]
    log_step(paste("Layer", layer),
             "Any missing value on a variable or scale admitting none",
             n0, nrow(current))
    layer <- layer + 1L
  }

  ## --- Layer 4: proportional, per tolerant scale ---------------------------
  tolerant <- Filter(function(s) s$max_missing > 0, SCALE_DEFS)
  if (length(tolerant) > 0) {
    n0 <- nrow(current)
    over <- vapply(tolerant, function(s)
      rowMeans(is.na(current[, s$items, drop = FALSE])) > s$max_missing,
      logical(nrow(current)))
    keep <- rowSums(matrix(over, nrow = nrow(current))) == 0
    excluded_ids[[paste0("Layer ", layer, " (per-scale tolerance)")]] <-
      current[[ID_VAR]][!keep]
    current <- current[keep, , drop = FALSE]
    log_step(paste("Layer", layer),
             paste0("Exceeded the per-scale tolerance on ",
                    paste(vapply(tolerant, function(s) s$target, character(1)),
                          collapse = " or ")),
             n0, nrow(current))
  }

  ## --- Proration ------------------------------------------------------------
  scored        <- compute_scale_scores(current, SCALE_DEFS, prorate = TRUE)
  data_analytic <- scored[, carry_cols, drop = FALSE]

  n_prorated <- sum(vapply(SCALE_DEFS, function(s)
    sum(rowSums(is.na(current[, s$items, drop = FALSE])) > 0), integer(1)))
  ok("Layered exclusions retained ", nrow(data_analytic), " of ",
     nrow(data_eligible), " eligible cases.")
  note(n_prorated, " scale score(s) were prorated from partially answered ",
       "scales. Under listwise deletion these cases would have been lost.")

} else if (MISSING_METHOD == "imputation") {

  need_pkg("mice", "multiple imputation (section 4.3, route C)")

  mi_vars  <- unique(c(missing_scope, MI_AUXILIARY_VARS))
  mi_frame <- data_eligible[, mi_vars, drop = FALSE]

  # Types are set before imputation, not after: mice chooses its elementary
  # imputation model from the class of each column (predictive mean matching
  # for numeric, logistic for binary, polytomous for nominal, proportional
  # odds for ordered). Handing it numeric codes for a nominal variable gets
  # you a linear model on arbitrary integers.
  for (v in CLUSTER_VARS) {
    if (v$name %in% direct_vars)
      mi_frame[[v$name]] <- coerce_gower_column(mi_frame[[v$name]], v)
  }

  note("Imputing ", MI_N_IMPUTATIONS, " datasets over ", ncol(mi_frame),
       " columns. This can take several minutes.")
  imp <- mice::mice(mi_frame, m = MI_N_IMPUTATIONS, printFlag = FALSE,
                    seed = SEED)

  if (length(imp$loggedEvents) > 0 && nrow(imp$loggedEvents) > 0) {
    warn("mice logged ", nrow(imp$loggedEvents), " event(s) -- usually ",
         "constant or collinear predictors that were dropped.")
    save_table(imp$loggedEvents, "imputation_logged_events.csv",
               "Events logged by mice during imputation")
  }

  imputed_frames <- lapply(seq_len(MI_N_IMPUTATIONS), function(i) {
    cd <- mice::complete(imp, i)
    cd[[ID_VAR]] <- data_eligible[[ID_VAR]]
    cd <- compute_scale_scores(cd, SCALE_DEFS, prorate = FALSE)
    cd[, unique(c(ID_VAR, cluster_var_names)), drop = FALSE]
  })

  data_analytic <- imputed_frames[[MI_REFERENCE]]
  if (length(external_vars) > 0) {
    data_analytic <- dplyr::left_join(
      data_analytic,
      data_eligible[, unique(c(ID_VAR, external_vars)), drop = FALSE],
      by = ID_VAR)
  }
  log_step("Layer 2", paste0("Multiple imputation (m = ", MI_N_IMPUTATIONS,
                             "); no cases deleted"),
           nrow(data_eligible), nrow(data_analytic))
  ok("Imputation complete. Reference dataset: imputation ", MI_REFERENCE,
     " of ", MI_N_IMPUTATIONS, ".")
}

#  4.6  Cross-check of recomputed scores ---------------------------------------
#  WHAT    Where a scale score of the same name already existed in the file,
#          compares it with the score recomputed here, on cases with no
#          missing item.
#  WHY     On a complete case, proration reduces to the ordinary score, so the
#          two figures MUST agree exactly. A non-zero difference means the
#          scoring rule in SCALE_DEFS is not the rule that produced the file:
#          a sum where the original took a mean, or -- far more often --
#          reverse-scoring applied a second time. This check is the reason to
#          recompute scores rather than trust them.
#  STATUS  OPTIONAL -- runs only when a stored score is available to compare

if (length(SCALE_DEFS) > 0 && MISSING_METHOD != "imputation") {
  step("4.6  Score cross-check")
  rows <- list()
  for (s in SCALE_DEFS) {
    stored_present <- s$target %in% names(data_eligible) &&
                      is.numeric(data_eligible[[s$target]])
    if (!stored_present) next
    recomputed <- compute_scale_scores(data_eligible, list(s), prorate = FALSE)
    complete_i <- rowSums(is.na(data_eligible[, s$items, drop = FALSE])) == 0
    d <- abs(recomputed[[s$target]][complete_i] -
             data_eligible[[s$target]][complete_i])
    rows[[s$target]] <- tibble::tibble(
      scale            = s$target,
      n_complete_cases = sum(complete_i),
      max_abs_diff     = round(max(d, na.rm = TRUE), 6),
      n_discrepant     = sum(d > 1e-8, na.rm = TRUE))
  }
  if (length(rows) > 0) {
    score_check <- dplyr::bind_rows(rows)
    print(as.data.frame(score_check), row.names = FALSE)
    save_table(score_check, "score_recomputation_check.csv",
               "Recomputed vs stored scale scores, complete cases only")
    if (any(score_check$n_discrepant > 0)) {
      warn("Recomputed scores differ from the stored ones on complete cases.")
      warn("Check: sum vs mean, and whether reverse-scoring is being applied ",
           "twice (see the note in 4.4). Do not proceed until this resolves.")
    } else {
      ok("Recomputed scores reproduce the stored scores exactly.")
    }
  } else {
    note("No stored scale score of the same name was available to compare.")
  }
}

#  4.7  Attrition log and excluded identifiers ---------------------------------
#  WHAT    Exports the cascade of exclusions, and the identifiers removed at
#          each step.
#  WHY     This table belongs in your supplementary materials whether or not a
#          journal asks for it: it is what lets a reader judge whom the
#          typology describes. Exporting the identifiers separately lets you
#          answer the reviewer question that follows -- "did the excluded
#          cases differ from the retained ones?" -- without rerunning
#          anything.
#  STATUS  REQUIRED

step("4.7  Attrition")
print(as.data.frame(attrition), row.names = FALSE)
save_table(attrition, "attrition_log.csv",
           "Cascade of exclusions from import to analytic sample")

for (nm in names(excluded_ids)) {
  if (length(excluded_ids[[nm]]) == 0) next
  fname <- paste0("excluded_ids_",
                  gsub("[^a-z0-9]+", "_", tolower(nm)), ".csv")
  save_table(tibble::tibble(id = excluded_ids[[nm]]), fname,
             paste("Identifiers excluded at:", nm))
}

#  4.8  Final completeness check -----------------------------------------------
#  WHAT    Refuses to continue if any clustering variable still contains NA.
#  WHY     This should be impossible by construction. If it fires, a threshold
#          or an item list was changed without updating the corresponding
#          layer, and the right response is to find that change -- not to add
#          an imputation step here to make the error go away.
#  STATUS  REQUIRED

residual <- colSums(is.na(data_analytic[, cluster_var_names, drop = FALSE]))
if (sum(residual) > 0) {
  print(residual[residual > 0])
  stop("\n  Missing values remain in the clustering variables after route '",
       MISSING_METHOD, "'.\n  This should not be possible. Trace the ",
       "configuration change that caused it.", call. = FALSE)
}

if (nrow(data_analytic) < 5 * max(K_RANGE)) {
  warn("Only ", nrow(data_analytic), " cases for up to ", max(K_RANGE),
       " clusters. Small samples make every index below unstable; consider ",
       "lowering the top of K_RANGE.")
}

# Which retained cases would ALSO have survived listwise deletion. Used in 7.8
# to show what the choice of route cost or bought.
complete_case_ids <- {
  s <- compute_scale_scores(data_eligible, SCALE_DEFS, prorate = FALSE)
  s[[ID_VAR]][stats::complete.cases(s[, cluster_var_names, drop = FALSE])]
}

ok("Analytic sample: ", nrow(data_analytic), " cases, complete on all ",
   length(cluster_var_names), " clustering variables.")
save_table(tibble::tibble(id = data_analytic[[ID_VAR]]),
           "analytic_sample_ids.csv",
           "Identifiers of the cases entering the cluster analysis")


#  5.  BUILDING THE DISSIMILARITY MATRIX =======================================
#
#  Everything the clustering knows about your cases is in the matrix built
#  here. The algorithm never sees a variable again. That makes this section,
#  not section 6, the place where most of the analytic damage gets done, and
#  the place where it is cheapest to prevent.

banner("5.  BUILDING THE DISSIMILARITY")

#  5.1  Coercion to declared measurement types ---------------------------------
#  WHAT    Converts each clustering variable to the R class that makes daisy()
#          treat it as the type you declared.
#  WHY     daisy() infers measurement type from the CLASS of a column, not
#          from its content. Left alone, it will treat a nominal variable
#          coded 1 to 5 as interval -- so category 1 sits closer to category 2
#          than to category 5, an ordering that does not exist -- and a 0/1
#          indicator as interval rather than binary. Neither produces a
#          warning. Declaring the types in CONFIG and coercing here makes the
#          treatment explicit and gives 5.5 something to verify it against.
#  CHOICE  the 'type' of every entry in CLUSTER_VARS (2.3)
#  STATUS  REQUIRED

step("5.1  Coercion")

gower_frame <- prepare_gower_frame(data_analytic, CLUSTER_VARS)
daisy_type  <- build_daisy_type(CLUSTER_VARS)

cat("  Classes passed to daisy():\n")
print(data.frame(variable = names(gower_frame),
                 declared = spec_types(CLUSTER_VARS),
                 r_class  = vapply(gower_frame, function(x) class(x)[1],
                                   character(1))),
      row.names = FALSE)

# Asymmetric binary variables must be coded 0/1 with 1 = present, because the
# asymmetric coefficient counts 1-1 matches and ignores 0-0 matches. A 1/2
# coding silently inverts the meaning.
for (v in CLUSTER_VARS) {
  if (v$type != "asymm_binary") next
  x <- gower_frame[[v$name]]
  if (!all(stats::na.omit(x) %in% c(0, 1)))
    stop("'", v$name, "' is asymmetric binary but is not coded 0/1 after ",
         "coercion.\n  Declare levels = c(<absent>, <present>) in ",
         "CLUSTER_VARS.", call. = FALSE)
  if (mean(x == 1) > 0.5)
    note("'", v$name, "': the 'present' category holds ",
         round(100 * mean(x == 1)), "% of cases. Asymmetric treatment is ",
         "intended for the rarer category; consider symm_binary.")
}
ok("All ", ncol(gower_frame), " variables coerced to their declared types.")

#  5.2  Redundancy and implicit weighting --------------------------------------
#  WHAT    Reports how much of the distance each variable carries, and flags
#          variables that may be measuring the same thing twice.
#  WHY     Gower averages the per-variable dissimilarities. A construct
#          represented by four variables therefore contributes four times as
#          much to every distance as a construct represented by one. Nobody
#          decides this; it happens as a by-product of how many columns each
#          construct happened to produce, and it is one of the most common
#          unreported decisions in applied cluster analysis (Hennig & Liao,
#          2013).
#
#          The classic case is a single-choice question expanded into one-hot
#          indicators: five binary columns from one question means that
#          question counts five times and everything else counts once. Enter
#          it as ONE nominal variable instead.
#  CHOICE  the 'weight' argument of cvar(), and 5.4
#  STATUS  REQUIRED (diagnostic only - nothing is changed automatically)

step("5.2  Redundancy and implicit weighting")

w <- spec_weights(CLUSTER_VARS)
weight_share <- tibble::tibble(
  variable     = spec_names(CLUSTER_VARS),
  type         = spec_types(CLUSTER_VARS),
  weight       = w,
  share_of_distance = round(w / sum(w), 4)
)
print(as.data.frame(weight_share), row.names = FALSE)
save_table(weight_share, "gower_weight_share.csv",
           "Share of the Gower distance carried by each variable")

# Pairwise association among clustering variables. Spearman for interval and
# ordinal pairs; Cramer's V wherever a categorical variable is involved.
cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (any(dim(tab) < 2)) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

assoc_rows <- list()
vn <- spec_names(CLUSTER_VARS); vt <- spec_types(CLUSTER_VARS)
for (i in seq_along(vn)) for (j in seq_along(vn)) {
  if (j <= i) next
  xi <- gower_frame[[vn[i]]]; xj <- gower_frame[[vn[j]]]
  both_quant <- vt[i] %in% c("interval", "ordinal") &&
                vt[j] %in% c("interval", "ordinal")
  if (both_quant) {
    val <- suppressWarnings(stats::cor(as.numeric(xi), as.numeric(xj),
                                       method = "spearman"))
    measure <- "Spearman rho"
  } else {
    val <- cramers_v(as.factor(xi), as.factor(xj))
    measure <- "Cramer's V"
  }
  assoc_rows[[length(assoc_rows) + 1]] <- tibble::tibble(
    variable_1 = vn[i], variable_2 = vn[j],
    measure = measure, value = round(val, 3))
}
variable_association <- dplyr::bind_rows(assoc_rows) %>%
  dplyr::arrange(dplyr::desc(abs(value)))
save_table(variable_association, "variable_association.csv",
           "Pairwise association among the clustering variables")

flagged <- dplyr::filter(variable_association, abs(value) >= 0.85)
if (nrow(flagged) > 0) {
  warn("Variable pairs with |association| >= .85:")
  print(as.data.frame(flagged), row.names = FALSE)
  warn("Consider whether these are one construct entered twice. If they are, ",
       "the construct is silently double-weighted.")
} else {
  ok("No pair of clustering variables exceeds |association| = .85.")
}

# One-hot detection: a set of binary variables of which exactly one is 1 in
# almost every row is a single categorical question in disguise.
bin_vars <- vn[vt %in% c("symm_binary", "asymm_binary")]
if (length(bin_vars) >= 3) {
  rs <- rowSums(gower_frame[, bin_vars, drop = FALSE])
  if (mean(rs == 1) > 0.95)
    warn("The ", length(bin_vars), " binary variables sum to exactly 1 in ",
         round(100 * mean(rs == 1)), "% of cases. They are almost certainly ",
         "one-hot indicators of a single question; enter that question once, ",
         "as a nominal variable.")
}

#  5.3  Observed range of the interval variables -------------------------------
#  WHAT    Reports, for each interval variable, how much of its range is due
#          to a handful of extreme values, and optionally winsorises.
#  WHY     Gower divides each interval difference by the variable's OBSERVED
#          range. One outlier therefore inflates the denominator and shrinks
#          every other difference on that variable towards zero: the variable
#          stops contributing without ever being removed. The column below to
#          watch is range_ratio -- the trimmed range as a fraction of the full
#          one. A value near 1 means the range is well supported. A value of
#          0.5 means half the scale is being spent on the tails.
#  CHOICE  WINSORIZE_INTERVAL, WINSORIZE_PROBS (2.7)
#  STATUS  REQUIRED as a diagnostic; the winsorising itself is OPTIONAL

step("5.3  Range of the interval variables")

int_vars <- vn[vt == "interval"]
if (length(int_vars) > 0) {
  range_diag <- purrr::map_dfr(int_vars, function(v) {
    x  <- gower_frame[[v]]
    tr <- stats::quantile(x, WINSORIZE_PROBS, na.rm = TRUE)
    tibble::tibble(
      variable      = v,
      min = signif(min(x), 5), max = signif(max(x), 5),
      full_range    = signif(diff(range(x)), 5),
      trimmed_range = signif(diff(tr), 5),
      range_ratio   = round(diff(tr) / diff(range(x)), 3),
      skewness      = round(mean((x - mean(x))^3) / stats::sd(x)^3, 2))
  })
  print(as.data.frame(range_diag), row.names = FALSE)
  save_table(range_diag, "interval_range_diagnostics.csv",
             "Observed vs trimmed range of each interval variable")

  thin <- dplyr::filter(range_diag, range_ratio < 0.7)
  if (nrow(thin) > 0)
    note("Range driven by extreme values (ratio < .70): ",
         paste(thin$variable, collapse = ", "),
         ". Consider WINSORIZE_INTERVAL = TRUE and compare the solutions.")

  if (isTRUE(WINSORIZE_INTERVAL)) {
    for (v in int_vars) {
      x  <- gower_frame[[v]]
      lim <- stats::quantile(x, WINSORIZE_PROBS, na.rm = TRUE)
      n_capped <- sum(x < lim[1] | x > lim[2])
      gower_frame[[v]] <- pmin(pmax(x, lim[1]), lim[2])
      note("Winsorised ", v, ": ", n_capped, " value(s) capped.")
    }
    warn("Winsorising changes the distance matrix. Report it as an analytic ",
         "decision, and ideally report the unwinsorised solution alongside.")
  }
}

#  5.4  Variable weights -------------------------------------------------------
#  WHAT    Passes the declared weights to daisy().
#  WHY     Equal weighting is a choice, not a neutral default (5.2). If one
#          construct is represented by three variables and another by one, and
#          you consider them equally important, then the weights that express
#          that belief are 1/3 each and 1, not 1 each.
#  CHOICE  the 'weight' argument of cvar() (2.3)
#  STATUS  OPTIONAL - no effect when all weights are 1

use_weights <- any(abs(w - 1) > .Machine$double.eps^0.5)
if (use_weights) {
  if (!"weights" %in% names(formals(cluster::daisy)))
    stop("Variable weights require a newer version of the 'cluster' package.",
         "\n  Either update it, or set all weights to 1.", call. = FALSE)
  note("Using unequal variable weights; see gower_weight_share.csv.")
}

#  5.5  The Gower dissimilarity matrix -----------------------------------------
#  WHAT    Computes the matrix, then audits the treatment daisy() actually
#          applied to each variable.
#  WHY     The audit is the point. The absence of a warning is not evidence
#          that a variable was treated as you intended; attr(d, "Types") is.
#          A mismatch here means the clustering is being computed in a space
#          you did not specify.
#
#          Gower's coefficient (Gower, 1971) handles each variable in its own
#          terms -- interval differences scaled by range, simple matching for
#          nominal, the rank-based extension of Podani (1999) for ordinal,
#          presence-matching for asymmetric binary -- and averages them. This
#          is why it can take mixed data at all, and also why the averaging
#          discussed in 5.2 matters so much.
#  STATUS  REQUIRED

step("5.5  Gower distance")

gower_args <- list(x = gower_frame, metric = "gower")
if (length(daisy_type) > 0) gower_args$type <- daisy_type
if (use_weights)            gower_args$weights <- w

gower_dist <- do.call(cluster::daisy, gower_args)

observed_codes <- attr(gower_dist, "Types")
names(observed_codes) <- names(gower_frame)
expected_codes <- expected_daisy_codes(CLUSTER_VARS)

type_check <- tibble::tibble(
  variable       = names(observed_codes),
  declared_type  = spec_types(CLUSTER_VARS),
  expected_code  = unname(expected_codes),
  applied_code   = unname(observed_codes),
  interpretation = dplyr::recode(unname(observed_codes),
    I = "interval, scaled by observed range",
    N = "nominal, simple matching",
    O = "ordinal, rank-based (Podani)",
    S = "symmetric binary",
    A = "asymmetric binary",
    T = "ratio, logarithmic",
    .default = "unrecognised"),
  agrees = unname(observed_codes) == unname(expected_codes))

print(as.data.frame(type_check[, c(1, 3, 4, 6)]), row.names = FALSE)
save_table(type_check, "gower_variable_types.csv",
           "Measurement treatment daisy() applied to each variable")

if (!all(type_check$agrees)) {
  bad <- dplyr::filter(type_check, !agrees)
  stop("\n  daisy() did not apply the declared measurement type to: ",
       paste(bad$variable, collapse = ", "),
       "\n  Expected ", paste(bad$expected_code, collapse = "/"),
       ", applied ", paste(bad$applied_code, collapse = "/"),
       ".\n  The clustering would run in a space you did not specify.",
       call. = FALSE)
}
ok("Every variable was treated as declared.")

saveRDS(gower_dist, file.path(OUTPUT_DIR, "gower_dist.rds"))
OUTPUT_LOG$files <- c(OUTPUT_LOG$files,
  paste("gower_dist.rds",
        "The dissimilarity matrix, so it need not be recomputed", sep = "\t"))

#  5.6  Sanity checks on the matrix --------------------------------------------
#  WHAT    Distribution of the pairwise dissimilarities, and the number of
#          cases with identical profiles.
#  WHY     Two pathologies are invisible later. If the dissimilarities are all
#          bunched near one value, there is nothing for the algorithm to
#          separate, and any partition it returns is arbitrary. If many cases
#          have identical profiles -- common when all the variables are
#          categorical with few levels -- then distance is highly tied, PAM's
#          medoid choice becomes partly arbitrary, and the effective sample
#          size for clustering is the number of DISTINCT profiles rather than
#          the number of cases.
#  STATUS  REQUIRED

step("5.6  Sanity checks")

dvec <- as.vector(gower_dist)
n_distinct_profiles <- nrow(unique(gower_frame))
dist_summary <- tibble::tibble(
  n_cases            = nrow(gower_frame),
  n_distinct_profiles = n_distinct_profiles,
  n_pairs            = length(dvec),
  min = round(min(dvec), 4), q1 = round(stats::quantile(dvec, .25), 4),
  median = round(stats::median(dvec), 4),
  q3 = round(stats::quantile(dvec, .75), 4), max = round(max(dvec), 4),
  sd = round(stats::sd(dvec), 4),
  pct_zero_distance = round(100 * mean(dvec == 0), 3))
print(as.data.frame(dist_summary), row.names = FALSE)
save_table(dist_summary, "gower_distance_summary.csv",
           "Distribution of the pairwise Gower dissimilarities")

if (dist_summary$sd < 0.05)
  warn("The dissimilarities are almost constant (SD = ", dist_summary$sd,
       "). There may be no structure to partition; see 6.11.")
if (n_distinct_profiles < 0.5 * nrow(gower_frame))
  warn("Only ", n_distinct_profiles, " distinct profiles among ",
       nrow(gower_frame), " cases. Distances are heavily tied; treat the ",
       "medoids and the silhouette widths with caution.")


#  6.  HOW MANY CLUSTERS? (MODEL SELECTION) ====================================
#
#  No index chooses k. Milligan and Cooper (1985) compared thirty of them on
#  data with known structure and none was reliably best; nothing since has
#  changed that. What the indices can do is narrow the field and, more
#  usefully, disagree with each other in ways that tell you something.
#
#  This section assembles seven kinds of evidence and refuses to rank them for
#  you:
#
#    6.2   internal validity indices across the whole range
#    6.5   nesting -- does a larger solution SPLIT a cluster or RESHUFFLE?
#    6.6   cluster sizes -- is every cluster large enough to mean anything?
#    6.7   substantive profiles -- is each cluster interpretable in its own
#          right, or an arbitrary cut along a continuum?
#    6.8   stability -- does the solution reproduce itself under resampling?
#    6.9   replication -- does it reproduce in an independent half?
#    6.10  robustness -- does a different algorithm find the same thing?
#    6.11  is there any structure here at all?
#
#  The last is the one most often skipped and the one Rapkin and Luke (1993)
#  are most insistent about. PAM will return k tidy clusters from pure noise
#  and report a respectable silhouette width for them. Section 6.11 builds the
#  comparison that makes the observed value mean something.

banner("6.  MODEL SELECTION")

#  6.1  Fit PAM across the range -----------------------------------------------
#  WHAT    One PAM solution per value of k in K_RANGE.
#  WHY     Everything downstream reads from these fits, so they are computed
#          once and reused rather than refitted per diagnostic.
#
#          On the algorithm: PAM chooses k actual observations as medoids and
#          assigns every case to the nearest, minimising total dissimilarity
#          to the medoid (Kaufman & Rousseeuw, 1990). Two consequences worth
#          knowing. Because a medoid is a real case rather than an average, it
#          can be looked up, described and quoted -- which matters for
#          interpretation (7.2). And because it minimises dissimilarity rather
#          than squared Euclidean distance, it is markedly less sensitive to
#          extreme cases than k-means.
#  CHOICE  K_RANGE (2.8)
#  STATUS  REQUIRED

step("6.1  Fitting PAM across K_RANGE")

pam_fits <- lapply(K_RANGE, function(k)
  cluster::pam(gower_dist, diss = TRUE, k = k))
names(pam_fits) <- as.character(K_RANGE)
ok("Fitted ", length(pam_fits), " solutions, k = ", min(K_RANGE), " to ",
   max(K_RANGE), ".")

#  6.2  Internal validity indices ----------------------------------------------
#  WHAT    A table of indices for every k.
#  WHY     They summarise geometric separation in the space built in section 5.
#          Read them as descriptions of that space, not as verdicts about
#          reality.
#
#     avg_silhouette   mean over cases of (b - a)/max(a, b), where a is the
#                      mean dissimilarity to one's own cluster and b to the
#                      nearest other. Ranges -1 to 1; higher is better
#                      separated. The single most useful index here, because
#                      it is defined directly on a dissimilarity matrix and
#                      assumes no geometry.
#     n_negative_sil   cases whose silhouette is below zero: they resemble a
#                      neighbouring cluster more than their own. A
#                      classification-quality measure, and a more concrete one
#                      than the average.
#     dunn             smallest between-cluster separation over largest
#                      within-cluster diameter. Very sensitive to single
#                      cases; use as corroboration only.
#     calinski_harab   between/within variance ratio. Assumes a Euclidean
#                      geometry that Gower does not provide, so it is reported
#                      for continuity with the literature and should not
#                      outvote the silhouette here.
#     within_dissim    the PAM objective after the swap phase -- the quantity
#                      being minimised. Falls monotonically with k, so it is
#                      read as an elbow, never as a maximum.
#     n_pairwise       pairs of clusters to compare afterwards. Not a fit
#                      index: a reminder that every extra cluster costs
#                      multiplicity in whatever comes next.
#  STATUS  REQUIRED

step("6.2  Internal validity indices")

k_indices <- purrr::map_dfr(K_RANGE, function(k) {
  p  <- pam_fits[[as.character(k)]]
  cs <- fpc::cluster.stats(gower_dist, p$clustering)
  sw <- cluster::silhouette(p$clustering, gower_dist)
  tibble::tibble(
    k                  = k,
    avg_silhouette     = round(p$silinfo$avg.width, 4),
    n_negative_sil     = sum(sw[, "sil_width"] < 0),
    pct_negative_sil   = round(100 * mean(sw[, "sil_width"] < 0), 2),
    dunn               = round(cs$dunn, 4),
    calinski_harabasz  = round(cs$ch, 1),
    within_dissim      = round(p$objective[["swap"]], 4),
    min_cluster_size   = min(cs$cluster.size),
    max_cluster_size   = max(cs$cluster.size),
    n_below_min_size   = sum(cs$cluster.size < MIN_CLUSTER_SIZE),
    n_pairwise         = choose(k, 2))
})

print(as.data.frame(k_indices), row.names = FALSE)
save_table(k_indices, "model_selection_indices.csv",
           "Internal validity indices for every k in K_RANGE")

best_sil <- k_indices$k[which.max(k_indices$avg_silhouette)]
note("Highest average silhouette at k = ", best_sil,
     ". This is one piece of evidence, not the answer.")

#  6.3  Index figures ----------------------------------------------------------
#  WHAT    The indices of 6.2, plotted against k.
#  WHY     The SHAPE of these curves carries information the numbers do not.
#          A silhouette curve with a clear peak is one situation; a flat or
#          non-monotonic curve across several values of k is a different one,
#          and it means the index is not identifying an optimum at all -- in
#          which case the later diagnostics, not this figure, decide.
#  STATUS  REQUIRED

idx_long <- k_indices %>%
  dplyr::select(k,
                `Average silhouette`            = avg_silhouette,
                `Dunn index`                    = dunn,
                `Calinski-Harabasz`             = calinski_harabasz,
                `Within-cluster dissimilarity`  = within_dissim,
                `Cases better matched elsewhere (%)` = pct_negative_sil,
                `Smallest cluster (n)`          = min_cluster_size) %>%
  tidyr::pivot_longer(-k, names_to = "index", values_to = "value")

p_idx <- ggplot(idx_long, aes(x = k, y = value)) +
  geom_line(colour = "grey40") +
  geom_point(size = 1.8) +
  facet_wrap(~ index, scales = "free_y") +
  scale_x_continuous(breaks = K_RANGE) +
  labs(title = "Model selection across the candidate range",
       subtitle = "PAM on Gower dissimilarity",
       x = "Number of clusters (k)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

save_figure(p_idx, "model_selection_indices.png",
            "Validity indices plotted against k", width = 9, height = 6)

#  6.4  Candidate solutions ----------------------------------------------------
#  WHAT    Selects the small set of solutions examined in depth below.
#  WHY     The diagnostics from 6.5 onward are expensive; running them on nine
#          solutions wastes time and buries the comparison that matters. The
#          default -- the three highest average silhouette widths -- is a
#          screening rule, not a decision rule.
#  CHOICE  K_CANDIDATES (2.8). Override it when theory points elsewhere: a
#          solution predicted by your framework belongs among the candidates
#          whether or not its silhouette is in the top three.
#  STATUS  REQUIRED

step("6.4  Candidate solutions")

k_candidates <- if (is.null(K_CANDIDATES)) {
  sort(k_indices$k[order(-k_indices$avg_silhouette)][seq_len(min(3, nrow(k_indices)))])
} else sort(K_CANDIDATES)

ok("Candidates: k = ", paste(k_candidates, collapse = ", "),
   if (is.null(K_CANDIDATES)) "  (top three by silhouette)" else
     "  (set in CONFIG)")

pam_candidates <- pam_fits[as.character(k_candidates)]

#  6.5  Nesting: does a larger solution split, or reshuffle? -------------------
#  WHAT    Cross-tabulates each pair of candidate solutions, and quantifies
#          the relation two ways.
#  WHY     This is the single most informative comparison between adjacent
#          solutions, and the least often reported.
#
#            A larger solution that SPLITS one cluster while leaving the
#            others intact is evidence of genuine internal heterogeneity in
#            that cluster: the extra cluster is telling you something.
#
#            A larger solution that RESHUFFLES members across several clusters
#            is evidence that the extra cluster is an artefact of forcing the
#            partition into more pieces. Prefer the smaller solution.
#
#          nesting_purity is the mean, over the clusters of the LARGER
#          solution, of the proportion of members drawn from a single cluster
#          of the smaller one. 1.00 is a perfect split; values near 1/k_small
#          indicate wholesale reorganisation. The adjusted Rand index measures
#          overall agreement, corrected for chance.
#  STATUS  REQUIRED

step("6.5  Nesting between candidate solutions")

nesting_purity <- function(small, large) {
  tab <- table(large, small)
  mean(apply(tab, 1, max) / rowSums(tab))
}

nesting_rows <- list()
for (i in seq_along(k_candidates)) for (j in seq_along(k_candidates)) {
  if (j <= i) next
  ks <- k_candidates[i]; kl <- k_candidates[j]
  cs <- pam_fits[[as.character(ks)]]$clustering
  cl <- pam_fits[[as.character(kl)]]$clustering
  tab <- table(smaller = cs, larger = cl)
  save_table(as.data.frame.matrix(tab) %>%
               tibble::rownames_to_column(paste0("k", ks, "_cluster")),
             paste0("nesting_k", ks, "_vs_k", kl, ".csv"),
             paste0("Cross-tabulation of the k = ", ks, " and k = ", kl,
                    " solutions"))
  nesting_rows[[length(nesting_rows) + 1]] <- tibble::tibble(
    comparison          = paste0("k", ks, " vs k", kl),
    adjusted_rand_index = round(adjusted_rand(cs, cl), 3),
    nesting_purity      = round(nesting_purity(cs, cl), 3),
    reading = dplyr::case_when(
      nesting_purity(cs, cl) >= 0.95 ~ "clean split of existing cluster(s)",
      nesting_purity(cs, cl) >= 0.80 ~ "mostly a split, with some reshuffling",
      TRUE                           ~ "substantial reshuffling"))
}
partition_agreement <- dplyr::bind_rows(nesting_rows)
print(as.data.frame(partition_agreement), row.names = FALSE)
save_table(partition_agreement, "nesting_summary.csv",
           "Agreement and nesting between the candidate solutions")

#  6.6  Cluster sizes ----------------------------------------------------------
#  WHAT    The size of every cluster in every candidate solution.
#  WHY     A cluster of eight cases is not a type; it is eight cases. It
#          cannot be described with any precision, cannot support a
#          comparison, and is usually the first thing to dissolve under
#          resampling in 6.8. Size is not a fit criterion, but a solution that
#          buys its silhouette with a micro-cluster has bought it on credit.
#  CHOICE  MIN_CLUSTER_SIZE (2.8) -- a flag, not a rule
#  STATUS  REQUIRED

cluster_sizes <- purrr::imap_dfr(pam_candidates, function(p, lbl)
  tibble::tibble(k = as.integer(lbl),
                 cluster = seq_along(p$clusinfo[, "size"]),
                 size = as.integer(p$clusinfo[, "size"]),
                 pct = round(100 * p$clusinfo[, "size"] / sum(p$clusinfo[, "size"]), 1)))
print(as.data.frame(tidyr::pivot_wider(cluster_sizes[, 1:3],
        names_from = k, values_from = size, names_prefix = "k=")),
      row.names = FALSE)
save_table(cluster_sizes, "cluster_sizes_candidates.csv",
           "Cluster sizes in each candidate solution")

small_clusters <- dplyr::filter(cluster_sizes, size < MIN_CLUSTER_SIZE)
if (nrow(small_clusters) > 0) {
  warn("Cluster(s) below MIN_CLUSTER_SIZE = ", MIN_CLUSTER_SIZE, ":")
  print(as.data.frame(small_clusters), row.names = FALSE)
}

#  6.7  Substantive profiles of the candidate solutions ------------------------
#  WHAT    Describes every cluster of every candidate on every clustering
#          variable, in a form built from the declared types: median and
#          interquartile range for interval variables, percentage per category
#          for the rest.
#  WHY     This is the decisive criterion and the only one that is not a
#          number. An extra cluster earns its place when it is interpretable
#          IN ITS OWN RIGHT. A cluster that differs from its neighbour only by
#          sitting slightly lower on one continuous variable, while sharing
#          every categorical characteristic, is a cut along a continuum, not a
#          type - however good the indices look.
#
#          Rapkin and Luke (1993) put this the other way round, and more
#          sharply: the researcher's obligation is to say what each cluster
#          MEANS, and a solution whose clusters cannot be described without
#          referring to their numbers has not been interpreted at all.
#  STATUS  REQUIRED

step("6.7  Cluster profiles")

profile_by_cluster <- function(data, spec, clustering) {
  frame <- prepare_gower_frame(data, spec)
  cl    <- factor(clustering)
  nmv   <- spec_names(spec); tyv <- spec_types(spec)
  ivars <- nmv[tyv == "interval"]
  cvars <- setdiff(nmv, ivars)

  interval_tbl <- if (length(ivars)) purrr::map_dfr(ivars, function(v) {
    x <- frame[[v]]
    tibble::tibble(cluster = cl, value = x) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(variable = v, n = dplyr::n(),
                       median = round(stats::median(value), 2),
                       q1 = round(stats::quantile(value, .25), 2),
                       q3 = round(stats::quantile(value, .75), 2),
                       mean = round(mean(value), 2),
                       sd = round(stats::sd(value), 2), .groups = "drop") %>%
      dplyr::relocate(variable)
  }) else tibble::tibble()

  categorical_tbl <- if (length(cvars)) purrr::map_dfr(cvars, function(v) {
    x <- frame[[v]]
    if (!is.factor(x)) {
      vv <- spec_get(spec, v)
      lb <- if (is.null(vv$labels)) c("0", "1") else as.character(vv$labels)
      if (length(lb) != 2) lb <- c("0", "1")
      x <- factor(x, levels = c(0, 1), labels = lb)
    }
    # 'level' is carried as character: different variables have different
    # level sets, and binding ordered with unordered factors would fail. The
    # original ordering is preserved in level_order so that figures and
    # tables can restore it.
    tibble::tibble(cluster = cl, level = x) %>%
      # .drop = FALSE keeps categories that no member of a cluster chose, so
      # they appear as 0% rather than as a hole in the table.
      dplyr::count(cluster, level, name = "n", .drop = FALSE) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(variable = v, pct = round(100 * n / sum(n), 1)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(level_order = as.integer(level),
                    level = as.character(level)) %>%
      dplyr::relocate(variable, cluster, level, level_order)
  }) else tibble::tibble()

  list(interval = interval_tbl, categorical = categorical_tbl)
}

# A compact, readable table: one row per cluster, one column per variable.
format_profile_wide <- function(prof) {
  parts <- list()
  if (nrow(prof$interval) > 0) {
    parts$int <- prof$interval %>%
      dplyr::mutate(value = sprintf("%.1f [%.1f-%.1f]", median, q1, q3)) %>%
      dplyr::select(cluster, variable, value) %>%
      tidyr::pivot_wider(names_from = variable, values_from = value)
  }
  if (nrow(prof$categorical) > 0) {
    parts$cat <- prof$categorical %>%
      dplyr::mutate(col = paste0(variable, ": ", level),
                    value = sprintf("%.1f%%", pct)) %>%
      dplyr::select(cluster, col, value) %>%
      tidyr::pivot_wider(names_from = col, values_from = value)
  }
  n_tbl <- if (nrow(prof$interval) > 0)
    dplyr::distinct(dplyr::select(prof$interval, cluster, n)) else
    prof$categorical %>% dplyr::group_by(cluster) %>%
      dplyr::summarise(n = sum(n) / dplyr::n_distinct(.data$variable),
                       .groups = "drop")
  Reduce(function(a, b) dplyr::left_join(a, b, by = "cluster"),
         c(list(n_tbl), unname(parts)))
}

for (kk in k_candidates) {
  prof <- profile_by_cluster(data_analytic, CLUSTER_VARS,
                             pam_fits[[as.character(kk)]]$clustering)
  wide <- format_profile_wide(prof)
  cat("\n  Cluster profiles, k = ", kk,
      "   (interval: median [Q1-Q3];  categorical: % of cluster)\n", sep = "")
  print(as.data.frame(wide), row.names = FALSE)
  save_table(wide, paste0("profiles_k", kk, ".csv"),
             paste0("Cluster profiles for the k = ", kk, " solution"))
}

#  6.8  Stability under resampling ---------------------------------------------
#  WHAT    Draws bootstrap samples, re-clusters each, and reports for every
#          original cluster the mean Jaccard similarity with its closest
#          counterpart in the bootstrap solution (Hennig, 2007).
#  WHY     A cluster that does not reappear when the sample is perturbed will
#          not reappear in a new sample either, whatever its silhouette. This
#          is usually the diagnostic that separates two solutions whose
#          indices are close: the extra cluster obtained by cutting a
#          continuum is markedly less stable than the clusters around it.
#
#          Conventional reading (Hennig, 2007):
#            below 0.60   unstable; the cluster should not be interpreted
#            0.60 - 0.75  a pattern is present, its membership uncertain
#            0.75 - 0.85  stable
#            above 0.85   highly stable
#
#          HOW MANY REPLICATIONS. B = 100 is convention rather than
#          derivation. With STABILITY_ADAPTIVE, blocks are added until the
#          mean Jaccard of EVERY cluster is estimated to within
#          STABILITY_TARGET_MCSE, so that a comparison between two solutions
#          rests on a difference larger than the noise in its own estimate.
#          The rule depends only on the spread of the bootstrap distribution,
#          never on the cluster means, so it cannot become a
#          stop-when-favourable rule. A solution whose weakest cluster has a
#          bimodal bootstrap distribution may never reach the target; it stops
#          at STABILITY_B_MAX with a warning, and that failure is itself
#          evidence against the solution.
#  CHOICE  RUN_STABILITY, STABILITY_* (2.11)
#  STATUS  OPTIONAL, strongly recommended

cluster_stability <- NULL
stability_summary <- NULL

if (isTRUE(RUN_STABILITY)) {
  step("6.8  Bootstrap stability")
  note("This is the slowest part of the run. Progress is printed per block.")

  run_block <- function(k, seed) {
    set.seed(seed)
    fpc::clusterboot(gower_dist,
                     B             = STABILITY_BLOCK,
                     distances     = TRUE,
                     bootmethod    = "boot",
                     clustermethod = fpc::claraCBI,
                     k             = k,
                     usepam        = TRUE,
                     diss          = TRUE,
                     count         = FALSE)
  }

  stability_for_k <- function(k) {
    pooled <- NULL
    block  <- 0L
    repeat {
      block  <- block + 1L
      cb     <- run_block(k, seed = SEED + 1000L * k + block)
      # Blocks pool column-wise: replications are independent and PAM on a
      # fixed dissimilarity matrix returns the same original partition every
      # time, so row i is the same cluster in every block.
      pooled <- cbind(pooled, cb$bootresult)
      B_done <- ncol(pooled)
      sd_k   <- apply(pooled, 1, stats::sd)
      mcse   <- sd_k / sqrt(B_done)
      cat(sprintf("     k = %d | B = %4d | max SD = %.3f | max MCSE = %.4f\n",
                  k, B_done, max(sd_k), max(mcse)))
      if (!isTRUE(STABILITY_ADAPTIVE) && B_done >= STABILITY_B_MIN) break
      if (B_done >= STABILITY_B_MIN && max(mcse) <= STABILITY_TARGET_MCSE) break
      if (B_done >= STABILITY_B_MAX) {
        warn("k = ", k, ": reached STABILITY_B_MAX with max MCSE = ",
             round(max(mcse), 4), " (target ", STABILITY_TARGET_MCSE,
             "). The estimate is less precise than requested, which is ",
             "itself informative about this solution.")
        break
      }
    }
    B_done <- ncol(pooled)
    tibble::tibble(
      k              = k,
      cluster        = seq_len(nrow(pooled)),
      B_used         = B_done,
      jaccard_mean   = round(rowMeans(pooled), 3),
      jaccard_sd     = round(apply(pooled, 1, stats::sd), 3),
      jaccard_mcse   = round(apply(pooled, 1, stats::sd) / sqrt(B_done), 4),
      prop_dissolved = round(rowMeans(pooled < 0.50), 3),
      prop_recovered = round(rowMeans(pooled > 0.75), 3))
  }

  cluster_stability <- purrr::map_dfr(k_candidates, stability_for_k)
  print(as.data.frame(cluster_stability), row.names = FALSE)

  # The weakest cluster decides whether a solution can be interpreted, so the
  # summary is built around the minimum rather than the average.
  stability_summary <- cluster_stability %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(
      B_used         = dplyr::first(B_used),
      jaccard_min    = min(jaccard_mean),
      jaccard_mean   = round(mean(jaccard_mean), 3),
      n_below_0.60   = sum(jaccard_mean < 0.60),
      n_below_0.75   = sum(jaccard_mean < 0.75),
      n_dissolving   = sum(prop_dissolved > 0),
      max_mcse       = max(jaccard_mcse), .groups = "drop")
  cat("\n")
  print(as.data.frame(stability_summary), row.names = FALSE)

  save_table(cluster_stability, "stability_bootstrap.csv",
             "Bootstrap Jaccard similarity, per cluster and candidate")
  save_table(stability_summary, "stability_summary.csv",
             "Bootstrap stability summarised per candidate solution")

  p_stab <- ggplot(cluster_stability,
                   aes(x = factor(cluster), y = jaccard_mean)) +
    geom_hline(yintercept = 0.75, linetype = 2, colour = "grey35") +
    geom_hline(yintercept = c(0.60, 0.85), linetype = 3, colour = "grey65") +
    geom_pointrange(aes(ymin = jaccard_mean - 1.96 * jaccard_mcse,
                        ymax = jaccard_mean + 1.96 * jaccard_mcse)) +
    facet_wrap(~ paste0("k = ", k), scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Bootstrap stability of each cluster",
         subtitle = paste("Dashed line: 0.75. Dotted: 0.60 and 0.85",
                          "(Hennig, 2007)"),
         x = "Cluster", y = "Mean Jaccard similarity") +
    theme_minimal(base_size = 11)
  save_figure(p_stab, "stability_bootstrap.png",
              "Bootstrap Jaccard similarity per cluster")
}

#  6.9  Replication across random halves ---------------------------------------
#  WHAT    Splits the sample at random, clusters one half, classifies the
#          other half to its nearest medoid, and compares that classification
#          with the clustering the second half produces on its own. Repeated
#          SPLIT_N times.
#  WHY     Bootstrap stability (6.8) perturbs the sample; this replicates the
#          whole procedure on data the solution has not seen. It answers a
#          different and slightly harder question -- would an independent
#          sample of this population yield the same typology? -- and it is the
#          form of validation Rapkin and Luke (1993) treat as the minimum
#          standard for a typology offered as a finding rather than as a
#          description.
#
#          Read the mean adjusted Rand index as agreement corrected for
#          chance. There is no threshold with a claim to authority; as a
#          working guide, above .70 is strong replication, .40 to .70 partial,
#          below .40 means the typology is sample-specific.
#
#          Note that this halves the sample, and small clusters suffer most.
#          A low value for a solution with a small cluster may reflect the
#          split rather than the solution.
#  CHOICE  RUN_SPLIT_HALF, SPLIT_N (2.11)
#  STATUS  OPTIONAL, recommended

split_half <- NULL

if (isTRUE(RUN_SPLIT_HALF)) {
  step("6.9  Split-half replication")
  dm <- as.matrix(gower_dist)
  n  <- nrow(dm)

  split_half_for_k <- function(k) {
    ari <- numeric(SPLIT_N)
    for (s in seq_len(SPLIT_N)) {
      set.seed(SEED + 5000L + 100L * k + s)
      idx <- sample.int(n)
      A <- idx[seq_len(floor(n / 2))]
      B <- idx[(floor(n / 2) + 1L):n]
      pa <- cluster::pam(stats::as.dist(dm[A, A]), diss = TRUE, k = k)
      pb <- cluster::pam(stats::as.dist(dm[B, B]), diss = TRUE, k = k)
      med_A <- A[pa$id.med]
      # Each case of half B is assigned to the nearest medoid of half A --
      # the rule PAM itself uses, applied out of sample.
      pred_B <- max.col(-dm[B, med_A, drop = FALSE], ties.method = "first")
      ari[s] <- adjusted_rand(pred_B, pb$clustering)
    }
    tibble::tibble(k = k, n_splits = SPLIT_N,
                   ari_mean   = round(mean(ari), 3),
                   ari_sd     = round(stats::sd(ari), 3),
                   ari_min    = round(min(ari), 3),
                   ari_max    = round(max(ari), 3))
  }

  split_half <- purrr::map_dfr(k_candidates, split_half_for_k)
  print(as.data.frame(split_half), row.names = FALSE)
  save_table(split_half, "replication_split_half.csv",
             "Adjusted Rand index between independent halves, per candidate")
}

#  6.10  Robustness to the choice of algorithm ---------------------------------
#  WHAT    Partitions the SAME Gower matrix with three other algorithms and
#          measures how far each departs from the PAM solution.
#  WHY     A typology that appears only under one algorithm is a property of
#          that algorithm. One that survives average-linkage, Ward and a fuzzy
#          method is a property of the data. This is cheap to run and almost
#          never reported.
#
#          What the comparison is NOT: a tournament. The other algorithms are
#          not competing to be chosen, and a low agreement does not mean PAM
#          is wrong. It means the partition is method-dependent, which is a
#          finding you are obliged to report rather than a problem to fix.
#
#          A caveat on Ward: it minimises an increase in squared Euclidean
#          error, a criterion that does not strictly apply to a Gower matrix.
#          It is included because it is ubiquitous in this literature, and
#          average linkage is reported alongside it precisely because it makes
#          no such assumption.
#  CHOICE  RUN_ALGO_ROBUSTNESS (2.11)
#  STATUS  OPTIONAL, recommended

algo_robustness <- NULL

if (isTRUE(RUN_ALGO_ROBUSTNESS)) {
  step("6.10  Robustness to the algorithm")

  algo_for_k <- function(k) {
    ref <- pam_fits[[as.character(k)]]$clustering
    out <- list()
    hc_a <- stats::hclust(gower_dist, method = "average")
    out$`Hierarchical (average linkage)` <- stats::cutree(hc_a, k)
    hc_w <- stats::hclust(gower_dist, method = "ward.D2")
    out$`Hierarchical (Ward D2)` <- stats::cutree(hc_w, k)
    fz <- tryCatch(cluster::fanny(gower_dist, k = k, diss = TRUE,
                                  memb.exp = 1.2, maxit = 1000),
                   error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fz)) out$`Fuzzy (fanny)` <- fz$clustering

    purrr::imap_dfr(out, function(cl, nm) tibble::tibble(
      k = k, method = nm,
      adjusted_rand_vs_pam = round(adjusted_rand(ref, cl), 3),
      avg_silhouette = round(mean(cluster::silhouette(cl, gower_dist)[, 3]), 3),
      min_cluster_size = min(table(cl))))
  }

  algo_robustness <- purrr::map_dfr(k_candidates, algo_for_k)
  print(as.data.frame(algo_robustness), row.names = FALSE)
  save_table(algo_robustness, "robustness_algorithms.csv",
             "Agreement between PAM and alternative algorithms on the same matrix")

  weak <- dplyr::filter(algo_robustness, adjusted_rand_vs_pam < 0.50)
  if (nrow(weak) > 0)
    note("Some alternatives agree only weakly with PAM (ARI < .50). The ",
         "partition is method-dependent; say so in the manuscript.")
}

#  6.11  Is there any structure here at all? -----------------------------------
#  WHAT    Permutes each clustering variable independently, which destroys the
#          associations between variables while leaving every marginal
#          distribution exactly as it was, then runs the whole clustering on
#          the permuted data. Repeated PERM_N times, for every k.
#  WHY     Because a clustering algorithm always returns clusters. PAM applied
#          to variables that carry no joint structure will still split the
#          cases into k compact groups and report a perfectly respectable
#          silhouette width for them -- and there is no way to tell from that
#          number alone. The permuted data provide the missing comparison:
#          they have your sample size, your variable types, your marginal
#          distributions and your missingness-free structure, and nothing
#          else. Whatever silhouette they produce is what your design yields
#          from noise.
#
#          Read it as: observed silhouette versus the permutation band. Well
#          above it, the structure is not an artefact of the design. Inside
#          it, the partition is indistinguishable from what noise produces,
#          and no amount of stability or profiling can rescue it -- stable
#          nonsense is still nonsense.
#
#          This is the operational answer to the objection Rapkin and Luke
#          (1993) raise against cluster analysis in the social sciences, and
#          the reason it is on by default.
#  CHOICE  RUN_PERMUTATION, PERM_N (2.11)
#  STATUS  OPTIONAL, strongly recommended

permutation_null <- NULL

if (isTRUE(RUN_PERMUTATION)) {
  step("6.11  Permutation reference distribution")
  note("Running ", PERM_N, " permutations over ", length(K_RANGE),
       " values of k.")

  perm_mat <- matrix(NA_real_, nrow = PERM_N, ncol = length(K_RANGE),
                     dimnames = list(NULL, as.character(K_RANGE)))
  for (r in seq_len(PERM_N)) {
    set.seed(SEED + 20000L + r)
    pf <- gower_frame
    for (j in seq_along(pf)) pf[[j]] <- sample(pf[[j]])
    pargs <- list(x = pf, metric = "gower")
    if (length(daisy_type) > 0) pargs$type <- daisy_type
    if (use_weights)            pargs$weights <- w
    dperm <- do.call(cluster::daisy, pargs)
    for (i in seq_along(K_RANGE))
      perm_mat[r, i] <- cluster::pam(dperm, diss = TRUE,
                                     k = K_RANGE[i])$silinfo$avg.width
    if (r %% max(1, floor(PERM_N / 5)) == 0)
      cat("     permutation ", r, " of ", PERM_N, "\n", sep = "")
  }

  obs_sil <- k_indices$avg_silhouette
  # Transposed twice so that the comparison recycles the observed value along
  # k rather than down the replications.
  exceed  <- colSums(t(t(perm_mat) >= obs_sil))

  permutation_null <- tibble::tibble(
    k                = K_RANGE,
    observed         = obs_sil,
    null_mean        = round(colMeans(perm_mat), 4),
    null_p95         = round(apply(perm_mat, 2, stats::quantile, .95), 4),
    null_max         = round(apply(perm_mat, 2, max), 4),
    exceeds_null_p95 = obs_sil > apply(perm_mat, 2, stats::quantile, .95),
    # Proportion of permutations reaching or exceeding the observed value:
    # a permutation p value for "no multivariate structure".
    p_permutation    = round((exceed + 1) / (PERM_N + 1), 4))
  print(as.data.frame(permutation_null), row.names = FALSE)
  save_table(permutation_null, "permutation_reference.csv",
             "Observed silhouette against a no-structure reference distribution")

  perm_long <- tibble::tibble(
    k = rep(K_RANGE, each = PERM_N), value = as.vector(perm_mat))
  p_perm <- ggplot(perm_long, aes(x = factor(k), y = value)) +
    geom_boxplot(outlier.size = 0.6, fill = "grey92", colour = "grey45") +
    geom_point(data = tibble::tibble(k = factor(K_RANGE),
                                     value = k_indices$avg_silhouette),
               colour = "firebrick", size = 2.4) +
    geom_line(data = tibble::tibble(k = factor(K_RANGE),
                                    value = k_indices$avg_silhouette),
              aes(group = 1), colour = "firebrick") +
    labs(title = "Observed silhouette against a no-structure reference",
         subtitle = paste0("Boxes: ", PERM_N, " permutations preserving each ",
                           "variable's marginal distribution.  Red: observed."),
         x = "Number of clusters (k)", y = "Average silhouette width") +
    theme_minimal(base_size = 11)
  save_figure(p_perm, "permutation_reference.png",
              "Observed silhouette vs the permutation reference distribution")

  if (!any(permutation_null$exceeds_null_p95))
    warn("No value of k produces a silhouette above the 95th percentile of ",
         "the no-structure reference. The data may not contain cluster ",
         "structure that this variable set can detect. Read 6.11 again ",
         "before continuing.")
}

#  6.12  Decision matrix -------------------------------------------------------
#  WHAT    Assembles every criterion into one table, one row per candidate,
#          and stops the run if K_FINAL has not been set.
#  WHY     The criteria pull in different directions more often than not, and
#          the table is where you see that. It deliberately contains no
#          overall score and no recommendation: a weighted composite would
#          hide precisely the disagreements you need to reason about, and
#          would put the decision back in the software's hands.
#
#          When you write your rationale, say which criteria you weighed most
#          and why. "The indices favoured seven but the seventh cluster was
#          unstable and interpretable only as a cut along one continuum" is a
#          defensible sentence. "k was chosen by silhouette" is not, when the
#          silhouette curve was flat.
#  CHOICE  K_FINAL, K_FINAL_RATIONALE (2.9)
#  STATUS  REQUIRED

step("6.12  Decision matrix")

decision <- k_indices %>% dplyr::filter(k %in% k_candidates) %>%
  dplyr::select(k, avg_silhouette, pct_negative_sil, dunn,
                calinski_harabasz, min_cluster_size, n_below_min_size,
                n_pairwise)
if (!is.null(stability_summary))
  decision <- dplyr::left_join(decision,
    dplyr::select(stability_summary, k, jaccard_min, n_below_0.75), by = "k")
if (!is.null(split_half))
  decision <- dplyr::left_join(decision,
    dplyr::select(split_half, k, replication_ari = ari_mean), by = "k")
if (!is.null(permutation_null))
  decision <- dplyr::left_join(decision,
    dplyr::select(permutation_null, k, exceeds_null_p95, p_permutation),
    by = "k")
if (!is.null(algo_robustness))
  decision <- dplyr::left_join(decision,
    algo_robustness %>% dplyr::group_by(k) %>%
      dplyr::summarise(algo_agreement_min = min(adjusted_rand_vs_pam),
                       .groups = "drop"), by = "k")
if (nrow(partition_agreement) > 0)
  decision <- dplyr::left_join(decision,
    dplyr::select(cluster_sizes %>% dplyr::group_by(k) %>%
                    dplyr::summarise(smallest_pct = min(pct), .groups = "drop"),
                  k, smallest_pct), by = "k")

cat("\n")
print(as.data.frame(decision), row.names = FALSE)
save_table(decision, "decision_matrix.csv",
           "All model-selection criteria, one row per candidate solution")

if (is.null(K_FINAL)) {
  banner("PASS 1 COMPLETE")
  cat("
  The model-selection diagnostics are written to ", dir_tables, "/ and
  ", dir_figures, "/. Read decision_matrix.csv first, then the profile
  tables for each candidate.

  To continue:
    1. decide how many clusters to retain;
    2. set K_FINAL in section 2.9;
    3. write your reasons in K_FINAL_RATIONALE -- they are reproduced in
       the reporting checklist and are part of the deposit;
    4. run the script again.

  The stop below is intentional. It is not an error in the pipeline.
", sep = "")
  stop("Pass 1 complete: set K_FINAL in section 2.9, then run again.",
       call. = FALSE)
}

ok("K_FINAL = ", K_FINAL, ". Continuing to the retained solution.")


#  7.  THE RETAINED SOLUTION ===================================================

banner(paste0("7.  RETAINED SOLUTION (k = ", K_FINAL, ")"))

#  7.1  Final partition --------------------------------------------------------
#  WHAT    Takes the PAM solution at K_FINAL. Under multiple imputation with
#          MI_MODE = "consensus", derives the partition from the co-membership
#          of the m imputed solutions instead.
#  CHOICE  K_FINAL (2.9), MI_MODE (2.6)
#  STATUS  REQUIRED

pam_final        <- pam_fits[[as.character(K_FINAL)]]
final_clustering <- pam_final$clustering
mi_partitions    <- NULL

if (MISSING_METHOD == "imputation") {
  step("7.1  Clustering each imputed dataset")
  mi_partitions <- lapply(imputed_frames, function(df) {
    fr <- prepare_gower_frame(df, CLUSTER_VARS)
    ar <- list(x = fr, metric = "gower")
    if (length(daisy_type) > 0) ar$type <- daisy_type
    if (use_weights)            ar$weights <- w
    cluster::pam(do.call(cluster::daisy, ar), diss = TRUE,
                 k = K_FINAL)$clustering
  })

  if (MI_MODE == "consensus") {
    co <- Reduce(`+`, lapply(mi_partitions,
                             function(p) outer(p, p, "==") * 1)) /
          length(mi_partitions)
    consensus_dist   <- stats::as.dist(1 - co)
    pam_consensus    <- cluster::pam(consensus_dist, diss = TRUE, k = K_FINAL)
    final_clustering <- pam_consensus$clustering
    note("Partition derived from the co-membership of ",
         length(mi_partitions), " imputed solutions.")
    warn("In consensus mode the partition is not a PAM solution on the Gower ",
         "matrix: it has no medoids in the original space, and section 7.2 ",
         "reports the medoids of the reference imputation instead.")
  }
}

cluster_n <- table(final_clustering)
cat("  Cluster sizes:\n")
print(data.frame(cluster = names(cluster_n),
                 n = as.integer(cluster_n),
                 pct = round(100 * as.integer(cluster_n) / sum(cluster_n), 1)),
      row.names = FALSE)

#  7.2  Medoids ----------------------------------------------------------------
#  WHAT    Prints the actual cases chosen as cluster centres.
#  WHY     This is what k-medoids gives you that k-means does not. A centroid
#          is an average that may correspond to no one: a mean neighbourhood
#          type of 2.4 describes nobody. A medoid is a real case, with a real
#          identifier, that you can look up, describe, and -- in a mixed-
#          methods design -- go back to. When you write up the typology, the
#          medoid is the most honest thing to call a representative case.
#  STATUS  REQUIRED

step("7.2  Medoids")
medoid_rows <- data_analytic[pam_final$id.med,
                             c(ID_VAR, spec_names(CLUSTER_VARS)),
                             drop = FALSE]
medoid_rows <- cbind(cluster = seq_len(K_FINAL), medoid_rows)
print(as.data.frame(medoid_rows), row.names = FALSE)
save_table(medoid_rows, "medoids.csv",
           "The representative case chosen for each cluster")

#  7.3  Profiles of the retained solution --------------------------------------
#  WHAT    Describes each cluster on each clustering variable, and identifies
#          what makes each cluster distinctive.
#  WHY     Description, not inference. And here is the rule that gets broken
#          most often in published cluster analyses:
#
#          DO NOT TEST WHETHER THE CLUSTERS DIFFER ON THE VARIABLES THAT BUILT
#          THEM. The algorithm maximised exactly that separation. The test
#          will be significant, the effect size will be large, and both
#          numbers are artefacts of the procedure. An F test here is not weak
#          evidence; it is no evidence.
#
#          The pipeline therefore reports these differences descriptively and
#          runs no test on them. Tests belong in 7.7, on variables that did
#          not help build the clusters.
#  STATUS  REQUIRED

step("7.3  Cluster profiles")

profile_final <- profile_by_cluster(data_analytic, CLUSTER_VARS,
                                    final_clustering)
profile_wide  <- format_profile_wide(profile_final)
print(as.data.frame(profile_wide), row.names = FALSE)

save_table(profile_wide, "final_profiles_wide.csv",
           "Cluster profiles, formatted for reading")
if (nrow(profile_final$interval) > 0)
  save_table(profile_final$interval, "final_profiles_interval.csv",
             "Interval variables by cluster: median, IQR, mean, SD")
if (nrow(profile_final$categorical) > 0)
  save_table(profile_final$categorical, "final_profiles_categorical.csv",
             "Categorical variables by cluster: counts and percentages")

# What makes each cluster distinctive: the largest departures from the sample
# as a whole. For interval variables the departure is expressed in robust
# standard units (median absolute deviation), for categorical ones in
# percentage points. This is the raw material for naming the clusters.
distinctive <- local({
  frame <- prepare_gower_frame(data_analytic, CLUSTER_VARS)
  rows  <- list()
  for (v in spec_names(CLUSTER_VARS)) {
    x  <- frame[[v]]
    ty <- spec_get(CLUSTER_VARS, v)$type
    if (ty == "interval") {
      s <- stats::mad(x); if (s == 0) s <- stats::sd(x); if (s == 0) next
      for (cl in sort(unique(final_clustering))) {
        d <- (stats::median(x[final_clustering == cl]) - stats::median(x)) / s
        rows[[length(rows) + 1]] <- tibble::tibble(
          cluster = cl, variable = v, feature = "median",
          departure = round(d, 2), unit = "robust SD")
      }
    } else {
      xf <- if (is.factor(x)) x else factor(x)
      for (cl in sort(unique(final_clustering))) {
        for (lv in levels(xf)) {
          d <- 100 * (mean(xf[final_clustering == cl] == lv) - mean(xf == lv))
          rows[[length(rows) + 1]] <- tibble::tibble(
            cluster = cl, variable = v, feature = as.character(lv),
            departure = round(d, 1), unit = "percentage points")
        }
      }
    }
  }
  dplyr::bind_rows(rows) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(abs(departure)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>% dplyr::ungroup()
})

cat("\n  Most distinctive features of each cluster ",
    "(departure from the whole sample):\n", sep = "")
print(as.data.frame(distinctive), row.names = FALSE)
save_table(distinctive, "final_distinctive_features.csv",
           "The four largest departures from the sample, per cluster")

#  7.4  Silhouette analysis of the retained solution ---------------------------
#  WHAT    Per-case and per-cluster silhouette widths, and the cases that are
#          better matched to a neighbouring cluster than to their own.
#  WHY     The average silhouette hides where the classification is weak. A
#          solution with one poorly separated cluster and five clean ones is a
#          different object from one that is uniformly mediocre, and the
#          write-up should say which it is. The negative-silhouette cases are
#          exported because they are the ones a reviewer will ask about.
#  STATUS  REQUIRED

step("7.4  Silhouette analysis")

sil <- cluster::silhouette(final_clustering, gower_dist)
sil_by_cluster <- tibble::tibble(
  cluster    = as.integer(sil[, "cluster"]),
  neighbor   = as.integer(sil[, "neighbor"]),
  sil_width  = as.numeric(sil[, "sil_width"])) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n = dplyr::n(),
                   avg_silhouette = round(mean(sil_width), 3),
                   min_silhouette = round(min(sil_width), 3),
                   n_negative = sum(sil_width < 0),
                   pct_negative = round(100 * mean(sil_width < 0), 1),
                   modal_neighbour = as.integer(names(which.max(table(neighbor)))),
                   .groups = "drop")
print(as.data.frame(sil_by_cluster), row.names = FALSE)
save_table(sil_by_cluster, "final_silhouette_by_cluster.csv",
           "Silhouette width summarised per cluster")
ok("Overall average silhouette width: ",
   round(mean(sil[, "sil_width"]), 3))

#  7.5  Labels -----------------------------------------------------------------
#  WHAT    Attaches substantive names, and exports the evidence for them.
#  WHY     Naming is interpretation, and it is where a typology acquires the
#          meaning it will carry through the rest of the literature. Two
#          cautions. A label is a claim, so it should be traceable to the
#          profile: use 7.3's distinctive-features table as the justification.
#          And a label can drift -- cluster numbering is arbitrary and changes
#          if the data or the seed change, so a label hard-coded against a
#          cluster NUMBER will silently attach to the wrong group after any
#          edit. The verification file exported here exists so that drift is
#          detectable.
#  CHOICE  CLUSTER_LABELS (2.9)
#  STATUS  OPTIONAL

cluster_label <- if (is.null(CLUSTER_LABELS))
  paste("Cluster", seq_len(K_FINAL)) else CLUSTER_LABELS

if (!is.null(CLUSTER_LABELS)) {
  if (anyDuplicated(CLUSTER_LABELS))
    stop("CLUSTER_LABELS contains duplicates.", call. = FALSE)
  label_evidence <- distinctive %>%
    dplyr::mutate(label = cluster_label[cluster]) %>%
    dplyr::relocate(cluster, label)
  save_table(label_evidence, "final_label_verification.csv",
             "Each label beside the features it is supposed to describe")
  note("Check final_label_verification.csv: every label should be readable ",
       "off the features beside it. If it is not, the labels have drifted.")
}

#  7.6  Cluster assignments ----------------------------------------------------
#  WHAT    One row per case: identifier, cluster, label, silhouette width,
#          and the nearest competing cluster.
#  WHY     This is the file every downstream analysis reads, and the file a
#          reader of your deposit needs in order to reproduce anything you did
#          with the typology. Exporting the silhouette width alongside the
#          assignment lets a later analysis exclude or flag ambiguous cases
#          without recomputing the distance matrix.
#  STATUS  REQUIRED

assignments <- tibble::tibble(
  id             = data_analytic[[ID_VAR]],
  cluster        = as.integer(final_clustering),
  cluster_label  = cluster_label[final_clustering],
  silhouette     = round(as.numeric(sil[, "sil_width"]), 4),
  nearest_other  = as.integer(sil[, "neighbor"]),
  is_ambiguous   = as.numeric(sil[, "sil_width"]) < 0)
names(assignments)[1] <- ID_VAR
save_table(assignments, "cluster_assignments.csv",
           "Cluster membership per case, with silhouette width")
ok("Assignments written for ", nrow(assignments), " cases.")

#  7.7  External validation ----------------------------------------------------
#  WHAT    Tests whether the clusters differ on variables that did NOT enter
#          the clustering.
#  WHY     This is where validity comes from. Everything before this section
#          establishes that the partition is internally coherent, stable and
#          not an artefact of noise -- necessary conditions, none of which
#          shows that the typology corresponds to anything outside itself. A
#          typology earns its name by predicting something it was not built
#          from (Rapkin & Luke, 1993).
#
#          HOW TO READ THE NUMBERS. Effect sizes first, p values a distant
#          second, for two reasons that are not pedantry:
#
#          (1) The clusters were derived from this same sample, so the groups
#              are not independent of the data the test uses. The p value does
#              not have its usual long-run interpretation, and it is not
#              adjusted for the search over k that produced the groups.
#          (2) Omega squared is reported rather than eta squared because eta
#              squared rises mechanically with the number of groups, which
#              would flatter a six-cluster solution against a three-cluster
#              one for no substantive reason.
#
#          Treat a large effect on an external variable as the evidence, and
#          the p value as a formality that a reviewer will expect to see.
#  CHOICE  EXTERNAL_VARS (2.10)
#  STATUS  OPTIONAL -- and the most important thing you can add

external_results <- NULL

if (length(EXTERNAL_VARS) > 0) {
  step("7.7  External validation")

  omega_squared <- function(y, g) {
    fit <- stats::aov(y ~ g)
    s   <- summary(fit)[[1]]
    ss_b <- s[1, "Sum Sq"]; ss_w <- s[2, "Sum Sq"]
    df_b <- s[1, "Df"];     ms_w <- s[2, "Mean Sq"]
    (ss_b - df_b * ms_w) / (ss_b + ss_w + ms_w)
  }

  rows <- list()
  for (v in EXTERNAL_VARS) {
    x  <- data_analytic[[v$name]]
    g  <- factor(final_clustering)
    ok_i <- !is.na(x)
    if (sum(ok_i) < 3 * K_FINAL) {
      warn("Too few complete observations on '", v$name, "'; skipped.")
      next
    }
    xi <- x[ok_i]; gi <- droplevels(g[ok_i])

    if (v$type == "interval") {
      wt <- stats::oneway.test(xi ~ gi, var.equal = FALSE)   # Welch
      kw <- stats::kruskal.test(xi, gi)
      eps <- (unname(kw$statistic) - nlevels(gi) + 1) /
             (length(xi) - nlevels(gi))
      rows[[length(rows) + 1]] <- tibble::tibble(
        variable = v$name, type = "interval", n = length(xi),
        test = "Welch one-way ANOVA",
        statistic = round(unname(wt$statistic), 3),
        df1 = round(unname(wt$parameter[1]), 1),
        df2 = round(unname(wt$parameter[2]), 1),
        p_value = signif(wt$p.value, 3),
        effect_size = round(omega_squared(xi, gi), 4),
        effect_name = "omega squared",
        secondary = paste0("Kruskal-Wallis epsilon squared = ",
                           round(eps, 4)))
    } else {
      xf  <- if (is.factor(xi)) xi else factor(xi)
      tab <- table(gi, xf)
      cs  <- suppressWarnings(stats::chisq.test(tab))
      vcr <- sqrt(unname(cs$statistic) /
                  (sum(tab) * (min(dim(tab)) - 1)))
      rows[[length(rows) + 1]] <- tibble::tibble(
        variable = v$name, type = v$type, n = sum(tab),
        test = "Pearson chi-square",
        statistic = round(unname(cs$statistic), 3),
        df1 = unname(cs$parameter), df2 = NA_real_,
        p_value = signif(cs$p.value, 3),
        effect_size = round(vcr, 4), effect_name = "Cramer's V",
        secondary = if (any(cs$expected < 5))
          "warning: expected count below 5 in at least one cell" else "")
      save_table(as.data.frame.matrix(tab) %>%
                   tibble::rownames_to_column("cluster"),
                 paste0("external_crosstab_", v$name, ".csv"),
                 paste0("Cluster by ", v$name))
    }
  }

  external_results <- dplyr::bind_rows(rows)
  print(as.data.frame(external_results), row.names = FALSE)
  save_table(external_results, "external_validation.csv",
             "Cluster differences on variables not used to build the clusters")

  # Multiplicity across the external variables: with several outcomes, some
  # will reach conventional significance by chance. Benjamini-Hochberg
  # controls the false discovery rate; Bonferroni is reported beside it as the
  # stricter reference. Applied once, to this family only.
  if (nrow(external_results) > 1) {
    external_results <- external_results %>%
      dplyr::mutate(p_BH = signif(stats::p.adjust(p_value, "BH"), 3),
                    p_bonferroni = signif(stats::p.adjust(p_value, "bonferroni"), 3))
    save_table(external_results, "external_validation.csv",
               "Cluster differences on external variables, with adjusted p values")
  }

  if (nrow(external_results) > 0) {
    top <- which.max(external_results$effect_size)
    note("Largest external effect: ", external_results$variable[top], " (",
         external_results$effect_name[top], " = ",
         external_results$effect_size[top], ").")
    if (max(external_results$effect_size) < 0.01)
      warn("No external variable is meaningfully related to the clusters. The ",
           "partition may be internally coherent without corresponding to ",
           "anything outside itself. Say so, or reconsider the variable set.")
  }
} else {
  note("No external validation variables declared. The solution rests on ",
       "internal evidence alone; see the header of 7.7 for why that is a ",
       "weak position to publish from.")
}

#  7.8  What the missing-data route cost or bought -----------------------------
#  WHAT    Re-clusters the subset of cases that would ALSO have survived
#          listwise deletion, and compares that partition with the retained
#          one on the same cases.
#  WHY     Routes B and C exist to retain cases that route A would delete. The
#          honest question is whether those extra cases changed the typology
#          or merely enlarged it. A high agreement means the missing-data
#          decision was inconsequential and can be reported in one sentence; a
#          low one means it was consequential and must be reported as an
#          analytic choice with an effect.
#  STATUS  OPTIONAL - skipped when no case was retained by proration or
#          imputation

route_sensitivity <- NULL
idx_cc <- which(data_analytic[[ID_VAR]] %in% complete_case_ids)

if (length(idx_cc) < nrow(data_analytic) && length(idx_cc) > 5 * K_FINAL) {
  step("7.8  Sensitivity to the missing-data route")
  dm_cc   <- stats::as.dist(as.matrix(gower_dist)[idx_cc, idx_cc])
  pam_cc  <- cluster::pam(dm_cc, diss = TRUE, k = K_FINAL)
  ari_cc  <- adjusted_rand(pam_cc$clustering, final_clustering[idx_cc])
  route_sensitivity <- tibble::tibble(
    n_analytic            = nrow(data_analytic),
    n_complete_case       = length(idx_cc),
    n_retained_by_route   = nrow(data_analytic) - length(idx_cc),
    adjusted_rand_index   = round(ari_cc, 3),
    reading = dplyr::case_when(
      ari_cc >= 0.90 ~ "the route did not change the typology",
      ari_cc >= 0.70 ~ "broadly the same typology, some cases reassigned",
      TRUE           ~ "the route changed the typology; report it as a choice"))
  print(as.data.frame(route_sensitivity), row.names = FALSE)
  save_table(route_sensitivity, "sensitivity_missing_data_route.csv",
             "Retained solution vs the same analysis on complete cases only")
} else {
  note("Every retained case is a complete case; no route sensitivity to test.")
}

#  7.9  How much does the imputation drive the typology? -----------------------
#  WHAT    Compares the partitions obtained from each imputed dataset.
#  WHY     See 4.3, route C. Partitions cannot be pooled, but their AGREEMENT
#          can be measured, and that agreement is the thing a reader needs in
#          order to know how much of your typology came out of the imputation
#          model. It is the answer to the question multiple imputation makes
#          unavoidable and most applications leave unasked.
#  STATUS  OPTIONAL - runs only under MISSING_METHOD = "imputation"

mi_sensitivity <- NULL

if (!is.null(mi_partitions) && length(mi_partitions) > 1) {
  step("7.9  Sensitivity to the imputation")
  m <- length(mi_partitions)
  pairs <- utils::combn(m, 2)
  aris <- apply(pairs, 2, function(ij)
    adjusted_rand(mi_partitions[[ij[1]]], mi_partitions[[ij[2]]]))

  co <- Reduce(`+`, lapply(mi_partitions,
                           function(p) outer(p, p, "==") * 1)) / m
  # For each case, how consistently it is placed with the same others.
  case_consistency <- apply(co, 1, function(r) mean(pmax(r, 1 - r)))

  mi_sensitivity <- tibble::tibble(
    m_imputations       = m,
    ari_mean            = round(mean(aris), 3),
    ari_min             = round(min(aris), 3),
    ari_max             = round(max(aris), 3),
    case_consistency_mean = round(mean(case_consistency), 3),
    pct_cases_below_0.90  = round(100 * mean(case_consistency < 0.90), 1),
    reading = dplyr::case_when(
      mean(aris) >= 0.90 ~ "the imputation has little influence on the typology",
      mean(aris) >= 0.70 ~ "moderate influence; report it",
      TRUE               ~ "the typology depends substantially on the imputation"))
  print(as.data.frame(mi_sensitivity), row.names = FALSE)
  save_table(mi_sensitivity, "sensitivity_imputation.csv",
             "Agreement between the partitions of the m imputed datasets")
  save_table(tibble::tibble(id = data_analytic[[ID_VAR]],
                            consistency = round(case_consistency, 3)),
             "imputation_case_consistency.csv",
             "Per-case stability of cluster membership across imputations")
}


#  8.  FIGURES =================================================================
#
#  A caution before the figures. None of them is evidence about the number of
#  clusters. A two-dimensional picture of a six-variable mixed-type space is a
#  projection, and projections lose things: clusters that are cleanly
#  separated in the full space can overlap on the page, and clusters that
#  overlap can look separated. Use these to COMMUNICATE a solution you decided
#  on in section 6, never to decide one. Section 8.1 therefore reports how
#  much of the distance structure the projection actually preserves, so the
#  reader can calibrate what they are looking at.

banner("8.  FIGURES")

plot_data <- tibble::tibble(
  cluster = factor(final_clustering,
                   labels = paste0(seq_len(K_FINAL), ". ", cluster_label)),
  silhouette = as.numeric(sil[, "sil_width"]))

#  8.1  Multidimensional scaling -----------------------------------------------
#  WHAT    Projects the Gower matrix into two and three dimensions, and
#          quantifies how much of the pairwise structure survives.
#  WHY     Gower dissimilarities are not generally Euclidean, so classical MDS
#          produces negative eigenvalues -- dimensions with no geometric
#          interpretation. The proportion of the eigenvalue mass they carry is
#          reported below: when it is large, the space genuinely cannot be
#          drawn, and the picture should be presented with that caveat or not
#          at all.
#  CHOICE  RUN_MDS (2.11)
#  STATUS  OPTIONAL

if (isTRUE(RUN_MDS)) {
  step("8.1  Multidimensional scaling")

  mds <- stats::cmdscale(gower_dist, k = 3, eig = TRUE)
  eig <- mds$eig
  pos <- pmax(eig, 0)
  mds_fit <- tibble::tibble(
    dimension          = 1:3,
    eigenvalue         = round(eig[1:3], 4),
    variance_explained = round(pos[1:3] / sum(pos), 4),
    cumulative         = round(cumsum(pos[1:3]) / sum(pos), 4))
  mds_fit$negative_eigenvalue_mass <-
    round(sum(abs(eig[eig < 0])) / sum(abs(eig)), 4)
  print(as.data.frame(mds_fit), row.names = FALSE)
  save_table(mds_fit, "mds_goodness_of_fit.csv",
             "Variance recovered by the MDS projection")

  if (mds_fit$cumulative[2] < 0.50)
    note("The two-dimensional projection recovers only ",
         round(100 * mds_fit$cumulative[2]), "% of the distance structure. ",
         "Present it with that figure attached, and prefer 8.2.")

  mds_df <- cbind(plot_data,
                  as.data.frame(mds$points[, 1:3]) %>%
                    stats::setNames(c("Dim1", "Dim2", "Dim3")))
  med_df <- mds_df[pam_final$id.med, , drop = FALSE]

  p_mds <- ggplot(mds_df, aes(Dim1, Dim2, colour = cluster)) +
    geom_point(alpha = 0.55, size = 1.5) +
    stat_ellipse(level = 0.68, linewidth = 0.4, show.legend = FALSE) +
    geom_point(data = med_df, shape = 21, size = 3.4, stroke = 1.1,
               fill = "white", show.legend = FALSE) +
    labs(title = paste0("Cluster solution projected by MDS (k = ", K_FINAL, ")"),
         subtitle = sprintf(
           "Two dimensions recover %.0f%% of the pairwise structure. Rings: medoids.",
           100 * mds_fit$cumulative[2]),
         x = sprintf("Dimension 1 (%.0f%%)", 100 * mds_fit$variance_explained[1]),
         y = sprintf("Dimension 2 (%.0f%%)", 100 * mds_fit$variance_explained[2]),
         colour = NULL) +
    theme_minimal(base_size = 11)
  save_figure(p_mds, "mds_2d.png",
              "Cluster solution in two MDS dimensions", width = 8, height = 5.5)

#  8.2  Pairwise projections of three dimensions -------------------------------
#  A matrix of two-dimensional views is preferred to a rotating or static
#  three-dimensional scatter: a printed 3-D scatter has no depth cues, so
#  the reader cannot tell near from far, and apparent overlap is an artefact
#  of the viewing angle.

pair_grid <- dplyr::bind_rows(
    transform(mds_df, panel = "Dim 1 x Dim 2", x = mds_df$Dim1, y = mds_df$Dim2),
    transform(mds_df, panel = "Dim 1 x Dim 3", x = mds_df$Dim1, y = mds_df$Dim3),
    transform(mds_df, panel = "Dim 2 x Dim 3", x = mds_df$Dim2, y = mds_df$Dim3))

p_pairs <- ggplot(pair_grid, aes(x, y, colour = cluster)) +
    geom_point(alpha = 0.5, size = 1.1) +
    facet_wrap(~ panel, scales = "free") +
    labs(title = "Pairwise projections of three MDS dimensions",
         subtitle = sprintf("Three dimensions recover %.0f%% of the structure",
                            100 * mds_fit$cumulative[3]),
         x = NULL, y = NULL, colour = NULL) +
    theme_minimal(base_size = 11)
  save_figure(p_pairs, "mds_pairwise_projections.png",
              "Three MDS dimensions shown as pairwise projections",
              width = 10, height = 4)
}

#  8.3  Silhouette plot --------------------------------------------------------
#  WHAT    Every case's silhouette width, sorted within cluster.
#  WHY     It shows the shape of the classification: whether a cluster has a
#          solid core with a ragged edge, or is weakly held together
#          throughout. The average alone cannot distinguish the two.
#  STATUS  REQUIRED

sil_df <- tibble::tibble(
  cluster = factor(sil[, "cluster"],
                   labels = paste0(seq_len(K_FINAL), ". ", cluster_label)),
  width   = as.numeric(sil[, "sil_width"])) %>%
  dplyr::arrange(cluster, dplyr::desc(width)) %>%
  dplyr::mutate(index = dplyr::row_number())

p_sil <- ggplot(sil_df, aes(x = index, y = width, fill = cluster)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = mean(sil_df$width), linetype = 2, colour = "grey30") +
  coord_flip() +
  labs(title = paste0("Silhouette widths (k = ", K_FINAL, ")"),
       subtitle = sprintf("Dashed line: overall mean = %.3f.  Bars below zero: cases better matched elsewhere.",
                          mean(sil_df$width)),
       x = NULL, y = "Silhouette width", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
save_figure(p_sil, "final_silhouette.png",
            "Silhouette width of every case, by cluster", width = 8, height = 6)

#  8.4  Profile figures --------------------------------------------------------
#  WHAT    Two panels: the interval variables as a percentage of their
#          observed range, and the categorical composition of each cluster.
#  WHY     Interval variables are shown on a common 0-100 scale because their
#          raw units are incomparable, and because that is precisely the
#          rescaling Gower applied when computing the distances. The figure
#          therefore shows the clusters in the same terms the algorithm saw
#          them.
#  STATUS  REQUIRED

if (nrow(profile_final$interval) > 0) {
  frame_i <- prepare_gower_frame(data_analytic, CLUSTER_VARS)
  heat <- profile_final$interval %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pct_of_range = {
      x <- frame_i[[variable]]
      round(100 * (median - min(x)) / diff(range(x)), 1)
    }) %>% dplyr::ungroup()

  p_heat <- ggplot(heat, aes(x = variable, y = factor(cluster),
                             fill = pct_of_range)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct_of_range)), size = 3.1,
              colour = "grey15") +
    scale_fill_gradient2(midpoint = 50, low = "#4575b4", mid = "#f7f7f7",
                         high = "#d73027", limits = c(0, 100)) +
    scale_y_discrete(labels = paste0(seq_len(K_FINAL), ". ", cluster_label)) +
    labs(title = "Cluster profile: interval variables",
         subtitle = "Cluster median as a percentage of the variable's observed range",
         x = NULL, y = NULL, fill = "% of range") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          panel.grid = element_blank())
  save_figure(p_heat, "final_profile_interval.png",
              "Cluster medians on the interval variables, as % of range",
              width = 8, height = 4.5)
}

if (nrow(profile_final$categorical) > 0) {
  cat_plot_data <- profile_final$categorical %>%
    dplyr::arrange(variable, level_order) %>%
    dplyr::mutate(level = factor(level, levels = unique(level)))
  p_cat <- ggplot(cat_plot_data,
                  aes(x = factor(cluster), y = pct, fill = level)) +
    geom_col(width = 0.8) +
    facet_wrap(~ variable, scales = "free_x") +
    scale_x_discrete(labels = seq_len(K_FINAL)) +
    labs(title = "Cluster profile: categorical variables",
         subtitle = "Composition of each cluster, in percent",
         x = "Cluster", y = "% of cluster", fill = NULL) +
    theme_minimal(base_size = 11)
  save_figure(p_cat, "final_profile_categorical.png",
              "Categorical composition of each cluster",
              width = 9, height = 5)
}

#  8.5  Distance to own medoid -------------------------------------------------
#  WHAT    How far each case sits from the medoid of its own cluster.
#  WHY     A compactness diagnostic the silhouette does not give: a cluster
#          can have good silhouette widths (it is far from the others) while
#          being internally diffuse (its members are far from each other). A
#          long tail here means the cluster is held together by the absence of
#          alternatives rather than by resemblance.
#  STATUS  REQUIRED

dm_full <- as.matrix(gower_dist)
dist_med <- vapply(seq_len(nrow(dm_full)),
                   function(i) dm_full[i, pam_final$id.med[final_clustering[i]]],
                   numeric(1))
p_dist <- ggplot(tibble::tibble(cluster = plot_data$cluster, d = dist_med),
                 aes(x = cluster, y = d, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.7, show.legend = FALSE) +
  labs(title = "Dissimilarity to the cluster medoid",
       subtitle = "A long upper tail indicates an internally diffuse cluster",
       x = NULL, y = "Gower dissimilarity to own medoid") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
save_figure(p_dist, "final_distance_to_medoid.png",
            "Distance from each case to its cluster medoid", width = 8, height = 5)

save_table(tibble::tibble(cluster = as.integer(final_clustering),
                          distance_to_medoid = round(dist_med, 4)) %>%
             dplyr::group_by(cluster) %>%
             dplyr::summarise(n = dplyr::n(),
                              mean = round(mean(distance_to_medoid), 4),
                              median = round(stats::median(distance_to_medoid), 4),
                              max = round(max(distance_to_medoid), 4),
                              .groups = "drop"),
           "final_compactness.csv",
           "Dissimilarity to the medoid, summarised per cluster")


#  9.  REPORTING ===============================================================

banner("9.  REPORTING")

#  9.1  Reporting checklist ----------------------------------------------------
#  WHAT    Writes a markdown file recording every decision the run made, with
#          the value actually used, plus the headline results.
#  WHY     Cluster analysis is a sequence of choices, and a report that
#          describes only the final solution is not reproducible even when the
#          code is public -- the reader cannot tell which choices were made
#          deliberately and which were defaults. This file is generated from
#          the configuration and the results, so it cannot drift out of step
#          with what was actually run, and it is written in an order that maps
#          onto a methods section.
#
#          It is a starting point, not a finished methods section. Edit it.
#  STATUS  REQUIRED

step("9.1  Reporting checklist")

fmt_num <- function(x, d = 3) {
  if (length(x) == 0 || all(is.na(x))) "not computed"
  else formatC(x, format = "f", digits = d)
}

checklist <- c(
"# Reporting checklist",
"",
paste0("Generated automatically by pam_gower_pipeline.R on ", Sys.Date(), "."),
"Every value below is the value the run actually used. Edit the prose,",
"not the numbers.",
"",
"## 1. Sample",
"",
paste0("- Source: `", DATA_PATH, "`"),
paste0("- Cases imported: ", nrow(data_raw)),
paste0("- Eligibility criterion: ",
       if (is.null(ELIGIBILITY_FILTER)) "none applied"
       else paste0("`", paste(deparse(ELIGIBILITY_FILTER), collapse = " "), "`")),
paste0("- Cases analysed: ", nrow(data_analytic)),
paste0("- Attrition table: `tables/attrition_log.csv`"),
"",
"## 2. Variables entering the clustering",
"",
paste0("- ", spec_names(CLUSTER_VARS), " (", spec_types(CLUSTER_VARS),
       ", weight ", spec_weights(CLUSTER_VARS), ")",
       ifelse(nzchar(vapply(CLUSTER_VARS, function(v) v$note, character(1))),
              paste0(" -- ", gsub("\\s+", " ",
                     vapply(CLUSTER_VARS, function(v) v$note, character(1)))),
              "")),
"",
paste0("Variables were entered ",
       if (any(abs(spec_weights(CLUSTER_VARS) - 1) > 1e-8))
         "with the unequal weights listed above"
       else "with equal weight"),
paste0("(share of the distance carried by each: `tables/gower_weight_share.csv`)."),
"",
"## 3. Missing data",
"",
paste0("- Method: ", MISSING_METHOD),
paste0("- Overall missingness in the eligible sample: ",
       overall_missing$pct_missing_cells, "% of cells; ",
       overall_missing$pct_complete_cases, "% of cases complete"),
paste0("- MCAR diagnostic: ",
       if (is.null(mcar_result)) "not computed"
       else paste0("Little's test, chi-square = ",
                   fmt_num(mcar_result$statistic, 2), ", df = ",
                   mcar_result$df, ", p = ", fmt_num(mcar_result$p.value))),
paste0("- Sensitivity to the route: ",
       if (is.null(route_sensitivity)) "not applicable"
       else paste0("ARI = ", route_sensitivity$adjusted_rand_index, " against ",
                   "the complete-case analysis (", route_sensitivity$reading, ")")),
paste0("- Sensitivity to the imputation: ",
       if (is.null(mi_sensitivity)) "not applicable"
       else paste0("mean ARI = ", mi_sensitivity$ari_mean, " across ",
                   mi_sensitivity$m_imputations, " imputations")),
"",
"## 4. Dissimilarity",
"",
"- Gower's general coefficient (Gower, 1971), computed with `cluster::daisy()`.",
paste0("- Measurement treatment applied: ",
       paste(paste0(type_check$variable, "=", type_check$applied_code),
             collapse = ", ")),
paste0("- Interval variables winsorised: ",
       if (isTRUE(WINSORIZE_INTERVAL))
         paste0("yes, at the ", WINSORIZE_PROBS[1] * 100, "th and ",
                WINSORIZE_PROBS[2] * 100, "th percentiles")
       else "no"),
paste0("- Distinct profiles: ", dist_summary$n_distinct_profiles, " among ",
       dist_summary$n_cases, " cases"),
"",
"## 5. Choosing the number of clusters",
"",
paste0("- Algorithm: partitioning around medoids (Kaufman & Rousseeuw, 1990)."),
paste0("- Range examined: k = ", min(K_RANGE), " to ", max(K_RANGE), "."),
paste0("- Candidates examined in detail: k = ",
       paste(k_candidates, collapse = ", "), "."),
paste0("- Highest average silhouette: k = ", best_sil, "."),
paste0("- Bootstrap stability: ",
       if (is.null(stability_summary)) "not computed"
       else paste0("Jaccard (Hennig, 2007), B up to ",
                   max(stability_summary$B_used),
                   "; weakest cluster at the retained k = ",
                   fmt_num(stability_summary$jaccard_min[
                     stability_summary$k == K_FINAL]))),
paste0("- Split-half replication: ",
       if (is.null(split_half)) "not computed"
       else paste0("mean ARI = ",
                   fmt_num(split_half$ari_mean[split_half$k == K_FINAL]),
                   " over ", SPLIT_N, " random splits")),
paste0("- Structure vs a no-structure reference: ",
       if (is.null(permutation_null)) "not computed"
       else paste0("permutation p = ",
                   fmt_num(permutation_null$p_permutation[
                     permutation_null$k == K_FINAL], 4),
                   " over ", PERM_N, " permutations")),
paste0("- Robustness to the algorithm: ",
       if (is.null(algo_robustness)) "not computed"
       else paste0("ARI with alternatives at the retained k: ",
                   paste(algo_robustness$adjusted_rand_vs_pam[
                     algo_robustness$k == K_FINAL], collapse = ", "))),
"",
paste0("**Retained: k = ", K_FINAL, ".**"),
"",
"Rationale as recorded in the configuration:",
"",
paste0("> ", trimws(unlist(strsplit(trimws(K_FINAL_RATIONALE), "\n")))),
"",
"## 6. The solution",
"",
paste0("- Cluster sizes: ",
       paste(paste0(names(cluster_n), " (n = ", as.integer(cluster_n), ")"),
             collapse = ", ")),
paste0("- Average silhouette width: ", fmt_num(mean(sil[, "sil_width"]))),
paste0("- Cases better matched to another cluster: ",
       sum(sil[, "sil_width"] < 0), " (",
       fmt_num(100 * mean(sil[, "sil_width"] < 0), 1), "%)"),
paste0("- Labels: ", paste(cluster_label, collapse = "; ")),
"",
"## 7. External validation",
"",
if (is.null(external_results)) "- No external validation was performed." else
  paste0("- ", external_results$variable, ": ", external_results$test,
         ", ", external_results$effect_name, " = ", external_results$effect_size,
         ", p = ", external_results$p_value),
"",
"## 8. What a reader still needs from you",
"",
"- Why these variables and not others.",
"- What each cluster means, in words, without reference to its number.",
"- Which decisions above you would expect to change the solution if reversed.",
"- Where the data and this script can be obtained.",
"",
"## 9. Software",
"",
paste0("- ", R.version.string),
paste0("- cluster ", utils::packageVersion("cluster"),
       ", fpc ", utils::packageVersion("fpc")),
paste0("- Random seed: ", SEED),
"",
"Full session information: `session_info.txt`.")

writeLines(checklist, file.path(OUTPUT_DIR, "reporting_checklist.md"))
OUTPUT_LOG$files <- c(OUTPUT_LOG$files,
  paste("reporting_checklist.md",
        "Every decision the run made, with the value used", sep = "\t"))
ok("Reporting checklist written to ", OUTPUT_DIR, "/reporting_checklist.md")

#  9.2  Output inventory -------------------------------------------------------
#  WHAT    Lists every file the run produced, with a one-line description.
#  WHY     So that the deposit is navigable by someone who did not write it,
#          and so that you can tell at a glance whether a file in the folder
#          came from this run or an earlier one.
#  STATUS  REQUIRED

inventory <- tibble::tibble(entry = OUTPUT_LOG$files) %>%
  tidyr::separate(entry, into = c("file", "description"),
                  sep = "\t", fill = "right") %>%
  dplyr::arrange(file)
save_table(inventory, "output_inventory.csv",
           "This inventory")
cat("\n  ", nrow(inventory), " files written to ", OUTPUT_DIR, "/\n", sep = "")
print(as.data.frame(inventory), row.names = FALSE)

#  9.3  Session information ----------------------------------------------------
#  WHAT    R and package versions.
#  WHY     Results can shift between package versions -- clustering algorithms
#          change their tie-breaking, indices get corrected. Without this file
#          "reproducible" means "reproducible on an unspecified machine".
#  STATUS  REQUIRED

writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(OUTPUT_DIR, "session_info.txt"))

banner("RUN COMPLETE")
cat("
  Retained solution : k = ", K_FINAL, "
  Cases analysed    : ", nrow(data_analytic), "
  Outputs           : ", OUTPUT_DIR, "/
  Start here        : ", OUTPUT_DIR, "/reporting_checklist.md
", sep = "")
