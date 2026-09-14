#  make_example_data.R =========================================================
#  Synthetic mixed-type dataset with known cluster structure
#
#  Purpose. To let anyone run pam_gower_pipeline.R end to end, and see what
#  every diagnostic looks like, before adapting the pipeline to their own
#  data. The structure is known, so the diagnostics can be checked against the
#  right answer: five groups were generated, and section 6 should recover
#  something close to five.
#
#  The data are invented. They resemble a community survey in shape only, and
#  no substantive conclusion should be drawn from them.
#
#  What the dataset exercises:
#    - all five measurement types the pipeline supports
#    - a scale score computed from items, so proration has something to do
#    - eligibility screening (non-consenting, unfinished and under-age cases)
#    - three kinds of missingness: scattered item omissions, a variable with
#      a heavier rate, and a handful of near-abandoned questionnaires
#    - two external variables, related to the groups but not used to build them
#
#  Used by the pipeline when DATA_PATH is "example". It can also be run on its
#  own to write the file to disk:
#
#      source("make_example_data.R")
#      d <- make_example_data()
#      saveRDS(d, "data/example_data.rds")
#
make_example_data <- function(n = 900, seed = 2026) {

  set.seed(seed)

  # --- Five latent groups, with deliberately different profiles -------------
  # Group 5 is small and only partly distinct: it is there so that the
  # small-cluster and stability diagnostics have something to react to.
  #  1  suburban, high efficacy, well supported, civically engaged
  #  2  urban, moderate on everything, engaged
  #  3  urban, low efficacy, strained, disengaged, high service use
  #  4  rural, moderate-low efficacy, moderately strained
  #  5  inner-suburban, small and only partly distinct from 2
  group_prob <- c(0.24, 0.22, 0.20, 0.24, 0.10)
  g <- sample.int(5, n, replace = TRUE, prob = group_prob)

  # Collective efficacy: eight items on a 1-5 scale.
  eff_mean <- c(4.3, 3.5, 2.3, 3.1, 3.7)[g]
  eff_items <- vapply(seq_len(8), function(i) {
    x <- round(stats::rnorm(n, eff_mean, 0.6))
    pmin(pmax(x, 1), 5)
  }, numeric(n))
  colnames(eff_items) <- paste0("eff_", 1:8)

  # Perceived social support, 0-100.
  social_support <- pmin(pmax(round(stats::rnorm(
    n, c(82, 62, 38, 55, 70)[g], 10)), 0), 100)

  # Economic strain, ordered 1-5. Each group has its own distribution rather
  # than its own mean, so the ordinal treatment has something to work with.
  strain_probs <- list(
    c(.55, .28, .12, .04, .01), c(.18, .34, .30, .13, .05),
    c(.02, .08, .22, .36, .32), c(.08, .20, .32, .27, .13),
    c(.30, .32, .22, .11, .05))
  econ_strain <- vapply(g, function(k)
    sample.int(5, 1, prob = strain_probs[[k]]), integer(1))

  # Neighbourhood type, unordered 1-4. Groups 2 and 3 share a dominant
  # category on purpose: the nominal variable alone cannot separate them, so
  # the other variables have to do the work.
  nb_probs <- list(
    c(.03, .03, .92, .02), c(.92, .03, .03, .02), c(.92, .03, .03, .02),
    c(.02, .03, .03, .92), c(.03, .92, .03, .02))
  neighborhood_type <- vapply(g, function(k)
    sample.int(4, 1, prob = nb_probs[[k]]), integer(1))

  # Civic association membership: symmetric binary.
  civic_member <- stats::rbinom(n, 1, c(.95, .92, .04, .07, .90)[g])

  # Specialised service use: asymmetric binary, deliberately rare overall.
  service_use <- stats::rbinom(n, 1, c(.02, .03, .60, .06, .03)[g])

  # --- External variables: related to the groups, not used to cluster -------
  distress <- pmin(pmax(round(stats::rnorm(
    n, c(8, 14, 25, 18, 13)[g], 5.5)), 0), 40)
  emp_probs <- list(
    c(.86, .05, .09), c(.72, .12, .16), c(.41, .32, .27),
    c(.60, .19, .21), c(.74, .10, .16))
  employment_status <- vapply(g, function(k)
    sample.int(3, 1, prob = emp_probs[[k]]), integer(1))

  # --- Administrative variables, for the eligibility filter -----------------
  age      <- pmin(pmax(round(stats::rnorm(n, 42, 15)), 15), 89)
  consent  <- stats::rbinom(n, 1, 0.99)
  finished <- stats::rbinom(n, 1, 0.96)

  d <- data.frame(
    id = sprintf("P%04d", seq_len(n)),
    consent = consent, finished = finished, age = age,
    eff_items,
    social_support = social_support,
    econ_strain = econ_strain,
    neighborhood_type = neighborhood_type,
    civic_member = civic_member,
    service_use = service_use,
    distress = distress,
    employment_status = employment_status,
    true_group = g,                    # the answer, for checking only
    stringsAsFactors = FALSE)

  # --- Missingness ----------------------------------------------------------
  # (a) scattered item omissions: most affected cases miss a single item, so
  #     proration has a clear job to do and listwise deletion would be
  #     needlessly costly.
  for (j in paste0("eff_", 1:8)) {
    idx <- sample.int(n, round(0.012 * n))
    d[idx, j] <- NA
  }
  # (b) one variable with a heavier rate.
  d[sample.int(n, round(0.05 * n)), "social_support"] <- NA
  d[sample.int(n, round(0.03 * n)), "econ_strain"]    <- NA
  d[sample.int(n, round(0.04 * n)), "distress"]       <- NA

  # (c) a handful of near-abandoned questionnaires, which the overall
  #     missingness layer should remove.
  aband <- sample.int(n, 18)
  d[aband, c(paste0("eff_", 1:8), "social_support", "econ_strain",
             "neighborhood_type", "civic_member", "service_use")] <- NA

  # A stored scale score, computed the way a source file usually would be, so
  # that the pipeline's score cross-check (section 4.6) has something to
  # verify itself against.
  d$efficacy_total <- rowSums(d[, paste0("eff_", 1:8)])

  d
}
