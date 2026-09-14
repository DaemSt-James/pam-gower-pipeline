# Changelog

All notable changes to this pipeline are recorded here. Versions follow
[semantic versioning](https://semver.org): the major number changes when a
result could change, the minor number when a feature is added compatibly, the
patch number for fixes and documentation.

## [1.0.0] — 2026-09-13

First public release.

### Added
- `pam_gower_pipeline.R`: configuration-driven pipeline for PAM clustering of
  mixed-type data on a Gower dissimilarity.
- Five declared measurement types (interval, ordinal, nominal, symmetric
  binary, asymmetric binary), with an audit that stops the run if `daisy()`
  applies a treatment other than the one declared.
- Three missing-data routes: listwise deletion, layered exclusion with
  proration, and multiple imputation with a partition-agreement sensitivity
  analysis.
- Multi-criteria model selection: internal indices, nesting and adjusted Rand
  index between candidate solutions, cluster sizes, type-aware profiles,
  adaptive bootstrap stability with a Monte Carlo precision stopping rule,
  split-half replication, robustness across four algorithms, and a permutation
  reference distribution for the absence of multivariate structure.
- Two-pass workflow: the run halts after model selection until `K_FINAL` and a
  written rationale are supplied.
- External validation on variables excluded from the clustering.
- Automatically generated reporting checklist and output inventory.
- `make_example_data.R`: synthetic mixed-type dataset with five known groups.
- `docs/decision_guide.md`: whether this method suits your data and question.
