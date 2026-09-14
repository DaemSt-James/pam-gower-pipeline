# Reporting checklist

Generated automatically by pam_gower_pipeline.R on 2026-09-13.
Every value below is the value the run actually used. Edit the prose,
not the numbers.

## 1. Sample

- Source: `example`
- Cases imported: 900
- Eligibility criterion: `finished == 1 & consent == 1 & age >= 18`
- Cases analysed: 731
- Attrition table: `tables/attrition_log.csv`

## 2. Variables entering the clustering

- efficacy_total (interval, weight 1) -- Collective efficacy, sum of 8 items (see SCALE_DEFS).
- social_support (interval, weight 1) -- Perceived social support, 0-100.
- econ_strain (ordinal, weight 1) -- Ordered severity; spacing between categories not assumed.
- neighborhood_type (nominal, weight 1) -- Unordered. Entered once, not as one-hot indicators.
- civic_member (symm_binary, weight 1) -- Member of a civic association. Symmetric: non-membership is itself a position, so two non-members count as similar.
- service_use (asymm_binary, weight 1) -- Used a specialised service in the past year. Asymmetric: the vast majority did not, and shared non-use carries little information about resemblance.

Variables were entered with equal weight
(share of the distance carried by each: `tables/gower_weight_share.csv`).

## 3. Missing data

- Method: layers
- Overall missingness in the eligible sample: 3.306% of cells; 81.4% of cases complete
- MCAR diagnostic: Little's test, chi-square = 191.29, df = 186, p = 0.380
- Sensitivity to the route: ARI = 1 against the complete-case analysis (the route did not change the typology)
- Sensitivity to the imputation: not applicable

## 4. Dissimilarity

- Gower's general coefficient (Gower, 1971), computed with `cluster::daisy()`.
- Measurement treatment applied: efficacy_total=I, social_support=I, econ_strain=O, neighborhood_type=N, civic_member=S, service_use=A
- Interval variables winsorised: no
- Distinct profiles: 687 among 731 cases

## 5. Choosing the number of clusters

- Algorithm: partitioning around medoids (Kaufman & Rousseeuw, 1990).
- Range examined: k = 2 to 10.
- Candidates examined in detail: k = 5, 6, 7.
- Highest average silhouette: k = 5.
- Bootstrap stability: Jaccard (Hennig, 2007), B up to 800; weakest cluster at the retained k = 0.907
- Split-half replication: mean ARI = 0.908 over 25 random splits
- Structure vs a no-structure reference: permutation p = 0.0196 over 50 permutations
- Robustness to the algorithm: ARI with alternatives at the retained k: 0.767, 0.844, 0.94

**Retained: k = 5.**

Rationale as recorded in the configuration:

> Example dataset; five clusters retained.
> 
> The average silhouette is effectively tied between k = 5 (.580) and k = 6
> (.579), so it does not decide between them, and on stability k = 6 is if
> anything the stronger of the two. The decision rests on the profiles in
> 6.7. The sixth cluster is obtained by splitting the low-efficacy urban
> cluster into service users and non-users; the two halves are otherwise
> indistinguishable -- collective efficacy, social support, economic strain,
> neighbourhood type and civic membership are all within sampling noise of
> each other. That is a subdivision on one indicator, not a distinct profile,
> so k = 5 is retained. Both solutions clear the permutation reference in
> 6.11 and replicate across halves in 6.9.
> 
> (Replace this text with your own reasoning. It is reproduced verbatim in
> the reporting checklist and is part of what you deposit.)

## 6. The solution

- Cluster sizes: 1 (n = 162), 2 (n = 184), 3 (n = 150), 4 (n = 165), 5 (n = 70)
- Average silhouette width: 0.580
- Cases better matched to another cluster: 12 (1.6%)
- Labels: Cluster 1; Cluster 2; Cluster 3; Cluster 4; Cluster 5

## 7. External validation

- distress: Welch one-way ANOVA, omega squared = 0.4617, p = 8.22e-73
- employment_status: Pearson chi-square, Cramer's V = 0.2333, p = 5.96e-14

## 8. What a reader still needs from you

- Why these variables and not others.
- What each cluster means, in words, without reference to its number.
- Which decisions above you would expect to change the solution if reversed.
- Where the data and this script can be obtained.

## 9. Software

- R version 4.4.2 (2024-10-31 ucrt)
- cluster 2.1.6, fpc 2.2.13
- Random seed: 2026

Full session information: `session_info.txt`.
