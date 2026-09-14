# Reporting checklist — blank template

The pipeline generates a filled-in version of this checklist at
`outputs/reporting_checklist.md` after every complete run. This blank copy is
for two other uses: planning an analysis before you run it, and reviewing
someone else's cluster analysis.

Each item is a decision that changes the solution. A report that omits it is
not reproducible even when the code is public, because the reader cannot tell
which choices were deliberate and which were defaults.

---

## 1. Sample

- [ ] Source of the data, and how it was collected
- [ ] Cases before any exclusion
- [ ] Eligibility criterion, stated as a rule rather than a description
- [ ] Cases analysed
- [ ] Attrition table, from import to analytic sample, one row per criterion
- [ ] Whether the excluded cases differ systematically from the retained ones

## 2. Variables entering the clustering

- [ ] Each variable, with its measurement type
- [ ] **Why each variable is there** — the theoretical basis for including it,
      and for excluding anything obvious you left out
- [ ] Whether any construct is represented by more than one variable, and what
      that does to its weight in the distance
- [ ] Weights, if unequal, and the reasoning behind them
- [ ] Confirmation that the software applied the intended measurement type

## 3. Missing data

- [ ] Extent, reported as both cells missing and cases incomplete
- [ ] Pattern: spread across cases, or concentrated in a few?
- [ ] Evidence about the mechanism, and the limits of that evidence
- [ ] The route taken, and why
- [ ] For proration: the tolerance, and how many scores were prorated
- [ ] For imputation: the model, the auxiliaries, the number of imputations,
      and how the partitions were combined or compared
- [ ] Sensitivity: does the typology change under listwise deletion?

## 4. Dissimilarity

- [ ] The coefficient used, with a citation
- [ ] How each measurement type is handled by it
- [ ] Any transformation applied to the variables first, and its effect
- [ ] How many distinct profiles the data contain, if the variables are largely
      categorical

## 5. Choosing the number of clusters

- [ ] The algorithm, with a citation
- [ ] The range of k examined, and why it stops where it does
- [ ] Every index computed — including the ones that disagreed with your choice
- [ ] Stability under resampling, with the number of replications and how that
      number was fixed
- [ ] Replication in independent subsamples
- [ ] Evidence that the structure is distinguishable from no structure
- [ ] Robustness to the choice of algorithm
- [ ] **The retained k, and the reasoning in prose.** Not "k was chosen by
      silhouette" when the silhouette curve was flat; which criteria you
      weighed most, and what you traded away

## 6. The solution

- [ ] Cluster sizes
- [ ] A profile of each cluster on each clustering variable
- [ ] Classification quality: average silhouette, and cases better matched to
      another cluster
- [ ] The labels, and the features each label is meant to describe
- [ ] Confirmation that no inferential test was run on the variables that built
      the clusters

## 7. External validation

- [ ] The variables used, and confirmation that they did not enter the
      clustering
- [ ] Effect sizes, with the p values secondary
- [ ] Correction for multiplicity across the external variables
- [ ] An honest statement of what the p values mean here, given that the groups
      were derived from the same sample

## 8. What a reader still needs

- [ ] What each cluster means, in words, without reference to its number
- [ ] Which decisions above you would expect to change the solution if reversed
- [ ] Where the data, the code and the outputs can be obtained
- [ ] What the typology is for — what question it was built to answer

## 9. Software

- [ ] R version and the version of every package that touched the analysis
- [ ] The random seed
- [ ] A session-information file deposited with the code
