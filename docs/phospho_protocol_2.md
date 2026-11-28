# Phosphoproteomics test, second kinase-activity estimators (written and committed BEFORE they were run)

Date written: 2026-11-29. Same data and same test as `docs/phospho_protocol.md` (whose primary criterion was NOT met with the
mean-logFC estimator); only the kinase-activity estimator changes. This is a second look at the same data, so the result cannot count
as an independent replication and no p-value is corrected for having tried several estimators.

Added input: protein-kinase list from kinase.com (Manning), via OmniPath annotations, downloaded 2026-11-29 (889 genes). Reason: the first
run produced "kinases" that are not protein kinases (PPP2CA, FGF2, TNF, ...).

Estimators (per drug, per replicate; the site background = all quantified sites in that replicate; a kinase needs >= 3 quantified substrate
sites in the replicate and must be defined in all 3 replicates; only kinases in the kinase.com list):
- B (primary), KSEA z-score: z = (mean substrate logFC - mean of all sites) * sqrt(m) / SD of all sites, m = number of quantified substrate sites.
- C (secondary), mean substrate logFC with the kinase filter only (isolates the effect of removing non-kinases).
- D (secondary), rank-based: Wilcoxon rank-sum z of the substrate-site logFCs versus all other quantified sites (normal approximation).

Everything else exactly as in protocol 1: sanity gate (already passed), observed set = 15 kinases (in the SIGNOR PKN, reachable, not a root,
not a drug target node) with the largest |t| over 3 replicates, state = sign of the replicate score, hypotheses, null of 2,000 reach-matched
random hypotheses (seed 20261129), L = 3 primary, L = 2, 4 and unlimited secondary. Primary criterion: estimator B, trametinib, L = 3:
share of random pairs with count >= MEK(-1) < 0.05 AND share for MEK(+1) >= 0.05. A cross-drug statement requires the criterion in >= 2 of 3
drugs at L = 3 with estimator B. Estimators C and D are descriptive. If B fails, the conclusion "mechanism recovery not established with
kinase-substrate scores" stands and no further estimator is tried on these data.
