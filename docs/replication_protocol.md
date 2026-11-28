# Replication protocol (written and committed BEFORE the replication data were analysed)

Date written: 2026-11-28. Dataset chosen before looking at any expression value: GEO GSE325686, A375 melanoma cells
(BRAF V600E), 36 h, DMSO vs 1 uM vemurafenib, EGR1 wild-type (2 + 2 replicates) and EGR1 knockout (2 + 2 replicates);
raw read counts (`GSE325686_geneID_symbol_readcount.xls.gz`). Known drug action: BRAF inhibition (BRAF = -1).
Why this dataset: different drug and target than the first study (GSE329052, MEK inhibitor), an acute treatment, DMSO controls.
Limitation known in advance: only 2 replicates per group; "units" are replicate pairs (vemurafenib rep k minus DMSO rep k
within a genotype), 4 units in total.

Fixed analysis (no changes after seeing results; any deviation will be reported as such):
1. Expression: log2(CPM + 1) from raw counts, symbols as given (rows whose symbol starts with "MSTRG" are dropped).
2. Sanity gate: the mean paired log2 change of DUSP6, SPRY4, ETV4, ETV5 and DUSP4 under vemurafenib must be negative for at
   least 4 of the 5 genes; otherwise the dataset is reported as unusable and nothing else is analysed.
3. TF scores: both variants of the first study, `lfc` (mean of sign x paired log2 change over CollecTRI targets, >= 8 targets)
   and `z` (z-scored expression). Primary variant: `lfc`.
4. Observed TFs: the 15 TFs (in the SIGNOR PKN, reachable, not structural roots) with the largest |t| of the 4 unit changes
   (one-sample t); state = sign of the unit's change.
5. Hypothesis: roots = {BRAF}, state -1. Explanation of a TF = a path from BRAF with the right sign (SignSAT semantics; walk-parity BFS
   with depth limit, validated against SAT).
6. Primary test: maximum path length L = 3 (chosen from the first study as "kinase -> kinase -> TF"; BRAF -> MAP2K -> MAPK -> TF
   has 3 edges). Statistic: summed over the 4 units, number of observed TFs explained. Null: 2,000 random single nodes with a reach
   (number of the observed TFs reachable without length limit) within a factor 2 of BRAF's, state -1. Primary success criterion:
   share of random nodes with a count >= BRAF -1 is < 0.05 in the `lfc` variant AND the share for BRAF +1 (wrong direction) is >= 0.05.
7. Secondary, reported without claims: L = 2 and L = 4, the `z` variant, MEK (MAP2K1/MAP2K2 = -1) as the downstream hypothesis,
   rank of BRAF among all single-node hypotheses, EGR1 WT and KO reported separately.
8. Multiplicity: one primary test; everything else is descriptive.
