# Phosphoproteomics test protocol (written and committed BEFORE any phosphosite value was analysed)

Date written: 2026-11-28. Only the file structure of the supplementary tables was inspected (column names, 5 example rows, a few
Table S2 rows shown by the loader); no analysis was run.

Data: Sci Rep 2019 (doi 10.1038/s41598-019-47245-x), PRIDE PXD013923, A375 melanoma cells, 30 min treatment with dabrafenib (BRAF),
trametinib (MEK) or SCH772984 (ERK), SILAC ratio versus DMSO, 3 biological replicates per drug, Supplementary Table S1
(`41598_2019_47245_MOESM2_ESM.xlsx`, obtained from the Europe PMC supplementaryFiles service). In the table a missing value is stored
as 0; a replicate value of exactly 0 is therefore treated as "not quantified".
Kinase-substrate pairs: OmniPath enzyme-substrate (phosphorylation, human), downloaded 2026-11-28.
PKN: the SIGNOR PKN used before.

Known drug actions: dabrafenib BRAF = -1; trametinib MAP2K1 = MAP2K2 = -1; SCH772984 MAPK1 = MAPK3 = -1.

Fixed analysis (no changes after seeing results; deviations are reported as such):
1. Sites: gene = first token of "Gene Names", residue = AA, position = Position. Matched to OmniPath substrates by (gene symbol, residue,
   position). The share of table sites matched is reported.
2. Sanity gate: the average logFC of the MAPK1 (T185, Y187) and MAPK3 (T202, Y204) activation-loop sites present in the table must be
   negative (median < 0) for all three drugs, with at least 3 of the 4 sites quantified; otherwise the data are reported as unusable.
3. Kinase score per replicate: mean of the replicate values of its quantified substrate sites; defined only if >= 3 sites are quantified
   in that replicate. A kinase is observed if defined in all 3 replicates of that drug (the kinase score is a change versus DMSO).
4. Observed kinases: in the SIGNOR PKN, reachable from the structural roots, not a root, and not one of the drug's own target nodes
   (BRAF, MAP2K1/2, MAPK1/3, to avoid the trivial case); the K = 15 (or all, if fewer) with the largest |t| (one-sample t of the 3 replicate
   scores); state per replicate = sign of the score.
5. Explanation: a path from the hypothesis nodes (fixed state -1) with the right sign, maximum length L (walk-parity BFS, validated against SAT).
6. Primary test: trametinib, L = 3, statistic = number of observed kinases explained, summed over the 3 replicates. Null: 2,000 random
   node pairs (state -1; each node with >= 2 outgoing edges, not an observed kinase, not MAP2K1/2) whose mean unlimited reach of the
   observed kinases is within a factor 2 of MEK's. Success: share of random pairs with count >= MEK(-1) is < 0.05 AND the share for
   MEK(+1) is >= 0.05.
7. Secondary, reported without claims: dabrafenib (single-node null, 2,000 nodes) and SCH772984 (pair null); L = 2, 4 and unlimited; rank
   of the hypothesis among all single nodes. A claim of consistency across drugs requires the primary criterion to hold in at least 2 of 3
   drugs at L = 3; with three drugs the p-values are not corrected and are descriptive.
8. Known circularity, stated in advance: OmniPath kinase-substrate pairs include literature ERK/MEK substrates and SIGNOR is among its
   sources, so kinase scores and the PKN are not independent of prior knowledge about this pathway.
