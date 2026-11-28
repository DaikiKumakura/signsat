# Validation on real drug-perturbation data (2026-11-22, one dataset)

Dataset: GEO GSE329052, patient-derived colorectal organoids from 4 patients, control vs trametinib (MEK inhibitor)
vs A1155463 (BCL-xL inhibitor), normalised counts. Sanity check passed: paired log2 change under the MEK inhibitor
DUSP6 -3.9, SPRY4 -2.2, ETV4 -6.2, ETV5 -4.0, DUSP4 -3.1 (MAPK output genes fall as expected).
Pipeline (`validation/gse329052_organoids.py`, `gse329052_k_sweep.py`): CollecTRI-weighted TF scores from z-scored
expression, per-patient change versus the patient's own control, observed TFs = top-k by paired |t| among 236 TFs that
are in the SIGNOR PKN and reachable, state = sign of the patient's change, SignSAT on SIGNOR (human proteins),
`classify=False`. Query: `target={"MAP2K1": -1, "MAP2K2": -1}` (the known drug action), count of observed TFs with a
consistent explanation through MEK. Only this dataset was tried (one earlier candidate, GSE296000, was rejected before
use because its MAPK output genes did not fall under trametinib, so it could not serve as ground truth).

## Results (numbers copied from the runs; `validation/*.csv`, `validation/ksweep_k5_k10_console.txt`)
| k observed TFs | real data, 4 patients | same TFs with the patient's signs shuffled across TFs (80 samples) |
|---|---|---|
| 5 | 4 SAT | 80 SAT |
| 10 | 4 SAT | 80 SAT |
| 15 | 4 SAT | 80 SAT |
| 20 | 4 SAT | 80 SAT |
| 30 | 4 UNSAT | 80 UNSAT |
| 40 | 4 UNSAT | 80 UNSAT |
At k = 30 and 40 every patient has the same minimal conflict, MLXIPL -1 together with NRF1 -1: in SIGNOR both are
reached only from PRKAA1 (MLXIPL with sign -1, NRF1 with +1), so the two TFs must have opposite states. This is the
model working as specified, and real data and sign-shuffled data behave identically, i.e. the all-TF consistency
test did not separate the real signal from noise in this dataset.

MEK query (60 s limit; "incomplete" = lower bound):
* k = 5: MEK -1 explains 5/5 TFs in all 4 patients. MEK +1 (wrong direction) also explains 5/5 in 2 patients and at
  least 3/5 in the other 2 (incomplete). The direction of the target therefore did not discriminate.
* k = 10: all 16 queries reached the 60 s limit (MEK -1: 2-5 of 10 proven; MEK +1: 2-4 of 10 proven).
* BCL-xL inhibitor (negative-control drug; 11 observed TFs, 4 SAT): MEK queries did not finish in 20 s.
* MEK inhibitor with 40 TFs (first script): 4 of 4 UNSAT.

## What this does and does not show
Shown: the pipeline runs end to end on a real perturbation dataset; the known target can be queried; real MAPK-axis
TFs can be routed through MEK. Not shown: that SignSAT recovers the true route better than chance. With this PKN
density the explanation through MEK exists with either state of MEK and for sign-shuffled data as well, so neither
all-TF consistency nor the target direction separated truth from noise here. One dataset, four patients, and a
TF-score method chosen by me; a different TF-activity method, a smaller (kinase-level) PKN restriction, or a
different dataset may behave differently. A fair test needs a null that matches the question (for example random
targets compared with MEK at equal path counts), which was not completed because queries exceeded the time limit.

# Follow-up (2026-11-27): hypothesis test, saturation and path length
All scripts: `validation/gse329052_hypothesis_test.py`, `gse329052_depth_test.py`, `benchmarks/parity_ambiguity.py`.
New package features used: `infer(root_states=...)` (impose a drug action on root nodes) and `infer(max_path_length=L)`
(explanations of at most L edges; encoded exactly with distance pruning plus a totalizer cardinality constraint; 80 random
cases agree with brute force, `tests/test_maxlen.py`).

**Hypothesis test.** Hypothesis h = intervened root nodes with imposed states (MAP2K1 = MAP2K2 = -1 for trametinib). For a patient,
n_explained(h) = number of the top-15 observed TFs (paired |t|) for which a path from h gives the observed sign. Null = 2,000
random node pairs (both -1) with a reach comparable to MEK's; also every single node as a hypothesis. Per-TF explainability was
computed by walk-parity BFS, which matched the SAT check on 120 of 120 (TF, state, hypothesis) cases without a length limit and
168 of 168 with limits L = 2, 3, 4 (`tests`-style check in the console notes, MEK and random hypotheses).

1. **No length limit: complete saturation.** Every hypothesis explained all 15 TFs (MEK -1 = MEK +1 = random median = 15 of 15
   in all 4 patients, both TF-score variants, both drugs; MEK ranked 707 of 3,965 single-node hypotheses, i.e. a tie). The
   cause is structural (`benchmarks/parity_ambiguity.py`, `docs/parity_ambiguity.csv`): on SIGNOR 99.9% of reachable
   (source, node) pairs are reachable with BOTH signs (98.3% even for direct phosphorylation/dephosphorylation edges only; 4.5% for
   transcriptional-regulation edges alone). Feedback and crosstalk make any sign possible, so sign consistency over unrestricted
   simple paths carries almost no information on signalling PKNs.
2. **Short paths restore discrimination (exploratory).** Means over 4 patients, TF-score variant `lfc` / `z`:
| L (max edges) | MEKi data: MEK -1 explained (of 15) | MEKi: MEK +1 | random-pair median | share of random pairs >= MEK -1 | BCL-xL-inhibitor data: MEK -1 |
|---|---|---|---|---|---|
| 2 | 5 / 6 | 0 / 1 | 0 | 0.016 / 0.015 | 0 / 2 |
| 3 | 10 / 10 | 6 / 8 | 1 / 2 | 0.015 / 0.019 | 5 / 6 |
| 4 | 13 / 12 | 13 / 10 | 5 / 5 | 0.076 / 0.099 | 11 / 10 |
| 5 | 15 / 14 | 15 / 13 | 11 / 10 | 0.149 / 0.172 | 12 / 12 |
| >= 8 | 15 / 15 | 15 / 15 | 15 / 15 | 0.82-0.94 | 14 / 14 |
   At L = 2-3 the known drug action explains clearly more of the MEK-inhibitor TF changes than reach-matched random pairs
   (share about 0.015-0.02) and the wrong direction does not (share 0.07-0.35); the negative-control drug is explained less at
   L = 2 (0 / 2 TFs). The effect disappears by L = 5. Caveats: one dataset, four patients, L chosen after looking at all values
   (8 length-by-variant comparisons, no multiplicity correction), nulls matched on unrestricted reach only.

**Consequence for the design.** The parameter-free formulation (no path-length limit) is nearly uninformative on signalling PKNs;
`max_path_length` is therefore offered as an explicit, default-off option, and the PKN-ambiguity statistic should be checked before
trusting an explanation result. Short paths mean the canonical cascade (MEK -> ERK -> TF within 2-3 edges) can be tested, which is
consistent with the biology but is also what a plain shortest-path analysis would find.

# Replication attempt (2026-11-28): A375 melanoma, vemurafenib (GSE325686) - NOT replicated
Protocol fixed and committed before analysis: `docs/replication_protocol.md` (commit 7b5994c). Script
`validation/replication_a375_vemurafenib.py` followed it without deviation; console output in
`validation/replication_a375_console.txt`, table `docs/replication_a375_results.csv`.

* Sanity gate passed: paired log2 change under vemurafenib DUSP6 -5.0, SPRY4 -7.2, ETV4 -7.2, ETV5 -5.7, DUSP4 -0.4.
* **Primary criterion (BRAF = -1, L = 3, `lfc`, 4 units, 15 TFs each): NOT met.** BRAF -1 explained 8 of 60 TF observations;
  the wrong direction (BRAF +1) also explained 8 of 60; 22.6% of 2,000 reach-matched random nodes explained at least as many
  (criterion: below 5% for BRAF -1 and at least 5% for BRAF +1). No directional specificity at all.
* Secondary (descriptive): L = 2: BRAF explains 0 of 60 (BRAF -> MAP2K -> MAPK -> TF needs 3 edges); L = 4: 32 of 60 for both
  directions (share of random nodes >= 6.4%); `z` variant L = 3: 4 of 60 (BRAF -1 and +1), share 0.42. MEK = -1 explains 20 of 60 at
  L = 3 (no direction control was run for it). WT and KO strata gave identical counts (the units had the same TF signs).
* Likely reasons (not tested): the strongest TF changes after 36 h of BRAF inhibition are cell-cycle regulators (POU2F2, KDM5A, E2F1,
  FOXM1, ...), which lie more than three edges below BRAF, unlike the immediate MAPK-output TFs that dominated the MEK-inhibitor
  organoids; only 2 replicates per group; selection by |t| of 4 units.

**Conclusion after both datasets.** The short-path discrimination seen in the MEK-inhibitor organoids (GSE329052) was not reproduced
in this dataset under a pre-registered test. With one positive and one negative dataset, the claim that `max_path_length` makes
explanations discriminate is withdrawn as established; it remains an unreplicated, exploratory observation. What stands:
unrestricted signed paths are saturated on SIGNOR (99.9% ambiguous pairs), so sign consistency alone does not identify a mechanism.

## Phosphoproteomics test (PXD013923, A375, 30 min; protocol `docs/phospho_protocol.md` committed before analysis, bef62f0)

Script `validation/phospho_a375.py`, output `validation/phospho_a375_console.txt`, `validation/phospho_a375_results.csv`.
Data: Table S1 of Sci Rep 2019 (doi 10.1038/s41598-019-47245-x), 7,827 sites, SILAC ratios versus DMSO, 3 replicates per drug.
Kinase scores from OmniPath enzyme-substrate pairs (1,347 of 7,827 sites matched, 4,262 pairs); 168-179 kinases defined per drug;
15 observed kinases per drug; null of 2,000 random reach-matched hypotheses.

Facts (from the output):
- Sanity gate passed: the MAPK1/MAPK3 activation-loop sites fall by 2.0 to 3.1 log2 units under all three drugs.
- Primary criterion (trametinib, L = 3): MEK(-1) explains 12 of 45 observed kinase-states (null median 0, max 24); share of random pairs >= 12 is
  0.0595 (not < 0.05); MEK(+1) 3, share 0.4565. **Criterion not met** (by a narrow margin).
- Dabrafenib, L = 3: BRAF(-1) and BRAF(+1) both explain 6; p = 0.1115 for both. No discrimination.
- SCH772984, L = 3: ERK(-1) 27, p = 0.001, but ERK(+1) 24, p = 0.0035: the wrong sign is also significant. At L = 2 (not primary) ERK(-1) 15, p = 0.0005
  versus ERK(+1) 6, p = 0.0855.
- Unlimited path length: every hypothesis explains almost everything (36-39 of 45; null median equals the hypothesis), p = 1.0: the saturation seen
  with TF data is the same here.
- Criterion met at L = 3 in 0 of 3 drugs, so no cross-drug claim.

Interpretation (not established): with kinase-substrate-derived scores, SignSAT does not recover the drug target at the pre-specified depth.
The only sign-specific signal is the secondary SCH772984 L = 2 result, which was not a pre-registered success and was not repeated. Caveats: the
observed "kinases" include non-kinases present in the enzyme-substrate table (e.g. PPP2CA, FGF2, TNF), scores come from few sites, ties between +1/-1
paths are common at L = 3, and prior knowledge (OmniPath/SIGNOR) is not independent of the MAPK pathway. Not done: phosphosite-level nodes,
other kinase-activity estimators, other PKNs.
