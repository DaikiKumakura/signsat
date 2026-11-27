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
