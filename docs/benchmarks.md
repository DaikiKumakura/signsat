# Measurements on real PKNs (2026-11-20)

Hardware/software: this Windows PC, Python 3.12.4, python-sat 1.9.dev24, single process. **All TF states are
simulated** (consistent random propagation from the structural roots, or random signs); nothing here measures
biological validity. Scripts: `benchmarks/pkn_stats.py`, `real_pkn.py`, `real_query.py`, `random_states.py`
(CSV outputs of the runs quoted below are copied to `docs/*.csv`; random_states.csv was transcribed from the console because the k=32 run did not finish).

## PKN structure (`pkn_stats.csv`)
| PKN | edges (candidate, signed) | nodes | root nodes (source SCCs) | largest SCC | CollecTRI TFs present (of 1,057) |
|---|---|---|---|---|---|
| OmniPath (directed, signed; 13,917 unsigned rows unusable) | 73,197 | 7,154 | 1,914 | 2,211 | 639 |
| SIGNOR 2026-11-20 (human proteins only; 5,314 rows unusable: complexes/families/chemicals/unknown sign) | 6,864 | 3,965 | 1,225 | 625 | 562 |
Every node of both PKNs is reachable from the structural roots. The number of structural roots (1,225 and 1,914,
several complexes among the OmniPath ones) shows that "structural source SCC" is a very permissive definition of an
upstream origin (risk 1 of the proposal).

## Satisfiability check only (`infer(..., classify=False, solver="auto", time_limit=60)`; `docs/benchmark_feasibility_auto.csv`)
Observed TFs k are drawn from CollecTRI regulators; states are consistent by construction, so every sample is SAT.
| PKN | k = 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| SIGNOR | 3/3 SAT, 0.5 s | 3/3, 0.7 s | 3/3, 1.3 s | 3/3, 2.7 s | 3/3, 5.5 s |
| OmniPath | 3/3, 3.0 s | 3/3, 5.3 s | 3/3, 9.9 s | 3/3, 19.5 s (max 41 s) | 0/3 (INCOMPLETE at 60 s) |
(median seconds; 3 samples each). The SAT formulas have up to 3.2 million variables / 9 million clauses (OmniPath, k=16).

## Solver variance (same formula)
OmniPath, k=2, seed 0: glucose4 INCOMPLETE at 60 s; maplechrono SAT 3.0 s; cadical153 SAT 4.8 s; minisat22
INCOMPLETE. SIGNOR k=2 classification, 60 s: glucose4 proved 6 of 3,111 candidate edges possible, maplechrono 488,
cadical153 1,396. No solver wins everywhere, hence `solver="auto"` (maplechrono 10 s -> glucose4 10 s ->
maplechrono 30 s -> cadical153; CaDiCaL cannot be interrupted through PySAT, so its time limit is soft) and
classification is done with cadical153.

## Exhaustive classification (`classify=True`)
Does not finish in 60 s on SIGNOR even for k = 2 (1,396 of 3,111 candidate edges proven possible with cadical153;
necessity not reached). Reported as INCOMPLETE with the proven edges labelled `Possible (partial)`. It finishes on small
or sparse synthetic PKNs (500 nodes / 1,000 edges / 15 TFs: 0.85 s) and on the 7-edge and random tiny tests.

## target / focus queries on SIGNOR (k = 4, 60 s limit, solver auto at the time of the run)
Of the first 6 queries, 1 was answered in 0.06 s (target not part of any possible route), 1 completed in 0.6 s
(`NO_CONNECTING_PATH`), 4 reached the 60 s limit (`complete=False`: 0 or 3 proven edges). Structural pruning
(target must reach the TF, edges must lie below the target) was added after this run and has not been re-benchmarked
on the whole grid; one remaining slow query was traced to the first SAT call (solver variance again).

## Random TF states (a null; `random_states.csv`, SIGNOR, k observed TFs with random signs, 20 samples each)
k=4: 19 SAT / 1 UNSAT; k=8: 20 SAT; k=16: 18 SAT / 2 UNSAT (k=32 not finished). Even random signs are almost always
jointly explainable on this PKN, so all-TF consistency alone discriminates very little at these sizes
(risk 2 of the proposal: it is not "most samples become UNSAT", it is the opposite).
