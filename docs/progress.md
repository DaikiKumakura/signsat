# Progress log

## 2026-11-20 - first implementation (P0 + core of P1)
Scope chosen by the author: P0 + core of P1 (PKN normalisation, TF sign encoding, SAT encoding with connectivity,
Possible/Necessary classification, UNSAT/unmapped diagnostics, brute-force agreement tests, per-sample network
output, target/focus queries, save/load). Environment: venv + pip (no uv on this PC); python-sat 1.9.dev24,
numpy 2.5.5, pandas 3.1.1, scipy 1.18.1, pytest 9.2.0, Python 3.12.4. Nothing is committed and no repository was created.

Done: modules `pkn.py graph.py encoding.py inference.py query.py results.py io.py`; 291 tests pass
(`python -m pytest -q`): the proposal's 7-edge example, 120 random inference cases and 120 random query cases
against brute force, I/O and normalisation tests, lazy-vs-order encoding agreement.

### Decisions made while implementing (not in the proposal)
* Cycle exclusion: binary order encoding per non-trivial SCC of a TF's region + 2-cycle clauses + lazy blocking as
  a safety net. The first version relied on lazy blocking only; models contained dozens to hundreds of spurious
  disjoint cycles, so the order encoding became the default (the lazy mode stays available as
  `order_encoding=False`; both give identical results on the tested cases).
* Solver default Glucose4. On a 200-node / 800-edge synthetic formula Glucose4 stalled with the order encoding,
  while CaDiCaL153 solved it in 1.8 s and MapleChrono in 0.5 s; on the 100-node case Glucose4 was fastest. No
  solver is best everywhere; `infer(..., solver=)` exposes the PySAT name.
* Per-TF edge classes were dropped (cost); only "is the explanation of a TF unique" is kept.
* `target` states are imposed for every target gene; qualification needs one target gene on the path (segment from
  that gene on); `focus` needs one focus gene at/after it. Reference semantics are in `tests/bruteforce.py`.

### Measurements (synthetic random PKNs, NOT biological; `benchmarks/scale.py`, `/tmp` one-off runs)
| nodes / edges / TFs | solver, cycle handling | result |
|---|---|---|
| 100 / 400 / 8 (dense) | glucose4, order | SAT, 3.5 s, 314 possible edges, 10 MB peak (Python) |
| 100 / 400 / 8 | glucose4, lazy only | 1.3 s |
| 100 / 400 / 8 | maplechrono, order | 13.6 s |
| 200 / 800 / 10 (dense) | glucose4 | did not finish in several minutes |
| 500 / 1000 / 15 (sparse) | glucose4, order / lazy | 0.85 s / 0.31 s |
| 1000 / 2500 / 20 | glucose4, order and lazy | did not finish in 100 s |
The generator is crude (random edges, a few roots); the numbers only show that the exhaustive per-edge
classification has a steep scaling wall in this implementation. They say nothing about SIGNOR/OmniPath yet.

### Open (next steps, in the order of risk in the proposal)
1. Test on a real PKN (SIGNOR / OmniPath): TF mapping rate, reachable fraction, UNSAT fraction, size of possible
   networks, runtime. Needs a data download (ask the author first).
2. If runtime is the blocker: cone-of-influence reduction per TF, a better acyclicity encoding, solver choice per
   instance, parallel samples. Do not add top-k / length caps (the proposal forbids them).
3. Loaders for SIGNOR / OmniPath, GraphML for whole networks, query result caching.
4. Decide licence and repository (the author's call); nothing is committed.

## 2026-11-20 (later) - real-PKN feasibility, performance work, local repository
Permission: the author allowed all downloads and said to proceed freely. Fetched OmniPath, CollecTRI (OmniPath) and
SIGNOR release dumps into `data/raw/` (git-ignored; commands and licence note in `docs/data.md`).

Changes made because of what the real-data runs showed:
* Classification rewritten as set cover / dual shrink loops (calls ~ models needed, not edges).
* `classify=False` (satisfiability only), `time_limit=` for `infer` and `subnetwork` with partial-but-labelled results,
  `solver="auto"` portfolio (no single solver is robust), structural pruning in queries (exact necessary conditions on
  the PKN), faster model processing (numpy truth table), per-TF edge classes already dropped earlier.
* Readers for OmniPath / SIGNOR / CollecTRI (`signsat.sources`), MIT licence chosen (rsemverse precedent; PySAT is MIT).
* Tests: 294 pass (added time-limit, classify=False, auto-solver tests; set-cover rewrite still matches brute force).

Findings that matter for the design (details `docs/benchmarks.md`): satisfiability scales to the real PKNs; exhaustive
classification and many queries do not within 60 s; random TF signs are almost always satisfiable on SIGNOR, so the
all-TF consistency condition discriminates little; structural roots are very numerous.

### Open
1. Real TF activities with a known ground truth (e.g. perturbation data with known targets) to test whether the
   method recovers true routes at all; not attempted (needs a dataset choice).
2. Possible fixes for scale if exhaustive output is required: a better acyclicity encoding, per-TF cone reduction,
   parallel solvers (process-level racing), or reporting classification only for queried regions.
3. Re-benchmark queries on the full grid after the structural pruning.
4. Repository: local git only; no remote created.

## 2026-11-22/27 - public repository and first real perturbation dataset
* Public repository created and pushed: https://github.com/DaikiKumakura/signsat (author's go-ahead: "やれ").
* `subnetwork(..., edges=False)` (explained TFs only); queries now use cadical153 (maplechrono needed >89 s on a
  query that cadical answered in 16 s).
* Validation on GSE329052 (MEK inhibitor, 4 patients): `docs/validation.md`. Result: pipeline works; no discrimination
  between real and sign-shuffled TF states, and the wrong target direction is also explained. Open: a better null,
  other TF-activity methods, a kinase-level PKN restriction, more datasets.

## 2026-11-27 - hypothesis test, saturation, max_path_length
Author: "すべて許可するから、全部進めて". Added `root_states`, `max_path_length` (exact; 80 brute-force tests), `edges=False`.
Findings (docs/validation.md follow-up): unrestricted signed paths on SIGNOR/OmniPath-scale PKNs are saturated (any hypothesis explains
any reachable TF with either sign), so the parameter-free design is almost uninformative there; with paths of 2-3 edges the known MEK action
beats reach-matched random targets on the MEK-inhibitor organoid data (exploratory, one dataset). Not done: a second dataset, a kinase-level
PKN restriction beyond the mechanism filters in parity_ambiguity.py, multiplicity-controlled choice of L, per-patient joint SAT at short L.
