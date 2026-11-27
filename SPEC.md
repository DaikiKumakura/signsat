# SignSAT specification (as implemented in 0.1.0.dev0)

Source of the design: the proposal pasted by the author on 2026-11-20 ("SignSAT: 詳細設計案"). This file fixes the
definitions that the code implements, records every deviation from the proposal, and states what has and has not
been verified. Where this file and the proposal differ, this file describes the code.

## 1. Inference problem

Inputs: a signed PKN (directed edges with sign +-1) and TF states per sample (+1 up, -1 down, NA unknown).
Continuous scores are converted by sign only (>0 -> +1, <0 -> -1, 0 or NA -> unknown); no threshold. Already
binarised input goes through `tf_states` (1 up, 0 down, NA unknown) so that a score of 0 is never confused with
state 0.

For one sample with observed TF set F:

* **Roots.** Nodes of the strongly connected components of the PKN that receive no edge from another component
  (or a user list `roots=`). Their states are free.
* **Explanation path of f in F:** a simple directed path from a root to f (empty if f is itself a root).
  Several candidate paths per TF are all considered; none is chosen in advance.
* **Explanation H:** the union of one path per observed TF such that some assignment z of +-1 to all nodes
  satisfies z_f = y_f for every observed f and z_v = sign(u,v) * z_u on every edge (u,v) of H. A node shared by
  several paths has one state.
* **Possible** edge: in at least one explanation. **Necessary**: in every explanation. **Alternative**:
  possible but not necessary. **Not inferred**: in none (not listed in the edge table).
* Node `possible_states` = states that node takes over the explanations that contain it. Node status:
  Observed (observed TF), Necessary (in every explanation), Alternative.
* TF `explanation_status`: Unique (all explanations use the same path for the TF) or Multiple alternatives.

Sample status: `SAT`, `UNSAT` (no explanation; a minimal conflicting set of TF observations is reported,
deletion-based), `UNREACHABLE` (an observed TF has no path from any root), `NO_OBSERVATIONS`, `INCOMPLETE`
(solver conflict budget exhausted). Observed TFs that are not in the PKN are dropped and reported
(`TF_NOT_IN_PKN`); a sample is never silently thinned in any other way.

`Necessary` is necessity under this PKN and this model, not experimental necessity. The model assumes a single
sufficient path per TF; it does not infer AND-type regulation.

## 2. PKN normalisation (fixed rules)

Gene symbols stripped and upper-cased; identical (source, target, sign) rows merged with all evidence kept
(`;`-joined); a pair with both signs stays as two candidate edges; unknown sign or missing endpoint are not used and
are reported; self loops are kept in the table but never used. Every other column is carried through as evidence.
No edge is dropped for low confidence.

## 3. SAT encoding (`encoding.py`)

Variables per sample: node state z_v; p[f,k] (TF f's path uses edge k); x_k <-> OR_f p[f,k]; u_v (node lies on a
path through an edge). For every observed f: in/out-degree <= 1 on every node of f's region (nodes that are
reachable from a root and reach f), f has no used out-edge and has an in-edge unless it is a root, non-root nodes
satisfy in <-> out, root nodes satisfy in -> out. These constraints force one path plus possibly disjoint cycles.
**Cycles are excluded by an order encoding**: inside every non-trivial SCC of f's region each node gets a binary
position and a used edge must increase it strictly; 2-cycles are excluded directly. A lazy cycle-blocking loop
remains as a safety net (every blocking clause is implied by "simple path", so no explanation is lost).
Observations are passed as assumptions so that an UNSAT core names the conflicting TFs.

Classification (`infer(classify=True)`, default): Possible edges are found by set cover (repeatedly ask for a model
that makes at least one still-uncovered edge true; the final UNSAT call proves the rest impossible), Necessary edges by
the dual (ask for a model that falsifies at least one remaining candidate). The number of SAT calls is the number of
models needed, not the number of edges. The same loops give the node states and node necessity. `classify=False`
only checks satisfiability (status, conflict diagnosis, one witness path per TF); queries do not need the
classification. Same TF pattern in several samples -> computed once.

Budgets: `infer(time_limit=)` and `subnetwork(time_limit=)` stop the search; what is proven so far is returned
and marked (`INCOMPLETE`, `Possible (partial)`, `Subnetwork.complete = False`), never presented as final.
Interruptible solvers stop exactly at the limit; CaDiCaL (PySAT cannot interrupt it) checks the limit between calls.
`solver="auto"` (default) runs the satisfiability check with maplechrono, glucose4, maplechrono, then cadical153 and
classifies with cadical153; results do not depend on the solver, only the time does.

### Deviations from the proposal
1. The proposal describes "order variables" for cycle removal in general terms; this is a concrete binary
   order encoding restricted to non-trivial SCCs (section 7.2 of the proposal).
2. Per-TF edge classes (Possible/Necessary per TF) are not computed; only whether the explanation of each TF is
   unique (one SAT call per TF). Per-TF edge sets are recomputed inside queries when needed.
3. The SAT model is not serialised; `save()` stores the normalised PKN and TF states and the constraints are
   regenerated on load (the proposal allows this).
4. `Result.diagnostics` also contains PKN rows that were not used (`PKN_*` codes) and TF columns absent from the
   PKN.

## 4. Queries (`query.py`)

`result.subnetwork(sample, target={gene: +-1}, focus=[genes])`, answered on the sample's explanations, never on
the union graph:

* `target`: every target gene is set to its state (an intervention); a path qualifies if it passes through at
  least one target gene; the returned part is the segment from that gene on.
* `focus` only: a path qualifies if it passes through a focus gene; the whole path root -> TF is returned.
* both: the path passes through a target gene and, at or after it, a focus gene (the same gene may be both);
  the segment from the target gene on is returned.
* neither: the whole possible network of the sample.
* Returned: edges (with PKN evidence), nodes and their possible states, explained/unexplained observed TFs,
  matched targets/foci, one witness path per edge. Empty results carry a reason: `GENE_NOT_IN_PKN`,
  `NOT_IN_SAMPLE_NETWORK`, `NO_CONNECTING_PATH`, `DIRECTION_CONFLICT` (feasible once the target states are
  dropped), `INFERENCE_UNAVAILABLE`.

The proposal's example table (section 8 and 11.5) is a unit test; its S1 query with `target={"A": -1}` is reported
as `DIRECTION_CONFLICT` (A must be +1 in S1), which the proposal's table leaves as "not applicable".

## 5. What has been verified (0.1.0.dev0)

* `tests/test_example.py`: the 7-edge example of the proposal (edge classes, node states, TF explanations, all
  query rows), UNSAT / UNREACHABLE / NO_OBSERVATIONS diagnostics, a cycle that must not appear.
* `tests/test_bruteforce.py`: 120 random tiny PKNs (4-7 nodes) against an independent brute-force enumeration
  (all simple paths, all path combinations, all node-state assignments): status, Possible and Necessary edges,
  node states and node necessity; 120 random target/focus queries (edges and explained TFs, empty reasons).
  Of the 240 cases 195 were SAT, 45 UNSAT, 155 had alternatives, 82 needed cycle handling.
* `tests/test_io_pkn.py`: normalisation rules, input rules, same-pattern caching, save/load round trip,
  query after load, GraphML output.
* The lazy-only and the order encodings give identical results on 40 random cases.
Verified by the same author who wrote the encoding; the brute force is independent code but not an independent
author.

## 6. Known limits (measured; full tables in `docs/benchmarks.md`)

* **Scale (real PKNs, simulated TF states).** The satisfiability check works: SIGNOR (6,864 edges) up to 32 TFs in
  about 5 s, OmniPath (73,197 edges) up to 16 TFs in about 20 s (3 of 15 OmniPath runs, all with 32 TFs, hit the 60 s
  limit). **Exhaustive classification does not finish in 60 s even for 2 TFs on SIGNOR**; queries often hit the
  limit. On a sparse synthetic PKN (500 nodes / 1,000 edges / 15 TFs) it took 0.85 s, on dense synthetic ones
  (200 nodes / 800 edges) it was too slow.
* The solver matters enormously (the same formula: 3 s with one solver, over 60 s with another).
* **Discrimination.** With random TF signs on SIGNOR 56 of 60 samples (k = 4, 8, 16) were jointly satisfiable, so
  "all TFs consistent" is a weak test at this PKN density; the structural root set (1,225 nodes in SIGNOR, 1,914 in
  OmniPath) is a very permissive definition of an upstream origin.
* Options added after the real-data work: `infer(root_states={gene: +-1})` imposes states on intervened nodes, `infer(max_path_length=L)`
  limits explanation paths to L edges (default off). Without a limit, 99.9% of reachable (source, node) pairs on SIGNOR are reachable
  with both signs, so path signs hardly discriminate (`docs/validation.md`, `benchmarks/parity_ambiguity.py`).
* Biological validity (earlier note, superseded by the follow-up in validation.md): one real perturbation dataset was tried (`docs/validation.md`): real and sign-shuffled TF states
  gave identical SAT/UNSAT results (SAT up to 20 TFs, UNSAT from 30) and both target directions explained the TFs, so
  no evidence of recovery beyond chance. Queries use cadical153 (soft time limit).
* Complexes are not supported (SIGNOR complexes/families and OmniPath complex nodes are excluded by the readers).
