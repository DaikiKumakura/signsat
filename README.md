# SignSAT

A SAT-based tool for checking whether observed transcription-factor (TF) changes can be explained, all at once, by sign-consistent
paths in a signed prior-knowledge network (PKN), and for listing the possible and necessary edges of such explanations.
Status: **0.1.0.dev0, early development** (not released, not on PyPI).

**What it is:** an exact consistency and explainability checker for small, sparse or hand-restricted signed networks
(SAT/UNSAT per sample with the minimal conflicting observations, Possible / Necessary edges, queries on a drug-action hypothesis).

**What it is not (so far):** a method that identifies a drug's mechanism of action. On the signalling PKNs tested (SIGNOR,
OmniPath), unrestricted signed paths explain almost anything: 99.9% of reachable (source, node) pairs on SIGNOR are reachable with
both signs, so real and sign-shuffled TF states behave identically and every hypothesis explains every reachable TF. Recovery of a
known drug target beyond chance was seen once (MEK inhibitor, paths of at most 2-3 edges) and was **not replicated** in a
pre-registered test on a second dataset (BRAF inhibitor). See [Validation and limitations](#validation-and-limitations).

Where it is informative: networks or edge sets in which feedback and crosstalk are limited (for example transcription-factor
regulation edges only, 4.5% ambiguous pairs), hand-curated pathway models, and sanity checks of a proposed explanation (does a
given set of TF changes contradict the network?).

## Installation

```
python -m venv .venv
.venv/Scripts/python -m pip install -e . pytest       # Windows; python-sat, numpy, pandas, scipy are dependencies
.venv/Scripts/python -m pytest -q                      # 375 tests, about 1 minute
```

## Quick start

```python
import pandas as pd
import signsat as ss

pkn = pd.DataFrame(
    [("R", "A", 1), ("R", "B", 1), ("R", "C", 1), ("A", "TF1", 1), ("A", "TF2", -1), ("B", "TF2", 1), ("C", "TF2", 1)],
    columns=["source", "target", "sign"])
tf = pd.DataFrame({"TF1": [1.2, 0.8, -2.0], "TF2": [0.5, -1.1, None]}, index=["S1", "S2", "S3"])

res = ss.infer(tf_scores=tf, pkn=pkn)          # or tf_states=... with 1 / 0 / NA
# large PKNs: res = ss.infer(..., classify=False, time_limit=60)  # satisfiability only; queries still work
# options: root_states={"MAP2K1": -1} (impose a drug action), max_path_length=3 (default: no limit; benefit unproven)
net = res.network("S1")
net.edges                  # Necessary / Alternative edges with signs
net.node_states            # possible states of each node
res.tf_explanations        # which observed TF is explained, unique or with alternatives
res.diagnostics            # unmapped TFs, UNSAT conflicts, unused PKN rows

sub = res.subnetwork("S2", target={"A": 1}, focus=["TF2"])   # query on one sample's explanations
sub.edges, sub.empty_reason, sub.witness
res.save("out/"); ss.Result.load("out/")
```

For sample `S1` the printed tables are (`net.edges`, then `res.tf_explanations` for all samples):

```text
  sample edge_id source target  sign       status
0     S1      e1      R      A     1    Necessary
1     S1      e2      R      B     1  Alternative
2     S1      e3      R      C     1  Alternative
3     S1      e4      A    TF1     1    Necessary
4     S1      e6      B    TF2     1  Alternative
5     S1      e7      C    TF2     1  Alternative

  sample   tf  observed_state  explained     explanation_status
0     S1  TF1               1       True                 Unique
1     S1  TF2               1       True  Multiple alternatives
2     S2  TF1               1       True                 Unique
3     S2  TF2              -1       True                 Unique
4     S3  TF1              -1       True                 Unique
```

Reading it: in `S1` both TFs are active, so the edge `R -> A -> TF1` is Necessary (the only way to explain TF1),
while TF2 can be reached through `B` or through `C` (Alternative). The edge `A -> TF2` (sign -1) is not in the
network, because it would make TF2 inactive. `S3` has no observation for TF2, so only TF1 is explained.
`res.diagnostics` is empty here: no unmapped TFs and no conflicts.

Inputs: a samples x TFs table (`tf_scores`: sign only; `tf_states`: 1 / 0 / NA) and a PKN with columns `source`,
`target`, `sign` (+ optional evidence columns). Statuses per sample: SAT, UNSAT (with the minimal conflicting
observations), UNREACHABLE, NO_OBSERVATIONS, INCOMPLETE.

## Documentation

[Usage guide](https://github.com/DaikiKumakura/signsat/blob/main/README-GUIDE.md)

## Validation and limitations

- **Correctness (established):** results equal an independent brute-force enumeration on 120 random small graphs (inference), 120 random
  queries and 80 random cases with a path-length limit; the SAT encoding agrees with a lazy-cycle variant.
- **Scale (measured, simulated TF states):** the satisfiability check runs on SIGNOR (6.9k edges) and OmniPath (73k edges) for 16-32
  observed TFs in seconds to under a minute. Exhaustive Possible/Necessary classification does **not** finish in 60 s even for 2 TFs on
  SIGNOR, and many queries reach a 60 s limit (reported as incomplete, never as final).
- **Discrimination (negative):** with unrestricted paths the all-TF consistency test does not separate real from sign-shuffled TF
  states (MEK-inhibitor organoids, GSE329052: both satisfiable up to 20 TFs, both unsatisfiable from 30), and the wrong drug direction is
  "explained" as well.
- **Mechanism recovery (not established):** with `max_path_length` 2-3 the known MEK action explained more TF changes than reach-matched
  random targets (exploratory, one dataset). A pre-registered replication on BRAF-inhibitor-treated A375 cells (GSE325686) failed: BRAF -1 and
  BRAF +1 explained the same number of TFs and 22.6% of random nodes did as well or better.
- **Not tested:** other TF-activity methods, kinase-substrate (phospho) readouts, other PKNs, more datasets, any clinical or drug-response use.
- Details, numbers and the pre-registered protocol: [docs/validation.md](https://github.com/DaikiKumakura/signsat/blob/main/docs/validation.md),
  [docs/replication_protocol.md](https://github.com/DaikiKumakura/signsat/blob/main/docs/replication_protocol.md),
  [docs/benchmarks.md](https://github.com/DaikiKumakura/signsat/blob/main/docs/benchmarks.md).

## Citation

Kumakura, D. *signsat* [Software]. [GitHub repository](https://github.com/DaikiKumakura/signsat).

## License

This project is licensed under MIT. See [LICENSE](LICENSE) for the license terms.
