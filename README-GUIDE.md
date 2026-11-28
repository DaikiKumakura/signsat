# SignSAT guide

SignSAT checks whether observed TF changes can be explained at once by sign-consistent paths in a signed prior-knowledge network
(PKN) and lists the possible and necessary edges. Status: **0.1.0.dev0, early development** (not released, not on PyPI).

For every sample the package returns the set of causal networks in which one simple path from an upstream root to each observed TF exists,
all paths are sign-consistent with each other, and each observed TF has the observed state. Edges are classified as Necessary / Alternative
(possible, not necessary); by default nothing is pruned by top-k, path length, regularisation or score thresholds.

**Read the limits first:** on signalling PKNs such as SIGNOR and OmniPath, unrestricted signed paths explain almost anything (99.9% of
reachable pairs have both signs), so the result is not evidence for a drug's mechanism; see the README section "Validation and limitations"
and [SPEC.md](SPEC.md).

<a id="install-for-development"></a>

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
net = res.network("S1")
net.edges                  # Necessary / Alternative edges with signs
net.node_states            # possible states of each node
res.tf_explanations        # which observed TF is explained, unique or with alternatives
res.diagnostics            # unmapped TFs, UNSAT conflicts, unused PKN rows

sub = res.subnetwork("S2", target={"A": 1}, focus=["TF2"])   # query on one sample's explanations
sub.edges, sub.empty_reason, sub.witness
res.save("out/"); ss.Result.load("out/")
```

Inputs: a samples x TFs table (`tf_scores`: sign only; `tf_states`: 1 / 0 / NA) and a PKN with columns `source`,
`target`, `sign` (+ optional evidence columns). Statuses per sample: SAT, UNSAT (with the minimal conflicting
observations), UNREACHABLE, NO_OBSERVATIONS, INCOMPLETE.

## Documentation

[Specification](SPEC.md) · [Data sources](docs/data.md) · [Benchmarks](docs/benchmarks.md)

Real resources: `signsat.sources.omnipath_pkn`, `signor_pkn`, `collectri_tfs` read OmniPath / SIGNOR / CollecTRI dumps
(download commands in `docs/data.md`). `infer(solver="auto")` is the default; `time_limit=` stops the search and
returns what is proven (`INCOMPLETE`, `Possible (partial)`, `Subnetwork.complete=False`).

<a id="what-is-and-is-not-established"></a>

## Validation and limitations

Correctness on small graphs is checked against an independent brute-force enumeration (SPEC.md section 5).
Measured on real PKNs with simulated TF states (`docs/benchmarks.md`): the satisfiability check works at SIGNOR
(6.9k edges) and OmniPath (73k edges) scale, but **exhaustive Possible/Necessary classification does not finish in 60 s
even for 2 TFs on SIGNOR**, queries often hit a 60 s limit, and random TF signs are almost always satisfiable, so
all-TF consistency discriminates little. Biological validity (real TF activities, drug perturbation data) is not
tested. Visualisation and caching beyond identical TF patterns are not implemented.

## Citation

Kumakura, D. *signsat* [Software]. [GitHub repository](https://github.com/DaikiKumakura/signsat).

## License

This project is licensed under MIT. See [LICENSE](LICENSE) for the license terms.
