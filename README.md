# SignSAT

SAT-based reconstruction of sample-specific signed causal networks from TF activity and a signed prior-knowledge
network (PKN). Status: **0.1.0.dev0, early development** (not released, not on PyPI).

## Installation

```
python -m venv .venv
.venv/Scripts/python -m pip install -e . pytest       # Windows; python-sat, numpy, pandas, scipy are dependencies
.venv/Scripts/python -m pytest -q                      # 291 tests, about 40 s
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

## Documentation

[Usage guide](https://github.com/DaikiKumakura/signsat/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Compared small graphs with brute-force enumeration and benchmarked public signed networks with simulated TF states.
- **Results:** Satisfiability checks ran on SIGNOR and OmniPath at the tested scales.
- **Limitations:** Exhaustive classifications and some queries timed out. These experiments do not validate biological predictions.
- **Reproduction:** [docs/benchmarks.md](https://github.com/DaikiKumakura/signsat/blob/main/docs/benchmarks.md).

## Citation

Kumakura, D. *signsat* [Software]. [GitHub repository](https://github.com/DaikiKumakura/signsat).

## License

This project is licensed under MIT. See [LICENSE](LICENSE) for the license terms.
