# SignSAT - working notes for coding agents

Read first: `SPEC.md` (definitions and deviations) -> `docs/progress.md` (decisions, measurements, open items).
Package in `src/signsat/`, tests in `tests/` (`python -m pytest -q` in the venv `.venv`; recreate it with
`python -m venv .venv` and `pip install -e . pytest`). `tests/bruteforce.py` is the reference semantics: if a
change alters any result on the random tests, the change is wrong or the spec must be updated deliberately.
Do not add top-k, path-length limits, regularisation or score thresholds (design rule). Do not commit, publish or
create a repository until the author asks. No AI-generation signatures in commits.
