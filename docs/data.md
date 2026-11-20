# Public data used for the feasibility runs

Nothing is downloaded by the package; the files below were fetched on 2026-11-20 into `data/raw/` (git-ignored,
not redistributed). Check each resource's licence before redistributing derived files.

```
# OmniPath interactions (genesymbols, academic licence), CollecTRI regulators
curl -o data/raw/omnipath_interactions.tsv "https://omnipathdb.org/interactions?genesymbols=1&datasets=omnipath&fields=sources,references,curation_effort&license=academic&format=tsv"
curl -o data/raw/collectri.tsv "https://omnipathdb.org/interactions?genesymbols=1&datasets=collectri&fields=sources&format=tsv"
# SIGNOR latest release (all species, human proteins are selected by the reader)
curl -o data/raw/signor_release.tsv "https://signor.uniroma2.it/releases/getLatestRelease.php"
```
Observed on 2026-11-20 (files as downloaded): OmniPath 86,188 rows (73,690 directed rows with a sign candidate
after splitting; 13,917 directed rows have no sign and are reported), CollecTRI 54,580 rows (1,057 regulators),
SIGNOR 7,877 usable-or-reported rows after reading (5,314 reported as unusable: unsupported entity type or
unknown sign). Readers: `signsat.sources.omnipath_pkn`, `signor_pkn`, `collectri_tfs`.
