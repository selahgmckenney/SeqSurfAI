# Runtime Data

SeqSurf AI starts without bundled personal datasets.

Generated study knowledge bases are written to `ai_knowledge_base/`.
GEO expression caches are written to `_geo_imports/`.
Imported app-ready datasets are written to `geo_<accession>/`.

For a fast Compare Studies demo, run:

```r
source("scripts/create_compare_demo_datasets.R")
```

This creates local precomputed fixtures for `geo_gse60450` and `geo_gse16476`
from the app dataset contract. The generated folders stay local by default so
the GitHub repository does not need to store bulky runtime data.
