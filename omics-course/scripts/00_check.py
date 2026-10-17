#!/usr/bin/env python3
"""Read-only readiness inventory; standard library only."""
from pathlib import Path
import hashlib
import importlib.util
import platform
import shutil
import subprocess

COURSE = Path(__file__).resolve().parents[1]
REPO = COURSE.parent
print("Python:", platform.python_version())
for name in ("numpy", "pandas", "anndata", "scanpy", "squidpy"):
    print(f"Python {name}: {'available' if importlib.util.find_spec(name) else 'missing'}")
for name in ("airway_counts.csv", "airway_coldata.csv"):
    path = REPO / "data" / "counts" / name
    if path.exists():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        print(f"Input {path}: {path.stat().st_size} bytes; sha256={digest}")
    else:
        print(f"MISSING input: {path}")
for path in (COURSE / "data/pbmc3k/filtered_gene_bc_matrices/hg19/matrix.mtx",
             COURSE / "data/spatial/visium_hne.h5ad"):
    print(f"{path}: {'present (content not validated)' if path.exists() else 'not downloaded'}")
if shutil.which("Rscript"):
    library = str(COURSE / "env/R-library").replace("\\", "/")
    # Pass R code as one argument, without invoking a shell.
    import json
    code = f'if (dir.exists({json.dumps(library)})) .libPaths(c({json.dumps(library)}, .libPaths())); '
    code += 'cat(R.version.string,"\\n"); for(p in c("DESeq2","Seurat","limma","xcms","faahKO","GEOquery","KEGGREST","clusterProfiler","metafor")) cat(p,requireNamespace(p,quietly=TRUE),"\\n")'
    result = subprocess.run(["Rscript", "-e", code], check=False)
    if result.returncode:
        raise SystemExit(result.returncode)
else:
    print("Rscript: missing")
