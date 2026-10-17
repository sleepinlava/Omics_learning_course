#!/usr/bin/env python3
"""Squidpy spatial statistics on its preprocessed Visium H&E tutorial object."""
from pathlib import Path
import argparse
import datetime
import hashlib
import importlib.metadata
import json

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, default=100, help="Top-variance panel size; defines the test family")
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--rings", type=int, default=1)
    args = parser.parse_args()
    if args.genes < 3 or args.permutations < 1 or args.rings < 1:
        parser.error("Require genes>=3, permutations>=1, rings>=1")
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp
        import squidpy as sq
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("Install and activate the spatial environment; see 00_setup/README.md") from exc
    course = Path(__file__).resolve().parents[1]
    path = course / "data/spatial/visium_hne.h5ad"
    if not path.exists():
        raise SystemExit("Fetch data first: python omics-course/scripts/04_fetch_data.py spatial")
    out = course / "results/spatial"
    out.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(path)
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError("Nonunique observation/feature IDs; inspect input rather than silently renaming")
    coords = np.asarray(a.obsm["spatial"])
    if coords.shape != (a.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError("Expected finite 2D coordinates, one row per spot")
    if len(a.uns.get("spatial", {})) > 1:
        raise ValueError("This example requires a single section; multiple libraries need separate spatial graphs")
    x = a.X
    values = x.data if sp.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite expression values")
    # Do not normalize/log a second time: this is a preprocessed tutorial object.
    means = np.asarray(x.mean(axis=0)).ravel()
    second = np.asarray(x.multiply(x).mean(axis=0)).ravel() if sp.issparse(x) else np.mean(x * x, axis=0)
    variance = np.maximum(second - means * means, 0)
    candidates = np.flatnonzero(variance > 0)
    indices = candidates[np.argsort(variance[candidates])[::-1][:args.genes]]
    if len(indices) < 3:
        raise ValueError("Too few variable genes")
    genes = a.var_names[indices].tolist()
    # New mode-specific API (Squidpy 1.8); no expression-space graph substituted.
    sq.gr.spatial_neighbors_grid(a, n_neighs=6, n_rings=args.rings, set_diag=False)
    graph = a.obsp["spatial_connectivities"]
    degree = np.asarray((graph > 0).sum(axis=1)).ravel()
    edge_rows, edge_cols = graph.nonzero()
    euclidean_lengths = np.linalg.norm(coords[edge_rows] - coords[edge_cols], axis=1)
    if not len(euclidean_lengths):
        raise ValueError("Spatial graph has no edges")
    sq.gr.spatial_autocorr(a, genes=genes, mode="moran", n_perms=args.permutations,
                          seed=1234, n_jobs=1, use_raw=False, corr_method="fdr_bh")
    results = a.uns["moranI"]
    results.to_csv(out / "moran_panel.csv")
    pd.DataFrame({"gene": genes, "selection_variance": variance[indices]}).to_csv(out / "tested_panel.csv", index=False)
    pd.DataFrame(coords, index=a.obs_names, columns=["x", "y"]).to_csv(out / "coordinates.csv")
    report = {"n_spots": a.n_obs, "panel_size": len(genes), "directed_nonzero_edges": graph.nnz,
              "isolated_spots": int((degree == 0).sum()), "degree_median": float(np.median(degree)),
              "max_euclidean_edge_length_in_saved_coordinate_units": float(euclidean_lengths.max()),
              "rings": args.rings, "permutations": args.permutations,
              "layers": list(a.layers), "raw_slot_present": a.raw is not None,
              "matrix_state": "upstream preprocessed X; no raw-count claim and no repeated normalization",
              "test_family": "variance-selected panel only; one-sided default; BH within this panel",
              "limitations": "single-section spatial association; no donor-level comparison or causal claims"}
    (out / "graph_qc.json").write_text(json.dumps(report, indent=2) + "\n")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, gene in zip(axes, results.sort_values("I", ascending=False).index[:3]):
        col = a[:, gene].X
        expr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
        dots = ax.scatter(coords[:, 0], coords[:, 1], c=expr, s=6, cmap="viridis")
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(str(gene))
        ax.set_xlabel("Saved spatial x"); ax.set_ylabel("Saved spatial y")
        fig.colorbar(dots, ax=ax, label="Preprocessed expression")
    fig.savefig(out / "spatial_expression.png", dpi=160)
    plt.close(fig)
    a.write_h5ad(out / "analysis.h5ad", compression="gzip")
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(chunk)
    versions = {p: importlib.metadata.version(p) for p in ("squidpy", "anndata", "numpy", "scipy", "matplotlib")}
    versions.update({"time": datetime.datetime.now(datetime.timezone.utc).isoformat(), "input_sha256": sha.hexdigest()})
    (out / "versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    print("Spatial analysis completed:", out)

if __name__ == "__main__":
    main()
