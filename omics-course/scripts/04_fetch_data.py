#!/usr/bin/env python3
"""Fetch explicitly selected public tutorial data, retaining provenance."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import tarfile
import urllib.request

COURSE = Path(__file__).resolve().parents[1]
DATA = COURSE / "data"

def receipt(path, source):
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    log = DATA / "download_receipts.json"
    entries = json.loads(log.read_text()) if log.exists() else []
    entries.append({"path": str(path.relative_to(COURSE)), "source": source,
                    "verified_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    "bytes": path.stat().st_size, "sha256": sha.hexdigest(),
                    "note": "Local digest; not a publisher-provided checksum."})
    log.write_text(json.dumps(entries, ensure_ascii=False, indent=2) + "\n")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=["pbmc3k", "spatial"])
    args = parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    if args.dataset == "pbmc3k":
        url = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
        dest = DATA / "pbmc3k"
        dest.mkdir(exist_ok=True)
        archive = dest / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
        if not archive.exists():
            partial = archive.with_suffix(".partial")
            request = urllib.request.Request(url, headers={"User-Agent": "omics-course/1.0"})
            try:
                with urllib.request.urlopen(request, timeout=60) as source, partial.open("wb") as sink:
                    while block := source.read(1024 * 1024):
                        sink.write(block)
                # Validate before committing the download.
                with tarfile.open(partial, "r:gz") as handle:
                    if not handle.getmembers():
                        raise ValueError("Empty archive")
                os.replace(partial, archive)
            except Exception:
                partial.unlink(missing_ok=True)
                raise
        # Extract only expected plain files; never follow symlinks or arbitrary paths.
        required = [f"filtered_gene_bc_matrices/hg19/{x}" for x in ("matrix.mtx", "genes.tsv", "barcodes.tsv")]
        with tarfile.open(archive, "r:gz") as handle:
            members = {m.name.removeprefix("./"): m for m in handle.getmembers()}
            for name in required:
                member = members.get(name)
                if member is None or not member.isfile():
                    raise ValueError(f"Missing expected regular file: {name}")
                path = dest / name
                path.parent.mkdir(parents=True, exist_ok=True)
                if path.exists():
                    if path.stat().st_size != member.size:
                        raise ValueError(f"Existing file has an unexpected size; inspect it before retrying: {path}")
                    print("Keeping existing file:", path)
                    continue
                extracted = handle.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Cannot read {name}")
                partial_file = path.with_suffix(path.suffix + ".partial")
                with extracted, partial_file.open("wb") as target:
                    while block := extracted.read(1024 * 1024):
                        target.write(block)
                os.replace(partial_file, path)
        receipt(archive, url)
        for name in required:
            receipt(dest / name, url + " [archive member: " + name + "]")
        print("PBMC3k ready for QC:", dest)
    else:
        try:
            import squidpy as sq
        except ImportError as exc:
            raise SystemExit("Install env/spatial-requirements.txt in the course venv first.") from exc
        dest = DATA / "spatial"
        dest.mkdir(exist_ok=True)
        path = dest / "visium_hne.h5ad"
        adata = sq.datasets.visium_hne_adata(path=path)
        if "spatial" not in adata.obsm or adata.n_obs < 2:
            raise ValueError("Dataset has no usable spatial coordinates")
        receipt(path, "squidpy.datasets.visium_hne_adata; "
                "https://squidpy.readthedocs.io/en/stable/api/squidpy.datasets.visium_hne_adata.html")
        print("Preprocessed Visium object ready:", adata.shape)

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"Download failed; no successful analysis is claimed: {exc}") from exc
