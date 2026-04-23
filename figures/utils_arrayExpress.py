# Imports
from __future__ import annotations

from pathlib import Path
from typing import Optional, List, Dict, Any, Set

import anndata as ad
import pandas as pd
import requests
import scipy.sparse as sp


# -----------------------------
# BioStudies metadata helpers
# -----------------------------

def get_biostudies_metadata(accession: str) -> dict:
    """
    Retrieve BioStudies metadata JSON for an accession.

    Example accession:
        E-MTAB-16583
    """
    candidate_urls = [
        f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}",
        f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}/info",
    ]

    last_error = None
    for url in candidate_urls:
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last_error = e

    raise RuntimeError(
        f"Could not retrieve BioStudies metadata for accession {accession}. "
        f"Last error: {last_error}"
    )


def _collect_file_records(obj: Any, out: List[dict]) -> None:
    """
    Recursively collect file-like entries from BioStudies JSON.
    This is intentionally flexible because the exact JSON structure
    may differ between study types / releases.
    """
    if isinstance(obj, dict):
        # A file-like record often contains path / size / attributes
        keys = set(obj.keys())
        if "path" in keys or "fileName" in keys or "relPath" in keys:
            out.append(obj)

        for value in obj.values():
            _collect_file_records(value, out)

    elif isinstance(obj, list):
        for item in obj:
            _collect_file_records(item, out)


def _normalize_blacklist(
    accession: str,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
) -> Set[str]:
    """
    Build the effective blacklist.
    """
    effective_blacklist = list(blacklist or [])
    if auto_blacklist:
        effective_blacklist.extend([
            f"{accession}.sdrf.txt",
            f"{accession}.idf.txt",
            "HTO_library.txt",
            "sample_metadata.txt",
        ])

    return set(effective_blacklist)


def _is_blacklisted(path: str, blacklist: Set[str]) -> bool:
    """
    Returns True if either the full path or the basename is blacklisted.
    """
    if not blacklist:
        return False

    return path in blacklist or Path(path).name in blacklist


def list_accession_files(
    accession: str,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
    *,
    as_dataframe: bool = True,
) -> pd.DataFrame | List[dict]:
    """
    Query BioStudies metadata and return available file entries for a study.

    Returns a DataFrame by default for easy inspection during development.

    Optional: Provide filenames to blacklist (e.g. known metadata files);
    these will be filtered out of the results.
    """
    meta = get_biostudies_metadata(accession)
    raw_files: List[dict] = []
    _collect_file_records(meta, raw_files)

    effective_blacklist = _normalize_blacklist(
        accession=accession,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    cleaned = []
    seen = set()

    for rec in raw_files:
        path = rec.get("path") or rec.get("relPath") or rec.get("fileName")
        if not path:
            continue

        if _is_blacklisted(path, effective_blacklist):
            continue

        row = {
            "path": path,
            "size": rec.get("size"),
            "type": rec.get("type"),
            "md5": rec.get("md5"),
        }

        key = tuple(row.items())
        if key not in seen:
            seen.add(key)
            cleaned.append(row)

    cleaned = sorted(cleaned, key=lambda x: x["path"])

    if as_dataframe:
        return pd.DataFrame(cleaned)

    return cleaned


def print_accession_files(
    accession: str,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
    max_rows: Optional[int] = None,
) -> None:
    """
    Convenience function for development.
    """
    df = list_accession_files(
        accession=accession,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
        as_dataframe=True,
    )

    if df.empty:
        print(f"No file entries found for {accession}.")
        return

    print(f"Files available for {accession}:")
    to_show = df if max_rows is None else df.head(max_rows)
    print(to_show.to_string(index=False))

    if max_rows is not None and len(df) > max_rows:
        print(f"\n... showing first {max_rows} of {len(df)} files.")


# -----------------------------
# Download helpers
# -----------------------------

def get_ebi_ftp_base_url(accession: str) -> str:
    """
    Build the EBI FTP base URL for a BioStudies accession.

    Example
    -------
    E-MTAB-16583 ->
    https://ftp.ebi.ac.uk/pub/databases/biostudies/E-MTAB-/583/E-MTAB-16583/Files/
    """
    suffix = accession.split("-")[-1]
    bucket = suffix[-3:]
    prefix = accession.rsplit("-", 1)[0] + "-"
    return (
        f"https://ftp.ebi.ac.uk/pub/databases/biostudies/"
        f"{prefix}/{bucket}/{accession}/Files/"
    )


def _stream_download(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> None:
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)


def ensure_accession_downloaded(
    accession: str,
    download_dir: str | Path,
    *,
    files_df: Optional[pd.DataFrame] = None,
    force: bool = False,
) -> Path:
    """
    Download all files for an accession into download_dir/accession.

    If files_df is provided, it is used directly and no metadata API call is made.
    """
    study_dir = Path(download_dir).expanduser().resolve() / accession
    study_dir.mkdir(parents=True, exist_ok=True)

    if files_df is None:
        files_df = list_accession_files(
            accession=accession,
            blacklist=None,
            auto_blacklist=False,
            as_dataframe=True,
        )

    if files_df.empty:
        raise RuntimeError(f"No files found for accession {accession}")

    ftp_base = get_ebi_ftp_base_url(accession)

    for relpath in files_df["path"]:
        relpath = str(relpath)
        dest = study_dir / relpath
        dest.parent.mkdir(parents=True, exist_ok=True)

        if dest.exists() and dest.stat().st_size > 0 and not force:
            continue

        url = ftp_base + relpath
        _stream_download(url, dest)

        if not dest.exists() or dest.stat().st_size == 0:
            raise RuntimeError(f"Downloaded empty file from {url}")

    return study_dir


# -----------------------------
# Dataset inference helpers
# -----------------------------

_SUFFIX_MAP = {
    "features": ["_features.txt", "_features.tsv"],
    "log1p": ["_counts_log1p.txt", "_counts_log1p_woCC.txt",],
    "counts": ["_counts.txt"],
    "metadata": ["_metadata.txt", "_metadata.tsv"],
    "umap": ["_umap.txt"],
}


def _strip_known_suffix(filename: str) -> Optional[tuple[str, str]]:
    """
    Return (dataset_prefix, file_role) if the filename matches one of the expected suffixes.
    Otherwise return None.

    Example:
        biobank_tissue_epi_counts.txt -> ("biobank_tissue_epi", "counts")
        treatment_HD4246_SG_features.tsv -> ("treatment_HD4246_SG", "features")
    """
    for role, suffixes in _SUFFIX_MAP.items():
        for suffix in suffixes:
            if filename.endswith(suffix):
                return filename[: -len(suffix)], role
    return None


def _dataset_choice_label(prefix: str) -> str:
    """
    Human-friendly short label derived from the last underscore-separated token.
    Example:
        biobank_tissue_epi -> epi
        biobank_tissue_full -> full
    """
    return prefix.split("_")[-1]


def infer_study_config_from_files(
    accession_dir: str | Path,
    dataset: Optional[str] = None,
    *,
    on_multiple: str = "pick",
    verbose: bool = True,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
) -> Dict[str, str]:
    """
    Infer a dataset-specific file mapping from downloaded files.
    """
    accession_dir = Path(accession_dir)
    accession = accession_dir.name

    if not accession_dir.exists():
        raise FileNotFoundError(f"Directory does not exist: {accession_dir}")

    effective_blacklist = _normalize_blacklist(
        accession=accession,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    grouped: Dict[str, Dict[str, str]] = {}

    for p in accession_dir.rglob("*"):
        if not p.is_file():
            continue

        relpath = str(p.relative_to(accession_dir))

        if _is_blacklisted(relpath, effective_blacklist):
            continue

        parsed = _strip_known_suffix(p.name)
        if parsed is None:
            continue

        prefix, role = parsed

        grouped.setdefault(prefix, {})

        if role in grouped[prefix]:
            raise ValueError(
                f"Multiple files found for dataset '{prefix}' and role '{role}':\n"
                f" - {grouped[prefix][role]}\n"
                f" - {relpath}"
            )

        grouped[prefix][role] = relpath

    if not grouped:
        raise ValueError(
            f"No matching dataset files found in {accession_dir}. "
            f"Expected files ending in: {', '.join(_SUFFIX_MAP.values())}"
        )

    available_prefixes = sorted(grouped.keys())

    if verbose:
        if len(available_prefixes) == 1:
            p = available_prefixes[0]
            print(
                f"Inferred 1 dataset representation: "
                f"{_dataset_choice_label(p)} (full prefix: {p})"
            )
        else:
            print(
                f"Inferred {len(available_prefixes)} dataset representations:\n - "
                + "\n - ".join(
                    f"{_dataset_choice_label(p)} (full prefix: {p})"
                    for p in available_prefixes
                )
            )

    if dataset is None:
        if len(grouped) > 1:
            if on_multiple == "error":
                raise ValueError(
                    "Multiple datasets were found in the study directory. "
                    "Please choose one via the 'dataset' argument.\n"
                    "Available choices:\n - "
                    + "\n - ".join(
                        f"{_dataset_choice_label(p)} (full prefix: {p})"
                        for p in available_prefixes
                    )
                )
            elif on_multiple == "pick":
                selected_prefix = available_prefixes[0]
                if verbose:
                    print(
                        f"No dataset specified. Automatically selected: "
                        f"{_dataset_choice_label(selected_prefix)} "
                        f"(full prefix: {selected_prefix})"
                    )
            else:
                raise ValueError("on_multiple must be either 'pick' or 'error'")
        else:
            selected_prefix = available_prefixes[0]
            if verbose:
                print(
                    f"Selected dataset: {_dataset_choice_label(selected_prefix)} "
                    f"(full prefix: {selected_prefix})"
                )
    else:
        dataset = str(dataset)

        exact_matches = [prefix for prefix in grouped if prefix == dataset]
        short_matches = [prefix for prefix in grouped if _dataset_choice_label(prefix) == dataset]
        matches = exact_matches or short_matches

        if len(matches) == 0:
            raise ValueError(
                f"Dataset '{dataset}' was not found.\n"
                f"Available choices:\n - "
                + "\n - ".join(
                    f"{_dataset_choice_label(p)} (full prefix: {p})"
                    for p in available_prefixes
                )
            )

        if len(matches) > 1:
            raise ValueError(
                f"Dataset selector '{dataset}' is ambiguous. Matching prefixes:\n - "
                + "\n - ".join(sorted(matches))
            )

        selected_prefix = matches[0]
        if verbose:
            print(
                f"Selected dataset: {_dataset_choice_label(selected_prefix)} "
                f"(full prefix: {selected_prefix})"
            )

    selected = grouped[selected_prefix]

    required = ["features", "metadata", "counts"]
    missing = [k for k in required if k not in selected]
    if missing:
        raise ValueError(
            f"Dataset '{selected_prefix}' is missing required file(s): {missing}\n"
            f"Found roles: {sorted(selected.keys())}"
        )

    file_map = {
        "dataset_prefix": selected_prefix,
        "features": selected["features"],
        "metadata": selected["metadata"],
        "counts": selected["counts"],
    }

    if 'log1p' in selected:
        file_map["log1p"] = selected["log1p"]
    if "umap" in selected:
        file_map["umap"] = selected["umap"]

    return file_map

# -----------------------------
# Assembly helpers
# -----------------------------

def _format_relative_path(path: str | Path, base: str | Path) -> str:
    """
    Return path relative to base if possible, otherwise return the original path.
    """
    path = Path(path)
    base = Path(base)

    try:
        return str(path.resolve().relative_to(base.resolve()))
    except Exception:
        return str(path)
    

def get_expected_h5ad_path(
    accession_dir: str | Path,
    dataset_prefix: str,
) -> Path:
    """
    Build the expected output path for the assembled h5ad file.

    Examples
    --------
    accession_dir = ./data/public/E-MTAB-16583

    dataset_prefix = biobank_tissue_epi
    -> ./data/public/E-MTAB-16583/E-MTAB-16583_epi.h5ad

    dataset_prefix = E-MTAB-16583
    -> ./data/public/E-MTAB-16583/E-MTAB-16583.h5ad
    """
    accession_dir = Path(accession_dir)
    accession = accession_dir.name
    dataset_suffix = _dataset_choice_label(dataset_prefix)

    if dataset_suffix and dataset_suffix != accession:
        return accession_dir / f"{accession}_{dataset_suffix}.h5ad"

    return accession_dir / f"{accession}.h5ad"


def list_inferred_datasets(
    accession_dir: str | Path,
    *,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
) -> List[str]:
    """
    Return all dataset prefixes inferred from downloaded study files.
    """
    accession_dir = Path(accession_dir)
    accession = accession_dir.name

    effective_blacklist = _normalize_blacklist(
        accession=accession,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    prefixes = set()

    for p in accession_dir.rglob("*"):
        if not p.is_file():
            continue

        relpath = str(p.relative_to(accession_dir))
        if _is_blacklisted(relpath, effective_blacklist):
            continue

        parsed = _strip_known_suffix(p.name)
        if parsed is None:
            continue

        prefix, _role = parsed
        prefixes.add(prefix)

    return sorted(prefixes)


def assemble_anndata_from_inferred_files(
    accession_dir: str | Path,
    dataset: Optional[str] = None,
    *,
    file_map: Optional[Dict[str, str]] = None,
    make_sparse: bool = True,
    save_adata: bool = True,
    on_multiple: str = "pick",
    verbose: bool = True,
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
) -> ad.AnnData:
    """
    Assemble an AnnData object from inferred study files.
    """
    accession_dir = Path(accession_dir)

    if file_map is None:
        file_map = infer_study_config_from_files(
            accession_dir=accession_dir,
            dataset=dataset,
            on_multiple=on_multiple,
            verbose=verbose,
            blacklist=blacklist,
            auto_blacklist=auto_blacklist,
        )

    features = pd.read_csv(accession_dir / file_map["features"], sep="\t", index_col=0)
    obs = pd.read_csv(accession_dir / file_map["metadata"], sep="\t", index_col=0)
    counts = pd.read_csv(accession_dir / file_map["counts"], sep="\t", index_col=0)
    log1p = None
    if 'log1p' in file_map:
        log1p = pd.read_csv(accession_dir / file_map["log1p"], sep="\t", index_col=0)
    umap = None
    if "umap" in file_map:
        umap = pd.read_csv(accession_dir / file_map["umap"], sep="\t", index_col=0)

    # Assumption: matrices are features x cells -> transpose to cells x features
    counts = counts.T
    log1p = log1p.T if log1p is not None else None
    umap = umap.T if umap is not None else None

    if log1p is not None:
        if not counts.index.equals(log1p.index):
            raise ValueError(
                "Raw counts and log1p matrices do not have matching obs axes after transpose."
            )

        if not counts.columns.equals(log1p.columns):
            raise ValueError(
                "Raw counts and log1p matrices do not have matching var axes after transpose."
            )
    if umap is not None:
        if not counts.index.equals(umap.index):
            raise ValueError(
                "Raw counts and umap matrices do not have matching obs axes after transpose."
            )

    missing_obs = counts.index.difference(obs.index)
    if len(missing_obs) > 0:
        raise ValueError(
            f"Metadata is missing {len(missing_obs)} observation IDs. "
            f"First few: {list(missing_obs[:5])}"
        )
    obs = obs.loc[counts.index].copy()

    # if features.index.equals(counts.columns):
    #     var = features.copy()
    #     var.index = var.index.astype(str)
    # else:
    #     first_col = features.iloc[:, 0].astype(str)
    #     if pd.Index(first_col).equals(counts.columns):
    #         var = features.copy()
    #         var.index = pd.Index(first_col, name="gene")
    #     else:
    #         raise ValueError(
    #             "Could not align features file to matrix columns.\n"
    #             f"First matrix columns: {list(counts.columns[:5])}\n"
    #             f"First features index: {list(features.index[:5])}\n"
    #             f"First features first column: {list(first_col[:5])}"
    #         )
    
    if features.index.equals(counts.columns):
        var = features.copy()
        var.index = var.index.astype(str)

    elif features.shape[1] > 0:
        first_col = features.iloc[:, 0].astype(str)
        if pd.Index(first_col).equals(counts.columns):
            var = features.copy()
            var.index = pd.Index(first_col, name="gene")
        else:
            raise ValueError(
                "Could not align features file to matrix columns.\n"
                f"First matrix columns: {list(counts.columns[:5])}\n"
                f"First features index: {list(features.index[:5])}\n"
                f"First features first column: {list(first_col[:5])}"
            )

    else:
        raise ValueError(
            "Could not align features file to matrix columns. "
            "The features table has no data columns after reading with index_col=0.\n"
            f"First matrix columns: {list(counts.columns[:5])}\n"
            f"First features index: {list(features.index[:5])}"
        )


    counts.columns = counts.columns.astype(str)
    counts_X = counts.to_numpy()
    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)

    if log1p is None:
        print("      Warning: No log1p file found for this dataset. Expression matrix X will contain raw counts.")
        X = counts_X
        log1p_X = None
    else:
        log1p.columns = log1p.columns.astype(str)
        log1p_X = log1p.to_numpy()
        X = log1p_X
        if make_sparse:
            log1p_X = sp.csr_matrix(log1p_X)

    if make_sparse:
        X = sp.csr_matrix(X)
        counts_X = sp.csr_matrix(counts_X)

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
    )

    adata.layers["counts"] = counts_X
    adata.uns["dataset_prefix"] = file_map["dataset_prefix"]
    adata.uns["source_files"] = file_map

    if log1p is not None:
        adata.layers["log1p"] = log1p_X
    if umap is not None:
        missing_umap = adata.obs_names.difference(umap.index.astype(str))
        if len(missing_umap) > 0:
            raise ValueError(
                f"UMAP file is missing {len(missing_umap)} observation IDs. "
                f"First few: {list(missing_umap[:5])}"
            )

        umap = umap.copy()
        umap.index = umap.index.astype(str)
        umap = umap.loc[adata.obs_names]

        adata.obsm["X_umap"] = umap.to_numpy()

    if save_adata:
        out_path = get_expected_h5ad_path(
            accession_dir=accession_dir,
            dataset_prefix=file_map["dataset_prefix"],
        )
        adata.write_h5ad(out_path)
        adata.uns["saved_h5ad"] = str(out_path)

    return adata


def get_public_study_adata(
    accession: str,
    download_dir: str | Path,
    dataset: Optional[str] = None,
    *,
    force_download: bool = False,
    force_rebuild: bool = False,
    save_adata: bool = True,
    make_sparse: bool = True,
    on_multiple: str = "pick",
    blacklist: Optional[List[str]] = None,
    auto_blacklist: bool = True,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Download public ArrayExpress/BioStudies files for one accession and
    return one AnnData object in memory.

    Behavior
    --------
    - If multiple datasets are available and save_adata=True, all datasets are
      assembled and saved to disk.
    - Only one dataset is loaded/returned in memory.
    """
    download_dir_rel = Path(download_dir)
    download_dir = Path(download_dir).expanduser().resolve()
    study_dir = download_dir / accession
    study_dir_rel = _format_relative_path(study_dir, download_dir)

    if verbose:
        print(f"[setup] Download directory: {download_dir_rel}")

    already_present = study_dir.exists() and any(study_dir.iterdir())

    if verbose:
        if already_present:
            msg = "re-downloading files." if force_download else "skipping download."
            print(f"[1/3] Found existing files for {accession} in: {download_dir_rel / study_dir_rel}")
            print(f"      force_download={force_download} -> {msg}")
        else:
            print(f"[1/3] No local files found for {accession}.")
            print("      Downloading accession files.")

    accession_dir = ensure_accession_downloaded(
        accession=accession,
        download_dir=download_dir,
        force=force_download,
    )

    if verbose:
        print(f"[2/3] Inferring dataset structure in: {download_dir_rel / study_dir_rel}")

    dataset_prefixes = list_inferred_datasets(
        accession_dir=accession_dir,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    if not dataset_prefixes:
        raise ValueError(f"No datasets could be inferred in {download_dir_rel / accession_dir}")

    if verbose:
        if len(dataset_prefixes) == 1:
            p = dataset_prefixes[0]
            print(
                f"      Inferred 1 dataset representation: "
                f"{_dataset_choice_label(p)} (full prefix: {p})"
            )
        else:
            print(
                f"      Inferred {len(dataset_prefixes)} dataset representations:\n      - "
                + "\n      - ".join(
                    f"{_dataset_choice_label(p)} (full prefix: {p})"
                    for p in dataset_prefixes
                )
            )

    # Resolve which dataset to load into memory
    if dataset is None:
        if len(dataset_prefixes) > 1:
            if on_multiple == "error":
                raise ValueError(
                    "Multiple datasets were found in the study directory. "
                    "Please choose one via the 'dataset' argument.\n"
                    "Available choices:\n - "
                    + "\n - ".join(
                        f"{_dataset_choice_label(p)} (full prefix: {p})"
                        for p in dataset_prefixes
                    )
                )
            elif on_multiple == "pick":
                selected_prefix = dataset_prefixes[0]
                if verbose:
                    print(
                        f"      No dataset specified. Automatically selected for in-memory loading: "
                        f"{_dataset_choice_label(selected_prefix)} "
                        f"(full prefix: {selected_prefix})"
                    )
            else:
                raise ValueError("on_multiple must be either 'pick' or 'error'")
        else:
            selected_prefix = dataset_prefixes[0]
            if verbose:
                print(
                    f"      Selected dataset for in-memory loading: "
                    f"{_dataset_choice_label(selected_prefix)} "
                    f"(full prefix: {selected_prefix})"
                )
    else:
        # Reuse your existing selector logic by calling inference for that dataset
        selected_file_map = infer_study_config_from_files(
            accession_dir=accession_dir,
            dataset=dataset,
            on_multiple="error",
            verbose=False,
            blacklist=blacklist,
            auto_blacklist=auto_blacklist,
        )
        selected_prefix = selected_file_map["dataset_prefix"]
        if verbose:
            print(
                f"      Selected dataset for in-memory loading: "
                f"{_dataset_choice_label(selected_prefix)} "
                f"(full prefix: {selected_prefix})"
            )

    # Save all datasets if requested
    if save_adata and len(dataset_prefixes) > 1:
        if verbose:
            print("[3/3] Multiple datasets detected.")
            print("      Saving all inferred datasets to disk; loading only one into memory.")

        for prefix in dataset_prefixes:
            out_path = get_expected_h5ad_path(
                accession_dir=accession_dir,
                dataset_prefix=prefix,
            )

            if out_path.exists() and not force_rebuild:
                if verbose:
                    print(
                        f"      Already exists, skipping rebuild: "
                        f"{_format_relative_path(out_path, download_dir)}"
                    )
                continue

            if verbose:
                print(
                    f"      Building dataset: {_dataset_choice_label(prefix)} "
                    f"-> {_format_relative_path(out_path, download_dir)}"
                )

            file_map = infer_study_config_from_files(
                accession_dir=accession_dir,
                dataset=prefix,
                on_multiple="error",
                verbose=False,
                blacklist=blacklist,
                auto_blacklist=auto_blacklist,
            )

            _ = assemble_anndata_from_inferred_files(
                accession_dir=accession_dir,
                dataset=prefix,
                file_map=file_map,
                make_sparse=make_sparse,
                save_adata=True,
                on_multiple="error",
                verbose=False,
                blacklist=blacklist,
                auto_blacklist=auto_blacklist,
            )

    else:
        if verbose:
            print("[3/3] Preparing selected dataset.")

    # Load selected dataset into memory
    selected_h5ad = get_expected_h5ad_path(
        accession_dir=accession_dir,
        dataset_prefix=selected_prefix,
    )

    if selected_h5ad.exists() and not force_rebuild:
        if verbose:
            print(
                f"      Loading selected dataset from disk: "
                f"{_format_relative_path(selected_h5ad, download_dir)}"
            )
        adata = ad.read_h5ad(selected_h5ad)
        adata.uns["loaded_from_h5ad"] = str(selected_h5ad)
        return adata

    if verbose:
        if selected_h5ad.exists() and force_rebuild:
            print("      force_rebuild=True -> rebuilding selected dataset from source files.")
        else:
            print("      No saved h5ad found for selected dataset; assembling from source files.")

    selected_file_map = infer_study_config_from_files(
        accession_dir=accession_dir,
        dataset=selected_prefix,
        on_multiple="error",
        verbose=False,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    adata = assemble_anndata_from_inferred_files(
        accession_dir=accession_dir,
        dataset=selected_prefix,
        file_map=selected_file_map,
        make_sparse=make_sparse,
        save_adata=save_adata,
        on_multiple="error",
        verbose=False,
        blacklist=blacklist,
        auto_blacklist=auto_blacklist,
    )

    if verbose:
        print(f"      Loaded dataset shape: {adata.shape}")
        if save_adata:
            saved_path = adata.uns.get("saved_h5ad", None)
            if saved_path is not None:
                print(
                    f"      Saved selected dataset to: "
                    f"{_format_relative_path(saved_path, download_dir)}"
                )

    return adata


