#!/usr/bin/env python3
"""
Downloads and processes mouse phenotype information from OpenTargets database,
converting the data from parquet format to CSV with simplified phenotype information.

Snippet of DataFrame output as CSV:

ensembl,mouse_phenotype
ENSG00000000003,{'behavior/neurological phenotype'}
ENSG00000000005,"{'muscle phenotype', 'skeleton phenotype'}"
ENSG00000000419,"{'behavior/neurological phenotype', 'mortality/aging', 'embryo phenotype'}"

Originally, we downloaded the OpenTargets metadata with the following wget command:

wget --recursive --no-parent --no-host-directories --cut-dirs=8 \
  --accept "*.parquet" \
   ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/mousePhenotypes/

We changed the download in the script below because only a single file records the mouse phenotype
information.
"""

from pathlib import Path

import pandas as pd
import requests
import typer


def download_file(url: str, output_path: Path, chunk_size: int = 8192) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        print(f"{output_path} already exists. Skipping download.")
        return

    print(f"Downloading {url} to {output_path}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with output_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                f.write(chunk)


def process_phenotypes(parquet_path: Path) -> pd.DataFrame:
    """

    Args:
        parquet_path: Path to the input parquet file

    Returns:
        Pandas DataFrame with processed phenotype data
    """
    mouse = pd.read_parquet(parquet_path)

    aggregated_mouse = (
        mouse.groupby("targetFromSourceId")
        .agg(
            {
                "modelPhenotypeClasses": lambda x: list(x),  # each element is a np.ndarray of dicts
                "modelPhenotypeLabel": lambda x: list(x),  # similarly collect these as lists
            }
        )
        .reset_index()
    )

    # After this step, aggregated_mouse["modelPhenotypeClasses"] is a list of arrays of
    # dictionaries. We need to flatten them so that the next explode step works easily.
    # We first explode at the top level, then handle individual elements if needed.

    # `explode` expects the column to be a list-like (which it is), but currently it's a
    # list of numpy arrays. Let's convert each np.ndarray into a Python list of dicts.
    def _convert_arrays_to_lists(arr_of_arrs):
        # arr_of_arrs is a list of numpy arrays (each array corresponding to a row in
        # the original data). We want to concatenate them into one list or keep them as
        # nested lists first since we want to explode. However, we must be careful: Each
        # grouped row now contains multiple original rows merged into a list. For each
        # grouped row, we have something like [array_of_dicts, array_of_dicts, ...].
        # It's probably best to flatten these so we have a single list of dicts per
        # grouped row. If we want to replicate the original Spark logic (which collected
        # all arrays into a single list), we should flatten all arrays into one list:
        flattened = []
        for arr in arr_of_arrs:
            flattened.extend(arr.tolist())
        return flattened

    aggregated_mouse["modelPhenotypeClasses"] = aggregated_mouse["modelPhenotypeClasses"].apply(
        _convert_arrays_to_lists
    )

    # Now each entry of aggregated_mouse["modelPhenotypeClasses"] is a *list of dicts*.

    # Explode `modelPhenotypeClasses` into individual rows.
    final_mouse = aggregated_mouse.explode("modelPhenotypeClasses")

    # Extract only the `label` from each row in `modelPhenotypeClasses`.
    # Each element is now a single dictionary.
    final_mouse["phenotype_label"] = final_mouse["modelPhenotypeClasses"].apply(
        lambda x: x["label"] if pd.notnull(x) else None
    )

    # Group by "targetFromSourceId" again and collect unique phenotype labels.
    final_mouse = (
        final_mouse.groupby("targetFromSourceId")["phenotype_label"]
        .agg(lambda vals: set(v for v in vals if v is not None))
        .reset_index()
        .rename(columns={"targetFromSourceId": "ensembl", "phenotype_label": "mouse_phenotype"})
    )

    return final_mouse


DOWNLOAD_DIR_OPTION = typer.Option(..., help="Directory to download files")
OUTPUT_OPTION = typer.Option(..., help="Path to the output file")


def main(
    download_dir: Path = DOWNLOAD_DIR_OPTION,
    output: Path = OUTPUT_OPTION,
):
    assert download_dir is not None
    assert output is not None

    root = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/"
    file = "mousePhenotypes/part-00000-41f2eb84-c5e3-4d67-b3b7-4bd87c1e23db-c000.snappy.parquet"
    url = f"{root}{file}"

    download_dir = Path(download_dir)
    download_dir.mkdir(exist_ok=True)
    parquet_path = download_dir / "mouse_phenotypes.parquet"

    download_file(url, parquet_path)
    df = process_phenotypes(parquet_path)

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)


if __name__ == "__main__":
    typer.run(main)
