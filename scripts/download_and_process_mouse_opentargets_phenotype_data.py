#!/usr/bin/env python3

import argparse
import os
import sys
import requests
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Process mouse phenotype data.')
    parser.add_argument('--download-dir', default='mousePhenotypes', help='Output directory for the downloaded file.')
    parser.add_argument('--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()
    
    download_dir = args.download_dir
    os.makedirs(download_dir, exist_ok=True)
    
    # wget --recursive --no-parent --no-host-directories --cut-dirs=8 \
    #     --accept "*.parquet" \
    #    ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/mousePhenotypes/

    url = 'https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/etl/parquet/mousePhenotypes/part-00000-41f2eb84-c5e3-4d67-b3b7-4bd87c1e23db-c000.snappy.parquet'
    local_filename = os.path.join(download_dir, 'mouse_phenotypes.parquet')

    if not os.path.exists(local_filename):
        print(f'Downloading {url} to {local_filename}')
        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
        except requests.exceptions.RequestException as e:
            print(f'Error downloading file: {e}')
            sys.exit(1)
    else:
        print(f'{local_filename} already exists. Skipping download.')

    # Initialize Spark session
    spark = SparkSession.builder.appName("LoadParquetFile").getOrCreate()

    mouse = spark.read.parquet(local_filename)

    # Aggregate modelPhenotypeClasses and modelPhenotypeLabel by targetFromSourceId
    aggregated_mouse = mouse.groupBy("targetFromSourceId").agg(
        F.collect_list("modelPhenotypeClasses").alias("modelPhenotypeClasses"),
        F.collect_list("modelPhenotypeLabel").alias("modelPhenotypeLabel"),
    )

    # Explode `modelPhenotypeClasses` to access each individual row
    final_mouse = aggregated_mouse.withColumn("exploded_phenotype", F.explode_outer("modelPhenotypeClasses"))

    # Extract only the `label` from each row in `modelPhenotypeClasses`
    final_mouse = final_mouse.withColumn(
        "phenotype_label",
        F.when(F.col("exploded_phenotype").isNotNull(), F.col("exploded_phenotype.label")).otherwise(None)
    )

    # Collect unique phenotype labels into a new list column `mouse_phenotype`
    final_mouse = final_mouse.groupBy("targetFromSourceId").agg(
        F.collect_set("phenotype_label").alias("mouse_phenotype")
    )

    final_mouse = final_mouse.withColumnRenamed("targetFromSourceId", "ensembl")

    final_mouse_pd = final_mouse.toPandas()
    final_mouse_pd.to_csv(args.output, index=False)

if __name__ == '__main__':
    main()
