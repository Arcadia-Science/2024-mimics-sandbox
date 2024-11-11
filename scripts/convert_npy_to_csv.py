import argparse

import numpy as np
import pandas as pd


def convert_npy_to_csv(npy_file_path, csv_file_path):
    # Load the .npy file
    embeddings = np.load(npy_file_path)

    # Convert to a pandas DataFrame
    df = pd.DataFrame(embeddings)

    # Save as a .csv file
    df.to_csv(csv_file_path, index=False)

    print(f"Converted {npy_file_path} to {csv_file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a .npy file to a .csv file.")
    parser.add_argument("--npy", type=str, required=True, help="Path to the input .npy file")
    parser.add_argument("--csv", type=str, required=True, help="Path to the output .csv file")

    args = parser.parse_args()

    convert_npy_to_csv(args.npy, args.csv)
