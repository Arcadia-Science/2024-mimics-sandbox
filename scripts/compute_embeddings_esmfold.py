#!/usr/bin/env python

import argparse
import torch
import torch.nn as nn
import esm
import pandas as pd
import numpy as np
from openfold.np import residue_constants
from torch.nn.utils.rnn import pad_sequence
from Bio import SeqIO
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute ESMFold embeddings for protein sequences.')
    parser.add_argument('--fasta', type=str, required=True, help='Path to input FASTA file containing protein sequences.')
    parser.add_argument('--output', type=str, required=True, help='Path to output CSV file to save embeddings.')
    parser.add_argument('--batch_size', type=int, default=1, help='Batch size for processing sequences.')
    parser.add_argument('--max_seq_length', type=int, default=1022, help='Maximum sequence length to process.')
    # Add arguments for ESMFold weight paths with default values
    parser.add_argument('--esm_s_combine_path', type=str, default='esm_s_combine_weights.pt', help='Path to esm_s_combine_weights.pt file.')
    parser.add_argument('--esm_s_mlp_state_dict_path', type=str, default='esm_s_mlp_state_dict.pt', help='Path to esm_s_mlp_state_dict.pt file.')
    parser.add_argument('--embedding_state_dict_path', type=str, default='embedding_state_dict.pt', help='Path to embedding_state_dict.pt file.')
    args = parser.parse_args()
    return args

def load_models(device, esm_s_combine_path, esm_s_mlp_state_dict_path, embedding_state_dict_path):
    # Load ESM2 model used in ESMFold
    esm_model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    esm_model = esm_model.eval().to(device)

    # Load the extracted weights from user-supplied paths
    esm_s_combine_weights = torch.load(esm_s_combine_path, map_location=device).to(device)
    esm_s_mlp_state_dict = torch.load(esm_s_mlp_state_dict_path, map_location=device)
    embedding_state_dict = torch.load(embedding_state_dict_path, map_location=device)

    # Reconstruct esm_s_combine as a learnable parameter (if needed)
    esm_s_combine = nn.Parameter(esm_s_combine_weights).to(device)

    # Reconstruct esm_s_mlp with correct dimensions
    embedding_dim = esm_model.embed_dim  # Should be 2560 for esm2_t36_3B_UR50D
    c_s = 1024  # From ESMFold's configuration

    esm_s_mlp = nn.Sequential(
        nn.LayerNorm(embedding_dim),
        nn.Linear(embedding_dim, c_s),
        nn.ReLU(),
        nn.Linear(c_s, c_s),
    ).to(device)
    esm_s_mlp.load_state_dict(esm_s_mlp_state_dict)

    # Reconstruct amino acid embedding layer
    n_tokens = 23  # Correct number of tokens in ESMFold's embedding layer
    pad_idx = 0

    embedding_layer = nn.Embedding(n_tokens, c_s, padding_idx=pad_idx).to(device)
    embedding_layer.load_state_dict(embedding_state_dict)

    return esm_model, alphabet, esm_s_combine, esm_s_mlp, embedding_layer

def seqs_to_aa_indices(sequences, aa_to_idx):
    aa_indices = []
    for seq in sequences:
        indices = [aa_to_idx.get(aa.upper(), residue_constants.restype_num) for aa in seq]
        aa_indices.append(indices)
    return aa_indices

def process_batch(data_batch, esm_model, alphabet, esm_s_combine, esm_s_mlp, embedding_layer, device):
    # Map amino acid letters to indices used in ESMFold
    aa_to_idx = residue_constants.restype_order_with_x  # List of amino acids
    aa_to_idx = {aa: idx for idx, aa in enumerate(aa_to_idx)}

    # Convert sequences to amino acid indices
    sequences_only = [seq for _, seq in data_batch]
    aa_indices = seqs_to_aa_indices(sequences_only, aa_to_idx)

    # Pad sequences to the same length
    aa_indices_padded = pad_sequence(
        [torch.tensor(seq) for seq in aa_indices], batch_first=True, padding_value=0
    ).to(device)

    # Process the sequences through the ESM model
    batch_converter = alphabet.get_batch_converter()
    batch_labels, batch_strs, batch_tokens = batch_converter(data_batch)
    batch_tokens = batch_tokens.to(device)

    # Get representations from all layers
    num_layers = esm_model.num_layers
    repr_layers = list(range(num_layers + 1))

    with torch.no_grad():
        results = esm_model(batch_tokens, repr_layers=repr_layers, return_contacts=False)

    # Stack embeddings from all layers
    embeddings_list = [results["representations"][i] for i in repr_layers]
    embeddings_stack = torch.stack(embeddings_list, dim=2)  # Shape: (B, L, num_layers+1, C)

    # Remove BOS and EOS tokens to align with amino acid indices
    embeddings_stack = embeddings_stack[:, 1:-1, :, :]  # Remove BOS and EOS

    # Adjust esm_s_combine_softmax shape for broadcasting
    esm_s_combine_softmax = esm_s_combine.softmax(dim=0).view(1, 1, -1, 1)  # Shape: (1, 1, num_layers+1, 1)

    # Perform the element-wise multiplication and sum over the layer dimension
    combined_embeddings = (esm_s_combine_softmax * embeddings_stack).sum(dim=2)  # Shape: (B, L, embedding_dim)

    # Process through esm_s_mlp
    processed_embeddings = esm_s_mlp(combined_embeddings)  # Shape: (B, L, c_s)

    # Get amino acid embeddings
    aa_embeddings = embedding_layer(aa_indices_padded)  # Shape: (B, L, c_s)

    # Add amino acid embeddings to processed embeddings
    final_embeddings = processed_embeddings + aa_embeddings  # Shape: (B, L, c_s)

    # Average over the sequence length to get fixed-size vectors
    protein_embeddings = final_embeddings.mean(dim=1)  # Shape: (B, c_s)

    # Convert to numpy arrays
    embeddings_np = protein_embeddings.detach().cpu().numpy()
    sequence_names = [name for name, _ in data_batch]

    return sequence_names, embeddings_np

def main():
    args = parse_arguments()

    # Set device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load models and weights using user-supplied paths
    esm_model, alphabet, esm_s_combine, esm_s_mlp, embedding_layer = load_models(
        device,
        args.esm_s_combine_path,
        args.esm_s_mlp_state_dict_path,
        args.embedding_state_dict_path
    )

    # Prepare to read sequences from FASTA file
    sequences = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        seq_str = str(record.seq).upper()
        # Filter out sequences that are too long
        if len(seq_str) > args.max_seq_length:
            print(f"Skipping sequence {record.id} of length {len(seq_str)} (exceeds max length).")
            continue
        sequences.append((record.id, seq_str))

    total_sequences = len(sequences)
    print(f"Total sequences to process: {total_sequences}")

    # Process sequences in batches
    batch_size = args.batch_size
    embeddings_list = []
    names_list = []

    for i in range(0, total_sequences, batch_size):
        data_batch = sequences[i:i+batch_size]
        print(f"Processing batch {i // batch_size + 1} / {((total_sequences - 1) // batch_size) + 1}")
        try:
            batch_names, batch_embeddings = process_batch(
                data_batch, esm_model, alphabet, esm_s_combine, esm_s_mlp, embedding_layer, device
            )
            embeddings_list.append(batch_embeddings)
            names_list.extend(batch_names)
        except Exception as e:
            print(f"Error processing batch starting at sequence {i}: {e}")
            continue

    # Concatenate all embeddings
    if embeddings_list:
        all_embeddings = np.vstack(embeddings_list)
    else:
        print("No embeddings were generated.")
        return

    # Create a DataFrame to store the embeddings and names
    embeddings_df = pd.DataFrame(all_embeddings)
    embeddings_df.insert(0, 'name', names_list)

    # Save to CSV
    embeddings_df.to_csv(args.output, index=False)
    print(f"Embeddings saved to {args.output}")

if __name__ == "__main__":
    main()
