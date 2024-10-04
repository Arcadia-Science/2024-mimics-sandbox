#!/usr/bin/env python3

import argparse
import requests

def get_uniprot_ids_from_proteome(proteome_id):
    """
    Retrieve a list of UniProt protein IDs associated with a given proteome ID.

    Parameters:
        proteome_id (str): The UniProt proteome identifier (e.g., 'UP000006826').

    Returns:
        list of str: A list of UniProt protein accession IDs.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"proteome:{proteome_id}",
        "fields": "accession",
        "format": "tsv",
        "size": "500",  # Maximum allowed size per page
    }

    uniprot_ids = []
    next_link = None

    while True:
        if next_link:
            response = requests.get(next_link)
        else:
            response = requests.get(base_url, params=params)

        response.raise_for_status()

        # The response content is a TSV with a header 'Entry'
        lines = response.text.strip().split('\n')
        if len(lines) <= 1:
            break  # No more entries

        # Skip the "Entry" header whenever it is present
        if lines[0] == 'Entry':
            data_lines = lines[1:]
        else:
            data_lines = lines

        uniprot_ids.extend([line.strip() for line in data_lines])

        # Check for a 'next' page in the response headers
        links = response.links
        if 'next' in links:
            next_link = links['next']['url']
        else:
            break

    return uniprot_ids

def main():
    parser = argparse.ArgumentParser(description='Retrieve UniProt IDs from a proteome ID.')
    parser.add_argument('--uniprot-proteome-id',
                        required=True,
                        help='UniProt proteome identifier (e.g., UP000006826)')
    parser.add_argument('--output', required=True, help='Output file to write the UniProt IDs')

    args = parser.parse_args()

    proteome_id = args.uniprot_proteome_id
    output_file = args.output

    uniprot_ids = get_uniprot_ids_from_proteome(proteome_id)

    with open(output_file, 'w') as f:
        for accession in uniprot_ids:
            f.write(f"{accession}\n")

    print(f"Retrieved {len(uniprot_ids)} UniProt IDs and saved to {output_file}")

if __name__ == '__main__':
    main()
