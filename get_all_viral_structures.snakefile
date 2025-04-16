"""
This snakefile creates standardized datasets of viral protein structures from Viro3D.
It identifies viral structures from viruses that infect humans and downloads them.
"""

from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs/")
INPUT_DIRPATH = Path("inputs")

# Read in host metadata
### we are only retrieving human data in this implementation but I left the structure to retrieve for multiple organisms
host_metadata = pd.read_csv("inputs/viral/host-information.csv", header=0).set_index(
    "organism", drop=False
)
HOST_ORGANISMS = host_metadata["organism"].unique().tolist()


rule download_kegg_virushostdb:
    """
    Downloaded 2024-12-09 version
    MD5sum 19ddaac3ffafc290793a69d367739b88
    This file includes information about which viruses (species, lineage, taxid) infect humans.
    We use the information in this file to determine which subsets of viro3d structures are
    from viruses that infect humans.
    """
    output:
        tsv=INPUT_DIRPATH / "viral" / "virushostdb.tsv",
    shell:
        """
        curl -JLo {output.tsv} https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
        """


rule download_vmr_metadata:
    """
    VMR_MSL38_v2 version downloaded. This is what Viro3D works with. 
    This is metadata from the International Committee on Taxonomy of Viruses. Entries include the virus name, isolate designation, 
    suggested abbreviation, GenBank accession number, genome composition, and host source.
    """
    output:
        xlsx=INPUT_DIRPATH / "viral" / "vmr_metadata.xlsx",
        tsv=INPUT_DIRPATH / "viral" / "vmr_metadata.tsv",
    conda:
        "envs/xlsx2csv.yml"
    shell:
        """
        curl -JLo {output.xlsx} https://ictv.global/sites/default/files/VMR/VMR_MSL38_v2.xlsx
        xlsx2csv {output.xlsx} --delimiter '\t' > {output.tsv}
        """


rule merge_virushostdb_into_vmr:
    """
    Merge all columns from virushostdb.tsv into vmr_metadata.tsv based on matching REFSEQ IDs.
    """
    input:
        virushostdb=INPUT_DIRPATH / "viral" / "virushostdb.tsv",
        vmr_metadata=INPUT_DIRPATH / "viral" / "vmr_metadata.tsv",
    output:
        merged_vmr=INPUT_DIRPATH / "viral" / "vmr_metadata_with_virushostdb.tsv",
    shell:
        """
        python scripts/merge_virushostdb_into_vmr.py {input.vmr_metadata} {input.virushostdb} {output.merged_vmr}
        """


rule extract_unique_virus_names:
    """
    Extract unique virus names from vmr_metadata_with_virushostdb.tsv.
    """
    input:
        merged_vmr=rules.merge_virushostdb_into_vmr.output.merged_vmr,
    output:
        unique_viruses_txt=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "unique_virus_names.txt",
    shell:
        """
        python scripts/extract_unique_viruses.py {input.merged_vmr} {output.unique_viruses_txt}
        """


rule fetch_viro3d_structures_metadata:
    """
    Fetch all viral structure metadata for filtered viruses from the Viro3D API.
    """
    input:
        unique_viruses=rules.extract_unique_virus_names.output.unique_viruses_txt,
    output:
        metadata_dir=directory(OUTPUT_DIRPATH / "metadata" / "viro3d_{host_organism}_metadata"),
    shell:
        """    
        mkdir -p {output.metadata_dir}
        # Read each virus name directly from the input file
        while read virus; do
            # Ensure safe filename (replace spaces & special characters)
            safe_virus=$(echo "$virus" | tr -cs '[:alnum:]_' '_')

            # Define the response file using the sanitized virus name
            response_file={output.metadata_dir}/"${{safe_virus}}.json"
        
            # Encode the virus name for the API call
            encoded_virus=$(python3 -c "import urllib.parse; print(urllib.parse.quote('''$virus'''))")
        
            # Make the API request and save the response
            curl -s -X 'GET' \
                "https://viro3d.cvr.gla.ac.uk/api/proteins/virus_name/?qualifier=$encoded_virus" \
                -H 'accept: application/json' > "$response_file"
        
            # Check if the response file is not empty and log the result
            if [[ -s "$response_file" ]]; then
                echo "Successfully fetched data for $virus" >> {output.metadata_dir}/debug.log
            else
                echo "No data fetched for $virus" >> {output.metadata_dir}/debug.log
            fi
        done < {input.unique_viruses}
        """


rule download_all_pdbs:
    """
    Download all PDB files for viruses listed in the Viro3D metadata.
    """
    input:
        metadata_dir=rules.fetch_viro3d_structures_metadata.output.metadata_dir,
    output:
        pdb_dir=directory(
            OUTPUT_DIRPATH / "viral" / "{host_organism}" / "viro3d_{host_organism}_pdbs"
        ),
        summary=OUTPUT_DIRPATH / "metadata" / "downloadedviro3d_{host_organism}_pdbs.txt",
        fails=OUTPUT_DIRPATH / "metadata" / "downloadedviro3d_{host_organism}_pdbs_fails.txt",
    shell:
        """
        python scripts/download_all_pdbs_viro3d.py {input.metadata_dir} {output.pdb_dir} {output.summary} {output.fails}
        """


rule download_fails:
    """
    Checks the failed downloads list and viro3d_all_pdbs folder for files with 22-byte sizes and redownloads them with the opposite prefix (CF or EF).
    Updates the summary file to add the new downloads.
    """
    input:
        pdb_dir=rules.download_all_pdbs.output.pdb_dir,
        summary=rules.download_all_pdbs.output.summary,
        fails=rules.download_all_pdbs.output.fails,
        metadata_dir=OUTPUT_DIRPATH / "metadata" / "viro3d_{host_organism}_metadata",
    output:
        new_summary=OUTPUT_DIRPATH
        / "metadata"
        / "downloadedviro3d_{host_organism}_pdbs_withfails.txt",
    shell:
        """
        python scripts/download_fails.py {input.pdb_dir} {input.summary} {input.metadata_dir} {input.fails} {output.new_summary}
        """


rule add_structure_file_column:
    """
    Adds a structure_file column to the summary file based on Virus Name, Chosen Method, and Record ID.
    """
    input:
        new_summary=rules.download_fails.output.new_summary,
    output:
        updated_summary=OUTPUT_DIRPATH
        / "metadata"
        / "downloadedviro3d_{host_organism}_pdbs_updated.txt",
    shell:
        """
        python scripts/add_structure_file.py {input.new_summary} {output.updated_summary}
        """


rule merge_metadata:
    """
    Merges summary file with metadata from JSON files.
    """
    input:
        updated_summary=rules.add_structure_file_column.output.updated_summary,
        json_dir=rules.fetch_viro3d_structures_metadata.output.metadata_dir,
    output:
        merged_metadata=OUTPUT_DIRPATH / "metadata" / "merged_viral_metadata_{host_organism}.tsv",
    shell:
        """
        python scripts/merge_metadata.py \
            --summary-file {input.updated_summary} \
            --json-dir {input.json_dir} \
            --output-file {output.merged_metadata}
        """


rule check_downloads:
    """
    Makes sure we have the same number of pdbs as record ids in the metadata.
    """
    input:
        pdb_dir=rules.download_all_pdbs.output.pdb_dir,
        merged_metadata=rules.merge_metadata.output.merged_metadata,
    output:
        done="logs/{host_organism}_compare_counts_done.txt",
    shell:
        """
        python scripts/check_download_number.py {input.merged_metadata} {input.pdb_dir}
        touch {output.done}
        """


rule all:
    default_target: True
    input:
        expand(rules.check_downloads.output.done, host_organism=HOST_ORGANISMS),
