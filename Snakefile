from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs")
INPUT_DIRPATH = Path("inputs")
# Indicate additional metadata fields to retrieve when fetching UniProt protein metadata.
# Options documented at https://www.uniprot.org/help/return_fields.
UNIPROT_ADDITIONAL_FIELDS = "ft_act_site,cc_activity_regulation,ft_binding,cc_catalytic_activity,cc_cofactor,ft_dna_bind,ec,cc_function,kinetics,cc_pathway,ph_dependence,redox_potential,rhea,ft_site,temp_dependence,fragment,organelle,mass,cc_rna_editing,reviewed,cc_interaction,cc_subunit,cc_developmental_stage,cc_induction,cc_tissue_specificity,go_id,cc_allergen,cc_biotechnology,cc_disruption_phenotype,cc_disease,ft_mutagen,cc_pharmaceutical,cc_toxic_dose,ft_intramem,cc_subcellular_location,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit,ft_coiled,ft_compbias,cc_domain,ft_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing,lit_pubmed_id"

# PROTEOME_IDS=['UP000008333', 'UP000002899']
proteome_metadata = pd.read_csv(
    "inputs/2024_mimics_uniprot_reference_proteomes_human_long_association.tsv", header=0, sep="\t"
)
PROTEOME_IDS = proteome_metadata["proteome_id"].unique().tolist()
print(PROTEOME_IDS)

###########################################################
## Download and preprocess proteins
###########################################################


rule download_proteincartography_scripts:
    """
    ProteinCartography (https://github.com/Arcadia-Science/ProteinCartography) contains many scripts
    for interacting with UniProt and AlphaFold. Because PC is not pip-installable, we download the
    scripts and environments we need in this workflow. We take this approach instead of making a
    copy of each script inside this repo. Technically, this should be broken out into many rules
    (one per file), but that seemed unnecessarily verbose so I confined it to one rule. 
    
    Note the envs that these scripts require (envs/plotting.yml, envs/web_apis.yml) need to already
    be present in the repo so they are duplicated.
    """
    output:
        # Create an empty file to use a pointer for this rule runnign successfully.
        # This will allow us to refer to the ProteinCartography scripts by name/filepath instead of
        # by snakemake output syntax.
        txt=touch("scripts/ProteinCartography_scripts_downloaded.txt"),
        api_utils="ProteinCartography/api_utils.py",
        artifact_generation_utils="ProteinCartography/tests/artifact_generation_utils.py",
        assess_pdbs="ProteinCartography/assess_pdbs.py",
        color_utils="ProteinCartography/color_utils.py",
        constants="ProteinCartography/constants.py",
        download_pdbs="ProteinCartography/download_pdbs.py",
        fetch_accession="ProteinCartography/fetch_accession.py",
        fetch_uniprot_metadata="ProteinCartography/fetch_uniprot_metadata.py",
        file_utils="ProteinCartography/file_utils.py",
        mocks="ProteinCartography/tests/mocks.py",
    params:
        commit="88160fcf098347a29124488f445ed1d9ad72bc12",
    shell:
        """
        curl -JLo {output.api_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/api_utils.py
        curl -JLo {output.artifact_generation_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/artifact_generation_utils.py
        curl -JLo {output.assess_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/assess_pdbs.py
        curl -JLo {output.color_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/color_utils.py
        curl -JLo {output.constants} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/constants.py
        curl -JLo {output.download_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/download_pdbs.py
        curl -JLo {output.fetch_accession} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_accession.py
        curl -JLo {output.fetch_uniprot_metadata} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_uniprot_metadata.py
        curl -JLo {output.file_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/files_utils.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py
        """


rule get_uniprot_ids_from_proteome:
    output:
        txt=OUTPUT_DIRPATH / "uniprot" / "proteomes" / "{proteome_id}_protein_identifiers.txt",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python scripts/get_uniprot_ids_from_proteome.py --uniprot-proteome-id {wildcards.proteome_id} --output {output.txt}
        """


rule fetch_uniprot_metadata:
    """
    Query Uniprot for the aggregated hits and download all metadata as a big ol' TSV.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.get_uniprot_ids_from_proteome.output.txt,
    output:
        tsv=OUTPUT_DIRPATH / "uniprot" / "proteomes" / "{proteome_id}_features.tsv",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py \
            --input {input.txt} \
            --output {output.tsv} \
            --additional-fields {UNIPROT_ADDITIONAL_FIELDS} 
        """


rule filter_to_proteins_that_contain_a_signal_peptide:
    """
    Parses the UniProt metadata to keep only proteins that are annotated as having a signal peptide.
    """
    input:
        tsv=rules.fetch_uniprot_metadata.output.tsv,
    output:
        tsv=OUTPUT_DIRPATH
        / "uniprot"
        / "proteomes"
        / "{proteome_id}_features_with_signal_peptides.tsv",
        txt=OUTPUT_DIRPATH
        / "uniprot"
        / "proteomes"
        / "{proteome_id}_protein_identifiers_with_signal_peptides.txt",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/filter_to_proteins_that_contain_a_signal_peptide.R \
            --input {input.tsv} \
            --output_tsv {output.tsv} \
            --output_txt {output.txt}
        """


#####################################################
## Retrieve human proteome informatio
#####################################################
"""
Some of these rules are duplicated from other places in the Snakefile.
We do this so that the human information is treated separately, which makes it simpler to do the
each-parasite-protein-versus-all-human-proteins comparison.
"""


rule download_uniprot_proteome_human:
    output:
        faa=INPUT_DIRPATH / "uniprot"/ "human"/ "UP000005640_9606.fasta.gz",
    shell:
        """
        curl -JLo {output.faa} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
        # curl -JLo {output.faa} https://rest.uniprot.org/uniprotkb/stream\?compressed\=true\&format\=fasta\&query\=%28%28proteome%3AUP000005640%29%29
        """


rule get_uniprot_ids_from_proteome_human:
    output:
        txt=INPUT_DIRPATH / "uniprot" / "human" / "UP000005640_protein_identifiers.txt",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python scripts/get_uniprot_ids_from_proteome.py --uniprot-proteome-id UP000005640 --output {output.txt}
        """


rule fetch_uniprot_metadata_human:
    """
    Query UniProt for human proteins and download all metadata as a TSV.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.get_uniprot_ids_from_proteome_human.output.txt,
    output:
        tsv=INPUT_DIRPATH / "uniprot" / "human" / "UP000005640_features.tsv",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py \
            --input {input.txt} \
            --output {output.tsv} \
            --additional-fields {UNIPROT_ADDITIONAL_FIELDS}
        """


#####################################################
## Download and compare structures with foldseek
#####################################################


checkpoint download_pdbs:
    """
    Download all PDB files from AlphaFold
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.filter_to_proteins_that_contain_a_signal_peptide.output.txt,
    output:
        protein_structures_dir=directory(
            OUTPUT_DIRPATH / "structures" / "alphafold_pdb_structures" / "{proteome_id}"
        ),
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py \
            --input {input.txt} \
            --output {output.protein_structures_dir} \
            --max-structures 100000
        """


def get_all_pdb_filepaths_per_proteome(wildcards):
    """
    Returns a list of all of the PDB files to use for the clustering analysis.

    This function references the `download_pdbs` checkpoint, triggering it to run,
    and then returns a list of all of the resulting downloaded PDB files
    """
    # note: referencing the `download_pdbs` checkpoint here is essential,
    # because this is what 'tells' snakemake to run the checkpoint
    pdb_dirpath = checkpoints.download_pdbs.get(**wildcards).output.protein_structures_dir
    pdb_filepaths = sorted(Path(pdb_dirpath).glob("*.pdb"))

    return pdb_filepaths


rule assess_pdbs:
    """
    Calculates the quality of all PDBs.
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
        get_all_pdb_filepaths_per_proteome,
        protein_structures_dir=rules.download_pdbs.output.protein_structures_dir,
    output:
        tsv=OUTPUT_DIRPATH / "structures" / "alphafold_pdb_quality" / "{proteome_id}.tsv",
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py \
            --input  {input.protein_structures_dir} \
            --output {output.tsv}
        """


rule download_human_proteome_structures_from_alphafold:
    output:
        tar=INPUT_DIRPATH / "uniprot" / "human"/ "UP000005640_9606_HUMAN_v4.tar",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        curl -JLo {output.tar} https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
        """


rule decompress_human_proteome_structures_from_alphafold:
    input:
        tar=rules.download_human_proteome_structures_from_alphafold.output.tar,
    output:
        human_protein_structures_dir=directory(INPUT_DIRPATH / "uniprot" / "human" / "structures"),
    shell:
        """
        tar xf {input.tar} --include "*.pdb.gz" -C {output}
        """


rule compare_each_parasite_pdb_against_human_pdb:
    """
    TER TODO: output like the foldseek server: foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3
    """
    input:
        human_protein_structures_dir=rules.decompress_human_proteome_structures_from_alphafold.output.human_protein_structures_dir,
        pdbs=get_all_pdb_filepaths_per_proteome,
    output:
        tsv=OUTPUT_DIRPATH / "structural_comparison" / "{proteome_id}.tsv",
    conda:
        "envs/foldseek.yml"
    shell:
        """
        foldseek easy-search \
            {input.pdbs} \
            {input.human_protein_structures_dir} \
            {output.tsv} \
            tmp_foldseek \
            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,lddtfull,prob,qcov,tcov,pident,bits,cigar,qseq,tseq,qstart,qend,tstart,tend,qaln,taln,qca,tca,u,t \
            --format-mode 4
        """

rule combine_foldseek_results_with_uniprot_metadata:
    input:
        foldseek_tsv=rules.compare_each_parasite_pdb_against_human_pdb.output.tsv,
        human_tsv=rules.fetch_uniprot_metadata_human.output.tsv,
        query_tsv=rules.fetch_uniprot_metadata.output.tsv
    output: 
        tsv=OUTPUT_DIRPATH / "structural_comparison" / "{proteome_id}_with_uniprot_metadata.tsv",
    conda: "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_foldseek_results_with_uniprot_metadata.R \
            --input_foldseek_results {input.foldseek_tsv} \
            --input_human_metadata {input.human_tsv} \
            --input_query_metadata {input.query_tsv} \
            --output {output}
        """
#####################################################
## Embed proteins in ESM2
#####################################################

# rule extract_protein_sequences_from_tsv:
#    """
#    The rule fetch_uniprot_metadata returns a TSV that includes the protein accession and sequence.
#    This rule extracts the sequence name and uniprot ID from the TSV file into a FASTA file.
#    """


rule combine_all_proteins:
    # input:
    output:
        faa="outputs/input_proteins/all_proteins.faa.gz",
    shell:
        """
        cat {input} > {output}
        """


rule filter_proteins_by_length:
    """
    ESM-2 has a maximum protein length of 1024, so we remove any proteins longer than this.
    If we kept these proteins, ESM would auto-truncate them at length 1024 which would likely
    impact our results/interpretation.
    """
    input:
        rules.combine_all_proteins.output.faa,
    output:
        faa="outputs/input_proteins/length_filtered_proteins.faa",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --max-len 1024 -o {output.faa} {input}
        """


rule plmutils_embed:
    """
    This rule embeds amino acid sequences into the embedding space of a protein language model.
    For now, plmutils only supports ESM-2.
    The parameter --layer-ind -1 means to extract the embedding from the last layer of the model.
    """
    input:
        rules.filter_proteins_by_length.output.faa,
    output:
        npy="outputs/plmutils/embedded_proteins.npy",
    conda:
        "envs/plmutils.yml"
    shell:
        """
        plmutils embed --model-name esm2_t48_3B_UR50D \
            --layer-ind -1 \
            --output-filepath {output.npy} \
            {input}
        """


rule convert_embeddings_to_csv:
    input:
        rules.plmutils_embed.output.npy,
    output:
        csv="outputs/plmutils/embedded_proteins.csv",
    conda:
        "envs/pandas.yml"
    shell:
        """
        python scripts/convert_npy_to_csv.py --npy {input} --csv {output.csv}
        """


rule extract_embedding_rownames:
    input:
        rules.filter_proteins_by_length.output.faa,
    output:
        txt="outputs/plmutils/embedded_protein_names.txt",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        seqkit seq --name --only-id {input} --out-file {output.txt}
        """


###################################################################
## Rule all
###################################################################


rule all:
    default_target: True
    input:
        expand(rules.fetch_uniprot_metadata.output.tsv, proteome_id=PROTEOME_IDS),
        expand(rules.assess_pdbs.output.tsv, proteome_id=PROTEOME_IDS),
        expand(rules.combine_foldseek_results_with_uniprot_metadata.output.tsv, proteome_id=PROTEOME_IDS),
        rules.fetch_uniprot_metadata_human.output.tsv,
        #"outputs/plmutils/embedded_proteins.csv",
        #"outputs/plmutils/embedded_protein_names.txt",
