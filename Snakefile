from pathlib import Path

OUTPUT_DIRPATH=Path("outputs/")
UNIPROT_ADDITIONAL_FIELDS="ft_act_site,cc_activity_regulation,ft_binding,cc_catalytic_activity,cc_cofactor,ft_dna_bind,ec,cc_function,kinetics,cc_pathway,ph_dependence,redox_potential,rhea,ft_site,temp_dependence,fragment,organelle,mass,cc_rna_editing,reviewed,cc_interaction,cc_subunit,cc_developmental_stage,cc_induction,cc_tissue_specificity,go_p,go_c,go,go_f,go_id,cc_allergen,cc_biotechnology,cc_disruption_phenotype,cc_disease,ft_mutagen,cc_pharmaceutical,cc_toxic_dose,ft_intramem,cc_subcellular_location,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit,ft_coiled,ft_compbias,cc_domain,ft_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing,lit_pubmed_id"

rule all:
    input:
        "outputs/plmutils/embedded_proteins.csv",
        "outputs/plmutils/embedded_protein_names.txt",

###########################################################
## Download and preprocess proteins
###########################################################

rule download_human_uniprot_proteome:
    output:
        faa="inputs/uniprot/human/UP000005640_9606.fasta.gz",
    shell:
        """
        curl -JLo {output.faa} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
        # curl -JLo {output.faa} https://rest.uniprot.org/uniprotkb/stream\?compressed\=true\&format\=fasta\&query\=%28%28proteome%3AUP000005640%29%29
        """

rule download_proteincartography_scripts:
    """
    ProteinCartography (https://github.com/Arcadia-Science/ProteinCartography) contains many scripts
    for interacting with UniProt and AlphaFold. Because PC is not pip-installable, we download the
    scripts and environments we need in this workflow. We take this approach instead of making a
    copy of each script inside this repo. Technically, this should be broken out into many rules
    (one per file), but that seemed unnecessarily verbose so I confined it to one rule. 
    """
    output:
        # Create an empty file to use a pointer for this rule runnign successfully.
        # This will allow us to refer to the ProteinCartography scripts by name/filepath instead of
        # by snakemake output syntax.
        txt=touch("scripts/ProteinCartography_scripts_downloaded.txt"),
        api_utils="ProteinCartography/api_utils.py",
        assess_pdbs="ProteinCartography/assess_pdbs.py",
        color_utils="ProteinCartography/color_utils.py",
        constants="ProteinCartography/constants.py",
        download_pdbs="ProteinCartography/download_pdbs.py",
        fetch_accession="ProteinCartography/fetch_accession.py",
        fetch_uniprot_metadata="ProteinCartography/fetch_uniprot_metadata.py",
        mocks="ProteinCartography/tests/mocks.py",
        envs_web_apis="envs/web_apis.yml",
        envs_plotting="envs/plotting.yml"
    params: commit="88160fcf098347a29124488f445ed1d9ad72bc12"
    shell:
        """
        curl -JLo {output.api_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/api_utils.py
        curl -JLo {output.assess_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/assess_pdbs.py
        curl -JLo {output.color_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/color_utils.py
        curl -JLo {output.constants} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/constants.py
        curl -JLo {output.download_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/download_pdbs.py
        curl -JLo {output.fetch_accession} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_accession.py
        curl -JLo {output.fetch_uniprot_metadata} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_uniprot_metadata.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py

        curl -JLo {output.envs_web_apis} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/envs/web_apis.yml
        curl -JLo {output.envs_plotting} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/envs/plotting.yml
        """

rule get_uniprot_ids_from_proteome:
    output: txt = OUTPUT_DIRPATH / "uniprot" / "proteomes" / "{proteome_id}.txt"
    conda: "envs/"
    shell:
        """
        scripts/get_uniprot_ids_from_proteome.py --uniprot-proteome-id {wildcards.proteome_id} --output {output.txt}
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



#####################################################
## Download and compare structures with foldseek
#####################################################

# TER TODO download all alphafolded human PDBs and put somehwere

checkpoint download_pdbs:
    """
    Download all PDB files from AlphaFold
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
        # TER TODO update input file
        rules.filter_aggregated_hits.output.filtered_aggregated_hits,
    output:
        protein_structures_dir=directory(OUTPUT_DIRPATH / "pdb_structures" / "{proteome_id}"),
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py \
            --input {input} \
            --output {output.protein_structures_dir} \
            --max-structures 100000
        """

# TER TODO update function so that it has a wildcard output to foldseek each PDB against all human PDBs
def get_pdb_filepaths(wildcards):
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
    Calculates the quality of all PDBs
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
        get_pdb_filepaths,
    # TER TODO: update output filepath
    # I think this will need the proteome_id wildcard bc that should be a part of the get_pdb_filepaths output
    output:
        pdb_features=PROTEIN_FEATURES_DIR / "pdb_features.tsv",
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py \
            --input {ANALYZED_PROTEIN_STRUCTURES_DIR} \
            --output {output.pdb_features}
        """


rule compare_each_parasite_pdb_against_human_pdb:
    input:
    output:
    conda: "envs/foldseek.yml"
    shell:
        """
        foldseek easy-search <i:PDB|mmCIF[.gz]> ... <i:PDB|mmCIF[.gz]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:alignmentFile> <tmpDir>
        """

#####################################################
## Embed proteins in ESM2
#####################################################

rule extract_protein_sequences_from_tsv:
    """
    The rule fetch_uniprot_metadata returns a TSV that includes the protein accession and sequence.
    This rule extracts the sequence name and uniprot ID from the TSV file into a FASTA file.
    """

rule combine_all_proteins:
    input:
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
