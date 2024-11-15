from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs")
INPUT_DIRPATH = Path("inputs")
# Indicate additional metadata fields to retrieve when fetching UniProt protein metadata.
# Options documented at https://www.uniprot.org/help/return_fields.
UNIPROT_ADDITIONAL_FIELDS = "ft_act_site,cc_activity_regulation,ft_binding,cc_catalytic_activity,cc_cofactor,ft_dna_bind,ec,cc_function,kinetics,cc_pathway,ph_dependence,redox_potential,rhea,ft_site,temp_dependence,fragment,organelle,mass,cc_rna_editing,reviewed,cc_interaction,cc_subunit,cc_developmental_stage,cc_induction,cc_tissue_specificity,go_id,cc_allergen,cc_biotechnology,cc_disruption_phenotype,cc_disease,ft_mutagen,cc_pharmaceutical,cc_toxic_dose,ft_intramem,cc_subcellular_location,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit,ft_coiled,ft_compbias,cc_domain,ft_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing,lit_pubmed_id"

proteome_metadata = pd.read_csv(
    INPUT_DIRPATH / "2024_mimics_uniprot_reference_proteomes_human_long_association.tsv",
    header=0,
    sep="\t",
)
PROTEOME_IDS = proteome_metadata["proteome_id"].unique().tolist()

###########################################################
## Download ProteinCartography scripts
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

    An alternative to this approach would be to use
    [Git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules). We may update to this
    approach in the future.
    """
    output:
        # Create an empty file to use a pointer for this rule running successfully.
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
        map_refseq_ids="ProteinCartography/map_refseq_ids.py",
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
        curl -JLo {output.file_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/file_utils.py
        curl -JLo {output.map_refseq_ids} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/map_refseq_ids.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py
        """


#####################################################
## Get human proteome information
#####################################################
"""
The human information downloaded in this section will be used to detect structural mimics.

Some of these rules are duplicated from other places in the Snakefile.
We do this so that the human information is treated separately, which makes it simpler to do the
each-parasite-protein-versus-all-human-proteins comparison.
"""


rule download_uniprot_proteome_human:
    output:
        faa=INPUT_DIRPATH / "uniprot" / "human" / "UP000005640_9606.fasta.gz",
    shell:
        """
        curl -JLo {output.faa} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
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


rule download_human_proteome_structures_from_alphafold:
    output:
        tar=INPUT_DIRPATH / "uniprot" / "human" / "UP000005640_9606_HUMAN_v4.tar",
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


###########################################################
## Eukaryotic structural comparisons
###########################################################


rule get_uniprot_ids_from_proteome_id:
    output:
        txt=OUTPUT_DIRPATH / "uniprot" / "proteomes" / "{proteome_id}_protein_identifiers.txt",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python scripts/get_uniprot_ids_from_proteome.py --uniprot-proteome-id {wildcards.proteome_id} --output {output.txt}
        """


rule fetch_uniprot_metadata_per_proteome:
    """
    Query Uniprot for the aggregated hits and download all metadata as a big ol' TSV.
    """
    input:
        txt1=rules.download_proteincartography_scripts.output.txt,
        txt2=rules.get_uniprot_ids_from_proteome_id.output.txt,
    output:
        tsv=OUTPUT_DIRPATH / "uniprot" / "proteomes" / "{proteome_id}_features.tsv",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py \
            --input {input.txt2} \
            --output {output.tsv} \
            --additional-fields {UNIPROT_ADDITIONAL_FIELDS} 
        """


rule filter_to_proteins_that_contain_a_signal_peptide:
    """
    Parses the UniProt metadata to keep only proteins that are annotated as having a signal peptide.
    """
    input:
        tsv=rules.fetch_uniprot_metadata_per_proteome.output.tsv,
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


rule download_pdbs:
    """
    Download all PDB files from AlphaFold.
    While this outputs many PDBs, we don't have any operations in Snakemake where the snakefile
    needs to be aware of all of the PDB accessions. Therefore, instead of treating this as a
    checkpoint and complicating the DAG, we'll only designate the directory as output.
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


rule assess_pdbs_per_proteome:
    """
    Calculates the quality of all PDBs.
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
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


rule compare_each_parasite_pdb_against_human_pdb:
    """
    TER TODO: output like the foldseek server: foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3
    For evalue filtering (-e), see
    [this issue](https://github.com/steineggerlab/foldseek/issues/167#issuecomment-1674254348),
    which recommends filtering with an evalue threshold of 0.01.
    """
    input:
        human_protein_structures_dir=rules.decompress_human_proteome_structures_from_alphafold.output.human_protein_structures_dir,
        protein_structures_dir=rules.download_pdbs.output.protein_structures_dir,
    output:
        tsv=OUTPUT_DIRPATH / "structural_comparison" / "foldseek_raw" / "{proteome_id}.tsv",
    conda:
        "envs/foldseek.yml"
    shell:
        """
        foldseek easy-search \
            {input.protein_structures_dir}/*pdb \
            {input.human_protein_structures_dir} \
            {output.tsv} \
            tmp_foldseek \
            -e 0.01 \
            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,lddtfull,prob,qcov,tcov,pident,bits,evalue,cigar,qseq,tseq,qstart,qend,tstart,tend,qaln,taln,qca,tca,u,t \
            --format-mode 4
        """


rule combine_foldseek_parasite_results_with_uniprot_metadata:
    input:
        foldseek_tsv=rules.compare_each_parasite_pdb_against_human_pdb.output.tsv,
        human_tsv=rules.fetch_uniprot_metadata_human.output.tsv,
        query_tsv=rules.fetch_uniprot_metadata_per_proteome.output.tsv,
    output:
        tsv=OUTPUT_DIRPATH
        / "structural_comparison"
        / "foldseek_with_uniprot_metadata"
        / "{proteome_id}_with_uniprot_metadata.tsv.gz",
        tsv_filtered=OUTPUT_DIRPATH
        / "structural_comparison"
        / "foldseek_with_uniprot_metadata"
        / "{proteome_id}_with_uniprot_metadata_filtered.tsv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_foldseek_results_with_uniprot_metadata.R \
            --input_foldseek_results {input.foldseek_tsv} \
            --input_human_metadata {input.human_tsv} \
            --input_query_metadata {input.query_tsv} \
            --output {output.tsv} \
            --output_filtered {output.tsv_filtered}
        """


###################################################################
## Rule all
###################################################################


rule all:
    default_target: True
    input:
        expand(rules.assess_pdbs_per_proteome.output.tsv, proteome_id=PROTEOME_IDS),
        expand(
            rules.combine_foldseek_parasite_results_with_uniprot_metadata.output.tsv,
            proteome_id=PROTEOME_IDS,
        ),
