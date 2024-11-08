from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs")
INPUT_DIRPATH = Path("inputs")
# Indicate additional metadata fields to retrieve when fetching UniProt protein metadata.
# Options documented at https://www.uniprot.org/help/return_fields.
UNIPROT_ADDITIONAL_FIELDS = "ft_act_site,cc_activity_regulation,ft_binding,cc_catalytic_activity,cc_cofactor,ft_dna_bind,ec,cc_function,kinetics,cc_pathway,ph_dependence,redox_potential,rhea,ft_site,temp_dependence,fragment,organelle,mass,cc_rna_editing,reviewed,cc_interaction,cc_subunit,cc_developmental_stage,cc_induction,cc_tissue_specificity,go_id,cc_allergen,cc_biotechnology,cc_disruption_phenotype,cc_disease,ft_mutagen,cc_pharmaceutical,cc_toxic_dose,ft_intramem,cc_subcellular_location,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit,ft_coiled,ft_compbias,cc_domain,ft_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing,lit_pubmed_id"

# Read in host metadata, which we'll use to link the organism name to identifiers like taxid and
# uniprot proteome id. Setting the index allows us to use the organism name as a look up to find 
# the correct value for other metadata.
host_metadata=pd.read_csv("inputs/viral/host-information.csv", header=0).set_index(["organism"], drop=False)
HOST_ORGANISMS=host_metadata["organism"].unique().tolist()

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
        curl -JLo {output.file_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/files_utils.py
        curl -JLo {output.map_refseq_ids} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/map_refseq_ids.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py
        """

##############################################################################
## Get host-specific viral structures
##############################################################################

rule download_kegg_virushostdb:
    """
    Downloaded 2024-10-07 version
    MD5sum 93a78881a8c2b9d01dcdddbb564276fd
    """
    output: tsv=INPUT_DIRPATH / "viral" / "virushostdb.tsv"
    shell:
        """
        curl -JLo {output.tsv} https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
        """

rule download_ncbi_taxdump_information:
    """
    Tax dump contains files to go from a NCBI taxonomy ID to full lineage or to name.
    """
    output:
        tar=INPUT_DIRPATH / "taxdump" / "taxdump.tar.gz",
        dmp=INPUT_DIRPATH / "taxdump" / "nodes.dmp"
    params: taxdump_dirpath=INPUT_DIRPATH / "taxdump"
    shell:
        """
        curl -JLo {output.tar} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && \
            tar xf {output.tar} -C {params.taxdump_dirpath}
        """

rule retrieve_ncbi_taxonomy_lineages_for_viral_taxids:
    """
    The NCBI taxonomy IDs don't match between KEGG virushostdb and Nomburg et al.
    Often, virushostdb uses a strain-level taxid while Nomburg uses a species.
    """
    input:
        tsv=rules.download_kegg_virushostdb.output.tsv,
        dmp=rules.download_ncbi_taxdump_information.output.dmp
    output: tsv=OUTPUT_DIRPATH / "viruses" / "viral_lineages.tsv"
    conda: "envs/taxonkit.yml"
    params: taxdump_dirpath=INPUT_DIRPATH / "taxdump"
    shell:
        """
        taxonkit lineage \
            --data-dir {params.taxdump_dirpath} \
            --show-lineage-taxids \
            --taxid-field 1 \
            --out-file {output.tsv} \
            {input.tsv}
        """

rule download_nomburg_supplementary_table:
    output: xlsx=INPUT_DIRPATH / "41586_2024_7809_MOESM4_ESM.xlsx"
    conda: "envs/web_apis.yml"
    shell:
        """
        curl -JLo {output.xlsx} https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07809-y/MediaObjects/41586_2024_7809_MOESM4_ESM.xlsx
        """

rule filter_nomburg_viruses_by_host:
    input:
        csv=INPUT_DIRPATH / "viral" / "host-information.csv",
        xlsx=rules.download_nomburg_supplementary_table.output.xlsx,
        tsv=rules.retrieve_ncbi_taxonomy_lineages_for_viral_taxids.output.tsv
    output:
        tsv=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "{host_organism}_structure_metadata.tsv",
        txt=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "{host_organism}_structure_paths.txt",
    conda: "envs/tidyverse.yml"
    shell:
        """
        Rscript filter_nomburg_viruses_by_host.R \
            --host_organism {wildcards.host_organism} \
            --host_metadata {input.csv} \
            --nomburg_metadata {input.xlsx} \
            --virushostdb {input.tsv} \
            --output_tsv {output.tsv} \
            --output_txt {output.txt}
        """


rule download_nomburg_eukaryotic_virus_structures:
    output: zipf=INPUT_DIRPATH / "viral" / "Nomburg_2023_structures.zip"
    shell:
        """
        curl -JLo {output} https://zenodo.org/records/10291581/files/Nomburg_2023_structures.zip?download=1
        """


rule decompress_viral_structures:
    """
    Only decompress structures that we want to compare against a given host.
    """
    input:
        zipf=rules.download_nomburg_eukaryotic_virus_structures.output.zipf,
        txt=rules.filter_nomburg_viruses_by_host.output.txt
    output: dest_dir=directory(INPUT_DIRPATH / "viral" / "{host_organism}" / "viral_structures")
    params: dest_dir = lambda wildcards: INPUT_DIRPATH / "viral" / wildcards.host_organism
    shell:
        """
        cat {input.txt} | tr '\n' '\0' | xargs -0 unzip {input.zipf} -d {params.dest_dir}
        """

#####################################################################
## Download host proteome structures
#####################################################################

rule get_uniprot_ids_from_proteome_id:
    output:
        txt=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "host_proteome_protein_identifiers.txt",
    params: uniprot_proteome_id=lambda wildcards: host_metadata.loc[wildcards.host_organism, "uniprot_proteome_id"],
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python scripts/get_uniprot_ids_from_proteome.py \
            --uniprot-proteome-id {params.uniprot_proteome_id} \
            --output {output.txt}
        """


rule fetch_uniprot_metadata_per_proteome:
    """
    Query Uniprot for the aggregated hits and download all metadata as a TSV.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.get_uniprot_ids_from_proteome_id.output.txt,
    output:
        tsv=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "host_proteome_protein_features.tsv",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/fetch_uniprot_metadata.py \
            --input {input.txt} \
            --output {output.tsv} \
            --additional-fields {UNIPROT_ADDITIONAL_FIELDS} 
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
        txt=rules.get_uniprot_ids_from_proteome_id.output.txt,
    output:
        protein_structures_dir=directory(
            OUTPUT_DIRPATH / "viral" / "{host_organism}" / "host_proteome_pdb_structures"
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
        tsv=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "host_proteome_pdb_structure_quality.tsv",
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py \
            --input  {input.protein_structures_dir} \
            --output {output.tsv}
        """


# TER TODO: figure out which are problematic columns and remove from output
#rule compare_each_viral_pdb_against_human_pdb:
#    """
#    TER TODO: output like the foldseek server: foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3
#    """
#    input:
#        human_protein_structures_dir=rules.decompress_human_proteome_structures_from_alphafold.output.human_protein_structures_dir,
#        pdbs=rules.decompress_viral_structures.output.dest_dir
#    output:
#        tsv=OUTPUT_DIRPATH / "viruses" / "foldseek_raw" / "nomburg_human_viruses.tsv",
#    conda:
#        "envs/foldseek.yml"
#    shell:
#        """
#        foldseek easy-search \
#            {input.pdbs} \
#            {input.human_protein_structures_dir} \
#            {output.tsv} \
#            tmp_foldseek \
#            -e 0.01 \
#            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,lddtfull,prob,qcov,tcov,pident,bits,evalue,cigar,qseq,tseq,qstart,qend,tstart,tend,qaln,taln,qca,tca,u,t \
#            --format-mode 4
#        """

rule all:
    default_target: True
    input:
        expand(rules.decompress_viral_structures.output.dest_dir, host_organism=HOST_ORGANISMS),
        expand(rules.assess_pdbs_per_proteome.output.tsv, host_organism=HOST_ORGANISMS)
