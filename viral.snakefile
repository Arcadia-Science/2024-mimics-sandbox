"""
This snakefile compares viral protein structures against proteomes of the virus hosts to identify
viral structural mimicry. The approach taken in this Snakefile takes advantage of two resources:
1. The Viro3D database of predicted virus protein structures
   ([preprint](https://doi.org/10.1101/2024.12.19.629443),
   [website](https://viro3d.cvr.gla.ac.uk/)) 
2. The [Kyoto Encyclopedia of Genes and Genomes Virus-Host DB](https://www.genome.jp/virushostdb/).  

The KEGG Virus-Host DB defines virus-host pairs for which we have evidence that the virus infects a
specific host. We use these relationships and the
[NCBI taxonomic identifiers](https://www.ncbi.nlm.nih.gov/taxonomy) in this resource to identify
viral protein structures folded by Nomburg et al. We then compare these structures against the
defined host proteomes. We chose to search for structural mimicry in viruses that infect humans,
mouse, rat, and three non-human primates. We think that if we identify the instances of mimicry of
the same gene or pathway in multiple host-viruses pairs, this may be a strong signal that the
protein can be used to change host physiology.

The approach in this snakefile is different than the strategy used in the eukaryote parasite mimicry
detection snakefile. As of Oct 2024, viral UniProt sequences are not a part of alphafold so can't be
handled in the same way as the eukaryotic parasites. If this changes, we could take the same
strategy as is used in the other snakefile. 
"""

from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs")
INPUT_DIRPATH = Path("inputs")
# Indicate additional metadata fields to retrieve when fetching UniProt protein metadata.
# Options documented at https://www.uniprot.org/help/return_fields.
UNIPROT_ADDITIONAL_FIELDS = "cc_alternative_products,ft_var_seq,ft_variant,cc_caution,ft_act_site,cc_activity_regulation,ft_binding,cc_catalytic_activity,cc_cofactor,ft_dna_bind,ec,cc_function,kinetics,cc_pathway,ph_dependence,redox_potential,rhea,ft_site,temp_dependence,fragment,organelle,mass,cc_rna_editing,reviewed,cc_interaction,cc_subunit,cc_developmental_stage,cc_induction,cc_tissue_specificity,go_id,cc_allergen,cc_biotechnology,cc_disruption_phenotype,cc_disease,ft_mutagen,cc_pharmaceutical,cc_toxic_dose,ft_intramem,cc_subcellular_location,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,cc_ptm,ft_propep,ft_signal,ft_transit,ft_coiled,ft_compbias,cc_domain,ft_domain,ft_motif,protein_families,ft_region,ft_repeat,ft_zn_fing,lit_pubmed_id"

# Read in host metadata, which we'll use to link the organism name to identifiers like taxid and
# uniprot proteome id. Setting the index allows us to use the organism name as a look up to find
# the correct value for other metadata.
host_metadata = pd.read_csv("inputs/viral/host-information.csv", header=0).set_index(
    ["organism"], drop=False
)
HOST_ORGANISMS = host_metadata["organism"].unique().tolist()
# limit to human for now
HOST_ORGANISMS = HOST_ORGANISMS[0]

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
        curl -JLo {output.file_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/file_utils.py
        curl -JLo {output.map_refseq_ids} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/map_refseq_ids.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py
        """


##############################################################################
## Get host-specific viral structures
##############################################################################

rule download_viro3d_virus_structures_that_infect_a_host:
    """
    We generated compressed files containing the structures that infect a specific host in the
    repository [Arcadia-Science/2024-mimic-benchmarking](https://github.com/Arcadia-Science/2024-mimic-benchmarking/blob/8eaeb9117b40910d0402cd9ace7a7436630d0514/getallviro.snakefile).  
    """
    output:
        zipf=INPUT_DIRPATH / "viral" / "viro3d_{host_organism}_pdbs.zip",
    shell:
        """
        curl -JLo {output} # placeholder rule for when these are available for download on zenodo 
        # we'll probably have to add the zenodo link to the host metadata CSV file and grab them
        # to make sure the correct set of structures is downloaded for each organism.
        """


rule decompress_viral_structures:
    input:
        zipf=rules.download_viro3d_virus_structures_that_infect_a_host.output.zipf,
    output:
        # TER TODO dest dir may change, test this out and adjust
        dest_dir=directory(OUTPUT_DIRPATH / "viral" / "{host_organism}" / "viral_structures"),
    shell:
        """
        unzip {input.zipf} -d {output.dest_dir}
        """

rule download_viro3d_virus_structure_metadata:
    output: tsv=INPUT_DIRPATH / "viral" / "{host_organism}" / "merged_viral_metadata.tsv"
    # TER TODO: add URL once it's public; metadata is now separated into different files
    # (emily's pipeline creates it by host), so this will probably be a url in the host org file
    # that we'll need to pull. It might make sense to change the name to have the host organism
    # in it instead of making that a directory. TBD based on if we expand this to other species.
    shell:
        """
        curl -JLo {output} #URL
        """


#####################################################################
## Download host proteome structures & metadata
#####################################################################


rule download_uniprot_proteome_canonical_sequence_ids:
    """
    The following rule uses the UniProt FASTA file that records one protein per gene for each of our
    reference proteomes as the source of UniProt Protein IDs. This allows us to use only the
    canonical host protein for comparisons, as these are more metadata-complete and facilitate a
    cleaner analysis.
    """
    output:
        txt=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "host_proteome_canonical_protein_ids.txt",
    conda:
        "envs/seqkit.yml"
    params:
        uniprot_proteome_id=lambda wildcards: host_metadata.loc[
            wildcards.host_organism, "uniprot_proteome_id"
        ],
        taxon_id=lambda wildcards: host_metadata.loc[wildcards.host_organism, "taxon_id"],
    shell:
        """
        curl -JL https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{params.uniprot_proteome_id}/{params.uniprot_proteome_id}_{params.taxon_id}.fasta.gz | \
            seqkit seq --only-id --name | cut -d'|' -f2 > {output.txt}
        """


rule fetch_uniprot_metadata_per_host_proteome:
    """
    Query UniProt using canonical protein IDs and download specified metadata as a TSV.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.download_uniprot_proteome_canonical_sequence_ids.output.txt,
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


rule download_host_pdbs:
    """
    Download all PDB files from AlphaFold.
    While this outputs many PDBs, we don't have any operations in Snakemake where the snakefile
    needs to be aware of all of the PDB accessions. Therefore, instead of treating this as a
    checkpoint and complicating the DAG, we'll only designate the directory as output.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.download_uniprot_proteome_canonical_sequence_ids.output.txt,
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


rule assess_pdbs_per_host_proteome:
    """
    Calculates the quality of all PDBs.
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
        protein_structures_dir=rules.download_host_pdbs.output.protein_structures_dir,
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


#####################################################################
## Download and format human-specific metadata
#####################################################################
"""
We're most interested in human phenotypes, and there is a lot of metadata available about human
genes/proteins. This section downloads and formats some human-specific resources. 
"""


rule download_human_expression_atlas_metadata:
    output:
        tsv=INPUT_DIRPATH / "human_proteinatlas.tsv.zip",
    shell:
        """
        curl -JLo {output} https://www.proteinatlas.org/download/proteinatlas.tsv.zip
        """


rule download_and_process_mouse_opentargets_phenotype_data:
    output:
        csv=OUTPUT_DIRPATH / "metadata" / "opentargets" / "mouse_phenotypes.csv",
    params:
        download_dir=INPUT_DIRPATH / "opentargets" / "mouse_phenotypes_raw",
    conda:
        "envs/pandas.yml"
    shell:
        """
        python scripts/download_and_process_mouse_opentargets_phenotype_data.py \
            --download-dir {params.download_dir} \
            --output {output.csv}
        """


rule combine_human_metadata:
    input:
        tsv1=expand(
            rules.fetch_uniprot_metadata_per_host_proteome.output.tsv, host_organism="human"
        ),
        tsv2=rules.download_human_expression_atlas_metadata.output.tsv,
        csv=rules.download_and_process_mouse_opentargets_phenotype_data.output.csv,
    output:
        csv1=OUTPUT_DIRPATH / "metadata" / "human_metadata_combined.csv",
        csv2=OUTPUT_DIRPATH / "metadata" / "human_metadata_combined_filtered.csv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_human_metadata.R \
            --input_uniprot {input.tsv1} \
            --input_protein_atlas {input.tsv2} \
            --input_mouse_opentargets {input.csv} \
            --output {output.csv1} \
            --output_filtered {output.csv2}
        """


#####################################################################
## Do structural comparisons between viral proteins and host proteins
#####################################################################


rule compare_each_viral_pdb_against_all_host_pdbs:
    """
    TER TODO: output like the foldseek server if useful to the team:
    foldseek easy-search example/d1asha_ example/ result.html tmp --format-mode 3

    Below, we use foldseek to compare viral structures against 
    For evalue filtering (-e), see
    [this issue](https://github.com/steineggerlab/foldseek/issues/167#issuecomment-1674254348),
    which recommends filtering with an evalue threshold of 0.01. This retains hits that likely have
    real structural homology while reducing the overall size of the initial results. In subsequent
    rules, we filter these results in different ways to pull out different kinds of mimics or mimics
    with different properties of interest.
    """
    input:
        protein_structures_dir=rules.download_host_pdbs.output.protein_structures_dir,
        pdbs=rules.decompress_viral_structures.output.dest_dir,
    output:
        tsv=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "viro3d_virus_matches.tsv",
    conda:
        "envs/foldseek.yml"
    shell:
        """
        foldseek easy-search \
            {input.pdbs} \
            {input.protein_structures_dir} \
            {output.tsv} \
            tmp_foldseek \
            -e 0.01 \
            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,prob,qcov,tcov,pident,bits,evalue,qstart,qend,tstart,tend \
            --format-mode 4
        """


rule combine_results_with_metadata_viral:
    input:
        foldseek_tsv=rules.compare_each_viral_pdb_against_all_host_pdbs.output.tsv,
        human_metadata_csv=rules.combine_human_metadata.output.csv1,
        host_metadata_tsv=rules.fetch_uniprot_metadata_per_host_proteome.output.tsv,
        host_lddt_tsv=rules.assess_pdbs_per_host_proteome.output.tsv,
        query_metadata_tsv=rules.download_viro3d_virus_structure_metadata.output.tsv,
    output:
        tsv=OUTPUT_DIRPATH / "viral" / "{host_organism}" / "viro3d_virus_matches_with_metadata.tsv",
    conda:
        "envs/tidyverse.yml"
    # TER TODO delete bit about query plddts in script
    # otherwise adjust for viro3d and potentially gtalign
    shell:
        """
        Rscript scripts/combine_results_with_metadata_viral.R \
            --host {wildcards.host_organism} \
            --input_results {input.foldseek_tsv} \
            --input_human_metadata {input.human_metadata_csv} \
            --input_host_metadata {input.host_metadata_tsv} \
            --input_host_lddt {input.host_lddt_tsv} \
            --input_query_metadata {input.query_metadata_tsv} 
            --output {output.tsv}
        """


#rule combine_all_foldseek_results:
#    input:
#        csvs=expand(
#            # TER TODO: figure out what needs to be combined and put the input files here
#        ),
#    output:
#        csv=OUTPUT_DIRPATH
#        / "viral"
#        # TER TODO: make a better output file name
#    conda:
#        "envs/csvtk.yml"
#    shell:
#        """
#        csvtk concat --out-file {output.csv}  {input.csvs}
#        """


rule all:
    default_target: True
    input:
        #rules.combine_all_foldseek_results.output.csv,
        expand(rules.combine_results_with_metadata_viral.output.tsv, host_organism = HOST_ORGANISMS),
        rules.combine_human_metadata.output.csv2,
