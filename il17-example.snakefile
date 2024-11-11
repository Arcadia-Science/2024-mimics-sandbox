rule all:
    input:
        "outputs/plmutils/embedded_proteins.csv",
        "outputs/plmutils/embedded_protein_names.txt",


rule download_human_uniprot_proteome:
    output:
        faa="inputs/uniprot/human/UP000005640_9606.fasta.gz",
    shell:
        """
        curl -JLo {output.faa} https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
        """


rule download_viral_IL17_mimics:
    """
    Note that not all of these viral proteins are IL17 mimics.
    This ends up being a good thing, as we can use them as controls.
    The viral proteins that do mimic human proteins should be more similar to human proteins (IL17)
    than the proteins that don't.
    """
    output:
        faa="inputs/uniprot/viral/IL-17_taxonomy_id_3A10239.fasta.gz",
    shell:
        """
        curl -JLo {output.faa} "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28IL-17%29+AND+%28taxonomy_id%3A10239%29%29"
        """


rule download_chelicerate_IL17_mimics:
    output:
        tsv="inputs/ticks-on-a-tree/OG0001324_targets.tsv",
        faa="inputs/ticks-on-a-tree/2024-06-24-all-chelicerate-noveltree-proteins.fasta",
    shell:
        """
        # place holder: s3://arcadia-pi-disco/datasheets/OG0001324_targets.tsv
        # place holder: s3://organism-selection-engines/noveltree/chelicerata-v1/orthogroup-analysis-v1/2024-06-24-all-chelicerate-noveltree-proteins.fasta
        """


rule extract_chelicerate_IL17_mimics:
    input:
        tsv=rules.download_chelicerate_IL17_mimics.output.tsv,
        faa=rules.download_chelicerate_IL17_mimics.output.faa,
    output:
        txt="outputs/ticks-on-a-tree/OG0001324_targets.txt",
        faa="outputs/ticks-on-a-tree/OG0001324_targets.faa.gz",
    conda:
        "envs/seqkit.yml"
    shell:
        """
        csvtk cut --tabs --fields gene_name {input.tsv} | csvtk del-header --out-tabs --out-file {output.txt}
        seqkit grep --pattern-file {output.txt} {input.faa} -o {output.faa} 
        """


rule combine_all_proteins:
    input:
        rules.download_human_uniprot_proteome.output.faa,
        rules.download_viral_IL17_mimics.output.faa,
        rules.extract_chelicerate_IL17_mimics.output.faa,
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
    Note we could probably help more proteins be included by cleaving signal peptides.
    Since IL17 isn't long enough for this to be a problem, we won't do that here.
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


#####################################################
## Embed proteins in ESM2
#####################################################


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
