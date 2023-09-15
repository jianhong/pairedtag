rule cellranger:
    input:
        r1 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R1_001.fastq.gz",
        r2 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R2_001.fastq.gz",
        ref= get_cellranger_ref
    output:
        out = "results/cellranger/{sample}_{unit}/outs",
    params:
        variousParams=config["tools"]["cellranger_count"]["variousParams"],
    log:
        "logs/cellranger/{sample}_{unit}.log",
    threads: config["computingResources"]["threads"]["high"],
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    singularity:
        "docker://quay.io/nf-core/cellranger:7.1.0"
    shell:
        """
        cellranger count --id={wildcards.sample}_{wildcards.unit} \
            --transcriptome={input.ref} \
            --fastqs=results/remap_barcode \
            --sample={wildcards.sample}_{wildcards.unit}_RNA \
            {params.variousParams}
        """
