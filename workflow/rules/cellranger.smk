rule cellranger:
    input:
        r1 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R1_001.fastq.gz",
        r2 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R2_001.fastq.gz",
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
        cellranger count --id={sample}_{unit} \
            --transcriptome={config[ref][reference_transcriptome]} \
            --fastqs=results/remap_barcode \
            --sample={sample}_{unit}_RNA \
            {params.variousParams}
        """
