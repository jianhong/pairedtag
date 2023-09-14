from multiprocessing import cpu_count

rule remap_barcode:
    input:
        unpack(get_fastq),
    output:
        r1 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R1_001.fastq.gz",
        r2 = "results/remap_barcode/{sample}_{unit}_RNA_S1_L001_R2_001.fastq.gz",
        d1 = "results/remap_barcode/{sample}_{unit}_DNA_S1_L001_R1_001.fastq.gz",
        d2 = "results/remap_barcode/{sample}_{unit}_DNA_S1_L001_R2_001.fastq.gz",
        ro = "results/remap_barcode/ori_{sample}_{unit}_RNA_S1_L001_R1_001.fastq.gz",
        do = "results/remap_barcode/ori_{sample}_{unit}_DNA_S1_L001_R2_001.fastq.gz",
        mp = "results/remap_barcode/{sample}_{unit}_barcodeMap.tsv.gz"
    log:
        "logs/remap_barcode/{sample}_{unit}.log",
    threads: lambda cores: max(2, cpu_count() - 1)
    conda:
        "../envs/remap_barcode.yaml"
    params:
        scr = workflow.source_path("../scripts/trim_R2.py")
    shell:
        """
        # swith R1 and R2 for RNAseq
        python {params.scr} \
            --cores {threads} \
            -i1 {input.r1} -i2 {input.r2} \
            -d1 {output.d1} -d2 {output.d2} \
            -r1 {output.r2} -r2 {output.r1} \
            --barcodeMap {output.mp} \
            --originBarcodePrefix ori_ \
            --rnaBarcode {config[params][rna_bar]} \
            --dnaBarcode {config[params][dna_bar]} \
            --cellranger \
            --barcodeWhitelist {config[params][whitelist]} \
            --pattern "{config[params][pattern]}"
        """
