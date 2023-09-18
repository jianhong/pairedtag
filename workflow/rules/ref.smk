rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v2.2.1/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    wrapper:
        "v2.2.1/bio/reference/ensembl-annotation"

rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v2.2.1/bio/samtools/faidx"

rule chrom_size:
    input:
        "resources/genome.fasta.fai",
    output:
        "resources/genome.sizes"
    log:
        "logs/genome-faidx.log",
    cache: True
    shell:
        "cut -f 1,2 genome.fasta.fai > genome.sizes"

rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v2.2.1/bio/bwa/index"

rule cellranger_index:
    input:
        fa  = "resources/genome.fasta",
        gtf = "resources/genome.gtf",
    output:
        cellranger_transcriptome = "resources/cellranger_transcriptome/"+config["ref"]["build"],
    params:
        mkgtfParams=config["tools"]["cellranger_count"]["mkgtfParams"],
        mkrefParams=config["tools"]["cellranger_count"]["mkrefParams"],
    log:
        "logs/cellranger_mkref",
    threads: config["computingResources"]["threads"]["high"],
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    cache: True
    singularity:
        "docker://quay.io/nf-core/cellranger:7.1.0"
    shell:
        """
        cellranger mkgtf \
            {input.gtf} \
            genome.filtered.gtf \
            {params.mkgtfParams}

        cellranger mkref \
            --genome={output.cellranger_transcriptome} \
            --fasta={input.fa} \
            --genes=genome.filtered.gtf \
            {params.mkrefParams}
        """
