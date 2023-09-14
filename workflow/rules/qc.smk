rule fastqc_R1:
    input:
        get_fq1,
    output:
        html="results/qc/fastqc/{sample}-{unit}_R1_fastqc.html",
        zip="results/qc/fastqc/{sample}-{unit}_R1_fastqc.zip",
    log:
        "logs/fastqc/{sample}-{unit}_R1.log"
    wrapper:
        "v2.6.0/bio/fastqc"
rule fastqc_R2:
    input:
        get_fq2,
    output:
        html="results/qc/fastqc/{sample}-{unit}_R2_fastqc.html",
        zip="results/qc/fastqc/{sample}-{unit}_R2_fastqc.zip",
    log:
        "logs/fastqc/{sample}-{unit}_R2.log"
    wrapper:
        "v2.6.0/bio/fastqc"
rule multiqc:
    input:
        expand(
            "results/qc/fastqc/{unit.sample_name}-{unit.unit_name}_R{reads}_fastqc.zip",
            unit=units.itertuples(),
            reads=[1,2]
        ),
    output:
        "results/qc/multiqc_report.html",
    params:
        extra="--config "+workflow.source_path("../resources/multiqc_config.yml"),
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    wrapper:
        "v2.2.1/bio/multiqc"
