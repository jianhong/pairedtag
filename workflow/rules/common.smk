import glob

import pandas as pd
# from snakemake.remote import FTP
from snakemake.utils import validate

# ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_fq1(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    assert len(fastqs) > 0
    return fastqs

def get_fq2(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()
    assert len(fastqs) > 0
    return fastqs

def get_final_output():
    final_output = expand(
        "results/remap_barcode/{sample}_{unit}_{seqtype}_S1_L001_R{reads}_001.fastq.gz",
        sample = samples["sample_name"],
        unit = units["unit_name"],
        seqtype=["DNA", "RNA"],
        reads=[1,2]
    )
    final_output.extend(
        expand("results/remap_barcode/ori_{sample}_{unit}_RNA_S1_L001_R1_001.fastq.gz",
        sample = samples["sample_name"],
        unit = units["unit_name"])
    )
    final_output.extend(
        expand("results/remap_barcode/ori_{sample}_{unit}_DNA_S1_L001_R2_001.fastq.gz",
        sample = samples["sample_name"],
        unit = units["unit_name"])
    )
    final_output.extend(
        expand("results/remap_barcode/{sample}_{unit}_barcodeMap.tsv.gz",
        sample = samples["sample_name"],
        unit = units["unit_name"])
    )
    if not config["ref"]["reference_transcriptome"]:
        final_output.append("resources/cellranger_transcriptome/"+config["ref"]["build"])
    final_output.append("results/qc/multiqc_report.html")
    return final_output