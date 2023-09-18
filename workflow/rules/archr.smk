rule archr:
    input:
        expand("results/snaptools/outs/fixed/{sample}_{unit}.fragments.tsv.gz",
                sample = samples["sample_name"],
                unit = units["unit_name"])
    output:
        archr = "results/archr/"
    log:
        "logs/archr.log",
    threads: config["computingResources"]["threads"]["high"],
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    conda:
        "../envs/archr.yaml"
    params:
        scr = workflow.source_path("../scripts/archr.R"),
        gn  = config["ref"]["ucscname"]
    shell:
        """
        # swith R1 and R2 for RNAseq
        Rscript {params.scr} \
            --cores {threads} \
            --input 'results/snaptools/outs/fixed' \
            --output 'results/archr' \
            --genome {params.gn}
        """
