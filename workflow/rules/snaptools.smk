rule snaptools:
    input:
        r1 = "results/remap_barcode/{sample}_{unit}_DNA_S1_L001_R1_001.fastq.gz",
        r2 = "results/remap_barcode/{sample}_{unit}_DNA_S1_L001_R2_001.fastq.gz",
        ref= config["ref"]["bwa_index"] if "bwa_index" in config["ref"].keys() else rules.bwa_index.output,
        gs = config["ref"]["chrom_size"] if "chrom_size" in config["ref"].keys() else "resources/genome.sizes",
    output:
        frag = "results/snaptools/outs/fixed/{sample}_{unit}.fragments.tsv.gz",
    params:
        ucscname      = config["ref"]["ucscname"],
    log:
        "logs/snaptools/{sample}_{unit}.log",
    threads: config["computingResources"]["threads"]["high"],
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    conda:
        "../envs/snaptools.yaml"
    shell:
        """
        mkdir -p results/snaptools/outs/bam
        mkdir -p results/snaptools/outs/fixed
        mkdir -p results/snaptools/tmp
        echo $(date +"%T")"align-paired-end..."
        snaptools align-paired-end \
        	--input-reference={input.ref} \
        	--input-fastq1={input.r1} \
        	--input-fastq2={input.r2} \
        	--output-bam=results/snaptools/outs/bam/{wildcards.sample}_{wildcards.unit}.bam \
        	--aligner=bwa \
        	--read-fastq-command=zcat \
        	--min-cov=0 \
        	--num-threads={threads} \
        	--if-sort=True \
        	--tmp-folder=results/snaptools/tmp/ \
        	--overwrite=TRUE

        echo $(date +"%T")"snap-pre..."
        snaptools snap-pre \
        	--input-file=results/snaptools/outs/bam/{wildcards.sample}_{wildcards.unit}.bam \
        	--output-snap=results/snaptools/outs/{wildcards.sample}_{wildcards.unit}.snap \
        	--genome-name={params.ucscname} \
        	--genome-size={input.gs} \
        	--min-mapq=10 \
        	--min-flen=0 \
        	--keep-chrm=TRUE \
        	--keep-single=False \
        	--keep-secondary=False \
        	--overwrite=True \
        	--max-num=10000000 \
        	--verbose=True

        echo $(date +"%T")"dump-fragments..."
        snaptools dump-fragment \
        	--snap-file=results/snaptools/outs/{wildcards.sample}_{wildcards.unit}.snap \
        	--output-file=results/snaptools/tmp/{wildcards.sample}_{wildcards.unit}.dump.frag.tsv.gz \
        	--buffer-size=10000 \
        	--tmp-folder=results/snaptools/tmp

        echo $(date +"%T")"sort fragments..."
        zcat results/snaptools/tmp/{wildcards.sample}_{wildcards.unit}.dump.frag.tsv.gz | \
        	perl -p -e "s/b\'//g" | \
        	perl -p -e "s/\'//g" | \
        	sort -V -k1,1 -k2,2n --stable  --parallel={threads} --temporary-directory=results/snaptools/tmp  -S 30G | \
            pbgzip -n {threads} -c > results/snaptools/outs/{wildcards.sample}_{wildcards.unit}.snaptools.frag.tsv.gz
        tabix -p bed results/snaptools/outs/{wildcards.sample}_{wildcards.unit}.snaptools.frag.tsv.gz

        echo $(date +"%T")"fix the fragments for ArchR..."
        zcat results/snaptools/outs/{wildcards.sample}_{wildcards.unit}.snaptools.frag.tsv.gz | \
            awk -v sample={wildcards.sample}_{wildcards.unit} 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,sample"_"$4, 1}}' | \
            pbgzip -n {threads} -c > results/snaptools/outs/fixed/{wildcards.sample}_{wildcards.unit}.fragments.tsv.gz
          tabix -p bed results/snaptools/outs/fixed/{wildcards.sample}_{wildcards.unit}.fragments.tsv.gz
        """
