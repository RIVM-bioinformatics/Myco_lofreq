#############################################################################
##### 			Alignment to reference genome  			#####
#############################################################################

rule bwa_alignment:
    input:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq")
    output:
        #bam=str(OUT / "bwa_alignment/{sample}_alignment.bam"),
        sorted=temp(str(OUT / "bwa_alignment/{sample}_sorted.bam")),
        markdup=str(OUT / "bwa_alignment/{sample}_markdup.bam")
    conda:
        "../../envs/bwa.yaml"
    benchmark:
        str(OUT / "log/benchmark/bwa_{sample}.txt")
    threads: config["threads"]["bwa"]
    resources:
        mem_mb=config["mem_mb"]["bwa"]
    params:
        ref_genome=config["bwa"]["ref_genome"]
    log:
        str(OUT / "log/bwa/bwa_alignment_{sample}.log")
    shell:
        """
bwa mem {params.ref_genome} {input} \
| samtools fixmate -m - - \
| samtools sort -O BAM \
| tee {output.sorted} \
| samtools markdup -r - {output.markdup}
    """

