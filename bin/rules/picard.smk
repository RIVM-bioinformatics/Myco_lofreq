#############################################################################
##### 			Alignment to reference genome  			#####
#############################################################################

rule picard:
    input:
        bam=str(OUT / "bwa_alignment/{sample}_alignment.bam")
    output:
        bam=str(OUT / "bwa_alignment/{sample}_alignment_markeddup.bam"),
        txt=str(OUT / "bwa_alignment/{sample}_alignment_markeddup.txt")
    conda:
        "../../envs/bwa.yaml"
    benchmark:
        str(OUT / "log/benchmark/picard_{sample}.txt")
    threads: config["threads"]["picard"]
    log:
        str(OUT / "log/picard/picard_markeddup_{sample}.log")
    shell:
        """
java -jar picard.jar MarkDuplicates \
      I={input} \
      O={output.bam} \
      M={output.txt}
    """
