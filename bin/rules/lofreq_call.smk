#############################################################################
##### 			Low-frequency variants analysis  		#####
#############################################################################

rule lofreq_call:
    input:
        bam=str(OUT / "bwa_alignment/{sample}_markdup.bam")
    output:
        viterbi=temp(str(OUT / "lofreq/{sample}_viterbi.bam")),
        aln=temp(str(OUT / "lofreq/{sample}_lofreq_aln.bam")),
        index=temp(str(OUT / "lofreq/{sample}_lofreq_aln.bai")),
        var=str(OUT / "lofreq/{sample}.vcf")
    conda:
        "../../envs/lofreq.yaml"
    benchmark:
        str(OUT / "log/benchmark/lofreq_{sample}.txt")
    threads: config["threads"]["lofreq"]
    resources:
        mem_mb=config["mem_mb"]["lofreq"]
    log:
        str(OUT / "log/lofreq/lofreq_{sample}.log")
    params:
        ref_genome=config["bwa"]["ref_genome"]
    shell:
        """
lofreq viterbi -f {params.ref_genome} {input} | \
samtools sort -O BAM | tee {output.viterbi} | \
lofreq indelqual -f {params.ref_genome} --dindel - | \
samtools view -b - > {output.aln} 2> {log}

samtools index -@ {threads} {output.aln} {output.index} 2> {log}

lofreq call-parallel --pp-threads {threads} --call-indels -f {params.ref_genome} -o {output.var} {output.aln} 2> {log}

        """

#lofreq alnqual -b -u - {params.ref_genome} | \ ####after indelqual
