"""
Myco_lofreq pipeline
Authors: Alejandra Hernandez Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: IDS - BPD - Bacteriology
Date: 22-09-2020


Snakemake rules (in order of execution):
    1 fastQC        # Asses quality of raw reads.
    2 trimmomatic   # Trim low quality reads and adapter sequences.
    3 fastQC        # Asses quality of trimmed reads.
    4 multiqc	    # MultiQC report of all individual samples
    5 bwa           # Alignment to reference genome.
    6 lofreq        # Detecting low frequency variants.

"""

#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

configfile: "config/pipeline_parameters.yaml"
configfile: "config/variables.yaml"

from pandas import *
import pathlib
#import pprint
import os
import yaml
import json


yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings

#################################################################################
##### Load samplesheet, load genus dict and define output directory         #####
#################################################################################

# SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = pathlib.Path(config["out"])

#@################################################################################
#@#### 				Processes                                    #####
#@################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################
include: "bin/rules/fastqc_raw_data.smk"
include: "bin/rules/trimmomatic.smk"
include: "bin/rules/fastqc_trimmed_data.smk"
include: "bin/rules/cat_unpaired.smk"
include: "bin/rules/multiqc.smk"
    #############################################################################
    #####    Alignment to Reference Genome                                  #####
    #############################################################################
include: "bin/rules/bwa_alignment.smk"
#include: "bin/rules/picard_markeddup.smk"

    #############################################################################
    ##### Low-frequency variants analysis    #####
    #############################################################################
include: "bin/rules/lofreq_call.smk"
include: "bin/rules/vcf_annotator.smk"


#@################################################################################
#@#### The `onstart` checker codeblock                                       #####
#@################################################################################

onstart:
    try:
        print("Checking if all specified files are accessible...")
        important_files = [ config["sample_sheet"],
                         'files/trimmomatic_0.36_adapters_lists/NexteraPE-PE.fa' ]
        for filename in important_files:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
    except FileNotFoundError as e:
        print("This file is not available or accessible: %s" % e)
        sys.exit(1)
    else:
        print("\tAll specified files are present!")
    shell("""
        mkdir -p {OUT}
        mkdir -p {OUT}/results
        echo -e "\nLogging pipeline settings..."
        echo -e "\tGenerating methodological hash (fingerprint)..."
        echo -e "This is the link to the code used for this analysis:\thttps://github.com/AleSR13/Myco_lofreq/tree/$(git log -n 1 --pretty=format:"%H")" > '{OUT}/results/log_git_myco_lofreq.txt'
        echo -e "This code with unique fingerprint $(git log -n1 --pretty=format:"%H") was committed by $(git log -n1 --pretty=format:"%an <%ae>") at $(git log -n1 --pretty=format:"%ad")" >> '{OUT}/results/log_git_myco_lofreq.txt'
        echo -e "\tGenerating full software list of current Conda environment (\"myco_lofreq_master\")..."
        conda list > '{OUT}/results/log_conda_myco_lofreq.txt'
        echo -e "\tGenerating config file log..."
        rm -f '{OUT}/results/log_config_myco_lofreq.txt'
        for file in config/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> '{OUT}/results/log_config_myco_lofreq.txt'
            cat ${{file}} >> '{OUT}/results/log_config_myco_lofreq.txt'
            echo -e "\n\n" >> '{OUT}/results/log_config_myco_lofreq.txt'
        done
    """)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

#onerror:
 #   shell("""""")


onsuccess:
    shell("""
        echo -e "\tGenerating HTML index of log files..."
        echo -e "\tGenerating Snakemake report..."
        snakemake --unlock --profile config --config out={OUT}
        snakemake --profile config --config out={OUT} --report '{OUT}/results/snakemake_report_myco_lofreq.html'
        echo -e "Finished"
    """)


#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all,
    cat_unpaired


rule all:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['R1', 'R2']),   
        expand(str(OUT / "trimmomatic/{sample}_{read}.fastq"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        str(OUT / "MultiQC/multiqc.html"),
        expand(str(OUT / "annotated_vcf/{sample}_annotated.vcf"), sample = SAMPLES)
