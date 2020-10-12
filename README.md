# Myco_lofreq pipeline

The main purpose of this pipeline is to call minority variants in _Mycobacterium_ samples. I takes paired-end raw fastq files as input (one should contain R1 and the other R2 on the name). The pipeline performs the following steps:

1. FastQC on the raw reads
2. Trimming with Trimmomatic
3. FastQC on the trimmed reads
4. Alignment to reference genome (using BWA)
5. Preparation to load in lofreq
6. Calling SNP and indels using LoFreq
7. Annotating the resulting vcf file

## Basic Usage

To run the pipeline you need to download this repository. Afterwards you need to open the command line in your Linux environment and enter the directory containing the code. For instance:

```cd /path/to/code/Myco_lofreq/```

Then you need to run the pipeline as follows:

```bash run_myco_lofreq_pipeline.sh -i <path/to/input_folder> -o <path/to/output_folder>```


