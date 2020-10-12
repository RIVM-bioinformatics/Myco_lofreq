# Myco_lofreq pipeline

The main purpose of this pipeline is to call minority variants in _Mycobacterium_ samples. I takes paired-end raw fastq files as input (one should contain R1 and the other R2 on the name). The pipeline performs the following steps:

1. FastQC on the raw reads
2. Trimming with Trimmomatic
3. FastQC on the trimmed reads
4. Alignment to reference genome (using BWA)
5. Preparation to load in lofreq
6. Calling minority variants using LoFreq

## Basic Usage

```bash run_myco_lofreq_pipeline.sh -i <input_dolder> -o <output_folder>```


