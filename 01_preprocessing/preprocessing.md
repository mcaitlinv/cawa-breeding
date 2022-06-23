# Pre-processing

Package versions
```
cutadapt 1.18
fastqc  0.11.9
trim-galore 0.6.7
fastp 0.22.0
multiqc 1.0.dev0  
bwa 0.7.17
picard  2.26.11
samtools  1.11
bedtools 2.30.0
```

The steps I used to go from raw fastqs to bams to be used in variant calling. I primarily used [Conda/Mamba](https://docs.conda.io/en/latest/) install and use software. Conda environments and loaded packages are listed with the overall commands.

Notes: I start all my scripts in a directory for the project that I use as the first command line argument.

Basic directory structure

```
cawa-breeding-lcwg  
-->number-units.txt
-->sample.txt
-->fastqs
---->a.fastq.R1.fq.gz
---->a.fastq.R2.fq.gz
-->results
---->bams
------>a.merged.mkDup.bam
---->qc
---->slurm-logs
---->trim
------>a.fastq.R1.trimmed.fq.gz
------>a.fastq.R2.trimmed.fq.gz
```

Workflow
1. Trimming adaptors and qc
2. Mapping and marking dups
3. Checking sample coverage

Each step of the workflow has 2 associated scripts, eg scripts 1a/1b are in the first step.

## 1) Trimming adaptors and quality control

Scripts: 
  1) [1a.conda.trimgalore.fastqc.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/1a.conda.trimgalore.fastqc.sbatch)
  2) [1b.conda.multiqc.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/1b.conda.multiqc.sbatch)

This dataset was generated with both HiSeq and NovaSeq platforms. NovaSeq platforms are 2-color chemistry and can have poly-G tails due to poor quality sequence that gets called as G, see [this discussion](https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/).

  - To check read quality I used `FASTQC` before any processing. Then, I trimmed adaptors using `TrimGalore`. `TrimGalore` has a library of adaptors built in, so adaptors are not specified. After trimming I used a sliding window cut using `fastp`.The sliding window trimmed the 3' end of the read after 4 bases in a row were below a mean quality of 20. This is to remove any potentially poor-quality bases at the end of the read, as in NovaSeq bams these may get called as G's with high confidence.

  - After doing individual quality control, I check the entire dataset together using `MULTIQC`. This is assigned to a separate step as I needed all samples completed to run `MULTIQC` and it only works through conda when installed by itself.

Notes: I use a script [line_assign.sh]() to assign sample ID, plate, flowcell ID, lane ID, and an index number from a reference file called [numbered-units.txt](). This is used to connect the generally pretty dense fastq names to the sample ID and let me assign each pair of fastqs to their own job in a slurm array.

## 2) Mapping and marking dups

Scripts: 
  1) [2a.conda.map.bwa.trim.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/2a.conda.map.bwa.trim.sbatch)
  2) [2b.conda.merge.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/2b.conda.merge.sbatch)

Canada Warbler reference genome can be found here. Before starting to map to the reference I indexed the genome. Note: the genome lives in a separate directory from where I run the analysis.

```bash
bwa index  /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta
```

 - I start off with mapping paired fqs to the reference using `bwa mem`. Then I use `samtools sort` to sort the sam and create a bam, `picard AddOrReplaceReadGroups` to add read groups, `samtools sort` to sort by name to add mate score tags, `samtools fixmate` to add mate score tags, `samtools sort` to sort by position, and `samtools markdup` to mark duplicates. 

  - After getting bams that have duplicates marked for each fastq pair, I merge the fastqs by samples. In the previous step I added read groups where `-RGSM` was a sample ID that corresponded to a single individual so I can merge multiple bams into a single sample bam. Some of my samples only have one bam, so I added a simple if/else statement to take care of that. Note: This step substitutes [samples.txt]() for [numnbered-units.txt]() when assigning jobs in the array because now I only want a single job for each sample, instead of a job for each pair of fastqs.

## 3) Getting coverage for each sample

Scripts: 
  1) [3a.conda.coverage.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/3a.conda.coverage.sbatch)
  2) [3b.conda.coverage.sbatch](https://github.com/mcaitlinv/cawa-breeding/blob/main/01_preprocessing/scripts/3b.conda.coverage.sbatch)


After getting the merged bams for each sample, I checked how coverage looked for each of the samples. Ideal target was >2X average. 

- I start by getting the coverage at all positions covered for a sample by using `bedtools`. I used an array where each segment was a sample.   

- Then I used a script that averages coverage across a sample and pastes it into a text file. 

Note that the individual coverage files get pretty big pretty fast, so an alternate strategy would be to calculate per-position coverage, get the average and then delete the individual coverage file in a loop, so one sample is calculated at a time. I did not do that because the trade-off is it takes forever. 

