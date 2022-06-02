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
1. Trimming adaptors and qc: 
2. Mapping and marking dups:
3. Checking sample coverage:

Each step of the workflow has 2 associated scripts, eg scripts 1a/1b are in the first step.

## 1) Trimming adaptors and quality control

This dataset was generated with both HiSeq and NovaSeq platforms. NovaSeq platforms are 2-color chemistry and can have poly-G tails due to poor quality sequence that gets called as G, see [this discussion](https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/).

  - To check read quality I used `FASTQC` before any processing. Then, I trimmed adaptors using `TrimGalore`. `TrimGalore` has a library of adaptors built in, so adaptors are not specified. After trimming I used a sliding window cut using `fastp`.The sliding window trimmed the 3' end of the read after 4 bases in a row were below a mean quality of 20. This is to remove any potentially poor-quality bases at the end of the read, as in NovaSeq bams these may get called as G's with high confidence.

  - After doing individual quality control, I check the entire dataset together using `MULTIQC`. This is assigned to a separate step as I needed all samples completed to run `MULTIQC` and it only works through conda when installed by itself.

Notes: I use a script [line_assign.sh]() to assign sample ID, plate, flowcell ID, lane ID, and an index number from a reference file called [numbered-units.txt](). This is used to connect the generally pretty dense fastq names to the sample ID and let me assign each pair of fastqs to their own job in a slurm array.

Trimming, cutting, and individual qc:

```bash
source ~/.bashrc
##This links sample ID, plate, and lane ID to the relevant fastq
##It also submits each line as a single job of an array
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID numbered-units.txt)
cd $1

##Activate conda
conda activate trimQC

##Set up directories to hold output
mkdir -p results/slurm_logs/trim/
mkdir -p results/trim/
mkdir -p results/qc/fqFASTQC
mkdir -p results/logs/qc/fqFASTQC
mkdir -p results/qc/trimFASTQC
mkdir -p results/logs/qc/trimFASTQC

T1=`basename $fq1 .fq.gz`_val_1.fq.gz
T2=`basename $fq2 .fq.gz`_val_2.fq.gz
O1=`basename $T1 _val_1.fq.gz`.trim.fq.gz
O2=`basename $T2 _val_2.fq.gz`.trim.fq.gz

##Initial QC step
fastqc -o results/qc/fqFASTQC --noextract fastqs/$fq1 fastqs/$fq2 > results/logs/qc/fqFASTQC/fastqc."$fq1".log 2>&1

###Trim low quality fragments and adapters
trim_galore -q 5 --paired --cores 8 --output_dir results/trim fastqs/$fq1 fastqs/$fq2

##Fastp for trimming the poly-g tails
/projects/caitlinv@colostate.edu/software/fastp --in1  results/trim/$T1 --out1 results/trim/$O1 \
	--in2 results/trim/$T2 --out2 results/trim/$O2  --thread 10 \
	-L -A --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20

##Post QC step
fastqc -o results/qc/trimFASTQC --noextract results/trim/$O1 results/trim/$O2 > results/logs/qc/trimFASTQC/fastqc."$O1".log 2>&1
```

Multi-sample qc:

```bash
source ~/.bashrc
cd $1

##Activate conda
conda activate multiqc

mkdir -p results/qc/multiqc

multiqc -f -v -o results/qc/multiqc results/qc/fqFASTQC results/qc/trimFASTQC
```

## 2) Mapping and marking dups

Canada Warbler reference genome can be found here. Before starting to map to the reference I indexed the genome. Note: the genome lives in a separate directory from where I run the analysis.

```bash
bwa index  /projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta
```

 - I start off with mapping paired fqs to the reference using `bwa mem`. Then I use `samtools sort` to sort the sam and create a bam, `picard AddOrReplaceReadGroups` to add read groups, `samtools sort` to sort by name to add mate score tags, `samtools fixmate` to add mate score tags, `samtools sort` to sort by position, and `samtools markdup` to mark duplicates. 

  - After getting bams that have duplicates marked for each fastq pair, I merge the fastqs by samples. In the previous step I added read groups where `-RGSM` was a sample ID that corresponded to a single individual so I can merge multiple bams into a single sample bam. Some of my samples only have one bam, so I added a simple if/else statement to take care of that. Note: This step substitutes [samples.txt]() for [numnbered-units.txt]() when assigning jobs in the array because now I only want a single job for each sample, instead of a job for each pair of fastqs.




Mapping and marking duplicates

```bash
source ~/.bashrc
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID numbered-units.txt)

cd $1
mkdir -p results/bams

##Variables: trimmed fq names, genome directory, new bam names
R1=`basename $fq1 .fq.gz`.trim.fq.gz
R2=`basename $fq2 .fq.gz`.trim.fq.gz
genomeDIR="/projects/caitlinv@colostate.edu/genomes/cardellina_canadensis_pseudohap_v1.fasta"
ID="$sample"."$library"."$flowcell"."$lane"

##Align each sample to genome. Note that genome reference must already be built through bwa
conda activate bwasam

##map paired, trimmed reads
bwa mem -t 8 $genomeDIR results/trim/$R1 results/trim/$R2 > results/bams/"$ID".sam

##Sort sam and create bam
samtools sort -o results/bams/"$ID".bam results/bams/"$ID".sam

##Add read groups, index RG bam
picard AddOrReplaceReadGroups -INPUT results/bams/"$ID".bam -RGID "$ID" -RGLB "$library" \
	-RGPL "$platform"."$flowcell"."$lane" -RGPU "$library"."$sample" -RGSM "$sample" \
	-OUTPUT results/bams/"$ID".RG.bam -VALIDATION_STRINGENCY SILENT
samtools index results/bams/"$ID".RG.bam

##Sort by name, not position
samtools sort -n -o results/bams/"$ID".namesort.bam results/bams/"$ID".RG.bam

##Add mate score tags for samtools markdup to select best reads
samtools fixmate -m results/bams/"$ID".namesort.bam results/bams/"$ID".fixm.bam

##Sort again by position
samtools sort -o results/bams/"$ID".fixm.sort.bam  results/bams/"$ID".fixm.bam

##Markdups
samtools markdup results/bams/"$ID".fixm.sort.bam results/bams/"$ID".mkdup.bam
samtools index results/bams/"$ID".mkdup.bam
```

Merging samples and marking duplicates again.

```bash
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc

cd $1
conda activate bwasam
##Get all bams of a single sample name
ls -d results/bams/*"$sample"*.mkdup.bam > results/bams/"$sample".txt

##Get number of bams for each sample
wcBam=$(wc -l < results/bams/"$sample".txt)

##If only 1 bam
##skip merge,change name, remove intermediate bams
##If more than 1 bam
##Merge bams, mark duplicates, remove intermediate bams
if expr "$wcBam" : '[1]'
then
echo "Single bam, skipping merge"
mv  results/bams/"$sample"*.mkdup.bam results/bams/"$sample".mergeMkDup.bam
mv  results/bams/"$sample"*.mkdup.bam.bai results/bams/"$sample".mergeMkDup.bam.bai
else
##merge all bams from file of bams
samtools merge --threads 8 -f results/bams/"$sample".merged.bam -b results/bams/"$sample".txt

##Removing PCR duplicates from BAM files
#sort by name, not position
samtools sort -n -o  results/bams/"$sample".namesort.bam  results/bams/"$sample".merged.bam
rm  results/bams/"$sample".merged.bam

#Add mate score tags for samtools markdup to select best reads
samtools fixmate -m  results/bams/"$sample".namesort.bam  results/bams/"$sample".fixm.bam
rm results/bams/"$sample".namesort.bam

#Sort again by posiition
samtools sort -o  results/bams/"$sample".fixm.sort.bam   results/bams/"$sample".fixm.bam
rm  results/bams/"$sample".fixm.bam

#Markdups
samtools markdup  results/bams/"$sample".fixm.sort.bam  results/bams/"$sample".mergeMkDup.bam
samtools index  results/bams/"$sample".mergeMkDup.bam
rm  results/bams/"$sample".fixm.sort.bam
fi
```

## 3) Getting coverage for each sample

After getting the merged bams for each sample, I checked how coverage looked for each of the samples. Ideal target was >2X average. 

- I start by getting the coverage at all positions covered for a sample by using `bedtools`. I used an array where each segment was a sample.   

- Then I used a script that averages coverage across a sample and pastes it into a text file. 

Note that the individual coverage files get pretty big pretty fast, so an alternate strategy would be to calculate per-position coverage, get the average and then delete the individual coverage file in a loop, so one sample is calculated at a time. I did not do that because the trade-off is it takes forever. 

Get individual coverage at each bp:

```bash
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID samples.txt)
source ~/.bashrc
cd $1

##Activate conda
conda activate bedtools

bedtools genomecov -d -ibam results/bams/"$sample".mergeMkDup.bam  > results/cov/"$sample".genomecov.txt
```
Get average coverage for all individuals:

```bash
source ~/.bashrc
cd $1

touch results/cov/covSummary.txt

for sample in `ls results/cov/*.genomecov.txt`
  	do
        awk '{if($3<500) {total+=$3; ++lines}} END {print FILENAME," ",total/lines}' \
        $sample >> results/cov/covSummary.txt
done
```
