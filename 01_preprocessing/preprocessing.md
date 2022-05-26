# Pre-processing

The steps I used to go from raw fastqs to bams to be used in variant calling.

Workflow
1. Trimming adaptors and remove low quality tails:
2. Mapping to reference:
3. Marking duplicates:
4. Merging samples with multiple bams and marking duplicates again:

## 1) Trimming adaptors and remove low quality tails

This dataset was generated with both HiSeq and NovaSeq platforms. NovaSeq platforms are 2-color chemistry and can have poly-G tails due to poor quality sequence that gets called as G, see [this discussion](https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/). 

In addition, I found that despite the sequencing centers processing to cut adaptors, I still had adaptor content. So I trimmed adaptors and used a sliding window cut on the 3' end of the read. The sliding window trimmed bases after 4 bases in a row were below 20.

Code overview
```
module purge
source ~/.bashrc
eval $(line_assign.sh $SLURM_ARRAY_TASK_ID numbered-units.txt)
cd $1

##Activate conda
conda activate trimQC

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

fastqc -o results/qc/fqFASTQC --noextract fastqs/$fq1 fastqs/$fq2 > results/logs/qc/fqFASTQC/fastqc."$fq1".log 2>&1



###Trim low quality fragments and adapters
trim_galore -q 5 --paired --cores 8 --output_dir results/trim fastqs/$fq1 fastqs/$fq2


##Fastp for trimming the poly-g tails
/projects/caitlinv@colostate.edu/software/fastp --in1  results/trim/$T1 --out1 results/trim/$O1 \
	--in2 results/trim/$T2 --out2 results/trim/$O2  --thread 10 \
	-L -A --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20

fastqc -o results/qc/trimFASTQC --noextract results/trim/$O1 results/trim/$O2 > results/logs/qc/trimFASTQC/fastqc."$O1".log 2>&1
```