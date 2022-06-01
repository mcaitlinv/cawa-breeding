# Pre-processing

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

Conda environment and package versions for trimming, cutting, and individual qc:

```
conda list -n trimQC

# packages in environment at /projects/caitlinv@colostate.edu/software/anaconda/envs/trimQC:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main
_openmp_mutex             4.5                       1_gnu
bz2file                   0.98             py37h06a4308_1
ca-certificates           2022.2.1             h06a4308_0
certifi                   2021.10.8        py37h06a4308_2
cutadapt                  1.18             py37h14c3975_1    bioconda
dbus                      1.13.18              hb2f20db_0
expat                     2.4.4                h295c915_0
fastqc                    0.11.9               hdfd78af_1    bioconda
font-ttf-dejavu-sans-mono 2.37                 hd3eb1b0_0
fontconfig                2.13.1               h6c09931_0
freetype                  2.11.0               h70c0345_0
glib                      2.69.1               h4ff587b_1
icu                       58.2                 he6710b0_3
ld_impl_linux-64          2.35.1               h7274673_9
libffi                    3.3                  he6710b0_2
libgcc-ng                 9.3.0               h5101ec6_17
libgomp                   9.3.0               h5101ec6_17
libpng                    1.6.37               hbc83047_0
libstdcxx-ng              9.3.0               hd4cf53a_17
libuuid                   1.0.3                h7f8727e_2
libxcb                    1.14                 h7b6447c_0
libxml2                   2.9.12               h03d6c58_0
ncurses                   6.3                  h7f8727e_2
openjdk                   11.0.13              h87a67e3_0
openssl                   1.1.1m               h7f8727e_0
pcre                      8.45                 h295c915_0
perl                      5.26.2               h14c3975_0
pigz                      2.6                  h27cfd23_0
pip                       21.2.2           py37h06a4308_0
python                    3.7.11               h12debd9_0
readline                  8.1.2                h7f8727e_1
setuptools                58.0.4           py37h06a4308_0
sqlite                    3.38.0               hc218d9a_0
tk                        8.6.11               h1ccaba5_0
trim-galore               0.6.7                hdfd78af_0    bioconda
wheel                     0.37.1             pyhd3eb1b0_0
xopen                     0.7.3                      py_0    bioconda
xz                        5.2.5                h7b6447c_0
zlib                      1.2.11               h7f8727e_4
```
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
Conda environment and package versions for multi-sample qc:

```
conda list -n multiqc
# packages in environment at /projects/caitlinv@colostate.edu/software/anaconda/envs/multiqc:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main
_openmp_mutex             4.5                       1_gnu
backports                 1.1                pyhd3eb1b0_0
backports.functools_lru_cache 1.6.4              pyhd3eb1b0_0
backports_abc             0.5                      py27_0
blas                      1.0                         mkl
bzip2                     1.0.8                h7b6447c_0
c-ares                    1.18.1               h7f8727e_0
ca-certificates           2021.10.26           h06a4308_2
certifi                   2020.6.20          pyhd3eb1b0_3
click                     7.1.2              pyhd3eb1b0_0
curl                      7.80.0               h7f8727e_0
cycler                    0.10.0                   py27_0
dbus                      1.13.18              hb2f20db_0
expat                     2.4.4                h295c915_0
fontconfig                2.13.1               h6c09931_0
freetype                  2.11.0               h70c0345_0
functools32               3.2.3.2          py36hc1f28f5_1
futures                   3.3.0                    py27_0
gettext                   0.21.0               hf68c758_0
git                       2.34.1          pl5262hc120c5b_0
glib                      2.69.1               h4ff587b_1
gst-plugins-base          1.14.0               h8213a91_2
gstreamer                 1.14.0               h28cd5cc_2
icu                       58.2                 he6710b0_3
intel-openmp              2022.0.1          h06a4308_3633
jinja2                    2.11.3             pyhd3eb1b0_0
jpeg                      9d                   h7f8727e_0
kiwisolver                1.1.0            py27he6710b0_0
krb5                      1.19.2               hac12032_0
libcurl                   7.80.0               h0b77cf5_0
libedit                   3.1.20210910         h7f8727e_0
libev                     4.33                 h7f8727e_1
libffi                    3.3                  he6710b0_2
libgcc-ng                 9.3.0               h5101ec6_17
libgfortran-ng            7.5.0               ha8ba4b0_17
libgfortran4              7.5.0               ha8ba4b0_17
libgomp                   9.3.0               h5101ec6_17
libnghttp2                1.46.0               hce63b2e_0
libpng                    1.6.37               hbc83047_0
libssh2                   1.9.0                h1ba5d50_1
libstdcxx-ng              9.3.0               hd4cf53a_17
libuuid                   1.0.3                h7f8727e_2
libxcb                    1.14                 h7b6447c_0
libxml2                   2.9.12               h03d6c58_0
markupsafe                1.1.1            py27h7b6447c_0
matplotlib                2.2.3            py27hb69df0a_0
mkl                       2020.2                      256
mkl-service               2.3.0            py27he904b0f_0
mkl_fft                   1.0.15           py27ha843d7b_0
mkl_random                1.1.0            py27hd6b4f25_0
multiqc                   1.0.dev0                 pypi_0    pypi
ncurses                   6.3                  h7f8727e_2
numpy                     1.16.6           py27hbc911f0_0
numpy-base                1.16.6           py27hde5b4d6_0
openssl                   1.1.1m               h7f8727e_0
pcre                      8.45                 h295c915_0
pcre2                     10.37                he7ceb23_1
perl                      5.26.2               h14c3975_0
pip                       19.3.1                   py27_0
pyparsing                 2.4.7              pyhd3eb1b0_0
pyqt                      5.9.2            py27h05f1152_2
python                    2.7.18               ha1903f6_2
python-dateutil           2.8.2              pyhd3eb1b0_0
pytz                      2021.3             pyhd3eb1b0_0
pyyaml                    5.2              py27h7b6447c_0
qt                        5.9.7                h5867ecd_1
readline                  8.1.2                h7f8727e_1
setuptools                44.0.0                   py27_0
simplejson                3.17.0           py27h7b6447c_0
singledispatch            3.7.0           pyhd3eb1b0_1001
sip                       4.19.8           py27hf484d3e_0
six                       1.16.0             pyhd3eb1b0_1
sqlite                    3.37.2               hc218d9a_0
subprocess32              3.5.4            py27h7b6447c_0
tk                        8.6.11               h1ccaba5_0
tornado                   5.1.1            py27h7b6447c_0
wheel                     0.37.1             pyhd3eb1b0_0
xz                        5.2.5                h7b6447c_0
yaml                      0.1.7                had09818_2
zlib                      1.2.11               h7f8727e_4
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

Conda environment

```
conda list -n bwasam
# packages in environment at /projects/caitlinv@colostate.edu/software/anaconda/envs/bwasam:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main
_openmp_mutex             4.5                       1_gnu
_r-mutex                  1.0.0               anacondar_1
bwa                       0.7.17               h5bf99c6_8    bioconda
bzip2                     1.0.8                h7b6447c_0
c-ares                    1.18.1               h7f8727e_0
ca-certificates           2022.2.1             h06a4308_0
cairo                     1.16.0               hf32fb01_1
dbus                      1.13.18              hb2f20db_0
expat                     2.4.4                h295c915_0
fontconfig                2.13.1               h6c09931_0
freetype                  2.11.0               h70c0345_0
fribidi                   1.0.10               h7b6447c_0
glib                      2.69.1               h4ff587b_1
graphite2                 1.3.14               h23475e2_0
harfbuzz                  2.8.1                h6f93f22_0
htslib                    1.11                 hd3b49d5_2    bioconda
icu                       58.2                 he6710b0_3
jpeg                      9d                   h7f8727e_0
krb5                      1.19.2               hac12032_0
libcurl                   7.80.0               h0b77cf5_0
libdeflate                1.7                  h27cfd23_5
libedit                   3.1.20210714         h7f8727e_0
libev                     4.33                 h7f8727e_1
libffi                    3.3                  he6710b0_2
libgcc                    7.2.0                h69d50b8_2
libgcc-ng                 9.3.0               h5101ec6_17
libgomp                   9.3.0               h5101ec6_17
libnghttp2                1.46.0               hce63b2e_0
libpng                    1.6.37               hbc83047_0
libssh2                   1.9.0                h1ba5d50_1
libstdcxx-ng              9.3.0               hd4cf53a_17
libtiff                   4.2.0                h85742a9_0
libuuid                   1.0.3                h7f8727e_2
libwebp-base              1.2.2                h7f8727e_0
libxcb                    1.14                 h7b6447c_0
libxml2                   2.9.12               h03d6c58_0
lz4-c                     1.9.3                h295c915_1
ncurses                   6.2                  he6710b0_1
openjdk                   11.0.13              h87a67e3_0
openssl                   1.1.1m               h7f8727e_0
pango                     1.45.3               hd140c19_0
pcre                      8.45                 h295c915_0
perl                      5.26.2               h14c3975_0
picard                    2.26.11              hdfd78af_0    bioconda
pixman                    0.40.0               h7f8727e_1
r-base                    3.2.2                         0
readline                  8.1                  h27cfd23_0
samtools                  1.11                 h6270b1f_0    bioconda
tk                        8.6.11               h1ccaba5_0
xz                        5.2.5                h7b6447c_0
zlib                      1.2.11               h7f8727e_4
zstd                      1.4.9                haebb681_0
```

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

  
