# Pre-processing

The steps I used to go from raw fastqs to bams to be used in variant calling. I primarily used [Conda/Mamba](https://docs.conda.io/en/latest/) install and use software. Conda environments and loaded packages are listed with the overall commands.

Notes: I start all my scripts in a directory for the project that I use as the first command line argument.

Workflow
1. Trimming adaptors and qc: 
2. Mapping to reference:
3. Marking duplicates:
4. Merging samples with multiple bams and marking duplicates again:

## 1) Trimming adaptors and quality control

This dataset was generated with both HiSeq and NovaSeq platforms. NovaSeq platforms are 2-color chemistry and can have poly-G tails due to poor quality sequence that gets called as G, see [this discussion](https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/).

  - To check read quality I used FASTQC before any processing. Then, I trimmed adaptors using TrimGalore. TrimGalore has a library of adaptors built in, so adaptors are not specified. After trimming I used a sliding window cut using Fastp.The sliding window trimmed the 3' end of the read after 4 bases in a row were below a mean quality of 20. This is to remove any potentially poor-quality bases at the end of the read, as in NovaSeq bams these may get called as G's with high confidence.

  - After doing individual quality control, I check the entire dataset together using MULTIQC. This is assigned to a separate step as I needed all samples completed to run MULTIQC and it only works through conda when installed by itself.

Notes: I use a [script]() to assign sample ID, plate, flowcell ID, lane ID, and an index number from a reference file called [numbered-units.txt](). This is used to connect the generally pretty dense fastq names to the sample ID and let me assign each pair of fastqs to their own job in a slurm array.

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

```{bash, eval=F}
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

```{bash}
source ~/.bashrc
cd $1

##Activate conda
conda activate multiqc

mkdir -p results/qc/multiqc

multiqc -f -v -o results/qc/multiqc results/qc/fqFASTQC results/qc/trimFASTQC
```

