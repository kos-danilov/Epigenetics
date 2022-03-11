#!/usr/bin/env bash
set -euo pipefail

#------------ [1] -------------------------------------
# Check conda env with tools for today
cat /mnt/meth/conda_envs/meth.yaml

# Activate conda env with required tools
conda activate meth

# ------------ [2] -------------------------------------
# Put data to home directory
mkdir -p ~/bs_seq

# ------------ [3] -------------------------------------
# Change working directory:
cd ~/bs_seq
mkdir reads_full
cd reads_full

fastq-dump --gzip SRR7449824
mv SRR7449824.fastq.gz DC552_HI48.fastq.gz

fastq-dump --gzip SRR7449827
mv SRR7449827.fastq.gz DC552_NI48.fastq.gz

fastq-dump --gzip SRR7449835
mv SRR7449835.fastq.gz DC555_HI48.fastq.gz

fastq-dump --gzip SRR7449839
mv SRR7449839.fastq.gz DC555_NI48.fastq.gz

# ------------ [4] -------------------------------------
# Reads QC for 4 samples:

cd ~/bs_seq
mkdir -p reads_qc/fastqc
mkdir -p reads_qc/fastp

# Find files:
find reads_full -name "*.fastq.gz"

# * find + FastQC:
find reads_full -name "*.fastq.gz" | xargs -I {} fastqc {} --outdir reads_qc/fastqc

# * find + FastP:
for F in $(find reads_full -name "*.fastq.gz"); do FN=$(basename $F); echo $FN; fastp --overrepresentation_analysis --thread 2 --in1 $F --html reads_qc/fastp/${FN/.fastq.gz/.fastp.html} --json reads_qc/fastp/${FN/.fastq.gz/.fastp.json}; done

# MultiQC
multiqc -f -o qc -n reads_qc reads_qc/fastqc reads_qc/fastp

# Check: ~/bs_seq/qc/reads_qc.html report

# ------------ [5] -------------------------------------
# Trim Reads using Trim Galore

mkdir -p ~/bs_seq/reads_trim
cd ~/bs_seq/reads_trim

# Single Sample:
# trim_galore --gzip --fastqc ../reads/DC552_NI48_chr16.fastq.gz

# All Samples
find ../reads_full -name "*.fastq.gz" | xargs -I {} trim_galore --gzip --fastqc {}

# MultiQC
cd ~/bs_seq
multiqc -f -o qc -n reads_trim reads_trim

mv SRR7449824_trimmed.fq.gz DC552_HI48_trimmed.fq.gz
mv SRR7449827_trimmed.fq.gz DC552_NI48_trimmed.fq.gz
mv SRR7449835_trimmed.fq.gz DC555_HI48_trimmed.fq.gz
mv SRR7449839_trimmed.fq.gz DC555_NI48_trimmed.fq.gz

# Check: ~/bs_seq/qc/reads_trim.html report

# ------------ [6] -------------------------------------
# Align to lambda genome: Build Bismark Indexes

mkdir -p ~/bs_seq/lambda

# Bismark expects `*.fa` files (not .fa.gz) in `indexes`, e.g. `indexes/hg38.fa`

# Download lambda genome:
mkdir -p ~/bs_seq/lambda/indexes
cd ~/bs_seq/lambda/indexes

# Download from NCBI: (https://www.ncbi.nlm.nih.gov/assembly/GCF_000840245.1/)
wget -O NC_001416.1.fa.gz 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz'
# In case of download issues, copy from:
cp /mnt/meth/genomes/NC_001416.1.fa.gz .

# Check content of FASTA file:
zcat NC_001416.1.fa.gz | head

# Build indexes:
cd ~/bs_seq/lambda
bismark_genome_preparation --verbose --bowtie2 --parallel 2 indexes 2>&1 | tee indexes/Bisulfite_Genome.log

# Inspect indexes
tree -h indexes

# ------------ [7] -------------------------------------
# Align to lambda genome: Align Reads

# Bismark expects reads to be in working directory, so symlink FASTQ files:
cd ~/bs_seq/lambda
mkdir -p reads_trim
cd reads_trim

for F in $(find ../../reads_trim -name "*.fq.gz"); do ln -s $F $(basename $F); done

# check files
tree

# Align reads
cd ~/bs_seq/lambda
mkdir -p bams

# All samples:
for F in $(find reads_trim -name "*.fq.gz"); do FN=$(basename $F); date; echo $FN; bismark --parallel 2 --output_dir bams --genome indexes --se $F >bams/${FN/.fq.gz/.log} 2>&1; done

# Check Conversion rate:
# Summary (multiqc doesn't show it, possible, but let's check manually)
find bams -name "*_SE_report.txt" | xargs -I {} bash -c 'echo {}; cat {} | grep "C methylated"'

# ------------ [8] -------------------------------------
# Align to hg38 (let's use only chr16)
cd ~/bs_seq

# Bismark Indexes:
mkdir -p ~/bs_seq/indexes
cd ~/bs_seq/indexes

# Download from UCSV:
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr16.fa.gz
# In case of download issues, copy from:
# cp /mnt/meth/genomes/chr16.fa.gz .

# Check chr content:
# * chr start
zcat chr16.fa.gz | head
# * some place far from start
zcat chr16.fa.gz | head -n 1000 | tail

# Build indexes: (bismark accepts *.fa, *.fasta, *.fa.gz, *.fasta.gz genomes)
cd ~/bs_seq
bismark_genome_preparation --verbose --bowtie2 --parallel 2 indexes 2>&1 | tee indexes/Bisulfite_Genome.log

# Align
cd ~/bs_seq
mkdir -p bams

# 1 sample example:
#bismark --parallel 2 --output_dir bams indexes reads_trim/DC552_TB48_chr16_trimmed.fq.gz > bams/DC552_TB48_chr16_trimmed.log 2>&1

# All samples:
for F in $(find reads_trim -name "*.fq.gz"); do FN=$(basename $F); date; echo $FN; bismark --nucleotide_coverage --parallel 2 --output_dir bams --genome indexes --se $F >bams/${FN/.fq.gz/.log} 2>&1; done

# Notice temp files like *_chr16_trimmed.fq.gz.temp.2_C_to_T.fastq
ls

# Summary (multiqc doesn't show it, possible to fix, but let's check manually)
find bams -name "*_SE_report.txt" | xargs -I {} bash -c 'echo {}; cat {} | grep "C methylated"'

# ------------ [9] -------------------------------------
# Bismark Deduplicate:
cd ~/bs_seq

# 1 sample example:
# deduplicate_bismark --bam --output_dir bams bams/DC552_TB48_chr16_trimmed_bismark_bt2.bam > bams/DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.bam.log 2>&1

# All samples:
for F in $(find bams -name "*_bt2.bam"); do echo $F; deduplicate_bismark --bam --output_dir bams $F >${F/.bam/.deduplicated.bam.log} 2>&1; done

# FastQC
mkdir -p bams_qc
find bams -name "*_bt2.deduplicated.bam" | xargs -I {} fastqc {} --outdir bams_qc

# MultiQC
multiqc -f -o qc -n bams bams bams_qc

# Check Report: ~/bs_seq/qc/bams.html

# ------------ [10] -------------------------------------
# Prepare BAMs for IGV (sort & index)

cd ~/bs_seq

# Sort & Index using Picard tools:
mkdir -p bams_sorted

# Single sample
# picard -Xmx4g SortSam --CREATE_INDEX true --SORT_ORDER coordinate --INPUT bams/DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.bam --OUTPUT bams_sorted/DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.sorted.bam

# All samples:
for F in $(find bams -name "*_bt2.deduplicated.bam"); do FN=$(basename $F); echo $FN; picard -Xmx4g SortSam --CREATE_INDEX true --SORT_ORDER coordinate --INPUT $F --OUTPUT bams_sorted/${FN/.bam/.sorted.bam} >"bams_sorted/${FN/.bam/.sorted.bam}.log" 2>&1; done

# check files:
tree bams_sorted

# Download files:
#/home/student/bs_seq/bams_sorted/DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated.sorted.bai
#/home/student/bs_seq/bams_sorted/DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated.sorted.bam

# * check: chr16:67,690,602-67,724,156
# * enable BS-Seq: CG coloring
# * check : chr16:67,720,273-67,720,309

# Another example: WGBS data
cd ~/bs_seq
ln -s /mnt/meth/bs_seq_examples
# * download '~/bs_seq_examples`
# * new 'hg18' session in IGV
# * open: chr11:8,973,794-8,988,048

# ------------ [11] -------------------------------------
#  Call methylation
cd ~/bs_seq

# Single sample
#bismark_methylation_extractor --gzip --comprehensive --ample_memory -o bams bams/DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.bam > bams/DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.bam.meth.log 2>&1

# All samples
for F in $(find bams -name "*_bt2.deduplicated.bam"); do echo $F; bismark_methylation_extractor --gzip --comprehensive --ample_memory -o bams $F >"$F.meth.log" 2>&1; done

# ------------ [12] -------------------------------------
# Bismark Methylation Results: Intermediate files
cd ~/bs_seq

# Check files for 1 sample:
find bams -name "*DC552_NI48*"
# * contexts splitting report
cat bams/DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated_splitting_report.txt
# * M-Bias
head bams/DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated.M-bias.txt
# * Methylation info: CpG context
# <read name> <strand> <chr> <position> <bismark methylation flag>
zcat bams/CpG_context_DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated.txt.gz | head
# * Check Bismark BAM file:
samtools view bams/DC552_NI48_chr16_trimmed_bismark_bt2.deduplicated.bam | head -n 1

# ------------ [13] -------------------------------------
# Export CpG methylome file
cd ~/bs_seq
mkdir meth_CpG

# Single sample
#bismark2bedGraph --dir meth_CpG --output DC552_TB48_chr16.bedGraph.gz  bams/CpG_context_DC552_TB48_chr16_trimmed_bismark_bt2.deduplicated.txt.gz

# All samples
for F in $(find bams -name "CpG*.txt.gz"); do FN=$(basename $F); echo $F; FN=${FN/CpG_context_/}; FN=${FN/_trimmed_bismark_bt2.deduplicated.txt.gz/.bedGraph.gz}; bismark2bedGraph --dir meth_CpG --output $FN $F > meth_CpG/"$FN.log" 2>&1; done

# Check files
ls meth_CpG

# Methylation info
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
zcat meth_CpG/DC552_NI48_chr16.bismark.cov.gz | head

# BedGraph with methylation info:
zcat meth_CpG/DC552_NI48_chr16.bedGraph.gz | head

# ------------ [14] -------------------------------------
# Bismark Methylation Results: QC
cd ~/bs_seq
mkdir meth_qc
cd bams
bismark2report --dir ../meth_qc

ls ~/bs_seq/meth_qc

# Open  ~/bs_seq/meth_qc/DC552_NI48_chr16_trimmed_bismark_bt2_SE_report.html
# in browser

# MultiQC
cd ~/bs_seq
multiqc -f -o qc -n meth bams bams_qc
# Open  ~/bs_seq/qc/meth.html
# in browser

# Final qc:
multiqc -f -o qc -n clean_final bams reads_trim reads_qc/fastp
# Open  ~/bs_seq/qc/clean_final.html
# in browser

# ------------ [15] -------------------------------------
# Methylome visualization

cd ~/bs_seq

# Download
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
# or copy from:
# cp /mnt/meth/genomes/hg38.chrom.sizes .

mv hg38.chrom.sizes indexes

# Create dir for BigWigs
mkdir meth_CpG_bw

# Convert "*.bismark.cov.gz to *.bw:
#   1. sort each *.bismark.cov.gz file:
for F in $(find meth_CpG -name "*.bismark.cov.gz"); do FN=$(basename $F); echo $FN; zcat $F | sort -k1,1 -k2,2n > meth_CpG_bw/${FN/.cov.gz/.sorted.cov}; done

#   2. Make BedGraph with MLevel info:
for F in $(find meth_CpG_bw -name "*.bismark.sorted.cov"); do echo $F; cat $F | awk '{{ print $1, $2-1, $3, $4 }}' > ${F/.sorted.cov/.sorted.mlevel.bedGraph}; done

#     BedGraph with MLevel info -> BigWig:
for F in $(find meth_CpG_bw -name "*.mlevel.bedGraph"); do echo $F; bedGraphToBigWig $F indexes/hg38.chrom.sizes ${F/.mlevel.bedGraph/.mlevel.bw}; done

#   3. Make BedGraph with C coverage info:
for F in $(find meth_CpG_bw -name "*.bismark.sorted.cov"); do echo $F; cat $F | awk '{{ print $1, $2-1, $3, $5+$6 }}' > ${F/.sorted.cov/.sorted.ccov.bedGraph}; done

#     BedGraph with C coverage info -> BigWig:
for F in $(find meth_CpG_bw -name "*.ccov.bedGraph"); do echo $F; bedGraphToBigWig $F indexes/hg38.chrom.sizes ${F/.ccov.bedGraph/.ccov.bw}; done

#   4. Cleanup
rm meth_CpG_bw/*.bedGraph
rm meth_CpG_bw/*.cov

# Check results:
tree -h meth_CpG_bw/

# Compare BW file size with BAM
ls -lah bams/*deduplicated*.bam

# Download `meth_CpG_bw` and open in JBR
# * navigate to: chr16:67648127-67745743
# * reorder tracks
#   Move Up/Down: Alt + Up/Down

# ------------ [16] -------------------------------------
# Combine +/- strand reads counts for same CpG

cd ~/bs_seq

# Single Sample
#coverage2cytosine --genome_folder indexes --merge_CpG --gzip --dir meth_CpG -o DC552_TB48_chr16 meth_CpG/DC552_TB48_chr16.bismark.cov.gz

# All Samples
for F in $(find meth_CpG -name "*.bismark.cov.gz"); do FN=$(basename $F); echo $FN; coverage2cytosine --genome_folder indexes --merge_CpG --gzip --dir meth_CpG -o ${FN/.bismark.cov.gz/} $F > meth_CpG/${FN/.bismark.cov.gz/_report.log} 2>&1; done

# Results
# Check files for 1 sample:
find meth_CpG -name "*DC552_NI48*"
zcat meth_CpG/DC552_NI48_chr16.CpG_report.txt.gz | awk '($2 >= 10400) && ($2 < 10450)'
zcat meth_CpG/DC552_NI48_chr16.CpG_report.merged_CpG_evidence.cov.gz | head

# ------------ [17] -------------------------------------
# Convert Bismark Methylome to Methpipe Tool Format
cd ~/bs_seq

# Single Sample
#zcat meth_CpG/DC5_TB48_chr16.CpG_report.merged_CpG_evidence.cov.gz  | awk -v OFS='\t' '{ print $1, $2, "+", "CpG", $4/100, $5+$6 }' > meth_CpG/DC552_TB48_chr16.merged_CpG.meth

# All Samples
for F in $(find meth_CpG -name "*.CpG_report.merged_CpG_evidence.cov.gz"); do echo $F; zcat $F | awk -v OFS='\t' '{ print $1, $2, "+", "CpG", $4/100, $5+$6 }' > ${F/.CpG_report.merged_CpG_evidence.cov.gz/.merged_CpG.meth}; done

# Check files:
cat meth_CpG/DC5_TB48_chr16.merged_CpG.meth | awk '($2 >= 10400) && ($2 < 10450)'


# ------------ [18] -------------------------------------
# Hyper/Hypo/Partially methylated Regions, partially methylated domains
cd ~/bs_seq

# Hypo methylated regions
mkdir meth_CpG_hmr
# Single sample
#/opt/methpipe-4.1.1/bin/hmr -o /mnt/meth/bs_seq/meth_CpG_hmr/DC5_TB48_chr16.hmr /mnt/meth/bs_seq/meth_CpG/DC5_TB48_chr16.merged_CpG.meth

# All Samples:
for F in $(find meth_CpG -name "*.meth"); do FN=$(basename $F); echo $FN; /opt/methpipe-4.1.1/bin/hmr -o meth_CpG_hmr/${FN/.meth/.hmr} $F; done

# Hyper methylated regions
mkdir meth_CpG_hypermr
# Single sample
#/opt/methpipe-4.1.1/bin/hypermr -o /mnt/meth/bs_seq/meth_CpG_hypermr/DC5_TB48_chr16.hypermr /mnt/meth/bs_seq/meth_CpG/DC5_TB48_chr16.merged_CpG.meth

# All Samples:
for F in $(find meth_CpG -name "*.meth"); do FN=$(basename $F); echo $FN; /opt/methpipe-4.1.1/bin/hypermr -o meth_CpG_hypermr/${FN/.meth/.hypermr} $F; done

# Partially methylated regions
mkdir meth_CpG_pmr
# Single sample
#/opt/methpipe-4.1.1/bin/hmr -partial -o /mnt/meth/bs_seq/meth_CpG_pmr/DC5_TB48_chr16.hypermr /mnt/meth/bs_seq/meth_CpG/DC5_TB48_chr16.merged_CpG.meth

# All Samples:
for F in $(find meth_CpG -name "*.meth"); do FN=$(basename $F); echo $FN; /opt/methpipe-4.1.1/bin/hmr -partial -o meth_CpG_pmr/${FN/.meth/.pmr} $F; done

# Partially methylated domains
mkdir meth_CpG_pmd
# Single sample
#/opt/methpipe-4.1.1/bin/pmd -partial -o /mnt/meth/bs_seq/meth_CpG_pmd/DC5_TB48_chr16.pmd /mnt/meth/bs_seq/meth_CpG/DC5_TB48_chr16.merged_CpG.meth

# All Samples:
for F in $(find meth_CpG -name "*.meth"); do  FN=$(basename $F); echo $FN; /opt/methpipe-4.1.1/bin/pmd -o meth_CpG_pmd/${FN/.meth/.pmd} $F; done

# Check in JBR
# E.g. check `chr16:68304582-68317062`

# Signal profile
cd ~/bs_seq
mkdir -p meth_CpG_profiles
# 
conda activate bio
 
# Merge all hypo methylated:
 cat  meth_CpG_hmr/DC555_HI48.merged_CpG.hmr meth_CpG_hmr/DC555_NI48.merged_CpG.hmr \
     meth_CpG_hmr/DC552_HI48.merged_CpG.hmr meth_CpG_hmr/DC552_NI48.merged_CpG.hmr \
     | sort -k1,1 -k2,2n > meth_CpG_profiles/all.hmr
 bedtools merge -i meth_CpG_profiles/all.hmr > meth_CpG_profiles/all_merged.hmr

computeMatrix scale-regions -S \
    meth_CpG_bw/DC552_HI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC552_NI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC555_HI48.bismark.sorted.mlevel.bw  \
    meth_CpG_bw/DC555_NI48.bismark.sorted.mlevel.bw \
    -R meth_CpG_profiles/all_merged.hmr \
    -a 3000 -b 3000 -out meth_CpG_profiles/hmr_matrix.mat.gz

#for dmr
computeMatrix scale-regions -S \
    meth_CpG_bw/DC552_HI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC552_NI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC555_HI48.bismark.sorted.mlevel.bw  \
    meth_CpG_bw/DC555_NI48.bismark.sorted.mlevel.bw \
    -R meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.down.bed \
    -a 3000 -b 3000 -out meth_CpG_profiles/hmr_matrix.mat.gz


plotProfile -m meth_CpG_profiles/hmr_matrix.mat.gz --outFileName meth_CpG_profiles/hmr_dmr.png \
   --plotTitle "Methylation profile"
 
plotProfile -m meth_CpG_profiles/hmr_matrix.mat.gz --outFileName meth_CpG_profiles/hmr_dmr2.png \
   --perGroup --kmeans 2 --plotTitle "Methylation profile"

plotProfile -m meth_CpG_profiles/hmr_matrix.mat.gz --outFileName meth_CpG_profiles/hmr_dmr_heatmap.png \
  --perGroup --kmeans 2 --plotType heatmap --yMin 0 --yMax 50 \
  --plotTitle "Methylation profile"

# Merge all hyper methylated:
cat  meth_CpG_hypermr/DC552_HI48.merged_CpG.hypermr meth_CpG_hypermr/DC552_NI48.merged_CpG.hypermr \
     meth_CpG_hypermr/DC555_HI48.merged_CpG.hypermr meth_CpG_hypermr/DC555_NI48.merged_CpG.hypermr \
     | sort -k1,1 -k2,2n > meth_CpG_profiles/all.hypermr
 bedtools merge -i meth_CpG_profiles/all.hypermr > meth_CpG_profiles/all_merged.hypermr

computeMatrix scale-regions -S \
    meth_CpG_bw/DC552_HI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC552_NI48.bismark.sorted.mlevel.bw \
    meth_CpG_bw/DC555_HI48.bismark.sorted.mlevel.bw  \
    meth_CpG_bw/DC555_NI48.bismark.sorted.mlevel.bw \
  -R meth_CpG_profiles/all_merged.hypermr \
  -a 3000 -b 3000 -out meth_CpG_profiles/hypermr_matrix.mat.gz

plotProfile -m meth_CpG_profiles/hypermr_matrix.mat.gz --outFileName meth_CpG_profiles/hypermr.png \
   --plotTitle "Methylation profile across hypermethylated regons"
   
plotProfile -m meth_CpG_profiles/hypermr_matrix.mat.gz --outFileName meth_CpG_profiles/hypermr2.png \
   --perGroup --kmeans 2 --plotTitle "Methylation profile across hypermethylated regons"
 
plotProfile -m meth_CpG_profiles/hypermr_matrix.mat.gz --outFileName meth_CpG_profiles/hypermr_heatmap.png \
   --perGroup --kmeans 2 --plotType heatmap --yMin 0 --yMax 50 \
   --plotTitle "Methylation profile across hypermethylated regons"

conda activate meth

# ------------ [19] -------------------------------------
# DMRs using Radmeth
cd ~/bs_seq

# Prepare input files: desing matrix, proportions table:
mkdir meth_CpG_dmrs

/opt/methpipe-4.1.1/bin/merge-methcounts -t meth_CpG/DC552_HI48.merged_CpG.meth meth_CpG/DC555_HI48.merged_CpG.meth meth_CpG/DC552_NI48.merged_CpG.meth meth_CpG/DC555_NI48.merged_CpG.meth > meth_CpG_dmrs/proportion_table_all.txt

# leave only CpG covered >= 5 reads in all samples (according to paper)
cat meth_CpG_dmrs/proportion_table_all.txt | awk '( NR == 1) || ($2 >= 5 && $4 >=5 && $6 >= 5 && $8 >= 5)' > meth_CpG_dmrs/proportion_table_all_5_cov.txt

# Design matrix
cat ~/bs_seq/design_matrix.tsv

# Radmeth

# regression
/opt/methpipe-4.1.1/bin/radmeth regression -v -factor infected -o meth_CpG_dmrs/dmc.bed /home/student/bs_seq/design_matrix.tsv meth_CpG_dmrs/proportion_table_all_5_cov.txt

# adjust p-vals
/opt/methpipe-4.1.1/bin/radmeth adjust -bins 1:200:1 -o meth_CpG_dmrs/dmc_adjusted.bed meth_CpG_dmrs/dmc.bed

# merge into dmrs
/opt/methpipe-4.1.1/bin/radmeth merge -p 0.05 meth_CpG_dmrs/dmc_adjusted.bed | awk -v OFS='\t' '{{ print $1,$2,$3,$4 NR, $5, $6 }}' > meth_CpG_dmrs/dmrs.bed

# filter dmrs
# vim /mnt/meth/scripts/split_dmrs.sh
cat /mnt/meth/scripts/split_dmrs.sh

/mnt/meth/scripts/split_dmrs.sh meth_CpG_dmrs/dmrs.bed meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3 0.025 3

# Check results:
find meth_CpG_dmrs -name "*.bed" | xargs -I {} wc -l {}

# Filter dmrs sets:
ls meth_CpG_dmrs/dmrs*.bed
# archive files:
tar -zcvf dmrs.tar.gz meth_CpG_dmrs/dmrs*.bed

