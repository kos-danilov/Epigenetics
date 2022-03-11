#!/bin/bash
set -euxo pipefail



echo "Activate bio"
source activate bio


# Put data to home directory
mkdir -p ~/chipseq
cd ~/chipseq

echo "Prepare data in fastq"
mkdir -p fastq
cd fastq

wget http://artyomovlab.wustl.edu/publications/supp_materials/4Oleg/2021_ITMO_epigenetics_practice/chipseq/14.fastq.gz
wget http://artyomovlab.wustl.edu/publications/supp_materials/4Oleg/2021_ITMO_epigenetics_practice/chipseq/15.fastq.gz

gunzip 14.fastq.gz
gunzip 15.fastq.gz

########################################################################################################################
echo "1. QC with fastq"
mkdir -p ~/chipseq/fastqc

for FILE in $(find ~/chipseq/fastq -name '*.f*q*')
do :
  FILE_NAME=${FILE##*/}
  NAME=${FILE_NAME%%.f*q} # file name without extension
  fastqc --outdir ~/chipseq/fastqc "${FILE}" 2>&1 | tee ${NAME}_fastqc.log
done

mv ~/chipseq/fastq/*_fastqc.log ~/chipseq/fastqc/


########################################################################################################################
echo "2. Reads QC report with multiqc"
multiqc -f -o ~/chipseq/fastqc ~/chipseq/fastqc


########################################################################################################################
echo "3. Alignment with bowtie"
GENOME=hg19
export BOWTIE_INDEXES=/mnt/chipseq/index/hg19/  # Configure index folder

for FILE in $(find ~/chipseq/fastq -name '*.f*q')
do :
  NAME=$(basename ${FILE%%.fastq})  # File name without extension
  ID=${NAME}_${GENOME}
  BAM_NAME="${ID}.bam"
  if [[ ! -f "${BAM_NAME}" ]]; then
    # Bowtie command line options used
    # -p/--threads <int> number of alignment threads to launch (default: 1)
    # -S/--sam           write hits in SAM format
    # -t/--time          print wall-clock time taken by search phases
    # -m <int>           suppress all alignments if > <int> exist (def: no limit)
    # -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
    # --best             hits guaranteed best stratum; ties broken by quality
    # --strata           hits in sub-optimal strata aren't reported (requires --best)
    ### Bowtie authors TIME MACHINE!
    # -x <ebwt>
    # The basename of the Bowtie, or Bowtie 2, index to be searched.
    # If a Bowtie and Bowtie 2 index are located in the same directory and share the same basename,
    # bowtie will use the Bowtie 2 index.
    bowtie -p 6 -St -m 1 -v 3 --best --strata ${GENOME} ${FILE} ${ID}.sam 2>&1 |\
      tee ${NAME}_bowtie.log
    samtools view -bS ${ID}.sam -o ${ID}_not_sorted.bam
    samtools sort ${ID}_not_sorted.bam -o ${BAM_NAME}
    # Cleanup
    rm ${ID}.sam ${ID}_not_sorted.bam
  fi
done

mkdir -p ~/chipseq/bam
mv ~/chipseq/fastq/*.bam ~/chipseq/bam/
mv ~/chipseq/fastq/*_bowtie.log ~/chipseq/bam/

########################################################################################################################
echo "4. Alignment QC with multiqc"
multiqc -f -o ~/chipseq/bams_qc ~/chipseq/bam/*_bowtie.log


########################################################################################################################
echo "5. Visualization with deeptools"
cd ~/chipseq

for FILE in $(find ~/chipseq/bam -name '*.bam')
do :
    NAME=$(basename ${FILE%%.bam})  # File name without extension
    BW=${NAME}.bw
    if [[ ! -f ${BW} ]]; then # Check that file is not created yet
        samtools index ${FILE}
        bamCoverage --bam ${FILE} -o ${BW} 2>&1 | tee ${NAME}_bw.log
    fi
done

mkdir ~/chipseq/bw
mv ~/chipseq/*.bw ~/chipseq/bw
mv ~/chipseq/*_bw.log ~/chipseq/bw



########################################################################################################################
echo "6. Prepare pileup coverage files for SICER"
for FILE in $(find ~/chipseq/bam -name '*.bam')
do :
  PILEUP=${FILE/.bam/_pileup.bed}
  if [[ ! -f ${PILEUP} ]]; then  # Not created yet
    bedtools bamtobed -i ${FILE} > ${PILEUP}
  fi
done


########################################################################################################################
echo "7. Peak calling with SICER"
WINDOW_SIZE=200
FRAGMENT_SIZE=150
GAP_SIZE=600
FDR=0.05
 # From MACS2 documentation:
 # The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
# Here are all precompiled parameters for effective genome size:
# hs: 2.7e9
# mm: 1.87e9
# ce: 9e7
# dm: 1.2e8"""
# See hg19 genome size at: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
# Effective genome fraction = (size / genome_length) * (1.0 * chromosomes_length / genome_length)
# (2.7e9 / 3,101,788,170) * (102,531,392 / 3,101,788,170) = 0.02877374213
# EFFECTIVE_GENOME_FRACTION=0.02877374213

# changes to effective genome fraction
EFFECTIVE_GENOME_FRACTION=0.8046990572332249

mv ~/chipseq/bam/15_hg19_pileup.bed ~/chipseq/bam/15_hg19_pileup_input.bed
INPUT_PILEUP_BED=~/chipseq/bam/15_hg19_pileup_input.bed

FILE=$(find ~/chipseq/bam -name '*_pileup.bed' | grep -v 'input')

for FILE in $(find ~/chipseq/bam -name '*_pileup.bed' | grep -v 'input')
do :
    NAME=$(basename ${FILE%%_pileup.bed}) # file name without extension
    ID="${NAME}-W${WINDOW_SIZE}-G${GAP_SIZE}-FDR${FDR}"

    # Naming example: OD_OD10_H3K27me3-W200-G0-FDR0.01-island.bed
    if [[ ! -f "${ID}-island.bed" ]]; then  # Not created yet
        FILE_BED=${NAME}.bed # It is used for results naming
        INPUT_BED=$(basename ${INPUT_PILEUP_BED/_pileup.bed/.bed})  # file name without extension

        # Create working folder
        SICER_FOLDER=~/chipseq/$(mktemp -d sicer.XXXXXX)
        SICER_OUT_FOLDER=${SICER_FOLDER}/out
        mkdir -p ${SICER_OUT_FOLDER}
        ln -sf ${FILE} ${SICER_FOLDER}/${FILE_BED}
        ln -sf ${INPUT_PILEUP_BED} ${SICER_FOLDER}/${INPUT_BED}

        cd ${SICER_FOLDER}
        # Usage:
        # SICER.sh    ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] \
        #   ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["FDR"]
        # SICER.sh [InputDir] [bed file] [control file] [OutputDir] [Species]
        # [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [FDR]
        # SICER-rb.sh ["InputDir"] ["bed file"]                  ["OutputDir"] ["species"] ["redundancy threshold"] \
        #   ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["E-value"]
        #
        # Defaults:
        #   redundancy threshold    = 1
        #   window size (bp)        = 200
        #   fragment size           = 150
        #   gap size (bp)           = 600
        SICER.sh ${SICER_FOLDER} ${FILE_BED} ${INPUT_BED} ${SICER_OUT_FOLDER} ${GENOME} 1 \
          ${WINDOW_SIZE} ${FRAGMENT_SIZE} ${EFFECTIVE_GENOME_FRACTION} ${GAP_SIZE} ${FDR} 2>&1 |\
          tee ~/chipseq/${ID}_sicer.log

        mv ${SICER_OUT_FOLDER}/* ~/chipseq/

        # Move back out of folder to be removed
        cd ~/chipseq/

        # Cleanup everything else
        rm -r ${SICER_FOLDER}
    fi
done

# Move results to sicer folder
mkdir ~/chipseq/sicer
ls ~/chipseq/* | grep "1-removed" | xargs -I {} mv {} ~/chipseq/sicer
ls ~/chipseq/* | grep "W${WINDOW_SIZE}" | xargs -I {} mv {} ~/chipseq/sicer


########################################################################################################################
echo "8. Peak calling with SPAN"
BIN=200
FDR=0.05
GAP=5

mv /home/student/chipseq/bam/15_hg19.bam /home/student/chipseq/bam/15_hg19_input.bam
INPUT=/home/student/chipseq/bam/15_hg19_input.bam


for FILE in $(find ~/chipseq/bam -name '*.bam' | sed 's#\./##g' | grep -v 'input')
do :
  $FILE
done
do :
    NAME=$(basename ${FILE%%.bam})  # File name without extension
    ID=${NAME}_${FDR}_${GAP}
    if [[ ! -f ${ID}.peak ]]; then # Not created yet
      java -Xmx16G -jar /mnt/chipseq/span-0.13.5244.jar analyze -t ${FILE} -c ${INPUT} \
        --chrom.sizes /mnt/chipseq/hg19.chrom.sizes \
        --bin ${BIN} --fdr ${FDR} --gap ${GAP} \
        --peaks ${ID}.peak \
        --threads 6 2>&1 | tee ${NAME}_span.log
    fi
done

# Move results to span folder
mkdir ~/chipseq/span
mv ~/chipseq/cache ~/chipseq/span/
mv ~/chipseq/logs ~/chipseq/span/
mv ~/chipseq/fit ~/chipseq/span/
mv ~/chipseq/*.peak ~/chipseq/span/
mv ~/chipseq/*span.log ~/chipseq/span/

# Command to get all SPAN models
tree ~/chipseq/span/ | grep "\.span"


########################################################################################################################
echo "9. Results exploration"

# Get most confident peaks
cat ~/chipseq/span/14_hg19_0.05_5.peak | sort -k9,9r

# Find out the bed file with biggest number of lines
ls ~/chipseq/sicer/*island.bed | xargs wc -l | grep -v total | sort -k1,1nr

# Sort BED file by chromosome and by start position
cat ~/chipseq/sicer/14_hg19-W200-G600-FDR0.05-island.bed | sort -k1,1 -k2,2n > sorted.bed

########################################################################################################################

echo "10. Motif analysis"

# Take 500 most significant peaks and take centres of them
sort -k 9,9nr ~/chipseq/span/14_hg19_0.05_5.peak | head -n 500 | sort -k1,1 -k2,2n |\
  awk -v OFS='\t' '{print($1, int(($3+$2)/2)-100, int(($3+$2)/2)+100)}' >\
  14_hg19_0.05_5_top500.peak

# Unzip fa file
gunzip -c /mnt/chipseq/hg19.fa.gz > ~/chipseq/hg19.fa

# Compute sequence for peaks
bedtools getfasta -fi ~/chipseq/hg19.fa -bed ~/chipseq/14_hg19_0.05_5_top500.peak >\
  ~/chipseq/14_hg19_0.05_5_top500.peak.fa

# Launch Homer
perl /opt/conda/envs/bio/share/homer/bin/findMotifsGenome.pl \
  ~/chipseq/14_hg19_0.05_5_top500.peak hg19 \
  ~/chipseq/14_hg19_0.05_5_top500.peak.motif -size 200 -mask

########################################################################################################################

echo "17. Prepare BED3 format peaks for ChIPpeakAnno"
cat ~/chipseq/span/14_hg19_0.05_5.peak | awk -v OFS='\t' '{N+=1; print $1,$2,$3,N}' >\
  ~/chipseq/span/14_hg19_0.05_5.peak3

#######################
# Switch to R console #
#######################


