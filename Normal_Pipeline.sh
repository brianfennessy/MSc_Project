#! /bin/bash

#$ -N normal_samples
#$ -q all.q
#$ -cwd

#$ -v PATH
#$ -v LD_LIBRARY_PATH
#$ -v PYTHONPATH
#$ -S /bin/bash

# Load library
# module load fastqc

OUTPUT="fastqc_N_results"
mkdir $OUTPUT

#Set directory with data
BASEDIR=/data4/bfennessy/project/raw_fastq/NORMAL_SAMPLES

for f in *.fastq;
do
        fastqc $f
done

mv $BASEDIR/*.html $OUTPUT
mv $BASEDIR/*.zip $OUTPUT

TRIMMOMATIC='/data4/bfennessy/project/raw_fastq/Trimmomatic-0.38/trimmomatic-0.38.jar'
ADAPTERS='/data4/bfennessy/project/raw_fastq/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa'

# Trim HN51
java -jar $TRIMMOMATIC PE -threads 8 -phred33 -trimlog HN51_N.trimlog \
HN51_S2_nml.BD1LYPACXX.lane_3_P1_I12.hg19.sequence.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P2_I12.hg19.sequence.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P1_I12_PAIRED.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P1_I12_UNPAIRED.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P2_I12_PAIRED.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P2_I12_UNPAIRED.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

# Trim HN60
java -jar $TRIMMOMATIC PE -threads 8 -phred33 -trimlog HN60_N.trimlog \
HN60_s2_normal.BD1LYPACXX.lane_7_P1_I16.hg19.sequence.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P2_I16.hg19.sequence.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P1_I16_PAIRED.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P1_I16_UNPAIRED.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P2_I16_PAIRED.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P2_I16_UNPAIRED.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

# Trim HN72AC
java -jar $TRIMMOMATIC PE -threads 8 -phred33 -trimlog HN72AC_N.trimlog \
HN72s2.AC254KACXX.lane_2_P1_I15.hg19.sequence.fastq \
HN72s2.AC254KACXX.lane_2_P2_I15.hg19.sequence.fastq \
HN72s2.AC254KACXX.lane_2_P1_I15_PAIRED.fastq \
HN72s2.AC254KACXX.lane_2_P1_I15_UNPAIRED.fastq \
HN72s2.AC254KACXX.lane_2_P2_I15_PAIRED.fastq \
HN72s2.AC254KACXX.lane_2_P2_I15_UNPARIED.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

# Trim HN72AH
java -jar $TRIMMOMATIC PE -threads 8 -phred33 -trimlog HN72AH_N.trimlog \
HN72s2.AH0LENADXX.lane_2_P1_I15.hg19.sequence.fastq \
HN72s2.AH0LENADXX.lane_2_P2_I15.hg19.sequence.fastq \
HN72s2.AH0LENADXX.lane_2_P1_I15_PAIRED.fastq \
HN72s2.AH0LENADXX.lane_2_P1_I15_UNPAIRED.fastq \
HN72s2.AH0LENADXX.lane_2_P2_I15_PAIRED.fastq \
HN72s2.AH0LENADXX.lane_2_P2_I15_UNPAIRED.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

# Set path to index and reference genome
INDEX='/data4/bfennessy/project/raw_fastq/hg19/hg19.fa.fai'
REFERENCE='/data4/bfennessy/project/raw_fastq/hg19/hg19.fa'

# Align HN51
bwa mem -M -t 8 -R '@RG\tID:BD1LYPACXX.3\tSM:HN51\tLB:I12' $REFERENCE HN51_S2_nml.BD1LYPACXX.lane_3_P1_I12_PAIRED.fastq \
HN51_S2_nml.BD1LYPACXX.lane_3_P2_I12_PAIRED.fastq | samtools view -Sb -> HN51_N.bam

# Align HN60
bwa mem -M -t 8 -R '@RG\tID:BD1LYPACXX.7\tSM:HN60\tLB:I16' $REFERENCE HN60_s2_normal.BD1LYPACXX.lane_7_P1_I16_PAIRED.fastq \
HN60_s2_normal.BD1LYPACXX.lane_7_P2_I16_PAIRED.fastq | samtools view -Sb -> HN60_N.bam

# Align HN72AC
bwa mem -M -t 8 -R '@RG\tID:AC254KACXX.2\tSM:HN72AC\tLB:I15' $REFERENCE HN72s2.AC254KACXX.lane_2_P1_I15_PAIRED.fastq \
HN72s2.AC254KACXX.lane_2_P2_I15_PAIRED.fastq | samtools view -Sb -> HN72AC_N.bam

# Align HN72AH
bwa mem -M -t 8 -R '@RG\tID:AH0LENADXX.2\tSM:HN72AH\tLB:I15' $REFERENCE HN72s2.AH0LENADXX.lane_2_P1_I15_PAIRED.fastq \
HN72s2.AH0LENADXX.lane_2_P2_I15_PAIRED.fastq | samtools view -Sb -> HN72AH_N.bam

# Set Picard
PICARD='java -Xms10g -Xmx20g -Djava.io.tmpdir=tmpdir -jar /data4/bfennessy/project/raw_fastq/picard.jar'

# Sort
$PICARD SortSam INPUT=HN51_N.bam OUTPUT=HN51_N.sorted.bam SORT_ORDER=coordinate
$PICARD SortSam INPUT=HN60_N.bam OUTPUT=HN60_N.sorted.bam SORT_ORDER=coordinate
$PICARD SortSam INPUT=HN72AC_N.bam OUTPUT=HN72AC_N.sorted.bam SORT_ORDER=coordinate
$PICARD SortSam INPUT=HN72AH_N.bam OUTPUT=HN72AH_N.sorted.bam SORT_ORDER=coordinate

# Mark Duplicates
$PICARD MarkDuplicates INPUT=HN51_N.sorted.bam OUTPUT=HN51_N.sorted.dedup.bam METRICS_FILE=HN51_N.metrics.txt
$PICARD MarkDuplicates INPUT=HN60_N.sorted.bam OUTPUT=HN60_N.sorted.dedup.bam METRICS_FILE=HN60_N.metrics.txt
$PICARD MarkDuplicates INPUT=HN72AC_N.sorted.bam INPUT=HN72AH_N.sorted.bam OUTPUT=HN72_N.sorted.dedup.bam METRICS_FILE=HN72_N.metrics.txt

# Build Index
$PICARD BuildBamIndex INPUT=HN51_N.sorted.dedup.bam
$PICARD BuildBamIndex INPUT=HN60_N.sorted.dedup.bam
$PICARD BuildBamIndex INPUT=HN72_N.sorted.dedup.bam

# Fix Read Group
$PICARD AddOrReplaceReadGroups I=HN51_N.sorted.dedup.bam O=HN51_N.sorted.dedup.RG.bam RGLB=I12 RGPL=illumina RGPU=unit2 RGSM=HN51_normal RGID=BD1LYPACXX.3
$PICARD AddOrReplaceReadGroups I=HN60_N.sorted.dedup.bam O=HN60_N.sorted.dedup.RG.bam RGLB=I16 RGPL=illumina RGPU=unit2 RGSM=HN60_normal RGID=BD1LYPACXX.7
$PICARD AddOrReplaceReadGroups I=HN72_N.sorted.dedup.bam O=HN72_N.sorted.dedup.RG.bam RGLB=I15 RGPL=illumina RGPU=unit2 RGSM=HN72_normal RGID=BD1LYPACXX.2

# Build Index
$PICARD BuildBamIndex INPUT=HN51_N.sorted.dedup.RG.bam
$PICARD BuildBamIndex INPUT=HN60_N.sorted.dedup.RG.bam
$PICARD BuildBamIndex INPUT=HN72_N.sorted.dedup.RG.bam

# Set path to index and reference genome
INDEX='/data4/bfennessy/project/raw_fastq/hg19/hg19.fa.fai'
REFERENCE='/data4/bfennessy/project/raw_fastq/hg19/hg19.fa'

# Provide path to GATK and VCF files
GATK='/data4/pilib/gatk/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar'
DBSNP='/data4/bfennessy/project/raw_fastq/dbsnp_138.hg19.vcf'
MILLS='/data4/bfennessy/project/raw_fastq/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'

# Base Recalibration
java -jar $GATK BaseRecalibrator --TMP_DIR=tmpdir -R $REFERENCE -I HN51_N.sorted.dedup.RG.bam --known-sites $DBSNP --known-sites $MILLS -O HN51_N.recal.table
java -jar $GATK BaseRecalibrator --TMP_DIR=tmpdir -R $REFERENCE -I HN60_N.sorted.dedup.RG.bam --known-sites $DBSNP --known-sites $MILLS -O HN60_N.recal.table
java -jar $GATK BaseRecalibrator --TMP_DIR=tmpdir -R $REFERENCE -I HN72_N.sorted.dedup.RG.bam --known-sites $DBSNP --known-sites $MILLS -O HN72_N.recal.table

# Apply Base Quality Score Recalibration
java -jar $GATK ApplyBQSR -R $REFERENCE -I HN51_N.sorted.dedup.RG.bam --bqsr-recal-file HN51_N.recal.table -O HN51_N.sorted.dedup.bqsr.bam
java -jar $GATK ApplyBQSR -R $REFERENCE -I HN60_N.sorted.dedup.RG.bam --bqsr-recal-file HN60_N.recal.table -O HN60_N.sorted.dedup.bqsr.bam
java -jar $GATK ApplyBQSR -R $REFERENCE -I HN72_N.sorted.dedup.RG.bam --bqsr-recal-file HN72_N.recal.table -O HN72_N.sorted.dedup.bqsr.bam

# Set path to target files
PRIMARY='/data4/bfennessy/project/raw_fastq/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed'
CAPTURE='/data4/bfennessy/project/raw_fastq/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed'
REFERENCE='/data4/bfennessy/project/raw_fastq/hg19/hg19.dict'

# Get interval list format
$PICARD BedToIntervalList I=$PRIMARY O=primary_targets_interval_list SD=$REFERENCE
$PICARD BedToIntervalList I=$CAPTURE O=capture_targets_interval_list SD=$REFERENCE

# Set path to files
PRIMARY='/data4/bfennessy/project/raw_fastq/primary_targets_interval_list'
CAPTURE='/data4/bfennessy/project/raw_fastq/capture_targets_interval_list'
REFERENCE='/data4/bfennessy/project/raw_fastq/hg19/hg19.fa'

# Collect hs metrics
$PICARD CollectHsMetrics I=HN51_N.sorted.dedup.RG.bam O=HN51_N_hs.metrics.txt R=$REFERENCE BAIT_INTERVALS=$CAPTURE TARGET_INTERVALS=$PRIMARY
$PICARD CollectHsMetrics I=HN60_N.sorted.dedup.RG.bam O=HN60_N_hs.metrics.txt R=$REFERENCE BAIT_INTERVALS=$CAPTURE TARGET_INTERVALS=$PRIMARY
$PICARD CollectHsMetrics I=HN72_N.sorted.dedup.RG.bam O=HN72_N_hs.metrics.txt R=$REFERENCE BAIT_INTERVALS=$CAPTURE TARGET_INTERVALS=$PRIMARY

