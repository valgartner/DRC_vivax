#!/bin/bash

set -e

if [ "$1" = '--dry-run' ] ; then
  dryrun=true
  shift
fi

if [ "$1" = '-h' ] ; then
  echo "Usage:"
  echo "   make-gvcf.sh <accessions-file> <ref>.fasta"
  echo "   make-gvcf.sh <accessions-file> (uses PVP01.fa)"
  exit 0
fi

if [ "$#" -lt 1 ] ; then
  echo "map.sh: only got $# arguments!  Need 1 arguments."
  echo "Usage:"
  echo "   make-gvcf.sh <accessions-file> <ref>.fasta"
  echo "   make-gvcf.sh <accessions-file> (uses PVP01.fa)"
  exit 1
fi

if [ "$#" -gt 2 ] ; then
  echo "map.sh: got $# arguments!"
  echo "Usage:"
  echo "   make-gvcf.sh <accessions-file> <ref>.fasta"
  echo "   make-gvcf.sh <accessions-file> (uses PVP01.fa)"
  exit 1
fi

AFILE="$1"

NACC=$(cat "${AFILE}" | wc -l)

REF="$2"
if [ -z "$REF" ] ; then
  REF=/path/to/PVP01.fa
fi

>&2 echo "$AFILE $REF"

EXE=sbatch
if [ "${dryrun}" = true ] ; then
  EXE=cat
fi

${EXE} << EOF
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-${NACC}%10
#SBATCH --output=make-gvcf%A-%a.out    # Standard output and error log
#SBATCH --mem=60G
#SBATCH -c16
#SBATCH --job-name=make-gvcf.sh


set -e
echo "parameters: $AFILE $REF"
ACC=\$(cat "${AFILE}" | sed -n \${SLURM_ARRAY_TASK_ID}p)
echo "task \${SLURM_ARRAY_TASK_ID} accession \${ACC}"
mkdir -p "\${ACC}"
cd \${ACC}
pwd

exec > slurm-\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out 2>&1

module load samtools
module load bwa
module load jdk/1.8.0_45-fasrc01
PICARD_JAR=/path/to/picard.jar
RUN_PICARD="java -jar -Xmx7g \${PICARD_JAR}"

#1. Get the FASTQ
R1=\${ACC}_1.fastq
R2=\${ACC}_2.fastq
R3=\${ACC}.fastq

if [ -e done ] ; then
  echo "Already done."
  exit 0
fi

if [ ! -e fastq-downloaded ] && [ ! -e bam-constructed ] && [ ! -e bam-deduped ] ; then
  echo command: fasterq-dump "\$ACC"
  fasterq-dump "\$ACC"
  if [ -e "\${R1}" ] ; then
     touch fastq-downloaded
     echo "done."
  else
     exit 1
  fi
fi

#2. Construct the BAM
BAM=\${ACC}.bam
DEDUP_BAM=\${ACC}.dedup.bam

if [ ! -e ${REF}.bwt ] ; then
  echo bwa index ${REF}
  bwa index ${REF}
  echo done
fi

if [ ! -e bam-constructed ] && [ ! -e bam-deduped ] ; then
  echo Making BAM...

  if [ ! -e "\${R1}" ] ; then 
    echo "Download finished, but can't find "\${R1}"?!  Quitting."
    exit 1
  fi

  if [ ! -e "\${R2}" ] ; then 
    echo "Download finished, but can't find "\${R2}"?!  Quitting."
    exit 1
  fi

  # Add a read group so other software doesn't crash.
  # The attributes here are entirely imaginary.
  RG="@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330"

  echo bwa mem -M -t 32 "${REF}" \${R1} \${R2} -R "\${RG}" \| samtools sort -@32 - -o "\${BAM}"
  bwa mem -M -t 32 "${REF}" \${R1} \${R2} -R "\${RG}"| samtools sort -@32 - -o "\${BAM}"
  echo samtools index \${BAM} -@32
  samtools index \${BAM} -@32
  echo "BAM constructed"
  touch bam-constructed
fi

if [ ! -e bam-deduped ] ; then

  #3. Mark Duplicates
  module load jdk/1.8.0_45-fasrc01
  PICARD_JAR=/path/to/picard.jar
  RUN_PICARD="java -jar -Xmx7g \${PICARD_JAR}"

  echo "Running MarkDuplicates..."
  echo \${RUN_PICARD} MarkDuplicates INPUT=\${BAM} OUTPUT=\${DEDUP_BAM} METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT
  \${RUN_PICARD} MarkDuplicates INPUT=\${BAM} OUTPUT=\${DEDUP_BAM} METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT
  echo "Done with MarkDuplicates."

  echo samtools index \${DEDUP_BAM} -@32
  samtools index \${DEDUP_BAM} -@32
  echo "Deduped BAM indexed"
  touch bam-deduped
fi

#remove these large files!
rm -rf \${ACC}.bam \{ACC}.bam.bai
rm -rf \${ACC}_1.fastq \${ACC}_2.fastq

#4. Get GVCF
module load jdk/1.8.0_45-fasrc01 #load java module again just in case this script is being rerun after the dedup bam is made
GATK3_JAR=/path/to/GenomeAnalysisTK.jar
RUN_GATK3="java -jar \${GATK3_JAR}"

# initial variant processing, using 4 cores. Note that we call full variants here, not the GVCF (which is what we want as the final output)
\${RUN_GATK3} -T HaplotypeCaller -R ${REF} -I \${DEDUP_BAM} -o \${ACC}-raw_variants.vcf -nct 4

# we split up the variant types since SNPs and indels need different processing
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants.vcf -selectType SNP -o \${ACC}-raw_snps.vcf
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants.vcf -selectType INDEL -o \${ACC}-raw_indels.vcf
rm -rf \${ACC}-raw_variants.vcf

#basic filters
# SNPs
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o \${ACC}-filtered_snps.vcf
# Indels
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o \${ACC}-filtered_indels.vcf


# Here is our first recalibration step. We use the hard filtered SNPs set as the training set. It's not an issue that this set isn't totally accurate because we're going to do more rounds of processing later.
\${RUN_GATK3} -T BaseRecalibrator -R ${REF} -I \${ACC}.dedup.bam -knownSites \${ACC}-filtered_snps.vcf -knownSites \${ACC}-filtered_indels.vcf -o \${ACC}-recal_data.table
\${RUN_GATK3} -T PrintReads -R ${REF} -I \${ACC}.dedup.bam -BQSR \${ACC}-recal_data.table -o \${ACC}-recal_reads.bam

# call raw variants round 2
# we call the second set here using the recalculated reads (variants2)
\${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${REF} -I \${ACC}-recal_reads.bam -o \${ACC}-raw_variants2.vcf
#cleanup
rm -rf *.idx \${ACC}-*.idx \${ACC}-raw_snps.vf \${ACC}-raw_indels.vcf \${ACC}-recal_reads.bam \${ACC}-recal_data.table \${ACC}-post_recal_data.table \${ACC}-filtered_snps.vcf \${ACC}-filtered_indels.vcf
# SNPS
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants2.vcf -selectType SNP -o \${ACC}-raw_snps.vcf
# Indels
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants2.vcf -selectType INDEL -o \${ACC}-raw_indels.vcf
rm -rf \${ACC}-raw_variants2.vcf

# filter round 2
# SNPs
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o \${ACC}-filtered_snps.vcf
# Indels
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o \${ACC}-filtered_indels.vcf

# really important to note that we use the dedup_reads file here as the input, not the printed reads from the last step. I have not tested what would happen if we used the recalculated reads here instead, but none of the guides I read seemed to do so.
\${RUN_GATK3} -T BaseRecalibrator -R ${REF} -I \${ACC}.dedup.bam  -knownSites \${ACC}-filtered_snps.vcf -knownSites \${ACC}-filtered_indels.vcf -o \${ACC}-recal_data.table
\${RUN_GATK3} -T PrintReads -R ${REF} -I \${ACC}.dedup.bam -BQSR \${ACC}-recal_data.table -o \${ACC}-recal_reads.bam


# call raw variants round 3
\${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${REF} -I \${ACC}-recal_reads.bam -o \${ACC}-raw_variants3.vcf
# SNPs
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants3.vcf -selectType SNP -o \${ACC}-raw_snps.vcf
# Indels
\${RUN_GATK3} -T SelectVariants -R ${REF} -V \${ACC}-raw_variants3.vcf -selectType INDEL -o \${ACC}-raw_indels.vcf
rm -rf \${ACC}-raw_variants3.vcf

# filter round 3
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o \${ACC}-filtered_snps.vcf
\${RUN_GATK3} -T VariantFiltration -R ${REF} -V \${ACC}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o \${ACC}-filtered_indels.vcf

# recalibrate round 3 
# really important to note that we use the dedup_reads file here as the input, not the printed reads from the last step. I have not tested what would happen if we used the recalculated reads here instead, but none of the guides I read seemed to do so.
\${RUN_GATK3} -T BaseRecalibrator -R ${REF} -I \${ACC}.dedup.bam -knownSites \${ACC}-filtered_snps.vcf -knownSites \${ACC}-filtered_indels.vcf -o \${ACC}-recal_data.table
\${RUN_GATK3} -T PrintReads -R ${REF} -I \${ACC}.dedup.bam -BQSR \${ACC}-recal_data.table -o \${ACC}-recal_reads.bam

# GVCF output step
\${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${REF} -I \${ACC}-recal_reads.bam --emitRefConfidence GVCF -ploidy 1 -o \${ACC}-haploid.g.vcf


#5. Clean up 

if [ -e "\${R1}" ] ; then
  echo -n "Removing reads \${R1} \${R2} \${R3}... "
  rm -f "\${R1}"
  rm -f "\${R2}"
  rm -f "\${R3}"
  echo done
fi
rm -f fastq-downloaded


if [ -e "\${BAM}" ] ; then
  echo -n "Removing \${BAM} ... "
  rm -f "\${BAM}"
  echo "done"
fi
rm -f "\${BAM}".bai
rm -f bam-constructed

if [ -e "\${DEDUP_BAM}" ] ; then
  echo -n "Removing \${DEDUP_BAM} ... "
  rm -f "\${DEDUP_BAM}"
  echo "done"
fi


rm -f ${ACC}.dedup.bam \${ACC}.dedup.bam.bai
rm -rf \${ACC}-raw_snps.vcf \${ACC}-raw_indels.vcf \${ACC}-recal_reads.bam \${ACC}-recal_data.table 
rm -rf \${ACC}-post_recal_data.table \${ACC}-filtered_snps.vcf \${ACC}-filtered_indels.vcf
rm -rf \${ACC}-filtered_indels.vcf.idx \${ACC}-filtered_snps.vcf.idx \${ACC}-raw_indels.vcf.idx 
rm -rf \${ACC}-raw_snps.vcf.idx \${ACC}-raw_variants3.vcf.idx \${ACC}-raw_variants2.vcf.idx metrics.txt \${ACC}-recal_reads.bai
rm -rf bam-deduped
echo "clean up large intermediate files: DONE!"

touch done
echo "finished!"
EOF
