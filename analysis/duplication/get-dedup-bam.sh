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
  REF=/data/wraycompute/malaria/reference/vivax/PVP01.fa
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
#SBATCH --output=dedup-bam%A-%a.out    # Standard output and error log
#SBATCH --mem=60G
#SBATCH -c16
#SBATCH --job-name=dedup-bam
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu


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
PICARD_JAR=/data/wraycompute/malaria/Applications/picard/picard.jar
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

  echo bwa mem -M -t 32 "${REF}" \${R1} \${R2} -R "\${RG}" \| samtools sort - -o "\${BAM}"
  bwa mem -M -t 32 "${REF}" \${R1} \${R2} -R "\${RG}"| samtools sort - -o "\${BAM}"
  echo samtools index \${BAM}
  samtools index \${BAM}
  echo "BAM constructed"
  touch bam-constructed
fi

if [ ! -e bam-deduped ] ; then

  #3. Mark Duplicates
  module load jdk/1.8.0_45-fasrc01
  PICARD_JAR=/data/wraycompute/malaria/Applications/picard/picard.jar
  RUN_PICARD="java -jar -Xmx7g \${PICARD_JAR}"

  echo "Running MarkDuplicates..."
  echo \${RUN_PICARD} MarkDuplicates INPUT=\${BAM} OUTPUT=\${DEDUP_BAM} METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT
  \${RUN_PICARD} MarkDuplicates INPUT=\${BAM} OUTPUT=\${DEDUP_BAM} METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT
  echo "Done with MarkDuplicates."

  echo samtools index \${DEDUP_BAM}
  samtools index \${DEDUP_BAM}
  echo "Deduped BAM indexed"
  touch bam-deduped
fi

#remove these large files!
rm -rf \${ACC}.bam \{ACC}.bam.bai
rm -rf \${ACC}_1.fastq \${ACC}_2.fastq

touch done
echo "finished!"
EOF
