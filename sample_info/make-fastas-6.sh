#!/bin/bash

set -e

if [ "$1" = '--dry-run' ] ; then
  dryrun=true
  shift
fi

if [ "$1" = '--keep' ] ; then
  keep=true
  shift
fi

if [ "$1" = '-h' ] ; then
  echo "Usage:"
  echo "   make-fastas-6.sh <accessions-file> <ref>.fasta"
  echo "   make-fastas-6.sh <accessions-file> (uses plasmo-combined.fasta)"
  exit 0
fi

if [ "$#" -lt 1 ] ; then
  echo "map.sh: only got $# arguments!  Need 1 arguments."
  echo "Usage:"
  echo "   make-fastas-6.sh <accessions-file> <ref>.fasta"
  echo "   make-fastas-6.sh <accessions-file> (uses plasmo-combined.fasta)"
  exit 1
fi

if [ "$#" -gt 2 ] ; then
  echo "map.sh: got $# arguments!"
  echo "Usage:"
  echo "   make-fastas-6.sh <accessions-file> <ref>.fasta"
  echo "   make-fastas-6.sh <accessions-file> (uses plasmo-combined.fasta)"
  exit 1
fi

AFILE="$1"

NACC=$(cat "${AFILE}" | wc -l)

REF="$2"
if [ -z "$REF" ] ; then
  REF=/data/wraycompute/malaria/reference/plasmo-combined.fasta
fi

>&2 echo "$AFILE $REF"

EXE=sbatch
if [ "${dryrun}" = true ] ; then
  EXE=cat
fi

REGIONS_FILE=/data/wraycompute/malaria/vivax/masks/regions-good.txt

${EXE} << EOF
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --array=1-${NACC}%75
#SBATCH --output=make-fastas-3%A-%a.out    # Standard output and error log
#SBATCH --mem=40G
#SBATCH -c16
#SBATCH --job-name=make-fastas-6.sh

set -e
echo "parameters: $AFILE $REF"
ACC=\$(cat "${AFILE}" | sed -n \${SLURM_ARRAY_TASK_ID}p)
echo "task \${SLURM_ARRAY_TASK_ID} accession \${ACC}"
mkdir -p "\${ACC}"
cd \${ACC}
pwd

exec > slurm-\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out 2>&1

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
BAM=combined.bam
DEDUP_BAM=combined.dedup.bam

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

  echo bwa mem -M -t 16 "${REF}" \${R1} \${R2} -R "\${RG}" \| samtools sort -@16 - -o "\${BAM}"
  bwa mem -M -t 16 "${REF}" \${R1} \${R2} -R "\${RG}"| samtools sort -@16 - -o "\${BAM}"
  echo samtools index \${BAM} -@16
  samtools index \${BAM} -@16
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

  echo samtools index \${DEDUP_BAM} -@16
  samtools index \${DEDUP_BAM} -@16
  echo "Deduped BAM indexed"
  touch bam-deduped
fi


#4. Get FASTA using get-consensus

# Activate the virtual-env
source ~/.virtualenvs/bio/bin/activate

for CHR in LT6356{12,26,27} Pf_M76611 Pf3D7_{14,API}_v3 PocGH01_API_v1 PKNH_API_v2 PKNH_MIT_v2 PmUG01_API_v1 PmUG01_MIT_v1 ; do
  if [ ! -e \${CHR}-1x.fasta ] ; then
    echo get-consensus.py \${DEDUP_BAM} ${REF} \${CHR} --chromosome \${CHR} --prefix \${ACC} --min-coverage 1 2 3 5
    get-consensus.py \${DEDUP_BAM} ${REF} \${CHR} --chromosome \${CHR} --prefix \${ACC} --min-coverage 1 2 3 5
  fi
done


if [ ! -e vcf-written-api ] ; then
  #5. Get VCF using octopus
  octopus -I \${DEDUP_BAM} -R ${REF} -T LT635626 -o api.poly3.vcf.gz --annotations AD -C polyclone  --max-clones 3 --threads 16 --sequence-error-model PCR
  octopus -I \${DEDUP_BAM} -R ${REF} -T LT635626 -o api.g.vcf.gz --annotations AD --refcall POSITIONAL --threads 16 --sequence-error-model PCR

  touch vcf-written-api
fi

if [ ! -e vcf-written-all-chrs ] ; then
  #5. Get VCF using octopus
  octopus -I \${DEDUP_BAM} -R ${REF} -t ${REGIONS_FILE} -o poly2.vcf.gz  --annotations AD -C polyclone  --max-clones 2 --threads 16 --sequence-error-model PCR

  touch vcf-written-all-chrs
fi

#if [ ! -e vcf-written-gatk3 ] ; then
#  gatk3 -Xmx38G -nct 16 -T HaplotypeCaller -R ${REF} -I combined.dedup.bam -ploidy 1 -o gatk3-p1.vcf.gz -ERC GVCF
#  gatk --java-options '-Xmx38G' HaplotypeCallerSpark -R ${REF} -I combined.dedup.bam -ploidy 1 -o gatk3-ploidy1.vcf.gz -ERC GVCF \
#     \$(sed 's/^/-L /g' $REGIONS_FILE)  --spark-master local[16]
#  gatk --java-options '-Xmx38G' HaplotypeCaller -R ${REF} -I combined.dedup.bam -ploidy 1 -O gatk4-ploidy1.vcf.gz -ERC GVCF \
#     \$(sed 's/^/-L /g' $REGIONS_FILE)  --native-pair-hmm-threads 16
#  touch vcf-written-gatk3
#fi


if [ "${keep}" = true ] ; then
  touch done
  echo "finished - not deleting artifacts!"
  exit
fi

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
rm -f "\${DEDUP_BAM}".bai
rm -f bam-deduped

touch done
echo "finished!"
EOF
