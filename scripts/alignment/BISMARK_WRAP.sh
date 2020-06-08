#!/usr/bin/env bash

for i in "$@"
do
  case $i in
    --bismark_path=*)
    BISMARK=${i#*=}
    ;;
    --seq_type=*)
    SEQ_TYPE=${i#*=}
    ;;
    --threads=*)
    THREADS=${i#*=}
    ;;
    --log_file=*)
    LOG=${i#*=}
    ;;
    --R1=*)
    R1=${i#*=}
    shift
    ;;
    --R1_unmapped=*)
    R1_unmapped=${i#*=}
    shift
    ;;
    --R2=*)
    R2=${i#*=}
    shift
    ;;
    --R2_unmapped=*)
    R2_unmapped=${i#*=}
    ;;
    --output=*)
    OUTPUT_FILE=${i#*=}
    ;;
    --methref=*)
    METHREF=${i#*=}
    ;;
    --output_dir=*)
    OUTDIR=${i#*=}
    ;;
    --sample_name=*)
    SAMPLE=${i#*=}
    ;;
  esac
done

if [[ $R2 == "" ]]; then
  R2_unmapped="${R1_trim/R1/R2}"
fi

if [[ $SEQ_TYPE == 2x* ]]
then
  $BISMARK --genome $METHREF --parralel $threads --quiet --unmapped -o $OUTDIR -1 $R1 -2 $R2
  mv $OUTDIR/$SAMPLE"_R1_trim_bismark_bt2_pe.bam" $OUTPUT_FILE
  mv $OUTDIR/$SAMPLE"R1_trim_bismark_bt2_PE_report.txt" $LOG

  mv $OUTDIR/$SAMPLE"_R1_trim.fastq_unmapped_reads_1.fq.gz" $UNMAPPED_R1
  mv $OUTDIR/$SAMPLE"_R1_trim.fastq_unmapped_reads_2.fq.gz" $UNMAPPED_R2
elif [[ $SEQ_TYPE == 1x* ]]
then
  $BISMARK --genome $METHREF --parralel $threads --quiet --unmapped -o $OUTDIR $R1
  mv $OUTDIR/$SAMPLE"_R1_trim_bismark_bt2.bam" $OUTPUT_FILE
  mv $OUTDIR/$SAMPLE"R1_trim_bismark_bt2_SE_report.txt" $LOG

  mv $OUTDIR/$SAMPLE"_R1_trim.fastq_unmapped_reads.fq.gz" $R1_unmapped
  echo "This was processed as a single-end dataset so this is just a fun placeholder file" > $R2_unmapped
else
  echo "Error. Library type not found for $SEQ_TYPE and $PREP_TYPE"
  exit 1
fi
