#!/usr/bin/env bash

for i in "$@"
do
  case $i in
    --bbduk_path=*)
    BBDUK=${i#*=}
    ;;
    --seq_type=*)
    SEQ_TYPE=${i#*=}
    ;;
    --prep_type=*)
    PREP_TYPE=${i#*=}
    ;;
    --adapter_file=*)
    ADAPTOR_SEQ=${i#*=}
    ;;
    --mem_mb=*)
    MEM=${i#*=}
    ;;
    --threads=*)
    THREADS=${i#*=}
    ;;
    --maq=*)
    MAQ=${i#*=}
    ;;
    --entropy=*)
    ENTROPY=${i#*=}
    ;;
    --log_file=*)
    LOG=${i#*=}
    ;;
    --R1=*)
    R1=${i#*=}
    shift
    ;;
    --R1_trim=*)
    R1_trim=${i#*=}
    shift
    ;;
    --R2=*)
    R2=${i#*=}
    shift
    ;;
    --R2_trim=*)
    R2_trim=${i#*=}
    ;;
  esac
done

if [[ $R2 == "" ]]; then
  R2_trim="${R1_trim/R1/R2}"
fi

if [[ $PREP_TYPE == "MEYER_SSLP" ]] && [[ $SEQ_TYPE == "2x75" ]]
then
  $BBDUK in1=$R1 in2=$R2 out1=$R1_trim out2=$R2_trim -Xmx1g -threads=$THREADS ref=$ADAPTOR_SEQ maq=$MAQ entropy=$ENTROPY tbo tpe &> $LOG
elif [[ $PREP_TYPE == "SRSLY_SSLP" ]] && [[ $SEQ_TYPE == "2x75" ]]
then
  $BBDUK in1=$R1 in2=$R2 out1=$R1_trim out2=$R2_trim -Xmx1g -threads=$THREADS ref=$ADAPTOR_SEQ maq=$MAQ entropy=$ENTROPY tbo tpe &> $LOG
elif [[ $PREP_TYPE == "SWIFT_ACCEL" ]] && [[ $SEQ_TYPE == 2x1* ]] #covers paired end of type 2x100 or 2x150bp
then
  $BBDUK in1=$R1 in2=$R2 out1=$R1.tmp.fastq out2=$R2.tmp.fastq -Xmx1g -threads=$THREADS ftr=74 &> $LOG
  $BBDUK in1=$R1.tmp.fastq in2=$R2.tmp.fastq out1=$R1_trim out2=$R2_trim -Xmx1g -threads=$THREADS ref=$ADAPTOR_SEQ maq=$MAQ entropy=$ENTROPY tbo tpe swift=t mink=11 ktrim=r &>> $LOG
  rm $R1.tmp.fastq $R2.tmp.fastq
elif [[ $PREP_TYPE == "SWIFT_ACCEL" ]] && [[ $SEQ_TYPE == "2x75" ]]
then
  $BBDUK in1=$R1 in2=$R2 out1=$R1_trim out2=$R2_trim -Xmx1g -threads=$THREADS ref=$ADAPTOR_SEQ maq=$MAQ entropy=$ENTROPY tbo tpe swift=t mink=11 ktrim=r &>> $LOG
elif [[ $PREP_TYPE == "SWIFT_ACCEL" ]] && [[ $SEQ_TYPE == "1x75" ]]
then
  $BBDUK in=$R1 out=$R1_trim ref=$ADAPTOR_SEQ maq=$MAQ entropy=$ENTROPY swift=t mink=11 ktrim=r &> $LOG
  echo "This was processed as a single-end dataset so this is just a fun placeholder file" > $R2_trim #the lazyest of ways ...
else
  echo "Error. Library type not found for $SEQ_TYPE and $PREP_TYPE"
  exit 1
fi
