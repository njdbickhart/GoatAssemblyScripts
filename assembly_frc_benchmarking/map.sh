#!/bin/bash
# Written by Sergey Koren

syst=`uname -s`
arch=`uname -m`
name=`uname -n`

if [ "$arch" = "x86_64" ] ; then
  arch="amd64"
fi
if [ "$arch" = "Power Macintosh" ] ; then
  arch="ppc"
fi

jobid=$SGE_TASK_ID
refid=$1

if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
   jobid=$2
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

NUM_JOBS=`cat mapping.fofn |wc -l |awk '{print $NF}'`
if [ $jobid -gt $NUM_JOBS ]; then
   echo "Invalid job id $jobid, max is $NUM_JOBS"
   exit
fi

# requirements: samtools 1.3, Lumpy 0.2.11 all in path

# find programs, assuming they are in path
LUMPY=`which lumpy`
LUMPY_SCRIPTS=`dirname $LUMPY |awk '{print "$1/../scripts/"}'`

reads=`cat mapping.fofn |head -n $jobid |tail -n 1 |awk '{print $1}'`
prefix=`echo $refid |sed s/.fasta//g`
readsPrefix=`basename $reads |sed s/.raw.fastq.gz//g |sed s/.fastq.gz//g |sed s/.fastq//g |sed s/.fq.gz//g |sed s/.fq//g`
prefix="${prefix}_${readsPrefix}"
type=`cat mapping.fofn |head -n $jobid |tail -n 1 |awk '{print $NF}'`
echo "Running and will map reads $reads to $refid will prefix to $prefix.bam"

INPUT=""
OPTIONS=""
if [ $type == "pacbio" ]; then
  INPUT="$reads"
  OPTIONS="-x pacbio "
elif [ $type == "illumina" ]; then
   INPUT="${reads}_1.fastq.gz ${reads}_2.fastq.gz"
elif [ $type == "illuminaR1" ]; then
   readsV2=`echo $reads |sed s/_R1_/_R2_/g`
   INPUT="${reads} ${readsV2}"
else
   echo "Unknown data type $type"
   exit
fi

mkdir -p map

if [ ! -e map/$prefix.splitters.bam ]; then
   if [ ! -e map/$prefix.unsorted.bam ]; then
      echo "Runing bwa mem -t 16 $OPTIONS index/$refid $INPUT 2> $prefix.err |samtools view -S -b - > map/$prefix.unsorted.bam"
      bwa mem -t 16 $OPTIONS index/$refid $INPUT 2> $prefix.err |samtools view -S -b - > map/$prefix.unsorted.bam
   fi

   if [ ! -e map/$prefix.splitters.unsorted.bam ]; then 
      echo "Extracting splitters/discordants"
      # extract
      samtools view -b -F 1294 map/$prefix.unsorted.bam > map/$prefix.discordants.unsorted.bam
      samtools view -h map/$prefix.unsorted.bam \
          | $LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin \
          | samtools view -Sb - \
          > map/$prefix.splitters.unsorted.bam
    fi
    if [ ! -e map/$prefix.splitters.bam ]; then
      # sort
      echo "Sorting"
      ulimit -Su 160000
      ulimit -all

      samtools sort -@ 16 -m 2G -o map/$prefix.bam -T `pwd`/$prefix.tmp map/$prefix.unsorted.bam
      samtools sort -@ 16 -m 2G -o map/$prefix.discordants.bam -T `pwd`/$prefix.tmp map/$prefix.discordants.unsorted.bam
      samtools sort -@ 16 -m 2G -o map/$prefix.splitters.bam -T `pwd`/$prefix.tmp map/$prefix.splitters.unsorted.bam
   fi
fi
