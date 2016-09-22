#!/bin/sh
# Script written by Sergey Koren

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
  jobid=$1
fi
if [ x$jobid = x ]; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line.
  exit 1
fi

file=`ls *.fasta |sort |head -n $jobid |tail -n 1`
if [ ! -e index/$file.bwt ]; then
   echo "Need to index $file"
   mkdir -p index
   bwa index $file -p index/$file
fi
