#!/bin/sh
# Script written by Sergey Koren

for file in `ls *.qv`; do
   prefix=`echo $file |sed s/.qv//g`
   lumpyFile="$prefix.lumpy.nop6.vcf"
   frcFile=`echo $prefix |awk -F "." '{STR=$1; for (i = 2; i <NF; i++) { STR=STR"."$i; } print STR"_"$NF}'`
   if [ ! -e $lumpyFile ]; then
      lumpyFile="$prefix.lumpy.vcf"
   fi
   QV=`cat $file`

   if [ -e $lumpyFile ]; then
      echo "QV: $QV, variants for $prefix"
      cat $lumpyFile |grep -v "#" |awk -F "\t" '{print $8}'|awk -F ";" '{print $1}'|sort |uniq -c
      echo "FRC features for $frcFile"
      cat $frcFile*_Features.txt |awk '{print $2}'|sort |uniq -c
   elif [ -e $prefix.1.lumpy.vcf ]; then
      echo "QV: $QV for $prefix, multiple samples:"
      for vcf in `ls $prefix.[0-9]*.lumpy.vcf`; do
         echo "Sample $vcf:"
         cat $vcf |grep -v "#" |awk -F "\t" '{print $8}'|awk -F ";" '{print $1}'|sort |uniq -c
      done
      for feat in `ls $prefix*_Features.txt`; do
         echo "Sample $feat:"
         cat $feat |awk '{print $2}'|sort |uniq -c
      done
   else
      echo "QV: $QV for $prefix"
   fi
done
