#!/bin/bash

#PBS -N cat
#PBS -L tasks=1:lprocs=1:memory=10gb
#PBS -l walltime=01:00:00

cd "$PBS_O_WORKDIR"


module load SAMtools

#for f in *.sorted.bam

#do
#	UNM=$(samtools view -c -f 4 $f)
#	M=$(samtools view -c -F 4 $f)
#	MQ25=$(samtools view -c -q 25 -F 4 $f)
#	BM=$(($M-$MQ25))
#	echo ${f%%.*} "mapped_bad_quality" $BM >> map.stat3.txt
#	echo ${f%%.*} "mapped_Q>25" $MQ25 >> map.stat3.txt
#	echo ${f%%.*} "unmapped" $UNM >> map.stat3.txt
#done


##### SET VARIABLE BEFORE SUBMITTING JOB !!!
#for f in *.sorted.bam; do qsub -v inputfile=$f -N MAPSTAT.$(echo ${f%%.*}) MAPSTAT.pbs; done

##### JOB
#REGEX="[^/]*$" 
SAMPLE=$(echo ${inputfile%%.*}) 

TR=$(samtools view -c $SAMPLE)
UNMR=$(samtools view -c -f 4 $SAMPLE)
MR=$(samtools view -c -F 4 $SAMPLE)
MR_MQ25=$(samtools view -c -q 25 -F 4 $SAMPLE)
MR_LQ=$(($M-$MQ25))
#echo ${f%%.*} "mapped_bad_quality" $BM >> map.stat3.txt
#echo ${f%%.*} "mapped_Q>25" $MQ25 >> map.stat3.txt
#echo ${f%%.*} "unmapped" $UNM >> map.stat3.txt

echo $SAMPLE $UNMR $MR $MR_MQ25 $MR_LQ > results.table.mapstat.txt
 

