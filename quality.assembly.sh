#!/bin/bash
#PBS -N cat
#PBS -L tasks=1:lprocs=1:memory=25gb
#PBS -l walltime=10:00:00
cd "$PBS_O_WORKDIR"

module load Python
module load SMALT
module load SAMtools

### EXTEND MC SEQUENCES 
#ls *minicircles.fasta > input.list
#while read line; do ./fasta_extend.py $(echo ${line%%.*}).extended.fasta $line 150; done < input.list

### INDEX
#ls *extended.fasta > input.list2
#while read line;do smalt index -k 5 -s 2 $line $line;done < input.list2 

### MAP
#while read line;do qsub -N $line quality.assembly.sh;done < input.list2
#ln -s ../../komics_BGI/*trimmed* .
SAMPLE=$(echo ${PBS_JOBNAME%%.*})
reads1=$SAMPLE.reads1.trimmed.fq.gz
reads2=$SAMPLE.reads2.trimmed.fq.gz
#smalt index -k 5 -s 2 $line $line
#smalt map -f sam -y 0.95 -o $SAMPLE.sam $PBS_JOBNAME $reads1 $reads2
./mapping_stats.sh $SAMPLE.sam 

#while read line;do ./mapping_stats.sh $line;done < input.list.sam
