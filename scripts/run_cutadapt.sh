#!/bin/bash
set -e # this makes the whole script exit on any error.
#fill these variables with the arguments given to this script
sample=$1
fq1=$2
fq2=$3
logDir=log # this is where all the files to keep track of progress will go.
mkdir -p log # make the directory where log files will go, if it doesn't exist already
echo running pipeline for $sample
if [ ! -e $logDir/$sample.cutadapt.done ] #run this code only if $logDir/$sample.fastqc.done is missing
then

        cutadapt --untrimmed-o untrimmed1.$1.fastq.gz --untrimmed-p untrimmed2.$1.fastq.gz -m 80 --too-short-o short1.$1.fastq.gz --too-short-p short2.$1.fastq.gz -M 80 --too-long-o long1.$1.fastq.gz --too-long-p long2.$1.fastq.gz -a ^TGCATTTTTTTCACATC...GGTTACGGCTGTT -A ^AACAGCCGTAACC...GATGTGAAAAAAATGCA -o trimmed.1.$1.fastq.gz -p trimmed.2.$1.fastq.gz /arc/project/st-cdeboer-1/nick/PlasmidLibraryAmplificationBias/data/$2 /arc/project/st-cdeboer-1/nick/PlasmidLibraryAmplificationBias/data/$3
        
        touch $logDir/$sample.cutadapt.done #create the file that we were looking for at the beginning of this if statement so that this same code is not run next time
else # $logDir/$sample.fastqc.done was not missing
        echo Already performed cutadapt of $sample
fi