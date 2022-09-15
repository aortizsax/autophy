#!/bin/sh
#$ -N nucreptrees      
#$ -cwd
#$ -o nucreptrees.screen.log
#$ -e nucreptrees.error.log
#$ -pe smp 16

module load beagle-lib/20140322-goolf-1.5.14
module load Java/1.8.0_5

echo A | java -jar ../../beast/lib/launcher.jar -seed 1111 -beagle_SSE ./uniprotnucreceptor_homo_beauti.xml

