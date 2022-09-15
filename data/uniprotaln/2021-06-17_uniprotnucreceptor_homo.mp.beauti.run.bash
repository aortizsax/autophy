#!/bin/sh
#$ -N nucreptrees      
#$ -cwd
#$ -o nucreptrees.screen.log
#$ -e nucreptrees.error.log
#$ -pe smp 16

module load beagle-lib/20140322-goolf-1.5.14
module load Java/1.8.0_5

BEAGLE_LIB_PATH="../../beagle-lib-master/beagle-lib-master/"

#model JTT+G+STRICTCLOCK+yule+exponBirthrate
echo A | srun java -Xmx10g -Djava.library.path=$BEAGLE_LIB_PATH -jar ../../beast/lib/launcher.jar -threads 16 -beagle_SSE -seed 777 ./2021-06-17_uniprotnucreceptor_homo_beauti.xml
