#!/bin/sh
#$ -N 2021-06-14_PF00005_seedtrees      
#$ -cwd
#$ -o 2021-06-14_PF00005_seed.screen.log
#$ -e 2021-06-14_PF00005_seed.error.log
#$ -pe smp 16

module load beagle-lib/20140322-goolf-1.5.14
module load Java/1.8.0_5

echo A | java -jar ../../beast/lib/launcher.jar -seed 1111 -beagle_SSE ./2021-06-14_PF00005_seed_fix_beauti.xml