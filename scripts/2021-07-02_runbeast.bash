#!/bin/bash
#$ -N runbeast      
#$ -cwd
#$ -o 2021-07-02_runbeast.screen.log
#$ -e 2021-07-01_PF00067.error.log
#$ -pe smp 2


echo "Positional Parameters"
echo '$0 = ' $0
echo '$1 = ' $1
echo '$2 = ' $2
echo '$3 = ' $3
