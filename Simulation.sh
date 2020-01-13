#!/bin/bash

#Running the main simulation
#sample size = 500, observation date = 1 (infection date)
qsub -N FullSimulation SimQsub.qsub 500 1

#Running the simulation using sampling date
#sample size = 500, observation date = 2 (sampling date)
qsub -N SampleDate SimQsub.qsub 500 2

#Running the genome length analysis
qsub -N GenomeQsub.qsub

#Running with different sample sizes and training proportions
sampleSize=(50 250 500)
for h in ${sampleSize[@]}; do
	qsub -N SampleSizeSim_${h} SampleSizeQsub.qsub $h
done

