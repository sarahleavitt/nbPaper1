#!/bin/bash

# Sample Size
sampleSize=(50 250 500)

#Running the full simulation
qsub -N FullSimulation SimQsub.qsub 500 1

#Running the simulation using sampling date
qsub -N SampleDate SimQsub.qsub 500 2

#Running the genome length analysis
qsub -N GenomeQsub.qsub

for h in ${sampleSize[@]}; do
	qsub -N SampleSizeSim_${h} SampleSizeQsub.qsub $h
done

