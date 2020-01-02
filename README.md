# nbPaper1

This directory contains the code to produce the results for Leavitt et al. 2020, 
published in ??? (LINK HERE). It contains code to run the simulations assessing 
the performance of the naive Bayes transmission method including estimating the 
reproductive number and the application to a TB outbreak in Hamburg, Germany. 
These programs use the nbTransmission R package 
(https://github.com/sarahleavitt/nbTransmission) as well as the simulation programs 
in the nbSimulation directory (https://github.com/sarahleavitt/nbSimulation).
 

## Simulation Programs

The following programs were used to run a main simulation to test the performance 
of the naive Bayes transmission method and then perform various sensitivity analyses. 

### SimRun.R

This program contains a function "simRun" which performs one iteration of a 
simulation which simulates an outbreak with pathogen genetics and covariates,
estimates the relative transmission probability using naive Bayes, assesses the 
performance, and estimates the reproductive number. 


### PerformSimulation.R

This program is a wrapper for "simRun" which runs the simulation coded in that 
function nSim times. It saves the raw pair-level results for one iteration 
and the effect estimates, reproductive number, and performance for all iterations. 
This is run by "SimQsub.qsub" with sample size and observation date as inputs.


### SampleSizeRun.R

This program has a function called "sampleSizeRun" which is analogous to "simRun" 
but is specifically designed to assess the affect of the sample size and the 
training dataset proportion on performance and so it varies the training proportion 
from 0.1 to 1.


### SampleSizePerform.R

This program is analogous to PerformSimulation.R but runs multiple iterations of 
"sampleSizeRun". This is run by "SampleSizeQsub.qsub" with sample size as an input.


### GenomeTesting.R

This program is used to run a sensitivity analysis of the pathogen genome length.
It simulates one outbreak and then simulates pathogen genomes of various lengths
multiple times. Then it calculates various summaries of the SNP distances across 
the outbreaks. This is run by "GenomeQsub.qsub" with no inputs.


### SimQsub.qsub

This a qsub file to run PerformSimulation.R as a batch job. It runs t parallel jobs
and therefore the total number of simulations is t*nSim times. It needs sample size
and observation date (1 for infection date, 2 for sampling date) as inputs and 
gets them from the shell script "Simulation.sh".


### SampleSizeQsub.qsub

This a qsub file to run SampleSizePerform.R as a batch job. It runs t parallel jobs
and therefore the total number of simulations is t*nSim times. It needs sample size
as an input and gets it from the shell script "Simulation.sh".


### GenomeQsub.qsub

This a qsub file to run GenomeTesting.R as a batch job, it does not have any inputs. 
It runs t parallel jobs and therefore the total number of simulations is t*nSim times.


### Simulation.sh

This is a shell script that runs SimQsub.qsub with infection date for the main 
simulation, SimQsub.qsub with sampling date for the date sensitivity analysis,
SampleSizeQsub.qsub with different sample sizes for the training proportion analysis, 
and GenomeQsub.qsub for the pathogen genome length sensitivity analysis.


### Paper1Results.R

This program reads in the results from all of the simulation programs and then 
calculates all results used in Leavitt et al. 2020 including values in the text,
tables and figures.


***

## Hamburg Analysis Programs

The following programs were used to prepare, analyze, and evaluate applying the 
naive Bayes transmission method to a TB outbreak in Hamburg Germany anaylzed by
Roetzer et al. 2013 (https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001387)
The full dataset used is not publically available. 


### HamburgPrep.R

This program reads in the demographic, contact, and pathogen genetic data, creates 
cleaned individual-level and pair-level datasets to be used for analysis.


### HamburgAnalysis.R

This program calculates relative transmission probabilities using naive Bayes, 
training with contact investigation, training with SNP distance, using various 
serial intervals, and randomly. It then estimates the monthly and over all average 
reproductive numbers with 95% confidence intervals. It also runs a sensitivity 
analysis for the small sample size correction.


### HamburgEvaluate.R

This program takes the results saved in HamburgAnalysis.R and creates all results 
used in Leavitt et al. 2020 including values in the text, tables and figures.


