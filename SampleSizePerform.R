#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program runs the full set of simulations to assess sample size and
# training proportion by running SampleSizeRun.R
################################################################################

rm(list = ls())
options(scipen=999)

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(pROC)
library(caret)
library(gtools)
library(devtools)
library(roxygen2)


#### Batch Mode ####

setwd("/project/sv-thesis/nbPaper1/")
#Getting sample size from the arguements
args <- commandArgs(trailingOnly = TRUE)
sampleSize <- as.numeric(args[])
#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#The number of simulations per split
nSim <- 10


#### Interactive Mode ####

# setwd("~/Boston University/Dissertation/nbPaper1")
# iTask <- 2
# sampleSize <- 50
# nSim <- 1


########### Getting access to functions in other programs ##############

load_all("../nbTransmission")
source("../nbSimulation/SimOutbreak.R")
source("../nbSimulation/SimulateOutbreakS.R")
source("../nbSimulation/SimCovariates.R")
source("../nbSimulation/PerformInterval.R")
source("../nbSimulation/PerformRandom.R")
source("../nbSimulation/SimEvaluate.R")
source("SampleSizeRun.R")



################## Setting the oubreak parameters #################

#print(sampleSize)

#Effective population size times generation time (1 day or 1/365 years)
neg <- 0.25
pi <- 1
#Sets reproductive number to off.r if off.p = 0.5
off.r <- 1.2
off.p <- 0.5
#From Yincheng's work
w.shape <- 1.05
w.scale <- 1 / (0.0014 * 365)
w.shift <- 0.25
ws.shape <- w.shape
ws.scale <- w.scale
ws.shift <- w.shift
time <- 20
#Do you want multiple outbreaks
multOutbreaks <- TRUE 

#Genotype parameters
rootseq <- NULL
length <- 300
#This gives the equivalent of 0.5 mutations/genome/year
rate <- 0.5 / length
covariates <- c("Y1", "Y2", "Y3", "Y4", "timeCat")
#Thresholds for SNP distance
thresholds <- c(2, 12)


############### Running full simulation ##############

repNumI <- NULL
performance <- NULL

#Running the simuliaton for the rest of the outbreaks
for (iteration in 1:nSim){
  
  #Setting the seed for each run so I can re-run specific iterations with errors
  set.seed(iTask * 1000 + iteration)
  
  #Running the simulation
  res <- sampleSizeRun()
  print("Finished calculating probabilities")
  
  #Combining results
  pTemp <- res[[1]] %>% mutate(runID = paste(iTask, iteration, sep = "_"))
  r0Temp <- res[[2]] %>% mutate(runID = paste(iTask, iteration, sep = "_"))
  
  #Adding to the existing dataframes
  performance <- bind_rows(performance, pTemp)
  repNumI <- bind_rows(repNumI, r0Temp)
  
  #Printing message that the run is finished
  print(paste0("Completed run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
}

#Saving dataframes with summary of results for all of the runs
saveRDS(performance, file=paste0("../Simulation_Results/performanceSS", sampleSize, "_", iTask, ".rds"))
saveRDS(repNumI, file=paste0("../Simulation_Results/repNumISS", sampleSize, "_", iTask, ".rds"))

