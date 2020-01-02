#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

##################################################################################
# This program tests the effect of changing the genome length including running
# a set of simulations and then reading in the results and making a figure
##################################################################################

rm(list = ls())
options(scipen=999)

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(glmnet)
library(FSelector)
library(caret)
library(gtools)

source("SimulateOutbreakS.R")

#### Batch Mode ####

setwd("/project/sv-thesis/nbPaper1/")
#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#The number of simulations per split
nSim <- 10

#### Interactive Mode ####

# setwd("~/Boston University/Dissertation/nbPaper1")
# iTask <- 1
# nSim <- 2


################## Setting the oubreak parameters #################

#Effective population size times generation time
neg <- 0.25
pi <- 1
#Sets reproductive number to off.r if off.p = 0.5
off.r <- 1.2
off.p <- 0.5
#From Yincheng's work
w.shape <- 1.05
w.scale <- 1 / (0.0014 * 365)
ws.shape <- w.shape
ws.scale <- w.scale
#Total sample size
sampleSize <- 200
#Do you want multiple outbreaks
multOutbreaks <- FALSE 
rootseq <- NULL
#Date of observation
observationDate = "infectionDate"
#Probability calculation parameters
covariates <- c("Y1", "Y2", "Y3", "Y4", "timeCat1")
pTraining <- 0.5
goldStandard <- "transmission"
scheme = "info"



################# Testing genotype length ###############

#Setting seed here so same outbreak is used each time
set.seed(15039)
simu <- simulateOutbreakS(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                          w.shape = w.shape, w.scale = w.scale,
                          ws.shape = ws.shape, ws.scale = ws.scale,
                          nSampled = sampleSize, dateStartOutbreak = 1980, dateT = 2000)
#Extract the phylogenetic tree and convert it to a phylo object
p <- phyloFromPTree(extractPTree(simu))


calcDistance <- function(length){
  #Simulating the sequences
  rate <- 0.5 / length
  seq <- simSeq(p, l = length, rate = rate, rootseq = rootseq)
  names(seq) <- 10000 + as.numeric(names(seq))
  #Making the sequences into a DNAbin object
  seqBin <- as.DNAbin(seq)
  snpDist <- dist.gene(seqBin)
  snpDistDf <- melt(as.matrix(snpDist), varnames = c("individualID.1", "individualID.2"))
  snpDistDf <- (snpDistDf
                %>% rename(snpDist = value)
                %>% filter(individualID.1 > individualID.2)
  )
  
  pLow <- round(100 * sum(snpDistDf$snpDist < 2) / nrow(snpDistDf), 0)
  pHigh <- round(100 * sum(snpDistDf$snpDist > 12) / nrow(snpDistDf), 0)
  results <- c("length" = length, round(summary(snpDistDf$snpDist), 0),
               "p<2" = pLow, "P>12" = pHigh)
  
  return(as.data.frame(results))
}

#Set seed for this run's set of genomes
set.seed(iTask * 1000 + 1)

results <- NULL
for (i in 1:nSim){
  lengths <- c(50, 100, 500, 1000, 5000, 10000, 100000, 4400000)
  #lengths <- c(50, 100, 500)
  for (l in lengths){
    rTemp <- calcDistance(l)
    rTemp$runID <- paste0(iTask, "_", i)
    results <- bind_rows(results, rTemp)
    print(paste0("Finished simulated length of ", l))
  }
  print(paste0("Finished run ", i))
}

saveRDS(results, paste0("../Simulation_Results/GenomeTest_", iTask, ".rds"))
