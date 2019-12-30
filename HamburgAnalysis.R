#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program estimates transmission probabilities and reproductive number
# for the Hamburg dataset using the NB transmission method.
# The data are cleaned in HamburgPrep.R
################################################################################

setwd("~/Boston University/Dissertation")
rm(list = ls())
set.seed(103020)

library(dplyr)
library(purrr)
library(tidyr)
library(devtools)

#Sourcing functions
load_all("../nbTransmission")
source("nbSimulation/PerformRandom.R")
source("nbSimulation/PerformInterval.R")


#Reading in datasets from HamburgPrep.R
hamInd <- readRDS("Datasets/HamburgInd.rds")
hamPair <- readRDS("Datasets/HamburgPair.rds")

hamPair <- hamPair %>% mutate(snpClose = ifelse(snpDist < 2, TRUE,
                                         ifelse(snpDist > 12, FALSE, NA)))

orderedHam <- hamPair  %>% filter(!is.na(IsolationDiff) & IsolationDiff >= 0)



############## Estimating Probabilities #############

covariates <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

resHam1 <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "snpClose",
                           covariates = covariates, label = "HamSNPs",
                           l = 1, n = 10, m = 1, nReps = 50)

resHam2 <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                           pairIDVar = "edgeID", goldStdVar = "SameGroup",
                           covariates = covariates, label = "HamContacts",
                           l = 1, n = 10, m = 1, nReps = 50)

resHamCov1 <- resHam1$probabilities %>% full_join(orderedHam, by = "edgeID")
resHamCov2 <- resHam2$probabilities %>% full_join(orderedHam, by = "edgeID")



#### Adding Comparison Methods ####

#For interval methods, not allowing intervals of less than 3 months
#Yicheng's estimate from Brazil
Int1 <- performInterval(orderedHam, genShape = 1.05, genScale = 2.0, shift = 0.25,
                        observationDiff = "IsolationDiff")
#Prior from Didelot MolBioEvol 2017 - applied to Hamburg
Int2 <- performInterval(orderedHam, genShape = 1.3, genScale = 3.3, shift = 0.25,
                        observationDiff = "IsolationDiff")
#Posterior from Didelot MolBiolEvol 2017 - estimated from Hamburg
Int3 <- performInterval(orderedHam, genShape = 0.54, genScale = 1.9, shift = 0.25,
                        observationDiff = "IsolationDiff")
#Randomly assigning probabilities
rand <- performRandom(orderedHam)
rComp <- orderedHam %>% full_join(bind_rows(Int1, Int2, Int3, rand), by = "edgeID")

resultsHam <- bind_rows(resHamCov1, resHamCov2, rComp)


#Saving results
saveRDS(resultsHam, "Datasets/HamburgResults_12.17.19.rds")



####################### Reproductive Number ###########################

rInitial <- estimateR(resultsHam %>% filter(label == "HamSNPs"),
                dateVar = "IsolationDate", indIDVar = "individualID",
                pVar = "pScaled", timeFrame = "months")
rt <- rInitial[[2]]

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.9 * totalTime)

#Plotting where to cut
ggplot(data = rt) +
  geom_histogram(aes(x = timeRank), alpha = 0.5, bins = 20) +
  geom_line(aes(x = timeRank, y = Rt)) +
  scale_y_continuous(name = "Rt") +
  geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7, col = "blue") +
  geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7, col = "blue")


#Calculating the reproductive number using the various methods
resList <- split(resultsHam, resultsHam$label)

rListFull <- purrr::map(resList, estimateR, dateVar = "IsolationDate",
                 indIDVar = "individualID", pVar = "pScaled",
                 timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                 bootSamples = 1000, alpha = 0.05)

#Adding label names and combining results for each method
repNumIH <- NULL
repNumMH <- NULL
monthR0H <- NULL

for(i in 1:length(rListFull)){
  rListFull[[i]]$RiDf$label <- names(rListFull)[i]
  rListFull[[i]]$RtDf$label <- names(rListFull)[i]
  rListFull[[i]]$RtAvg$label <- names(rListFull)[i]
  
  repNumIH <- bind_rows(repNumIH, rListFull[[i]]$RiDf)
  repNumMH <- bind_rows(repNumMH, rListFull[[i]]$RtDf)
  monthR0H <- bind_rows(monthR0H, rListFull[[i]]$RtAvg)
}
monthR0H


#Saving the confidence interval datasets
saveRDS(repNumMH, "Datasets/HamburgRtCI_12.17.19.rds")
saveRDS(monthR0H, "Datasets/HamburgR0CI_12.17.19.rds")



################### Sensitivity Analysis for Correction ####################

covariates <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV",
                "SubstanceAbuse", "Residence", "Milieu", "TimeCat")

monthR0H_c <- (monthR0H
               %>% filter(label %in% c("HamContacts", "HamSNPs"))
               %>% mutate(label = paste0(label, "_1"))
)

corrections <- c(0.001, 0.01, 0.1, 0.5)

for(l in corrections){
  
  resHam1_c <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                             pairIDVar = "edgeID", goldStdVar = "snpClose",
                             covariates = covariates, label = paste0("HamSNPs_", l),
                             l = l, n = 10, m = 1, nReps = 50)
  resHamCov1_c <- resHam1_c$probabilities %>% full_join(orderedHam, by = "edgeID")
  
  resHam2_c <- nbProbabilities(orderedPair = orderedHam, indIDVar = "individualID",
                             pairIDVar = "edgeID", goldStdVar = "SameGroup",
                             covariates = covariates, label = paste0("HamContacts_", l),
                             l = l, n = 10, m = 1, nReps = 50)
  resHamCov2_c <- resHam2_c$probabilities %>% full_join(orderedHam, by = "edgeID")
  
  resultsHam_c <- bind_rows(resHamCov1_c, resHamCov2_c)
  #Calculating the reproductive number using the various methods
  resList_c <- split(resultsHam_c, resultsHam_c$label)

  rListFull_c <- purrr::map(resList_c, estimateR, dateVar = "IsolationDate",
                            indIDVar = "individualID", pVar = "pScaled",
                            timeFrame = "months", rangeForAvg = c(monthCut1, monthCut2),
                            bootSamples = 1000, alpha = 0.05)
  
  #Adding label names and combining results for each method
  monthR0H_t <- NULL
  for(i in 1:length(rListFull_c)){
    rListFull_c[[i]]$RtAvg$label <- names(rListFull_c)[i]
    monthR0H_t <- bind_rows(monthR0H_t, rListFull_c[[i]]$RtAvg)
  }
  
  monthR0H_c <- bind_rows(monthR0H_c, monthR0H_t)
}

saveRDS(monthR0H_c, "Datasets/HamburgCorrection_12.17.19.rds")

