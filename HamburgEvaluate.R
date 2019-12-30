#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program creates the figures and tables to analyze the Hamburg data
# The data are cleaned in HamburgPrep.R and analyzed in HamburgAnalysis.R
################################################################################

setwd("~/Boston University/Dissertation")
#rm(list = ls())
options(scipen = 999)

library(dplyr)
library(tidyr)
library(tableone)
library(ggplot2)
library(igraph)
library(RColorBrewer)
library(pheatmap)

#Reading in cleaned datasets from HamburgPrep.R and results from HamburgAnalysis.R
hamInd <- readRDS("Datasets/HamburgInd.rds")
hamPair <- readRDS("Datasets/HamburgPair.rds")
resultsHam <- readRDS("Datasets/HamburgResults_12.17.19.rds")
R0CI <- readRDS("Datasets/HamburgR0CI_12.17.19.rds")
RtCI <- readRDS("Datasets/HamburgRtCI_12.17.19.rds")
R0Correction <- readRDS("Datasets/HamburgCorrection_12.17.19.rds")


#Making clear scenario labels
R0CI <- R0CI %>% mutate(Scenario = factor(label, levels = c("HamContacts",
                                                            "HamSNPs",
                                                            "GenInt0.5/1.9",
                                                            "GenInt1.1/2",
                                                            "GenInt1.3/3.3",
                                                            "Random"),
                                          labels =  c("Training: Contacts",
                                                      "Training: SNP Distance",
                                                      "Narrow Serial Interval",
                                                      "Medium Serial Interval",
                                                      "Wide Serial Interval",
                                                      "Random Probabilities")))

RtCI <- RtCI %>% mutate(Scenario = factor(label, levels = c("HamContacts",
                                                            "HamSNPs",
                                                            "GenInt0.5/1.9",
                                                            "GenInt1.1/2",
                                                            "GenInt1.3/3.3",
                                                            "Random"),
                                          labels =  c("Training: Contacts",
                                                      "Training: SNP Distance",
                                                      "Narrow Serial Interval",
                                                      "Medium Serial Interval",
                                                      "Wide Serial Interval",
                                                      "Random Probabilities")))

resultsHam <- resultsHam %>% mutate(Scenario = factor(label, levels = c("HamContacts",
                                                                        "HamSNPs",
                                                                        "GenInt0.5/1.9",
                                                                        "GenInt1.1/2",
                                                                        "GenInt1.3/3.3",
                                                                        "Random"),
                                                      labels =  c("Training: Contacts",
                                                                  "Training: SNP Distance",
                                                                  "Narrow Serial Interval",
                                                                  "Medium Serial Interval",
                                                                  "Wide Serial Interval",
                                                                  "Random Probabilities")))




######################## Description of Cohort ##########################

#Basic descriptive statistics
summary(hamInd$Age)
table(hamInd$Sex)
prop.table(table(hamInd$Sex))
#Number in each study (Hamburg and Schleswig-Holstein)
table(hamInd$Study)
prop.table(table(hamInd$Study))
#Number with contact tracing
sum(!is.na(hamInd$Group))
sum(!is.na(hamInd$Group)) / nrow(hamInd)

#SNP distance distribution
summary(hamPair$snpDist)



#### Table: Individual-Level Covariates ####
hamInd <- hamInd %>% mutate(NationalityC = ifelse(Nationality == "Germany", "Germany",
                                                  "Other"))
indCat <- c("Study", "NationalityC", "Sex", "AgeGroup", "SmearPos", "HIV",
            "SubstanceAbuse", "Residence", "Milieu")
covarInd <- CreateTableOne(vars = indCat, factorVars = indCat,
                           data = hamInd, test = FALSE)
covarInd <- as.data.frame(print(covarInd, showAllLevels = TRUE))



#### Table: Pair-level Covariates ####
pairTable <- (resultsHam
              %>% filter(label == "RandomHam")
              %>% replace_na(list(SameGroup = "missing", snpClose = "missing"))
)
pairCat <- c("Study", "Nationality", "Sex", "Age", "SmearPos", "HIV", 
             "SubstanceAbuse", "Residence", "Milieu",  
             "TimeCat", "snpClose", "SameGroup")
covarPair <- CreateTableOne(vars = pairCat, factorVars = pairCat,
                            data = pairTable, test = FALSE)

covarPair <- as.data.frame(print(covarPair, showAllLevels = TRUE))



#### Figure: Barplot of Case Counts ####

table(hamInd$IsolationYear)

#Plot of cases by month
ggplot(data = hamInd) +
  geom_bar(aes(x = IsolationYear)) +
  xlab("Isolation Year") +
  ylab("Number of Cases") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10)) +
  ggsave(file = "Figures/Hamburg_Years.tiff",
         width = 4, height = 3, units = "in", dpi = 300)




####################### Transmission Probabilities ########################


#### Figure: Heatmaps of Probabilties ####

nodes <- hamInd %>% select(individualID, IsolationDate)
colBreaks <- c(-0.01, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)

## SNP Distance ##
edgesSNPs <- (resultsHam
              %>% filter(label == "HamSNPs")
              %>% select(individualID.1, individualID.2, pScaled,
                         snpDist, snpClose, SameGroup)
              #Arrange so true edges and high probability edges are drawn first
              %>% arrange(pScaled)
)
netSNPs <- graph_from_data_frame(d = edgesSNPs, vertices = nodes, directed = T)
E(netSNPs)$pGroup <- cut(E(netSNPs)$pScaled, breaks = colBreaks, labels = 1:9)
#First get adjacency version of the network using pScaled
netSNPs.adj <- get.adjacency(netSNPs, attr = "pScaled", sparse = FALSE)

## Contact Tracing ##
edgesContact <- (resultsHam
                 %>% filter(label == "HamContacts")
                 %>% select(individualID.1, individualID.2, pScaled,
                            snpDist, snpClose, SameGroup)
                 %>% arrange(pScaled)
)
netContact <- graph_from_data_frame(d = edgesContact, vertices = nodes, directed = T)
E(netContact)$pGroup <- cut(E(netContact)$pScaled, breaks = colBreaks, labels = 1:9)
#First get adjacency version of the network using pScaled
netContact.adj <- get.adjacency(netContact, attr = "pScaled", sparse = FALSE)

## Serial Interval ##
edgesSI <- (resultsHam
            %>% filter(label == "GenInt1.1/2")
            %>% select(individualID.1, individualID.2, pScaled,
                       snpDist, snpClose, SameGroup)
            %>% arrange(pScaled)
)
netSI <- graph_from_data_frame(d = edgesSI, vertices = nodes, directed = T)
E(netSI)$pGroup <- cut(E(netSI)$pScaled, breaks = colBreaks, labels = 1:9)
#First get adjacency version of the network using pScaled
netSI.adj <- get.adjacency(netSI, attr = "pScaled", sparse = FALSE)

## Random ##
edgesRand <- (resultsHam
              %>% filter(label == "Random")
              %>% select(individualID.1, individualID.2, pScaled,
                         snpDist, snpClose, SameGroup)
              %>% arrange(pScaled)
)
netRand <- graph_from_data_frame(d = edgesRand, vertices = nodes, directed = T)
E(netRand)$pGroup <- cut(E(netRand)$pScaled, breaks = colBreaks, labels = 1:9)

#First get adjacency version of the network using pScaled
netRand.adj <- get.adjacency(netRand, attr = "pScaled", sparse = FALSE)



png("Figures/HeatmapRand.png", width = 4, height = 4,
    units = "in", res = 300)
par(mar = c(0, 0, 1, 0))
pheatmap(t(netRand.adj), cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         col=brewer.pal(9,"Blues"), breaks = colBreaks, legend = FALSE)
dev.off()

png("Figures/HeatmapSI.png", width = 4.4, height = 4,
    units = "in", res = 300)
par(mar = c(0, 0, 1, 0))
pheatmap(t(netSI.adj), cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         col=brewer.pal(9,"Blues"), breaks = colBreaks)
dev.off()

png("Figures/HeatmapSNPs.png", width = 4, height = 4,
    units = "in", res = 300)
par(mar = c(0, 0, 1, 0))
pheatmap(t(netSNPs.adj), cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         col=brewer.pal(9,"Blues"), breaks = colBreaks, legend = FALSE)
dev.off()

png("Figures/HeatmapContact.png", width = 4, height = 4,
    units = "in", res = 300)
par(mar = c(0, 0, 1, 0))
pheatmap(t(netContact.adj), cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, show_rownames = FALSE, show_colnames = FALSE,
         col=brewer.pal(9,"Blues"), breaks = colBreaks, legend = FALSE)
dev.off()





###################### Reproductive Number ########################

#### Figure: Plot of Monthly Rt with CIs ####

totalTime <- max(RtCI$timeRank) - min(RtCI$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.9 * totalTime)

ggplot(data = RtCI, aes(x = timeRank, y = Rt)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.1) +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Isolation Year", breaks = seq(0, 167, 12),
                     labels = seq(1997, 2010, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7) +
  geom_hline(data = R0CI, aes(yintercept = RtAvg), size = 0.7) +
  ggsave(file = "Figures/Hamburg_Rt.tiff",
         width = 8, height = 5, units = "in", dpi = 300)


#### Supplementary Table: Average Rt ####
R0CI

#### Figure: Plot of Average Rt CIs ####

ggplot(data = R0CI, aes(x = Scenario, y = RtAvg)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.3) +
  scale_y_continuous(name = "Average Effective Reproductive Number", 
                     limits = c(0.5, 1.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  ggsave(file = "Figures/Hamburg_R0.tiff",
         width = 5, height = 4, units = "in", dpi = 300)




################# Correction Sensitivity Analysis #####################

#### Supplementary Table: Average Rt by Correction ####
R0Correction %>% arrange(label)


