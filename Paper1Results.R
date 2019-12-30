#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program makes tables and figures for the simulation results of the 
# NB transmission method and reproductive number estimation.
################################################################################

rm(list = ls())
options(scipen = 999)
setwd("~/Boston University/Dissertation/Simulation_Results_11.22.19")

library(dplyr)
library(tidyr)
library(ggplot2)
library(plotROC)
library(tableone)
library(igraph)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(data.table)



####################### Full Simulation ########################

sampleSize <- 500
results <- readRDS(paste0("results", sampleSize, "_1ID.rds"))

#Initializing dataframes
performance <- NULL
coeff <- NULL
repNumI <- NULL

#Reading in the results
for (i in 1:20){
  pTemp <- readRDS(file=paste0("performance", sampleSize, "_", i, "ID.rds"))
  performance <- bind_rows(performance, pTemp)
  cTemp <- readRDS(file=paste0("coefficients", sampleSize, "_", i, "ID.rds"))
  coeff <- bind_rows(coeff, cTemp)
  r0Temp <- readRDS(file=paste0("repNumI", sampleSize, "_", i, "ID.rds"))
  repNumI <- bind_rows(repNumI, r0Temp)
}
methods <- unique(performance$label)


#Making clear scenario labels
performance <- (performance
                %>% mutate(Scenario = factor(label,
                                             levels = c("Truth_S1T0.6C1",
                                                        "SNPs_S1T0.6C1",
                                                        "GenInt1.1/2S1T0.6C1",
                                                        "GenInt1.3/3.3S1T0.6C1",
                                                        "GenInt0.5/1.9S1T0.6C1",
                                                        "RandomS1T0.6C1"),
                                             labels = c("Training: Truth",
                                                        "Training: SNP Distance",
                                                        "Correct Serial Interval",
                                                        "Wide Serial Interval",
                                                        "Narrow Serial Interval",
                                                        "Random Probabilities")))
)
results <- (results
            %>% mutate(Scenario = factor(label,
                                         levels = c("Truth_S1T0.6C1",
                                                    "SNPs_S1T0.6C1",
                                                    "GenInt1.1/2S1T0.6C1",
                                                    "GenInt1.3/3.3S1T0.6C1",
                                                    "GenInt0.5/1.9S1T0.6C1",
                                                    "RandomS1T0.6C1"),
                                         labels = c("Training: Truth",
                                                    "Training: SNP Distance",
                                                    "Correct Serial Interval",
                                                    "Wide Serial Interval",
                                                    "Narrow Serial Interval",
                                                    "Random Probabilities")))
)
repNumI <- (repNumI
            %>% unite(uniqueID, runID, label, sep = "_", remove = FALSE)
            %>% filter(!is.na(label))
            %>% mutate(Scenario = factor(label,
                                         levels = c("Truth_S1T0.6C1",
                                                    "SNPs_S1T0.6C1",
                                                    "GenInt1.1/2S1T0.6C1",
                                                    "GenInt1.3/3.3S1T0.6C1",
                                                    "GenInt0.5/1.9S1T0.6C1",
                                                    "RandomS1T0.6C1"),
                                         labels = c("Training: Truth",
                                                    "Training: SNP Distance",
                                                    "Correct Serial Interval",
                                                    "Wide Serial Interval",
                                                    "Narrow Serial Interval",
                                                    "Random Probabilities")))
)

#Creating a long version of performance measures
longData <- (performance
             %>% select(runID, label, Scenario, aucVal, pCorrect,
                        pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -Scenario, -label, -runID)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area Under the ROC", "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")))
)



##################### General Descriptive Statistics ####################

#Information about number of cases per run and number of outbreaks per run
sizes <- (performance
          %>% filter(!duplicated(runID))
          %>% select(runID, nCases, nOutbreaks)
)
summary(sizes$nCases)
summary(sizes$nOutbreaks)

#Find out information about the individual outbreak sizes
outbreakSize <- (repNumI
                 %>% filter(label == "RandomS1T0.6C1")
                 %>% group_by(runID, outbreakID)
                 %>% summarize(n = n())
)
summary(outbreakSize$n)



######################## Results from Run 1 ########################

#### Summary Information ####

#Separating the dataset by scenario
truth <- (results
          %>% filter(label == "Truth_S1T0.6C1")
          %>% mutate(pScaledT = pScaled)
)
snps <- (results
         %>% filter(label == "SNPs_S1T0.6C1")
         %>% mutate(pScaledS = pScaled)
)
int <- (results
        %>% filter(label == "GenInt1.1/2S1T0.6C1")
        %>% mutate(pScaledI = pScaled)
)

#Finding n(%) of pairs that have p < 0.005
table(truth$transmission, truth$pScaled < 0.005)
prop.table(table(truth$transmission, truth$pScaled < 0.005), 1)
table(snps$transmission, snps$pScaled < 0.005)
prop.table(table(snps$transmission, snps$pScaled < 0.005), 1)

#Finding the n(%) of true links that are given a higher probability
#in our method than the serial interval
comb <- (truth
         %>% select(edgeID, pScaledT, transmission)
         %>% full_join(select(snps, edgeID, pScaledS), by = "edgeID")
         %>% full_join(select(int, edgeID, pScaledI), by = "edgeID")
         %>% filter(transmission == TRUE)
)
sum(comb$pScaledT >= comb$pScaledI)
sum(comb$pScaledT >= comb$pScaledI) / nrow(comb)
sum(comb$pScaledS >= comb$pScaledI)
sum(comb$pScaledS >= comb$pScaledI) / nrow(comb)



#### Supplementary Figure: Barplot of Probabilities ####

#Creating a variable that indicates probability group
colBreaks <- c(-0.01, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
results <- results %>% mutate(pGroup = cut(pScaled, breaks = colBreaks,
                                           labels = c("<0.001", "0.001-0.005",
                                                      "0.005-0.01", "0.01-0.05",
                                                      "0.05-0.1", "0.1-0.25",
                                                      "0.25-0.5", "0.5-0.75", ">0.75")),
                              linked = ifelse(transmission == TRUE, "Linked Pairs",
                                              ifelse(transmission == FALSE, "Unlinked Pairs", NA)))

#Adding missing rows (probablity group not represented)
sum(results$label == "GenInt1.3/3.3S1T0.6C1" & results$transmission == TRUE)
newRows1 <- as.data.frame(list("Wide Serial Interval", "Linked Pairs", "<0.001", 0, 576, 0))
names(newRows1) <- c("Scenario", "linked", "pGroup", "n", "nLinked", "p")

sum(results$label == "Truth_S1T0.6C1" & results$transmission == FALSE)
newRows2 <- as.data.frame(list("Training: Truth", "Unlinked Pairs", ">0.75", 0, 173175, 0))
names(newRows2) <- c("Scenario", "linked", "pGroup", "n", "nLinked", "p")
newRows <- bind_rows(newRows1, newRows2)


counts <- (results
           %>% group_by(Scenario, linked)
           %>% summarize(nLinked = n())
           %>% left_join(results, by = c("Scenario", "linked"))
           %>% group_by(Scenario, linked, pGroup)
           %>% summarize(n = n(),
                         nLinked = first(nLinked),
                         p = 100 * n/nLinked)
           %>% bind_rows(newRows)
           %>% ungroup()
           %>% mutate(pGroup = factor(pGroup, levels = c("<0.001", "0.001-0.005",
                                                         "0.005-0.01", "0.01-0.05",
                                                         "0.05-0.1", "0.1-0.25",
                                                         "0.25-0.5", "0.5-0.75", ">0.75")),
                      Scenario = factor(Scenario, levels = c("Training: Truth", "Training: SNP Distance",
                                                   "Correct Serial Interval", "Wide Serial Interval",
                                                   "Narrow Serial Interval", "Random Probabilities")))
)

ggplot(data = counts) +
  geom_bar(aes(x = pGroup, y = p, fill = linked),
           position = position_dodge(), stat = "identity") +
  facet_wrap(~ Scenario) +
  labs(x = "Relative Transmission Probability", y = "Percentage of Pairs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(), legend.position = "bottom",
        legend.spacing.x = unit(0.25, "cm")) +
  ggsave(file = "../Figures/Probabilities.png",
         width = 7, height = 5, units = "in", dpi = 300)





#### Supplementary Figure: True Transmission Networks ####

#Getting the outbreakIDs for the first run
outbreakID1 <- (results
                %>% ungroup()
                %>% filter(!duplicated(individualID.1))
                %>% select(individualID = individualID.1, outbreakID = outbreakID.1)
)
outbreakID2 <- (results
                %>% ungroup()
                %>% filter(!duplicated(individualID.2))
                %>% select(individualID = individualID.2, outbreakID = outbreakID.2)
                %>% filter(!individualID %in% outbreakID1$individualID)
)
outbreakIDdf <- rbind(outbreakID1, outbreakID2)


#Extracting infection date for all other cases and creating nodes list
nodes <- (repNumI
          %>% filter(runID == "1_1")
          %>% group_by(individualID)
          %>% slice(1)
          %>% arrange(observationDate)
          %>% full_join(outbreakIDdf, by = "individualID")
          %>% select(individualID, outbreakID = outbreakID.y, observationDate)
)

#Colors for vertex based on outbreak
vColors <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),
             brewer.pal(12, "Set3"))

#Creating a variable that indicates probability group:
#1: <0.001, 2: 0.001-0.005, 3: 0.005-0.01, 4: 0.01-0.05,
#5: 0.05-0.1, 6: 0.1-0.25, 7: 0.25-0.5, 8: 0.5-0.75, 9: >0.75
#Breaks for colors (9 shades)
colBreaks <- c(-0.01, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)


## Random Network ##
edges.rand <- (results
               %>% ungroup()
               %>% filter(label == "RandomS1T0.6C1")
               %>% select(individualID.1, individualID.2, pScaled, transmission, sameOutbreak)
               #Arrange so true edges and high probability edges are drawn first
               %>% arrange(transmission, pScaled)
)
#Creating network
net.rand <- graph_from_data_frame(d = edges.rand, vertices = nodes, directed = T)
#Creating an indicator variable for transmission: 1=not linked, 2=linked
E(net.rand)$transID <- E(net.rand)$transmission + 1
#Creating a variable that indicates probability group:
E(net.rand)$pGroup <- cut(E(net.rand)$pScaled, breaks = colBreaks, labels = 1:9)
net.rand <- delete.edges(net.rand, E(net.rand)[transmission == FALSE])

## Serial Interval ##
edges.gen <- (results
              %>% ungroup()
              %>% filter(label == "GenInt1.1/2S1T0.6C1")
              %>% select(individualID.1, individualID.2, pScaled, transmission, sameOutbreak)
              #Arrange so true edges and high probability edges are drawn first
              %>% arrange(transmission, pScaled)
)
net.gen <- graph_from_data_frame(d = edges.gen, vertices = nodes, directed = T)
E(net.gen)$transID <- E(net.gen)$transmission + 1
E(net.gen)$pGroup <- cut(E(net.gen)$pScaled, breaks = colBreaks, labels = 1:9)
net.gen <- delete.edges(net.gen, E(net.gen)[transmission == FALSE])

## Training: SNP Distance Network ##
edges.snps <- (results
               %>% ungroup()
               %>% filter(label == "SNPs_S1T0.6C1")
               %>% select(individualID.1, individualID.2, pScaled, transmission, sameOutbreak)
               #Arrange so true edges and high probability edges are drawn first
               %>% arrange(transmission, pScaled)
)
net.snps <- graph_from_data_frame(d = edges.snps, vertices = nodes, directed = T)
E(net.snps)$transID <- E(net.snps)$transmission + 1
E(net.snps)$pGroup <- cut(E(net.snps)$pScaled, breaks = colBreaks, labels = 1:9)
net.snps <- delete.edges(net.snps, E(net.snps)[transmission == FALSE])

## Training: Truth Network ##
edges.truth <- (results
                %>% ungroup()
                %>% filter(label == "Truth_S1T0.6C1")
                %>% select(individualID.1, individualID.2, pScaled, transmission, sameOutbreak)
                #Arrange so true edges and high probability edges are drawn first
                %>% arrange(transmission, pScaled)
)
net.truth <- graph_from_data_frame(d = edges.truth, vertices = nodes, directed = T)
E(net.truth)$pGroup <- cut(E(net.truth)$pScaled, breaks = colBreaks, labels = 1:9)
net.truth <- delete.edges(net.truth, E(net.truth)[transmission == FALSE])

set.seed(20)
#Using Fruchterman-Reingold layout
l <- layout.fruchterman.reingold(net.truth)


png("../Figures/NetworkRaw.png", width = 8, height = 8,
     units = "in", res = 300)
par(mar = c(0, 0, 0.2, 0), mfrow = c(2, 2))
#Random Network
plot(net.rand, vertex.size = 2.5, vertex.label = NA,
     vertex.color = vColors[V(net.rand)$outbreakID],
     vertex.frame.color = vColors[V(net.rand)$outbreakID],
     edge.width = 2, layout = l, edge.arrow.size = 0.2,
     edge.color = brewer.pal(9,"Blues")[E(net.rand)$pGroup])
#Serial Interval
plot(net.gen, vertex.size = 2.5, vertex.label = NA,
     vertex.color = vColors[V(net.gen)$outbreakID],
     vertex.frame.color = vColors[V(net.gen)$outbreakID],
     edge.width = 2, layout = l, edge.arrow.size = 0.2,
     edge.color = brewer.pal(9,"Blues")[E(net.gen)$pGroup])
#SNPs Network with time
plot(net.snps, vertex.size = 2.5, vertex.label = NA,
     vertex.color = vColors[V(net.snps)$outbreakID],
     vertex.frame.color = vColors[V(net.snps)$outbreakID],
     edge.width = 2, layout = l, edge.arrow.size = 0.2,
     edge.color = brewer.pal(9,"Blues")[E(net.snps)$pGroup])
#Truth Network
plot(net.truth, vertex.size = 2.5, vertex.label = NA,
     vertex.color = vColors[V(net.truth)$outbreakID],
     vertex.frame.color = vColors[V(net.truth)$outbreakID],
     edge.width = 2, layout = l, edge.arrow.size = 0.2,
     edge.color = brewer.pal(9,"Blues")[E(net.truth)$pGroup])
dev.off()




################## Assessment of Performance Metrics #######################


#### Supplementary Table: Performance Metrics ####

performTable <- (longData
                 %>% group_by(Scenario, metric)
                 %>% summarize(mean = mean(value),
                               sd = sd(value))
                 %>% mutate(meanSD = paste0(100 * round(mean, 3),
                                            "% (", 100 * round(sd, 3), ")"))
                 %>% select(Scenario, metric, meanSD)
                 %>% spread(metric, meanSD)
)
as.data.frame(performTable)


#### Figure: Violin Plot of Performance Metrics ####

ggplot(data = longData %>% filter(!metric %in% c("Sensitivity", "Specificity")),
       aes(x = Scenario, y = value, fill = Scenario, color = Scenario)) +
  geom_violin(alpha = 0.75) +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  ggsave(file = "../Figures/Metrics.tiff",
         width = 8, height = 5, units = "in", dpi = 300)

## Legend instead of axis labels ##
ggplot(data = longData, aes(x = Scenario, y = value, fill = Scenario, color = Scenario)) +
  geom_violin(alpha = 0.75) +
  scale_x_discrete(name = "Scenario", breaks = NULL) +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.25, "cm")) +
  ggsave(file = "../Figures/Metrics.png",
         width = 8, height = 5, units = "in", dpi = 300)

## BLACK AND WHITE VERSION ##
ggplot(data = longData %>% filter(!metric %in% c("Sensitivity", "Specificity")),
       aes(x = Scenario, y = value)) +
  geom_violin(alpha = 0.75, fill = "grey") +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  ggsave(file = "../Figures/Metrics_BW.png",
         width = 8, height = 5, units = "in", dpi = 300)






################## Reproductive Number Estimation #################

off.r <- 1.2

#Calculating the average monthly Rt by run by month
repNumI <- as.data.table(repNumI)
repNumM <- repNumI[, .(Rt = mean(Ri, na.rm = TRUE),
                     month = first(month),
                     label = first(label),
                     Scenario = first(Scenario),
                     runID = first(runID)),
                  by=.(uniqueID, monthR)]

#Finding the first and last months for each outbreak
#Then adding that information back to the month data frame
totals <- repNumM[, .(totalTime = max(monthR, na.rm = TRUE) - min(monthR, na.rm = TRUE)),
                   by=runID]

repNumM2 <- (repNumM
             %>% full_join(totals, by = "runID")
             %>% mutate(monthCut1 = ceiling(0.1 * totalTime),
                        monthCut2 = ceiling(0.7 * totalTime))
)

#Calculating the average monthly Rt by run by month
repNumCutM <- repNumM2  %>% filter(monthR > monthCut1 & monthR < monthCut2)


#Calculating the average RtAvg for each method
monthR0 <- as.data.table(repNumCutM)[, .(R0 = mean(Rt, na.rm = TRUE),
                                         sdR0 = sd(Rt, na.rm = TRUE),
                                         label = first(label),
                                         Scenario = first(Scenario)),
                                     by=uniqueID]


#### Supplementary Table: Average Rt ####
monthSummaryR0 <- (monthR0
                   %>% group_by(Scenario)
                   %>% summarize(R0avg = mean(R0),
                                  R0sd = sd(R0))
)
monthSummaryR0


#### Figure: Violin Plot of Average Rt ####

ggplot(data = monthR0) +
  geom_violin(aes(x = Scenario, y = R0, fill = Scenario, color = Scenario),
              alpha = 0.5, draw_quantiles = 0.5) +
  scale_x_discrete(name = "Scenario") +
  scale_y_continuous(name = "Average Effective Reproductive Number") +
  geom_hline(aes(yintercept = off.r), linetype = 2, size = 1) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        legend.position = "none") +
  ggsave(file = "../Figures/RepNum.tiff",
         width = 5, height = 4, units = "in", dpi = 300)


## BLACK AND WHITE VERSION ##
ggplot(data = monthR0,
       aes(x = Scenario, y = R0)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, fill = "grey") +
  scale_x_discrete(name = "Scenario") +
  scale_y_continuous(name = "Average Effective Reproductive Number") +
  geom_hline(aes(yintercept = off.r), linetype = 2, size = 1) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        legend.position = "none") +
  ggsave(file = "../Figures/RepNum_BW.png",
         width = 5, height = 4, units = "in", dpi = 300)






###################### Training Proportion/Sample Size Sensitivity Analysis ######################


#Initializing dataframes with first run
performanceSS <- NULL
repNumISS <- NULL

sampleSizes <- c(50, 250, 500)

#Appending the rest of the runs (different sample sizes and files)
for (s in sampleSizes){
  for (i in 1:10){
    pTemp <- readRDS(file=paste0("performanceSS", s, "_", i, ".rds"))
    pTemp <- pTemp %>% mutate(runID2 = paste(runID, s, sep = "_"))
    performanceSS <- bind_rows(performanceSS, pTemp)
    r0Temp <- readRDS(file=paste0("repNumISS", s, "_", i, ".rds"))
    r0Temp <- r0Temp %>% mutate(runID2 = paste(runID, s, sep = "_"))
    repNumISS <- bind_rows(repNumISS, r0Temp)
  }
}


#Creating a categorical variable to group the training dataset sizes
performanceSS <- (performanceSS
                   %>% mutate(nCasesR = ifelse(nCases < 100, "< 100 Cases",
                                        ifelse(nCases >= 100 & nCases < 200, "100-199 Cases",
                                        ifelse(nCases >= 200 & nCases < 300, "200-299 Cases",
                                        ifelse(nCases >= 300 & nCases < 400, "300-399 Cases",
                                        ifelse(nCases >= 400 & nCases < 500, "400-499 Cases",
                                        ifelse(nCases >= 500 & nCases < 600, "500-599 Cases",
                                                 "600+ Cases")))))),
                              nCasesR = factor(nCasesR, levels = c("< 100 Cases", "100-199 Cases",
                                                                   "200-299 Cases", "300-399 Cases",
                                                                   "400-499 Cases", "500-599 Cases",
                                                                   "600+ Cases")))
                   %>% filter(!is.na(goldStandard))
)
repNumISS <- (repNumISS
               %>% unite(label2, goldStandard, trainingP, remove = FALSE)
               %>% unite(uniqueID, label2, runID2, sep = "_", remove = FALSE)
               %>% mutate(nCasesR = ifelse(nCases < 100, "< 100 Cases",
                                    ifelse(nCases >= 100 & nCases < 200, "100-199 Cases",
                                    ifelse(nCases >= 200 & nCases < 300, "200-299 Cases",
                                    ifelse(nCases >= 300 & nCases < 400, "300-399 Cases",
                                    ifelse(nCases >= 400 & nCases < 500, "400-499 Cases",
                                    ifelse(nCases >= 500 & nCases < 600, "500-599 Cases",
                                            "600+ Cases")))))),
                          nCasesR = factor(nCasesR, levels = c("< 100 Cases", "100-199 Cases",
                                                               "200-299 Cases", "300-399 Cases",
                                                               "400-499 Cases", "500-599 Cases",
                                                               "600+ Cases")))
               %>% filter(!is.na(goldStandard))
)

#Creating a long version of the performance dataset
longDataSS <- (performanceSS
             %>% select(label, trainingP, goldStandard, aucVal, nCases, nCasesR,
                        pCorrect, pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -label, -trainingP, -goldStandard, -nCases, -nCasesR)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area Under the ROC",
                                                   "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")))
)



#### Supplementary Figure: Boxplot of Metrics by Training Proportion ####

ggplot(data = longDataSS, aes(x = factor(trainingP), y = value,
                            fill = goldStandard, color = goldStandard)) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = "bottom") +
  labs(y = "Value", x = "Proportion of Cases in Training Dataset") +
  scale_fill_discrete(name = "Training Set") +
  scale_color_discrete(name = "Training Set") +
  ggsave(file = "../Figures/SampleSize_MetricsP.png",
         width = 8, height = 5, units = "in", dpi = 300)



#### Figure: Boxplots of Metrics by Training Proportion and Size ####

lines <- cbind.data.frame(metric = c("Area Under the ROC", "Proportion in Top 25%"),
                          lines = c(0.9, 0.9))

ggplot(data = longDataSS %>% filter(goldStandard == "SNP Distance",
                                  metric %in% c("Area Under the ROC", "Proportion in Top 25%")),
       aes(x = factor(trainingP), y = value, fill = trainingP, color = trainingP)) +
  facet_grid(metric ~ nCasesR, scales = "free_y") +
  theme_bw() +
  geom_boxplot(alpha = 0.5) +
  geom_hline(data = lines, aes(yintercept = lines), linetype = 2, size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = "none") +
  labs(y = "Value", x = "Proportion of Cases in Training Dataset") +
  ggsave(file = "../Figures/SampleSize_MetricsN.tiff",
         width = 9, height = 4, units = "in", dpi = 300)

## BLACK AND WHITE VERSION ##
ggplot(data = longDataSS %>% filter(goldStandard == "SNP Distance",
                                  metric %in% c("Area Under the ROC", "Proportion in Top 25%")),
       aes(x = factor(trainingP), y = value)) +
  facet_grid(metric ~ nCasesR, scales = "free_y") +
  theme_bw() +
  geom_boxplot(alpha = 0.5, fill = "grey") +
  geom_hline(aes(yintercept = 0.9), linetype = 2, size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = "none") +
  labs(y = "Value", x = "Proportion of Cases in Training Dataset") +
  ggsave(file = "../Figures/SampleSize_MetricsN_BW.png",
         width = 9, height = 4, units = "in", dpi = 300)




#### Reproductive Number ####

off.r <- 1.2

repNumMSS <- as.data.table(repNumISS)[, .(Rt = mean(Ri, na.rm = TRUE),
                                         runID2 = first(runID2), 
                                         label2 = first(label2),
                                         trainingP = first(trainingP),
                                         nCases = first(nCases),
                                         nCasesR = first(nCasesR),
                                         goldStandard = first(goldStandard)),
                                     by=.(uniqueID, monthR)]

#Finding the first and last months for each outbreak
#Then adding that information back to the month data frame
repNumMSS2 <- (repNumMSS
               %>% group_by(runID2)
               %>% summarize(totalTime = max(monthR) - min(monthR),
                             monthCut1 = ceiling(0.1 * totalTime),
                             monthCut2 = ceiling(0.7 * totalTime))
               %>% full_join(repNumMSS, by = "runID2")
)

#Calculating the average monthly Rt by run by month
repNumCutMSS <- repNumMSS2  %>% filter(monthR > monthCut1 & monthR <  monthCut2)

monthR0SS <- as.data.table(repNumCutMSS)[, .(R0 = mean(Rt, na.rm = TRUE),
                                         sdR0 = sd(Rt, na.rm = TRUE),
                                         runID2 = first(runID2),
                                         label2 = first(label2),
                                         trainingP = first(trainingP),
                                         nCases = first(nCases),
                                         nCasesR = first(nCasesR),
                                         goldStandard = first(goldStandard)),
                                     by=uniqueID]




#### Supplementary Figure: Violin Plot of Average Rt by Training Proportion ####

monthR0SS <- monthR0SS %>% mutate(nCasesR2 = ifelse(nCases < 250, "<250",
                                             ifelse(nCases < 500, "250-499",
                                                ">500")))

ggplot(data = monthR0SS, aes(x = factor(trainingP), y = R0, 
                           color = goldStandard, fill = goldStandard)) +
  theme_bw() +
  geom_violin(alpha = 0.25, draw_quantiles = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom") +
  labs(x = "Proportion of Cases in the Training Dataset",
       y = "Average Effective Reproductive Number") +
  geom_hline(aes(yintercept = off.r), linetype = 2, size = 1) +
  scale_fill_discrete(name = "Training") +
  scale_color_discrete(name = "Training") +
  scale_y_continuous(breaks = seq(0, 3.4, 0.2)) +
  ggsave(file = "../Figures/SampleSize_R0.png",
         width = 6, height = 4, units = "in", dpi = 300)






####################### Sampling Date Sensitivity Analysis ###########################

sampleSize <- 500
resultsSD <- readRDS(paste0("results", sampleSize, "_1SD.rds"))

#Initializing dataframes
performanceSD <- NULL
repNumISD <- NULL

#Reading in the results
for (i in 1:20){
  pTemp <- readRDS(file=paste0("performance", sampleSize, "_", i, "SD.rds"))
  performanceSD <- bind_rows(performanceSD, pTemp)
  r0Temp <- readRDS(file=paste0("repNumI", sampleSize, "_", i, "SD.rds"))
  repNumISD <- bind_rows(repNumISD, r0Temp)
}

performanceSD <- (performanceSD
                %>% mutate(Scenario = factor(label,
                                             levels = c("Truth_S1T0.6C1",
                                                        "SNPs_S1T0.6C1",
                                                        "GenInt1.1/2S1T0.6C1",
                                                        "GenInt1.3/3.3S1T0.6C1",
                                                        "GenInt0.5/1.9S1T0.6C1",
                                                        "RandomS1T0.6C1"),
                                             labels = c("Training: Truth",
                                                        "Training: SNP Distance",
                                                        "Correct Serial Interval",
                                                        "Wide Serial Interval",
                                                        "Narrow Serial Interval",
                                                        "Random Probabilities")))
)
repNumISD <- (repNumISD
            %>% unite(uniqueID, runID, label, sep = "_", remove = FALSE)
            %>% filter(!is.na(label))
            %>% mutate(Scenario = factor(label,
                                         levels = c("Truth_S1T0.6C1",
                                                    "SNPs_S1T0.6C1",
                                                    "GenInt1.1/2S1T0.6C1",
                                                    "GenInt1.3/3.3S1T0.6C1",
                                                    "GenInt0.5/1.9S1T0.6C1",
                                                    "RandomS1T0.6C1"),
                                         labels = c("Training: Truth",
                                                    "Training: SNP Distance",
                                                    "Correct Serial Interval",
                                                    "Wide Serial Interval",
                                                    "Narrow Serial Interval",
                                                    "Random Probabilities")))
)

#Creating a long version of performance measures
longDataSD <- (performanceSD
             %>% select(runID, label, Scenario, aucVal, pCorrect,
                        pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -Scenario, -label, -runID)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area Under the ROC", "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")),
                        dateVar = "Observation Date")
)

longData <- longData %>% mutate(dateVar = "Infection Date")
longDataAll <- bind_rows(longData, longDataSD)



#### Supplementary Figure: Violin Plot of Metrics by Observation Date ####

ggplot(data = longDataAll, aes(x = Scenario, y = value, fill = dateVar, color = dateVar)) +
  geom_violin(alpha = 0.75) +
  scale_x_discrete(name = "Scenario") +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.25, "cm")) +
  ggsave(file = "../Figures/Metrics_SampleDate.png",
         width = 8, height = 6, units = "in", dpi = 300)



#### Table of RtAvg by label ####

#Calculating the average monthly Rt by run by month
repNumISD <- as.data.table(repNumISD)
repNumMSD <- repNumISD[, .(Rt = mean(Ri, na.rm = TRUE),
                           month = first(month),
                           label = first(label),
                           Scenario = first(Scenario),
                           runID = first(runID)),
                       by=.(uniqueID, monthR)]

#Finding the first and last months for each outbreak
#Then adding that information back to the month data frame
totalsSD <- repNumMSD[, .(totalTime = max(monthR, na.rm = TRUE) - min(monthR, na.rm = TRUE)),
                  by=runID]

repNumMSD2 <- (repNumMSD
             %>% full_join(totalsSD, by = "runID")
             %>% mutate(monthCut1 = ceiling(0.1 * totalTime),
                        monthCut2 = ceiling(0.7 * totalTime))
)

#Calculating the average monthly Rt by run by month
repNumCutMSD <- repNumMSD2  %>% filter(monthR > monthCut1 & monthR < monthCut2)

#Calculating the average RtAvg for each method
monthR0SD <- as.data.table(repNumCutMSD)[, .(R0 = mean(Rt, na.rm = TRUE),
                                         sdR0 = sd(Rt, na.rm = TRUE),
                                         label = first(label),
                                         Scenario = first(Scenario)),
                                     by=uniqueID]
monthR0$dateVar <- "Infection Date"
monthR0SD$dateVar <- "Observation Date"

monthR0All <- bind_rows(monthR0, monthR0SD)

monthSummaryR0SD <- (monthR0SD
                   %>% group_by(Scenario)
                   %>% summarize(R0avg = mean(R0),
                                 R0sd = sd(R0))
)
monthSummaryR0SD

