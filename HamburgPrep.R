#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

#################################################################################
# This program reads in and formats the Hamburg data to prepare for analysis
#################################################################################

setwd("~/Boston University/Dissertation/nbPaper1")
rm(list = ls())

library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)
library(reshape2)
library(stringdist)


######################### Reading in the Original Data #############################


#### Reading in the demographic dataset ####

#Key: id for each case (includes month of isolation)
#Nationality: N.D. means not determined/foreingers
#Age: age at infection
#SubstanceAbuse: 1 = abuse of alcohol or other drugs
#HIV: 1 = HIV seropositive, 0 = HIV negative
#Sex: m = male, f = female
#Residence: 1 = permanent residence, 0 = homeless or resident at the bar
#Resistance: 1 = any antibiotic resistance, 0 = no resistance (all 0)
#Milieu: 1 = affiliation to alcohol-conusming milieu/street scene, 0 = no affiliation
#Study: HH = Hanseatic city, Hamburg, SH = Schleswig-Holstein
#Branch: HC = Hamburg clone, UN = Unsuccessful strains

dem <- read.table("../Datasets/HamburgDemographic_corrected.txt",
                  header = TRUE, stringsAsFactors = FALSE)


#### Creating a variable for contact tracing clusters ####

cont <- cbind.data.frame("Key" = c("5687/01", "4576/00", "6056/00", "1078/00",
                                   "4285/98", "11838/97", "10921/97",
                                   "8855/98", "9538/99", "6257/00",
                                   "749/98", "8349/99", "4127/02", "4217/02", "8496/99", "3015/08",
                                   "358/00", "4283/99", "5112/99", "5882/99", "7063/99", "7199/99", "10965/98",
                                   "4606/09", "446/07", "1880/09", "4414/09",
                                   "581/10", "2375/10",
                                   "9965/07", "10471/07"),
                         "Group" = c(rep("G1", 4), rep("G2", 3), rep("G3", 3), rep("G4", 6),
                                     rep("G5", 7), rep("G6", 4), rep("G7", 2), rep("G8", 2)),
                         stringsAsFactors = FALSE)



#### Reading in the genetic dataset ####

#SNPs: String with sequence from all locations that differed across the isolates (85 total)
#IsolationMY: month/year of isolate isolation

gen <- read.table("../Datasets/HamburgGenetics.txt", header = TRUE, stringsAsFactors = FALSE)
gen <- (gen
        %>% select(Key, IsolationMY, SNPs)
        %>% filter(Key != "Ref")
)



################## Creating an individual-level dataset ##################

hamInd <- (dem
            %>% full_join(gen, by = "Key")
            %>% full_join(cont, by = "Key")
            %>% mutate(IsolationMonth = substr(IsolationMY, 1, 2),
                       IsolationYear = ifelse(substr(IsolationMY, 4, 4) == 9,
                                              paste0("19", substr(IsolationMY, 4, 5)),
                                              paste0("20", substr(IsolationMY, 4, 5))),
                       IsolationDate = ymd(paste0(IsolationYear, IsolationMonth, "15")))
            %>% arrange(IsolationDate)
            %>% mutate(individualID = 1:nrow(.),
                       complete = complete.cases(.[, c("Nationality", "Age", "SubstanceAbuse",
                                                             "HIV", "SmearPos", "Residence",
                                                             "Milieu", "Study")]),
                       Nationality = ifelse(Nationality == "Germany/Turkey", "Turkey",
                                     ifelse(Nationality == "Ukraine/Poland", "Poland",
                                            Nationality)),
                       AgeGroup = ifelse(Age < 25, "<25",
                                  ifelse(Age >= 25 & Age < 35, "25-34",
                                  ifelse(Age >= 35 & Age < 45, "35-44",
                                  ifelse(Age >= 45 & Age < 55, "45-54",
                                  ifelse(Age >= 55 & Age < 65, "55-64",   
                                  ifelse(Age >= 65, "65+", NA)))))))
            %>% select(-IsolationMY)
)



################## Creating a pair-level dataset ####################

#Finding the SNP distance for all pairs
snpDistMatrix <- stringdistmatrix(hamInd$SNPs)
snpDistDf <- melt(as.matrix(snpDistMatrix), varnames = c("individualID.1", "individualID.2"))
snpDistDf <- (snpDistDf
              %>% rename(snpDist = value)
              %>% filter(individualID.1 != individualID.2)
)

#Defining threshold that denotes close genetic distance
hamPair <- (snpDistDf
             #Joining with the individual data to get the information about the infector
             #Making ID.1 the infector information and ID.2 the infectee information
             %>% full_join(hamInd, by = c("individualID.1" = "individualID"))
             %>% full_join(hamInd, by = c("individualID.2" = "individualID"),
                           suffix = c(".1", ".2"))
             #Creating time difference variables
             %>% mutate(IsolationDiff = as.numeric(difftime(IsolationDate.2, IsolationDate.1, units = "days")),
                        IsolationDiffY = IsolationDiff / 365)
             #Creating an edgeID
             %>% unite(edgeID, individualID.1, individualID.2, remove = FALSE)
             #Creating pair-level covariate values
             %>% mutate(Nationality = ifelse(Nationality.1 == "Germany" &
                                               Nationality.2 == "Germany", 1,
                                      ifelse(Nationality.1 != "Germany" &
                                               Nationality.2 != "Germany" &
                                               Nationality.1 == Nationality.2, 2,
                                      ifelse(Nationality.1 != Nationality.2 &
                                               (Nationality.1 == "Germany" |
                                                  Nationality.2 == "Germany"), 3, 4))),
                        Nationality = factor(Nationality, levels = c(3, 1, 4, 2),
                                             labels = c("Diff-Germany", "Same-Germany",
                                                        "Diff-Other","Same-Other")),
                        Age = ifelse(AgeGroup.1 == AgeGroup.2, 1, 0),
                        Age = factor(Age, levels = c(0, 1),
                                     labels = c("Different", "Same")),
                        SubstanceAbuse = ifelse(SubstanceAbuse.1 == 1 & SubstanceAbuse.2 == 1, 1,
                                         ifelse(SubstanceAbuse.1 == 0 & SubstanceAbuse.2 == 0, 2, 3)),
                        SubstanceAbuse = factor(SubstanceAbuse, levels = c(3, 2, 1),
                                                labels = c("Different", "Neither", "Both")),
                        HIV = ifelse(HIV.1 == 0, 1, 2),
                        HIV = factor(HIV, levels = c(1, 2), labels = c("InfectorHIV-", "InfectorHIV+")),
                        SmearPos = ifelse(SmearPos.1 == 0, 1, 2),
                        SmearPos = factor(SmearPos, levels = c(1, 2),
                                          labels = c("InfectorSmear-", "InfectorSmear+")),
                        Sex = ifelse(Sex.1 == "m" & Sex.2 == "m", 1,
                              ifelse(Sex.1 == "f" & Sex.2 == "f", 2,
                              ifelse(Sex.1 == "m" & Sex.2 == "f", 3, 4))),
                        Sex = factor(Sex, levels = c(2, 1, 3, 4),
                                     labels = c("f-f", "m-m", "m-f", "f-m")),
                        Residence = ifelse(Residence.1 == 1 & Residence.2 == 1, 1,
                                    ifelse(Residence.1 == 0 & Residence.2 == 0, 2, 3)),
                        Residence = factor(Residence, levels = c(3, 1, 2),
                                           labels = c("Different", "BothStable", "BothHomeless")),
                        Milieu = ifelse(Milieu.1 == 1 & Milieu.2 == 1, 1,
                                 ifelse(Milieu.1 == 0 & Milieu.2 == 0, 2, 3)),
                        Milieu = factor(Milieu, levels = c(3, 2, 1),
                                        labels = c("Different", "BothNot", "BothAssociated")),
                        Study = ifelse(Study.1 == Study.2, 1, 0),
                        Study = factor(Study, levels = c(0, 1), labels = c("Different", "Same")),
                        Branch = ifelse(Branch.1 == Branch.2, 1, 0),
                        Branch = factor(Branch, levels = c(0, 1), labels = c("Different", "Same")),
                        SameGroup = ifelse(!is.na(Group.1) & !is.na(Group.2) & Group.1 == Group.2, TRUE,
                                    ifelse(!is.na(Group.1) & !is.na(Group.2) & Group.1 != Group.2, FALSE, NA)),
                        TimeCat = ifelse(IsolationDiffY <= 1, 1,
                                  ifelse(IsolationDiffY > 1 & IsolationDiffY <= 2, 2,
                                  ifelse(IsolationDiffY > 2 & IsolationDiffY <= 3, 3,
                                  ifelse(IsolationDiffY > 3 & IsolationDiffY <= 4, 4, 5)))),
                        TimeCat = factor(TimeCat, levels = c(1, 2, 3, 4, 5),
                                         labels = c("<=1y", "1-2y", "2-3y", "3-4y", ">4y")),
                        
                        #Second set of covariates that are simply match not
                        Nationality.D = ifelse(Nationality.1 == Nationality.2, 1, 0),
                        Age.D = ifelse(AgeGroup.1 == AgeGroup.2, 1, 0),
                        SubstanceAbuse.D = ifelse(SubstanceAbuse.1 == SubstanceAbuse.2, 1, 0),
                        HIV.D = ifelse(HIV.1 == HIV.2, 1, 0),
                        SmearPos.D = ifelse(SmearPos.1 == SmearPos.2, 1, 0),
                        Sex.D = ifelse(Sex.1 == Sex.2, 1, 0),
                        Residence.D = ifelse(Residence.1 == Residence.2, 1, 0),
                        Milieu.D = ifelse(Milieu.1 == Milieu.2, 1, 0),
                        Study.D = ifelse(Study.1 == Study.2, 1, 0),
                        SameGroup = ifelse(!is.na(Group.1) & !is.na(Group.2) & Group.1 == Group.2, TRUE,
                                           ifelse(!is.na(Group.1) & !is.na(Group.2) & Group.1 != Group.2, FALSE, NA)))
             %>% select(edgeID, individualID.1, individualID.2, IsolationDate.1, IsolationDate.2, IsolationDiff,
                        IsolationDiffY, TimeCat, snpDist, SameGroup, Nationality, Age, SubstanceAbuse, HIV,
                        SmearPos, Sex, Residence, Milieu, Study, Branch, Nationality.D, Age.D, SubstanceAbuse.D,
                        HIV.D, SmearPos.D, Sex.D, Residence.D, Milieu.D, Study.D)
)


## Saving datasets ##
saveRDS(hamInd, "../Datasets/HamburgInd.rds")
saveRDS(hamPair, "../Datasets/HamburgPair.rds")




