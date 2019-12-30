#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program runs one iteration of a simulation to assess the effect of
# training proportion and sample size
################################################################################


sampleSizeRun <- function(){

  #Simulate outbreak  
  obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                     w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                     ws.scale = ws.scale, ws.shape = ws.shape, ws.shift = ws.shift,
                     sampleSize = sampleSize, time = time,
                     multOutbreaks = multOutbreaks,
                     length = length, rate = rate)
  indData <- obk[[1]]
  pairData <- obk[[2]]
  print(paste0("Simulated outbreak, n = ", nrow(indData)))
  
  #Simulating covariates
  covar <- simCovariates(indData, pairData)
  covarPair <- covar[[1]]
  covarInd <- covar[[2]]
  print("Simulated covariates")
  
  #Subseting to the pairs with the potential infector observed before the infectee
  covarOrderedPair <- (covarPair
                       %>% filter(observationDiff > 0)
                       %>% mutate(snpClose = ifelse(snpDist < thresholds[1], TRUE,
                                             ifelse(snpDist > thresholds[2], FALSE, NA)))
  )
  
  #Initializing dataframes to hold results
  rTemp <- NULL
  
  #Vector of proportions of cases to use in training dataset 
  trainingSizes <- seq(0.1, 1, by = 0.1)
  #trainingSizes <- 0.5
  
  for(pTraining in trainingSizes){
    
    label <- paste0("T", pTraining)
    
    #### Choosing pTraining cases for training set ####
    
    #Make sure that the training dataset has at least 3 true links (relevant for small sample sizes)
    nLinked <- 0
    nTries <- 0
    while(nLinked < 3){
      #Finding all pairs that can be included in the training dataset (no missing variables)
      trainingID <- (covarInd
                     %>% filter(complete == TRUE, !is.na(sampleDate))
                     %>% sample_frac(pTraining)
                     %>% pull(individualID)
      )
      
      covarOrderedPair <- covarOrderedPair %>% mutate(trainPair = ifelse(individualID.1 %in% trainingID &
                                                                           individualID.2 %in% trainingID,
                                                                         TRUE, FALSE))
      
      nLinkedSNPs <- sum(covarOrderedPair$trainPair == TRUE &
                           covarOrderedPair$snpClose == TRUE, na.rm = TRUE)
      nLinkedTruth <- sum(covarOrderedPair$trainPair == TRUE &
                            covarOrderedPair$transmission == TRUE, na.rm = TRUE)
      nLinked <- min(nLinkedSNPs, nLinkedTruth)
      
      nTries <- nTries + 1
    }
    if(nTries > 1){print(nTries)}
    
    
    
    #### Calculating probabilities ####
    
    #Creating new gold standard variables that take into account if the pair should be
    #in the training dataset (only relevant in simulations with full information)
    covarOrderedPair <- covarOrderedPair %>% mutate(snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                                                    transmissionGS = ifelse(trainPair == TRUE, transmission, NA))
    
    #Gold Std: True Transmission
    res1 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                              pairIDVar = "edgeID", goldStdVar = "transmissionGS",
                              covariates = covariates, label = paste0("Truth_", label),
                              n = 10, m = 1, nReps = 10)
    probs1 <- res1$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed transmission gold standard analysis")
    
    
    #Gold Std: SNP Distance
    res2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                              pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                              covariates = covariates, label = paste0("SNPs_", label),
                              n = 10, m = 1, nReps = 10)
    probs2 <- res2$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
    print("Completed SNP threshold gold standard analysis")
    
    rTemp <- bind_rows(rTemp, probs1, probs2)
    
    print(paste0("Completed analysis with pTraining = ", pTraining))
  }  
  

  
  ################# Assessing Performance ################
  
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  rTemp <- (rTemp
            %>% group_by(label, individualID.2)
            %>% arrange(desc(pScaled))
            %>% mutate(pRank = rank(desc(pScaled), ties.method = "min"))
  )
  
  ## Evaluating the performance ##
  pTemp <- (rTemp
            %>% group_by(label)
            %>% do(simEvaluate(.))
  )
  
  ## Calculating Ri ##
  
  r0Temp <- (rTemp
             %>% group_by(label)
             %>% do(estimateRi(., pVar = "pScaled", indIDVar = "individualID",
                               dateVar = "observationDate"))
             %>% ungroup()
  )
  
  #Adding the month ranks
  months <- format(seq(min(r0Temp$observationDate),
                       max(r0Temp$observationDate), by = "months"), "%Y-%m")
  monthsDf <- cbind.data.frame(month = months, monthR = 1:length(months),
                               stringsAsFactors = FALSE)
  
  #Adding the number of cases each case infected
  nInfected <- (covarPair
                %>% filter(transmission == TRUE)
                %>% group_by(individualID.1)
                %>% summarize(nInfected = n())
                %>% rename(individualID = individualID.1)
                %>% ungroup()
  )
  
  #Combining the year and month rankings with the raw data
  r0Temp <- (r0Temp
             %>% mutate(month = format(observationDate, "%Y-%m"))
             %>% left_join(monthsDf, by = "month")
             %>% left_join(select(covarInd, individualID, outbreakID), by = "individualID")
             %>% left_join(nInfected, by = "individualID")
             %>% replace_na(list(nInfected = 0))
  )
  
  #Extracting the training proportion and gold standard as separate variables
  pTemp2 <- (pTemp
             %>% mutate(trainingP = as.numeric(str_extract(str_extract(label, "T[:digit:]+\\.*[:digit:]*"),
                                                           "[:digit:]+\\.*[:digit:]*")),
                        trainingN = nCases * trainingP,
                        goldStandard = ifelse(grepl("Truth", label), "Truth", 
                                        ifelse(grepl("SNPs", label), "SNP Distance", NA)))
  )
  r0Temp2 <- (r0Temp
              %>% full_join(select(pTemp, label, nCases), by = "label")
              %>% mutate(trainingP = as.numeric(str_extract(str_extract(label, "T[:digit:]+\\.*[:digit:]*"),
                                                            "[:digit:]+\\.*[:digit:]*")),
                         trainingN = nCases * trainingP,
                         goldStandard = ifelse(grepl("Truth", label), "Truth", 
                                               ifelse(grepl("SNPs", label), "SNP Distance", NA)))
  )
  
  return(list(pTemp2, r0Temp2))
}


