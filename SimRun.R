#Sarah V. Leavitt
#Boston University Dissertation
#Paper 1

################################################################################
# This program runs one iteration of a simulation
################################################################################


simRun <- function(multOutbreaks, observationDate) {

  #Simulate outbreak  
  obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                     w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                     ws.scale = ws.scale, ws.shape = ws.shape, ws.shift = ws.shift,
                     sampleSize = sampleSize, time = time, multOutbreaks = multOutbreaks,
                     length = length, rate = rate)
  indData <- obk[[1]]
  pairData <- obk[[2]]
  print(paste0("Simulated outbreak, n = ", nrow(indData)))
  
  #Simulating covariates
  covar <- simCovariates(indData, pairData, observationDate = observationDate)
  covarPair <- covar[[1]]
  covarInd <- covar[[2]]
  print("Simulated covariates")

  
  #Initializing dataframes to hold results, and coefficients
  rTemp <- NULL
  cTemp <- NULL
  
  #List of sets of covariates to test
  covarList <- list(c("Y1", "Y2", "Y3", "Y4", "timeCat"))
  #covarList <- list(c("Y1", "Y2", "Y3", "Y4"),
  #                  c("Y1", "Y2", "Y3", "Y4", "timeCat"))
  
  #Vector of proportions of cases to use in training dataset 
  trainingSizes <- 0.6
  #Vector of proportion of cases to sample
  samplingProp <- 1
  
  ## Looping over sampling proportions ##
  for(pSampling in samplingProp){
    
    #Sampling pSampling percent of cases
    finalInd <- covarInd %>% sample_frac(pSampling)
    finalPair <- covarPair %>% filter(individualID.1 %in% finalInd$individualID &
                                        individualID.2 %in% finalInd$individualID)
    #Subseting to the pairs with the potential infector observed before the infectee
    covarOrderedPair <- finalPair %>% filter(observationDate.2 > observationDate.1)
    
    ## Looping over covariate sets ##
    for(covarI in 1:length(covarList)){
      
      ## Looping over proportion in the training dataset ##
      for(pTraining in trainingSizes){
        
        #Creating a label that identifies the set of inputs
        label <- paste0("S", pSampling, "T", pTraining, "C", covarI)
        
        
        #### Choosing pTraining cases for training set ####
        
        #Make sure that the training dataset has at least 3 true links (relevant for small sample sizes)
        nLinked <- 0
        nTries <- 0
        while(nLinked < 3){
          #Finding all pairs that can be included in the training dataset (no missing variables)
          trainingID <- (finalInd
                         %>% filter(complete == TRUE, !is.na(sampleDate))
                         %>% sample_frac(pTraining)
                         %>% pull(individualID)
          )
          
          covarOrderedPair <- covarOrderedPair %>% mutate(trainPair = ifelse(individualID.1 %in% trainingID &
                                                                               individualID.2 %in% trainingID,
                                                                             TRUE, FALSE))
          nLinked <- sum(covarOrderedPair$trainPair == TRUE & covarOrderedPair$transmission == TRUE)
          nTries <- nTries + 1
        }
        if(nTries > 1){print(nTries)}
        
        
        
        #### Estimating probabilities ####
        
        #Creating gold standard variables that take into account if the pair should be in the training dataset
        covarOrderedPair <- (covarOrderedPair
                             %>% mutate(snpClose = ifelse(snpDist < thresholds[1], TRUE,
                                                          ifelse(snpDist > thresholds[2], FALSE, NA)),
                                        snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                                        transmissionGS = ifelse(trainPair == TRUE, transmission, NA))
        )
        
        #Training: True Transmission
        res1 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                                  pairIDVar = "edgeID", goldStdVar = "transmissionGS",
                                  covariates = covarList[[covarI]], label = paste0("Truth_", label),
                                  n = 10, m = 1, nReps = 10)
        probs1 <- res1$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
        print("Completed transmission gold standard analysis")

        
        #Training: SNP Distance
        res2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                                  pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                                  covariates = covarList[[covarI]], label = paste0("SNPs_", label),
                                  n = 10, m = 1, nReps = 10)
        probs2 <- res2$probabilities %>% full_join(covarOrderedPair, by = "edgeID")
        print("Completed SNP threshold gold standard analysis")
        
        rTemp <- bind_rows(rTemp, probs1, probs2)
        cTemp <- bind_rows(cTemp, res1[[2]], res2[[2]])
      }
      print(paste0("Completed analysis with sampling proportion ", pSampling, 
                   ", covariate set ", covarI, ", and pTraining = ", pTraining))
    }

  
      
  #### Adding Comparison Methods ####
  
  #Using true generation interval
  Int1 <- performInterval(covarOrderedPair, genShape = w.shape, genScale = w.scale, shift = w.shift,
                          observationDiff = "observationDiff")
  #Using miss-specified generation interval (prior from Didelot et al. 2017)
  Int2 <- performInterval(covarOrderedPair, genShape = 1.3, genScale = 1/0.3, shift = w.shift,
                          observationDiff = "observationDiff")
  #Using another miss-specified generation interval (posterior from Didelot et al. 2017)
  Int3 <- performInterval(covarOrderedPair, genShape = 0.54, genScale = 1/0.54, shift = w.shift,
                          observationDiff = "observationDiff")
  #Randomly assigning probabilities
  rand <- performRandom(covarOrderedPair)
  rComp <- covarOrderedPair %>% full_join(bind_rows(Int1, Int2, Int3, rand), by = "edgeID")
  
  #Combining these tests with the previous methods
  rTemp <- bind_rows(rTemp, rComp)
  }
  
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  rTemp <- (rTemp
            %>% group_by(label, individualID.2)
            %>% arrange(desc(pScaled))
            %>% mutate(pRank = rank(desc(pScaled), ties.method = "min"))
            %>% ungroup()
  )
  
  
  ## Evaluating the performance ##
  pTemp <- (rTemp
            %>% group_by(label)
            %>% do(simEvaluate(., truthVar = "transmission"))
            %>% ungroup()
  )
  
  
  ## Calculating Ri ##
  
  #Calculating individual-level reproductive numbers
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
  
    
  return(list(rTemp, pTemp, r0Temp, cTemp))
} 



