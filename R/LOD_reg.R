### Packages

#library(stats)
#library(haven)
#library(dplyr)
#library(mvtnorm)
#library( HI ) # Used with ARMS

#################                     Function
LOD_regression <- function(dataset, frmla, d, var_LOD, invalid_ind,
                              resultsPath, nSamples=250, convergenceCriterion=0.001,
                              betaCensStartIndex, betaCensEndIndex){
  
  finalEstimates <- list()
    
    ######################################################################
    # Create required datasets                                           #
    ######################################################################
    data_remove_invalid <- dataset %>% filter_at(.vars=invalid_ind, all_vars(.==0))
    
    dataset_Y <- model.frame(frmla,data=data_remove_invalid)[1]
    dataset_X <- model.frame(frmla, data_remove_invalid)[-1]
    dataset_YX <- data.frame(cbind(dataset_Y,Intercept=1,dataset_X))
    
    subData <- dataset_YX[complete.cases(dataset_YX),]
    
    Data <- subData 
    sub2Data <- subData
    
    sqrt2LOd <- function(x){
      ifelse(x>0, x/sqrt(2), x*sqrt(2))
    }
    subd <- sqrt2LOd(d)
    
    for(i in 1:length(var_LOD)){
      Data[[var_LOD[i]]] <- ifelse(subData[[var_LOD[i]]]>d[i], 
                                   subData[[var_LOD[i]]], NA)
      
      sub2Data[[var_LOD[i]]] <- ifelse(subData[[var_LOD[i]]]>d[i], 
                                       subData[[var_LOD[i]]], subd[i])
    }
    
    # Create data for analysis
    ## Obs (Data), complete case (ccData), sub (subData), and sub sqrt2 (sub2Data) datasets
    ccData <- Data[complete.cases(Data),]
    
    ######################################################################
    # Enter convergence criterion, LOD, etc.                             #
    ######################################################################
    
    n <- dim(Data)[1]
    nObservations <- n
    
    #######################################################################
    # Perform Complete-Case, Substitution Analysis                        #
    #######################################################################
    
    ccModel <- glm( frmla,
                    family = gaussian(), data = ccData )
    summary( ccModel )
    
    ######################################################################
    # Perform all Substitution Analyses                                  #
    ######################################################################
    
    # Using LOD for sub
    ccSub_LOD <- glm( frmla,
                      family = gaussian(), data = subData )
    summary( ccSub_LOD )
    
    # Using LOD/sqrt(2) for sub
    ccSub_LODsqrt2 <- glm( frmla,
                           family = gaussian(), data = sub2Data )
    summary( ccSub_LODsqrt2 )
    
    ######################################################################
    # Obtain parameter estimates to be used as initial estimates in ARMS #
    ######################################################################
    
    # Extract Beta and residual variance
    BetaEstimatesCC <- as.numeric( summary( ccModel )$coefficients[,1])
    names(BetaEstimatesCC) <- rownames(summary( ccModel )$coefficients)
    BetaEstimatesSubsqrt2 <- as.numeric( summary( ccSub_LODsqrt2 )$coefficients[,1])
    names(BetaEstimatesSubsqrt2) <- rownames(summary( ccSub_LODsqrt2 )$coefficients)
    
    # Extract mean vector estimate and covariance matrix estimate of covariates
    cat_var <- names(ccModel$contrasts)
    length_unique <- function(x){length(unique(x))}
    unique_values <- apply(dataset_X, MARGIN=2,FUN=length_unique)
    binary_vars <- names(unique_values[unique_values<3])
    remove_vars <- c(var_LOD, cat_var, binary_vars)
    var_noLOD <- dataset_X %>% select(-one_of(remove_vars)) %>% names()
    
    xMeanCC <- apply( ccData[,c(var_noLOD, var_LOD)], 2, mean ) #-1 to eliminate outcome Y
    xVarCC <- apply( ccData[,c(var_noLOD, var_LOD)], 2, var )
    xCovCC <- cov( ccData[,c(var_noLOD, var_LOD)] )
    
    #######################################################################
    # Perform ARMS MLE Sampling                                           #
    #######################################################################
    
    ######################################################################
    # Make the log-likelihood function                                   #
    ######################################################################
    
    ### This is to sample the missing covariate ###
    
    LogDensity <- function( x ) {
      
      Ycur <- currentY
      currentBetas <- betaCurrent
      currentCovariates[covariateTargetIndex+1] <- x
      
      # Component 1 - (y-X'B)^2/sigma2
      xBeta <- sum( currentBetas[1:betaCensEndIndex]*currentCovariates[c("Intercept",var_noLOD,var_LOD)] )
      y_diff_XBeta <- (Ycur - xBeta)^2
      firstComponent <- y_diff_XBeta/sigma2Current
      
      # Component 3 - (x-mu)'Sigma^(-1)(x-mu)
      targetX <- currentCovariates[c(var_noLOD,var_LOD)]
      targetMeanX <- currentXMeans
      xMu <- targetX - targetMeanX
      thirdComponent <- t(unlist(xMu)) %*% solve(currentSigmaX) %*% unlist(xMu)
      
      finalResult <- -0.5*firstComponent - 0.5*thirdComponent
      
      return( finalResult )
      
    }
    
    ######################################################################
    # Perform the Sampling                                               #
    ######################################################################
    
    
    
    observedDataXY <- ccData
    observedData <- observedDataXY[,-1]
    missingcnt <- if(length(var_LOD)>1){apply(is.na(Data[,var_LOD]),1,sum)}else 
      if(length(var_LOD==1)){as.numeric(is.na(Data[,var_LOD]))}else{rep(0, dim(Data)[1])}
    censoredDataXY <- Data[missingcnt>0,]
    censoredData <- Data[missingcnt>0,-1]
    yObserved <- observedDataXY[,1]
    yCensored <- censoredDataXY[,1]
    nObserved <- length( yObserved )
    nCensored <- length( yCensored )
    
    subDataOrdered <- subData
    subDataCensored <- sub2Data[missingcnt>0,]
    
    subPts <- c( rep(0,length(var_noLOD)), subd )
    
    ### Make data structures for output ###
    betaEstimatesByIteration <- NA
    betaSEByIteration <- NA
    xMeansByIteration <- NA
    xSigmaByIteration <- NA
    
    ### Apply Initial Parameter Estimates ###
    currentXMeans <- xMeanCC
    currentXVars <- xVarCC
    betaCurrent <- if(ccModel$df.residual>0 & length(BetaEstimatesCC)==length(BetaEstimatesSubsqrt2)){BetaEstimatesCC}else{BetaEstimatesSubsqrt2}
    nBetas <- length(betaCurrent)
    sigma2Current <- if(ccModel$df.residual>0 & sigma(ccModel)^2>0.01){sigma(ccModel)^2}else{sigma(ccSub_LODsqrt2)^2}
    currentSigmaX <- if(abs(det(xCovCC))>0.0001){xCovCC}else{xCovSub2}
    
    betaConverged <- rep(FALSE,nBetas)
    converged <- FALSE
    
    iterationNumber <- 0
    
    while ( !converged ) {
      
      iterationNumber <- iterationNumber + 1
      
      print( paste( "Iteration", iterationNumber ) )
      
      finalSampledData <- NA
      
      ### E-Step ###
      
      for ( j in 1:nCensored ) {
        print( paste( "Censored Observation", j, "of", nCensored ) )
        
        # Count the number of censored variables in this row
        nCensIndiv <- 0
        censIndexes <- c()
        for ( k in betaCensStartIndex:betaCensEndIndex ) {
          if ( is.na( censoredData[j,k] ) ) {
            nCensIndiv <- nCensIndiv + 1
            censIndexes <- c(censIndexes,k) #Stores indexes of censored observations
          }
        }
        
        currentCovariates <- subDataCensored[j,-1] #USe the substituted data as starting values
        
        currentY <- yCensored[j]
        
        newDataFrame <- matrix( rep( censoredData[j,], nSamples ), ncol = length(censoredData[j,]), byrow = T )
        
        if ( nCensIndiv == 1 ) { #Sampling is simpler, can do all in one step
          
          covariateTargetIndex <- censIndexes[1]-1
          
          supportFunction <- function(x){
            ifelse(-100<x & x<d[covariateTargetIndex-length(var_noLOD)],1,0)
          }
          currentSample <- arms( subPts[covariateTargetIndex], LogDensity, supportFunction, nSamples )
          
          newDataFrame[,covariateTargetIndex+1] <- currentSample
          
        } else { #Have more than 1 censored covariate, use gibbs sampling
          
          newDataFrame <- matrix( NA, ncol = length(censoredData[j,]), nrow = nSamples )
          
          for ( sampleIndex in 1:nSamples ) {
            
            for ( censIndex in 1:nCensIndiv ) {
              
              covariateTargetIndex <- censIndexes[censIndex]-1
              supportFunction <- function(x){
                ifelse(-100<x & x<d[covariateTargetIndex-length(var_noLOD)],1,0)
              }
              currentSample <- arms( subPts[covariateTargetIndex], LogDensity, supportFunction, 1 )
              
              currentCovariates[covariateTargetIndex+1] <- currentSample
              
            }
            
            newDataFrame[sampleIndex,] <- unlist(currentCovariates)
            
          }
          
        }
        
        # Add this new data frame to the augmented dataset
        if ( j == 1 ) {
          allSampledData <- cbind( rep(currentY,nSamples), newDataFrame )
        } else {
          allSampledData <- rbind( allSampledData, cbind( rep(currentY,nSamples), newDataFrame ) )
        }
        
      }
      ### M-Step ###
      
      #Add the weights
      totalSamplesTaken <- length( allSampledData[,1] )
      allSampledDataFinal <-  as.data.frame(cbind(allSampledData, rep(1/nSamples,totalSamplesTaken) ))
      allObservedDataFinal <- cbind( as.numeric(yObserved), observedData, rep(1,nObserved) )
      names(allSampledDataFinal) <- names(allObservedDataFinal)
      
      allObservedDataFinal[] <- lapply(allObservedDataFinal, as.numeric)  
      allSampledDataFinal[] <- lapply(allSampledDataFinal, as.numeric) 
      
      finalMaximizationData <- rbind( allObservedDataFinal, allSampledDataFinal)
      for(i in 1:dim(finalMaximizationData)[2]){
        finalMaximizationData[,i] <- unlist(finalMaximizationData[,i])
      }
      names( finalMaximizationData ) <- c(names(Data), "DATAWEIGHTS")
      finalMaximizationData <- finalMaximizationData %>% mutate_at(.vars=cat_var, .funs = as.factor)
      #Fit the weighted GLM
      
      currentModel <- glm( frmla,
                           family = gaussian(), weights = DATAWEIGHTS, data = finalMaximizationData )
      
      sumARMS <- summary( currentModel )
      betaCurrent <- as.numeric( sumARMS$coefficients[,1] )
      names(betaCurrent) <- names(sumARMS$coefficients[,1])
      
      #Calc error variance
      X_indic <- model.matrix(currentModel)
      X_obs <- X_indic[finalMaximizationData$DATAWEIGHTS==1,]
      X_cen <- X_indic[finalMaximizationData$DATAWEIGHTS!=1,]
      
      yCensored_ARMs <- finalMaximizationData$ATTENTION_SCALE[finalMaximizationData$DATAWEIGHTS!=1]
      yObs_ARMs <- finalMaximizationData$ATTENTION_SCALE[finalMaximizationData$DATAWEIGHTS==1]
      
      Obs_error <- sum((yObserved-X_obs%*%betaCurrent)^2)
      Cens_error <- sum((yCensored_ARMs-X_cen%*%betaCurrent)^2)/nSamples
      sigma2Current <- (Obs_error+Cens_error)/(nObservations-1)
      
      ### Calculate X-Distribution Estimates ###
      
      # Weighted mean
      for ( cCovIndex in 3:(betaCensEndIndex)+1 ) {
        currentXMeans[cCovIndex-2] <- (sum(finalMaximizationData[,cCovIndex]*finalMaximizationData$DATAWEIGHTS, na.rm=TRUE))/nObservations
      }
      
      # Weighted variance
      for ( cCovIndex in 3:(betaCensEndIndex+1) ) {
        currentCovariate <- finalMaximizationData[,cCovIndex]
        currentMean <- currentXMeans[cCovIndex-2]
        currentSS <- finalMaximizationData$DATAWEIGHTS * (currentCovariate - currentMean) * (currentCovariate - currentMean)
        currentVar <- sum( currentSS, na.rm=TRUE ) / nObservations
        currentXVars[cCovIndex-2] <- currentVar
      }
      currentXVars <- round( currentXVars, 7 )
      
      xCensData <- finalMaximizationData[,c(var_noLOD,var_LOD)]
      xCensDataMatrix <- matrix( as.numeric( unlist( xCensData ) ), ncol = dim(xCensData)[2] )
      colnames(xCensDataMatrix) <- names(xCensData)
      currentSigmaX <- cov.wt( xCensDataMatrix[complete.cases(xCensDataMatrix),], 
                               finalMaximizationData$DATAWEIGHTS[complete.cases(xCensDataMatrix)] / nObservations)$cov
      
      ### Record the new estimates ###
      
      if ( iterationNumber == 1 ) {
        betaEstimatesByIteration <- betaCurrent
        xMeansByIteration <- currentXMeans
        xSigmaByIteration <- as.numeric( currentSigmaX )
      } else {
        betaEstimatesByIteration <- rbind( betaEstimatesByIteration, betaCurrent )
        xMeansByIteration <- rbind( xMeansByIteration, currentXMeans )
        xSigmaByIteration <- rbind( xSigmaByIteration, as.numeric( currentSigmaX ) ) 
      }
      ### Check for convergence of beta estimates ###
      
      if ( iterationNumber >= 20 ) {
        for ( betaIndex in 1:nBetas ) {
          if ( !betaConverged[betaIndex] ) { #If this beta hasn't converged yet
            mean1 <- mean( betaEstimatesByIteration[(iterationNumber-19):(iterationNumber-10),betaIndex] )
            mean2 <- mean( betaEstimatesByIteration[(iterationNumber-9):iterationNumber,betaIndex] )
            if ( abs(mean1-mean2) < convergenceCriterion ) {
              betaConverged[betaIndex] <- TRUE
              print( paste( "Beta", betaIndex, "Converged" ) )
            }
          }
        }
        converged <- TRUE
        for ( betaIndex in 1:nBetas ) {
          if ( !betaConverged[betaIndex] ) {
            converged <- FALSE
          }
        }
      }
      
      if ( converged ) { #Write out all the estimates 
        lastTenIterations <- betaEstimatesByIteration[(iterationNumber-9):iterationNumber,]
        finalEstimates <- apply( lastTenIterations, 2, mean )
        
      }
      
    }       
  
  save(finalEstimates, file=paste(resultsPath,"finalEstimates.Rdata", sep="/"))
  #results <- matrix(unlist(finalEstimates), nrow=length(finalEstimates), byrow=TRUE)
  #SEs <- apply(results, MARGIN = 2, FUN=sd)
  #write.table( SEs paste( resultsPath, "BetaSeEsts_Boot.txt", sep = "/" ),
  #sep = "\t", quote = F, row.names = F, col.names = F )

}
##############################################################
