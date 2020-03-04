#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:06:42 2020

@author: tom holding
"""

import numpy as np;
import pandas as pd;


#Calculates the exected uncertainty on an in situ data variable (i.e. DIC or AT)
#   using "nominal state of the art" uncertainties for in situ measurements
#   see: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
#Returns an iterable containing an uncertainty for each row in insituData
#insituData: dataframe containing the in situ data
#varoable: target variable from which the (e.g. 'DIC' or 'AT')
#settings: the global setting dictionary
def calc_insitu_uncertainty(insituData, variable, settings):
    if settings["useErrorRatios"]:
        insituUncertainty = (insituData[variable]*settings["insituErrorRatio"][variable])**2;
    else:
        insituUncertainty = [settings["insituError"][variable]**2]*len(insituData);
    return insituUncertainty


###Calculate weights based on normalised inverse size of combined insitu and algorithm uncertainty
def calc_weights(combinedUncertainty):
    weights = 1.0 / combinedUncertainty;
    weights /= weights.sum(); #normalise
    return weights;


#Calculates basic (single algorithm) metrics. Calculates:
#   mean model and reference output
#   sd (standard deviation)
#   rmsd (root mean squared difference)
#   mad (mean absolute difference)
#   r (correlation coefficient)
#Weighted (by estimated uncertainty) and unweighted versions are calculated.
#'Model output' refers to output estimated by the algorithm
#'Reference output' refers to the measured in situ values in the matchup database
#'w' prefix on variables refers to weighted versions
#Arguments as follows:
#   modelOutput: algorithm estimates
#   algorithm: the algorithm functor used to generate modelOutput
#   dataUsed: dataframe containing the matchup database subsetted for the rows used to calculate modelOutput
#   outputVariable: the common name of the output variable (i.e. 'DIC' or 'AT') (see global settings for definitions)
#   settings: the global settings dictionary
def calc_basic_metrics(modelOutput, algorithm, dataUsed, outputVariable, settings):
    ####Calculate basic statistics
    insituOutputUncertainties = calc_insitu_uncertainty(dataUsed, outputVariable, settings);
    
    #Uncertainty combined for in situ output variable and model output variable #Note: Called "squaredErrors" in Peter's 'printFit' function.
    #This is used to calculate weights based on normalised inverse size of combined insitu and algorithm uncertainty
    if np.any(algorithm.rmsd != None):
        combinedUncertainty = algorithm.rmsd + insituOutputUncertainties;
        weights = calc_weights(combinedUncertainty);
        hasWeights = True;
    else:
        weights = np.nan;
        hasWeights = False;
        print("*** No RMSD for", algorithm.__class__.__name__, "Weighted metrics will not be calculated.");
    
    #Calculate prediction errors: the difference between reference measurement and model output
    predictionErrors = modelOutput - dataUsed[outputVariable];
    bias = np.mean(predictionErrors);
    predictionErrorsSquared = predictionErrors**2; #Named 'squaredError' in Peter's 'printFit' function
    absPredictionError = abs(predictionErrors);
    
    #Metrics of ability of the algorithm to reproduce measured values
    predictionRmsd = np.sqrt(predictionErrorsSquared.mean()); #rmse
    predictionMad = absPredictionError.mean(); #mean absolute difference
    if hasWeights:
        wpredictionRmsd = np.sqrt((weights*predictionErrorsSquared).sum()); #weighted sum of the squared prediction error. Note sum not mean because the weights are normalised to sum to 1.
        wpredictionMad = (weights*absPredictionError).sum(); #weighted mean absolute difference. Note sum is used, not mean, as weights are already normalised to sum to 1
    else:
        wpredictionRmsd = np.nan;
        wpredictionMad = np.nan;
    
    #Calculates means of model and reference values (unweighted and weighted)
    meanModelOutput = modelOutput.mean();
    meanReferenceOutput = dataUsed[outputVariable].mean();
    if hasWeights:
        wMeanModelOutput = (weights*modelOutput).sum(); #weighted mean of the model/predicted output. Note sum is used, not mean, as weights are already normalised to sum to 1
        wMeanInSituOutput = (weights*dataUsed[outputVariable]).sum(); #weighted mean of the in situ measured output. Note sum is used, not mean, as weights are already normalised to sum to 1
    else:
        wMeanModelOutput = np.nan;
        wMeanInSituOutput = np.nan;
    
    #Calculating standard deviation
    #Using Peters method: TODO: repeat / check with other method
    if len(dataUsed) > 1: #Can't calculate SD from 1 data point
        covariance = np.ma.cov(dataUsed[outputVariable], y = modelOutput) # covariance matrix
        referenceOutputSD = covariance[0, 0] ** .5
        modelOutputSD = covariance[1, 1] ** .5
        r = covariance[0, 1] / (referenceOutputSD * modelOutputSD)
        
        #referenceOutputSD = np.std(dataUsed[outputVariable], ddof=1); #ddof=1 to calculate sample SD
        #modelOutputSD = np.std(modelOutput, ddof=1); #ddof=1 to calculate sample SD
        #r = np.sum (np.sum(dataUsed[outputVariable]-np.mean(dataUsed[outputVariable])) * np.sum(modelOutput-np.mean(modelOutput)) ) / referenceOutputSD*modelOutputSD;
        
        #weighted standard deviations and correlation coefficient
        if hasWeights:
            sumX2 = (weights * (dataUsed[outputVariable] - wMeanInSituOutput) ** 2).sum()
            sumY2 = (weights * (modelOutput - wMeanModelOutput) ** 2).sum()
            sumXY = ((weights * (dataUsed[outputVariable] - wMeanInSituOutput) * (modelOutput - wMeanModelOutput)).sum())
            try:
                referenceOutputWSD = sumX2 ** .5 #wsd = weighted standard deviation
            except:
                raise ValueError(sumX2)
            modelOutputWSD = sumY2 ** .5 #wsd = weighted standard deviation
        
            with np.errstate(divide='raise'): #Array divide by zero warning should be treated as full errors/exceptions
                wr = sumXY / (referenceOutputSD * modelOutputSD) #weighted r
        else:
            referenceOutputWSD = np.nan;
            modelOutputWSD = np.nan;
            wr = np.nan;
    else: #0 or 1 data point, so no standard deviations can be calculated
        modelOutputSD = np.nan;
        referenceOutputSD = np.nan;
        r = np.nan;
        modelOutputWSD = np.nan;
        referenceOutputWSD = np.nan;
        wr = np.nan;
        
    
    ####Calculate performance metrics (pairwise statistics)
    basicMetrics = {"model_output": modelOutput,
                    "model_output_mean": meanModelOutput,
                    "reference_output_mean": meanReferenceOutput,
                    "model_output_sd": modelOutputSD,
                    "reference_output_sd": referenceOutputSD,
                    "model_output_wsd": modelOutputWSD, #weighted standard deviation
                    "reference_output_wsd": referenceOutputWSD,
                    "rmsd": predictionRmsd,
                    "wrmsd": wpredictionRmsd,
                    "mad": predictionMad,
                    "wmad": wpredictionMad,
                    "r": r,
                    "wr": wr,
                    "weights": weights, #values used for weighted metrics
                    "n": len(dataUsed), #number of matchup rows
                    "bias": bias, #bias = mean difference: mean(model-reference)
                    };
    
    return basicMetrics;


#Calculates all the metrics for a set of algorithms. This includes:
#   basic metrics such as mean output, standard deviation, rmsd, mad and correlation coefficients
#   pairwise metrics: algorithm score, weighted score, intersect size, pairwise rmsd
#Returns a list of dictionaries containing the basic metrics for each given by algorithmFunctorList, and matrices for each pairwise metric
#Arguments as follows: 
#   algorithmFunctorList: list of algorithm functor instances used to calculate the model outputs
#   modelOutputList: list of predictions from the algorithms in algorithmFunctorList
#   matchupRowsUsedList: Defines the subset of the matchup database used by the algorithm
#   matchupData: pandas dataframe containing all the matchup data
#   settings: the global settings dictionary
def calc_all_metrics(algorithmFunctorList, modelOutputList, matchupRowsUsedList, matchupData, settings):
    #Sanity check:
    if len(algorithmFunctorList) != len(modelOutputList) != len(matchupRowsUsedList):
        raise ValueError("Number of algorithms used, number of model outputs provides and number of specified matchup subsets must match!");
    
    #Create some space for storing metrics (pairwise metrics are stored as matrices.)
    basicMetrics = []; #Stores simple non-paired metrics like mean, standard deviation, RMSD and MAD
#    nIntersectMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=int); #Symetrical
#    pairedScoreMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise score (see Land2019 Remote Sensing of Environment, section 2.3.2). Index order is [ith, jth]
#    pairedWScoreMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise weighted score. Index order is [ith, jth]
#    pairedRmsdMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise RMSD. Index order is [ith, jth]
#    pairedWRmsdMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise weighted RMSD. Index order is [ith, jth]
    nIntersectMatrix = np.ma.zeros((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=int); nIntersectMatrix.mask=False; #Symetrical
    pairedScoreMatrix = np.ma.empty((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); pairedScoreMatrix[:]=np.nan; pairedScoreMatrix.mask=False; #Pairwise score (see Land2019 Remote Sensing of Environment, section 2.3.2). Index order is [ith, jth]
    pairedWScoreMatrix = np.ma.empty((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); pairedWScoreMatrix[:]=np.nan; pairedWScoreMatrix.mask=False; #Pairwise weighted score. Index order is [ith, jth]
    pairedRmsdMatrix = np.ma.empty((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); pairedRmsdMatrix[:]=np.nan; pairedRmsdMatrix.mask=False; #Pairwise RMSD. Index order is [ith, jth]
    pairedWRmsdMatrix = np.ma.empty((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); pairedWRmsdMatrix[:]=np.nan; pairedWRmsdMatrix.mask=False; #Pairwise weighted RMSD. Index order is [ith, jth]
    
    #### Calculate pairwise metrics as described in section 2.3.2 of Land, P.E., Findlay, H.S., Shutler, J.D., Ashton, I.G., Holding, T., Grouazel, A., Girard-Ardhuin, F., Reul, N., Piolle, J.F., Chapron, B. and Quilfen, Y., 2019. Optimum satellite remote sensing of the marine carbonate system using empirical algorithms in the global ocean, the Greater Caribbean, the Amazon Plume and the Bay of Bengal. Remote Sensing of Environment, 235, p.111469.
    for i in range(len(algorithmFunctorList)):
        outputVariable = algorithmFunctorList[i].output_name();
        
        #Calculate means, standard deviations, rmse, mad, r
        basicMetrics.append(calc_basic_metrics(modelOutputList[i], algorithmFunctorList[i], matchupData.loc[matchupRowsUsedList[i]], algorithmFunctorList[i].output_name(), settings));
        
        #Calculate pairwise metrics
        for j in range(i): #for each pair of algorithms ((self, self) pair excluded)
            #Sanity check - can only compare algorithms which are predicting the same variable
            if outputVariable != algorithmFunctorList[j].output_name():
                raise ValueError("Cannot compare algorithms which estimate different output variables.");
            
            #Only use data points from the matchup dataset that have been predicted by both algorithms
            intersectingRows = set(matchupRowsUsedList[i]).intersection(set(matchupRowsUsedList[j]));
            if len(intersectingRows) == 0:
                continue;
            
            #subset weights for ith and jth algorithm to give just the weights corresponding to intersecting rows
            if (np.any(np.isfinite(basicMetrics[i]["weights"]))) & (np.any(np.isfinite(basicMetrics[j]["weights"]))):
                weightsi = basicMetrics[i]["weights"].loc[intersectingRows];
                weightsj = basicMetrics[j]["weights"].loc[intersectingRows];
                pairHasWeights = True;
            else:
                pairHasWeights = False; #Indicate to the rest of the function that this pair won't calculate weighted metrics
            
            #Calculate the 'error' (difference between model/predicted and measured reference)
            diffFromReferencei = modelOutputList[i].loc[intersectingRows] - matchupData[outputVariable].loc[intersectingRows]; #Differences between model and reference output for ith algorithm at only intersecting locations
            diffFromReferencej = modelOutputList[j].loc[intersectingRows] - matchupData[outputVariable].loc[intersectingRows]; #Differences between model and reference output for jth algorithm at only intersecting locations;
            #And again with weighted difference
            if pairHasWeights:
                wdiffFromReferencei = weightsi * diffFromReferencei;
                wdiffFromReferencej = weightsj * diffFromReferencej;
            
            #Calculate rmsd of current and ith models using only intersecting datapoints
            rmsdi = ((diffFromReferencei**2).mean())**0.5;
            rmsdj = ((diffFromReferencej**2).mean())**0.5;
            if pairHasWeights:
                wrmsdi = ((wdiffFromReferencei**2).sum() / weightsi.sum())**0.5;
                wrmsdj = ((wdiffFromReferencej**2).sum() / weightsj.sum())**0.5;
            
            #score models in the pair based on their rmse (weighted and unweighted)
            #model with the smaller rmsd gets a score of 1
            #model with larger RMSD get a score of otherRMSD/RMSD (always > 1)
            if rmsdj > rmsdi:
                scorei = 1.0;
                scorej = rmsdj/rmsdi;
            else:
                scorei = rmsdi/rmsdj;
                scorej = 1.0;
            #And again for weighted scores
            if pairHasWeights:
                if wrmsdj > wrmsdi:
                    wscorei = 1.0;
                    wscorej = wrmsdj/wrmsdi;
                else:
                    wscorei = wrmsdi/wrmsdj;
                    wscorej = 1.0;
            
            #Write to pairwise matrices
            nIntersectMatrix[i, j] = nIntersectMatrix[j, i] = len(intersectingRows); #Copy value symetrically in matrix
            pairedScoreMatrix[i, j] = scorei;
            pairedScoreMatrix[j, i] = scorej;
            pairedRmsdMatrix[i, j] = rmsdi;
            pairedRmsdMatrix[j, i] = rmsdj;
            
            if pairHasWeights:
                pairedWScoreMatrix[i, j] = wscorei;
                pairedWScoreMatrix[j, i] = wscorej;
                pairedWRmsdMatrix[i, j] = wrmsdi;
                pairedWRmsdMatrix[j, i] = wrmsdj;
    
    #TMP debugger issue.
    #globals().update(locals())
    
    #Calculate final scores and final RMSD (RMSDs)
    finalScores = pd.DataFrame();
    finalScores["algorithm"] = [algorithm.__class__.__name__ for algorithm in algorithmFunctorList];
    finalScoreArray = np.nanmean(pairedScoreMatrix, axis=1);
    finalScores["final_score"] = finalScoreArray;
    finalScores["algos_compared"] = [np.sum(np.isfinite(pairedScoreMatrix[i,:])) for i in range(0, len(algorithmFunctorList))];
    finalWScoreArray = np.nanmean(pairedWScoreMatrix, axis=1);
    finalScores["final_wscore"] = finalWScoreArray;
    finalScores["w_algos_compared"] = [np.sum(np.isfinite(pairedWScoreMatrix[i,:])) for i in range(0, len(algorithmFunctorList))];
    finalScores["n"] = [basicMetric["n"] for basicMetric in basicMetrics];
    
    #Calculate representative RMSD
    nArray = np.array([metrics["n"] for metrics in basicMetrics], dtype=float);
    nArray[nArray==0] = np.nan;
    
    iRMSDrep = np.nanargmin(finalScoreArray/nArray);
    RMSDrep = basicMetrics[iRMSDrep]["rmsd"];
    finalRMSDs = [RMSDrep*finalScoreArray[i]/finalScoreArray[iRMSDrep] for i in range(len(algorithmFunctorList))];
    finalScores["final_rmsd"] = finalRMSDs;
    
    wiRMSDrep = np.nanargmin(finalWScoreArray/nArray);
    wRMSDrep = basicMetrics[wiRMSDrep]["wrmsd"];
    wfinalRMSDs = [wRMSDrep*finalWScoreArray[i]/finalWScoreArray[wiRMSDrep] for i in range(len(algorithmFunctorList))];
    finalScores["final_wrmsd"] = wfinalRMSDs;
    
    
    
    return basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores;
    
    








































    
    
    
    