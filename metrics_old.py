#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:06:42 2020

@author: tom holding
"""

import numpy as np;
import pandas as pd;


#Calculates the expected type B uncertainty on an in situ data variable (i.e. DIC or AT)
#Based on the uncertainties suggested in Bockmon 2015: https://www.sciencedirect.com/science/article/pii/S0304420315000213
#   using "nominal state of the art" uncertainties for in situ measurements
#   see: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
#Returns an iterable containing an uncertainty for each row in insituData
#insituData: dataframe containing the in situ data
#varoable: target variable from which the (e.g. 'DIC' or 'AT')
#settings: the global setting dictionary
def calc_insitu_uncertainty(insituData, variable, settings):
    if settings["useErrorRatios"]:
        insituUncertainty = (insituData[variable]*settings["insituErrorRatio"][variable]);
    else:
        insituUncertainty = [settings["insituError"][variable]]*len(insituData);
    return insituUncertainty

#Calculates/gets the matchup database reference uncertainty for the output variable
#Uses the type A uncertainty calculated from AT / DIC (std dev) if available
#   otherwise the 'state-of-the-art' uncertainty from Bockmon et al 2015: https://www.sciencedirect.com/science/article/pii/S0304420315000213
#Arguments:
#   dataUsed: dataframe containing all the rows used to assess the algorithm
#   outputVar: name of the output variable (i.e. 'AT' or 'DIC')
#   settings: the global settings dictionary
def calc_reference_uncertainty(dataUsed, outputVar, settings):
    
    uncertainties = dataUsed[outputVar+"_err"]; #type A uncertainties
    
    missingTypeA = np.isfinite(uncertainties) == False; #Which rows are missing type A uncertainty
    
    #Fill missing type A values with the nominal / state-of-the-art type B error from Bockmon 2015 ( https://www.sciencedirect.com/science/article/pii/S0304420315000213 )
    try:
        if settings["useErrorRatios"]:
            outputVarValues = dataUsed[outputVar];
            uncertainties.loc[missingTypeA] = outputVarValues[missingTypeA]*(settings["insituErrorRatio"][outputVar])
        else:
            uncertainties.loc[missingTypeA] = settings["insituError"][outputVar];
    except KeyError: #E.g. for pH and pCO2w the settings["insituError"] values aren't currently provided. In this case nothing to be done but discard the data
        uncertainties.loc[missingTypeA] = np.nan;
    
    return uncertainties;


###Calculate weights based on normalised inverse size of combined insitu and algorithm uncertainty
# algorithmUncertainty - rmsd from the algorithms original paper
# matchupOutputUncertainty - uncertainty estimate in the matchup dataset (type A or type B as a fallback)
def calc_weights(algorithmUncertainty, matchupOutputVarUncertainty):
    combinedUncertainty = np.sqrt( algorithmUncertainty**2 + matchupOutputVarUncertainty**2 );
    weights = 1.0 / combinedUncertainty;
    weights /= np.nansum(weights); #normalise
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
def calc_basic_metrics(algorithmOutput, dataUsed, settings):
    ####Calculate basic statistics
    outputVariable = algorithmOutput["outputVar"];
    
    #Uncertainty combined for in situ output variable and model output variable #Note: Called "squaredErrors" in Peter's 'printFit' function.
    #This is used to calculate weights based on normalised inverse size of combined insitu and algorithm uncertainty
    if np.any(algorithmOutput["rmsd"] is not None): #cannot calculate weights for algorithms with no combined uncertainty
        #insituOutputUncertainties = calc_insitu_uncertainty(dataUsed, outputVariable, settings);
        #combinedUncertainty = np.sqrt(algorithm.rmsd**2 + insituOutputUncertainties**2); #Old Pathfinders method
        #calculate combined uncertainty
        if type(algorithmOutput["rmsd"]) is float:
            #If a single algorithm uncertainty is used, turn it into an array (one copied value for each model output)
            algorithmOutput["rmsd"] = np.array([algorithmOutput["rmsd"]]*len(algorithmOutput["modelOutput"]));
        
        matchupOutputVarUncertainty = calc_reference_uncertainty(dataUsed, outputVariable, settings);
        
        #Calculate weights
        weights = calc_weights(algorithmOutput["rmsd"], matchupOutputVarUncertainty);
        hasWeights = True;
    else:
        weights = np.nan;
        matchupOutputVarUncertainty = np.nan;
        hasWeights = False;
        print("*** No RMSD for", algorithmOutput["name"], "Weighted metrics will not be calculated.");
    
    #Calculate prediction errors: the difference between reference measurement and model output
    predictionErrors = algorithmOutput["modelOutput"] - dataUsed[outputVariable];
    bias = np.mean(predictionErrors);
    predictionErrorsSquared = predictionErrors**2; #Named 'squaredError' in Peter's 'printFit' function
    absPredictionError = abs(predictionErrors);
    
    #Metrics of ability of the algorithm to reproduce measured values
    predictionRmsd = np.sqrt(predictionErrorsSquared.mean()); #rmsd
    predictionMad = absPredictionError.mean(); #mean absolute difference
    if hasWeights:
        wpredictionRmsd = np.sqrt(np.nansum(weights*predictionErrorsSquared)); #weighted sum of the squared prediction error. Note sum not mean because the weights are normalised to sum to 1.
        wpredictionMad = (np.nansum(weights*absPredictionError)); #weighted mean absolute difference. Note sum is used, not mean, as weights are already normalised to sum to 1
        wpredictionBias = (np.nansum(weights*predictionErrors))/(np.nansum(weights)); #weighted mean absolute difference. Note sum is used, not mean, as weights are already normalised to sum to 1
    else:
        wpredictionRmsd = np.nan;
        wpredictionMad = np.nan;
        wpredictionBias = np.nan;
    

    
    #Calculates means of model and reference values (unweighted and weighted)
    meanModelOutput = algorithmOutput["modelOutput"].mean();
    meanReferenceOutput = dataUsed[outputVariable].mean();
    if hasWeights:
        wMeanModelOutput = np.nansum(weights*algorithmOutput["modelOutput"]); #weighted mean of the model/predicted output. Note sum is used, not mean, as weights are already normalised to sum to 1
        wMeanInSituOutput = np.nansum(weights*dataUsed[outputVariable]); #weighted mean of the in situ measured output. Note sum is used, not mean, as weights are already normalised to sum to 1
    else:
        wMeanModelOutput = np.nan;
        wMeanInSituOutput = np.nan;
    
    #Calculating standard deviation
    #Using Peters method: TODO: repeat / check with other method
    if len(dataUsed) > 1: #Can't calculate SD from 1 data point
        covariance = np.ma.cov(dataUsed[outputVariable], y = algorithmOutput["modelOutput"]) # covariance matrix
        referenceOutputSD = covariance[0, 0] ** .5
        modelOutputSD = covariance[1, 1] ** .5
        r = covariance[0, 1] / (referenceOutputSD * modelOutputSD)
        
        #referenceOutputSD = np.std(dataUsed[outputVariable], ddof=1); #ddof=1 to calculate sample SD
        #modelOutputSD = np.std(modelOutput, ddof=1); #ddof=1 to calculate sample SD
        #r = np.sum (np.sum(dataUsed[outputVariable]-np.mean(dataUsed[outputVariable])) * np.sum(modelOutput-np.mean(modelOutput)) ) / referenceOutputSD*modelOutputSD;
        
        #weighted standard deviations and correlation coefficient
        if hasWeights:
            sumX2 = np.nansum(weights * (dataUsed[outputVariable] - wMeanInSituOutput) ** 2)
            sumY2 = np.nansum(weights * (algorithmOutput["modelOutput"] - wMeanModelOutput) ** 2)
            sumXY = np.nansum(weights * (dataUsed[outputVariable] - wMeanInSituOutput) * (algorithmOutput["modelOutput"] - wMeanModelOutput))
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
    basicMetrics = {"model_output": algorithmOutput["modelOutput"],
                    "model_uncertainty": algorithmOutput["rmsd"],
                    "model_propagated_input_uncertainty": algorithmOutput["propagatedInputUncertainty"],
                    "model_combined_output_uncertainty": algorithmOutput["combinedUncertainty"],
                    "model_output_mean": meanModelOutput,
                    "reference_output_uncertainty": matchupOutputVarUncertainty,
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
                    "wbias": wpredictionBias, #bias = mean difference: mean(model-reference)
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
def calc_all_metrics(algorithmOutputList, matchupData, settings, outputVar):
    #Create some space for storing metrics (pairwise metrics are stored as matrices.)
    basicMetrics = []; #Stores simple non-paired metrics like mean, standard deviation, RMSD and MAD
#    nIntersectMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=int); #Symetrical
#    pairedScoreMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise score (see Land2019 Remote Sensing of Environment, section 2.3.2). Index order is [ith, jth]
#    pairedWScoreMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise weighted score. Index order is [ith, jth]
#    pairedRmsdMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise RMSD. Index order is [ith, jth]
#    pairedWRmsdMatrix = np.ma.masked_all((len(algorithmFunctorList), len(algorithmFunctorList)), dtype=float); #Pairwise weighted RMSD. Index order is [ith, jth]
    nIntersectMatrix = np.ma.zeros((len(algorithmOutputList), len(algorithmOutputList)), dtype=int); nIntersectMatrix.mask=False; #Symetrical
    pairedScoreMatrix = np.ma.empty((len(algorithmOutputList), len(algorithmOutputList)), dtype=float); pairedScoreMatrix[:]=np.nan; pairedScoreMatrix.mask=False; #Pairwise score (see Land2019 Remote Sensing of Environment, section 2.3.2). Index order is [ith, jth]
    pairedWScoreMatrix = np.ma.empty((len(algorithmOutputList), len(algorithmOutputList)), dtype=float); pairedWScoreMatrix[:]=np.nan; pairedWScoreMatrix.mask=False; #Pairwise weighted score. Index order is [ith, jth]
    pairedRmsdMatrix = np.ma.empty((len(algorithmOutputList), len(algorithmOutputList)), dtype=float); pairedRmsdMatrix[:]=np.nan; pairedRmsdMatrix.mask=False; #Pairwise RMSD. Index order is [ith, jth]
    pairedWRmsdMatrix = np.ma.empty((len(algorithmOutputList), len(algorithmOutputList)), dtype=float); pairedWRmsdMatrix[:]=np.nan; pairedWRmsdMatrix.mask=False; #Pairwise weighted RMSD. Index order is [ith, jth]
    
    #### Calculate pairwise metrics as described in section 2.3.2 of Land, P.E., Findlay, H.S., Shutler, J.D., Ashton, I.G., Holding, T., Grouazel, A., Girard-Ardhuin, F., Reul, N., Piolle, J.F., Chapron, B. and Quilfen, Y., 2019. Optimum satellite remote sensing of the marine carbonate system using empirical algorithms in the global ocean, the Greater Caribbean, the Amazon Plume and the Bay of Bengal. Remote Sensing of Environment, 235, p.111469.
    for i in range(len(algorithmOutputList)):
        algorithmOutput = algorithmOutputList[i];
        outputVariable = algorithmOutput["outputVar"];
        matchupDataUsed = matchupData.loc[algorithmOutput["dataUsedIndices"]];

        #Calculate means, standard deviations, rmse, mad, r
        basicMetrics.append(calc_basic_metrics(algorithmOutput, matchupDataUsed, settings));
        
        #Calculate pairwise metrics
        for j in range(i): #for each pair of algorithms ((self, self) pair excluded)
            #Sanity check - can only compare algorithms which are predicting the same variable
            if algorithmOutput["outputVar"] != algorithmOutputList[j]["outputVar"]:
                raise ValueError("Cannot compare algorithms which estimate different output variables.");
            
            #Only use data points from the matchup dataset that have been predicted by both algorithms
            intersectingRows = set(algorithmOutputList[i]["dataUsedIndices"]).intersection(set(algorithmOutputList[j]["dataUsedIndices"]));
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
            diffFromReferencei = algorithmOutputList[i]["modelOutput"].loc[intersectingRows] - matchupData[outputVariable].loc[intersectingRows]; #Differences between model and reference output for ith algorithm at only intersecting locations
            diffFromReferencej = algorithmOutputList[j]["modelOutput"].loc[intersectingRows] - matchupData[outputVariable].loc[intersectingRows]; #Differences between model and reference output for jth algorithm at only intersecting locations;
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
    finalScores["algorithm"] = [algorithmOutput["name"] for algorithmOutput in algorithmOutputList];
    finalScoreArray = np.nanmean(pairedScoreMatrix, axis=1);
    finalScores["final_score"] = finalScoreArray;
    finalScores["algos_compared"] = [np.sum(np.isfinite(pairedScoreMatrix[i,:])) for i in range(0, len(algorithmOutputList))];
    finalWScoreArray = np.nanmean(pairedWScoreMatrix, axis=1);
    finalScores["final_wscore"] = finalWScoreArray;
    #finalScores["algos_compared"] = [np.sum(np.isfinite(pairedScoreMatrix[i,:])) for i in range(0, len(algorithmFunctorList))];
    #finalScores["w_algos_compared"] = [np.sum(np.isfinite(pairedWScoreMatrix[i,:])) for i in range(0, len(algorithmFunctorList))];
    finalScores["n"] = [basicMetric["n"] for basicMetric in basicMetrics];
    finalScores["wbias"] = [basicMetric["wbias"] for basicMetric in basicMetrics];

    #Calculate representative RMSD
    nArray = np.array([metrics["n"] for metrics in basicMetrics], dtype=float);
    nArray[nArray==0] = np.nan;
    
    iRMSDrep = np.nanargmin(finalScoreArray/nArray);
    RMSDrep = basicMetrics[iRMSDrep]["rmsd"];
    finalRMSDs = [RMSDrep*finalScoreArray[i]/finalScoreArray[iRMSDrep] for i in range(len(algorithmOutputList))];
    finalScores["final_rmsd"] = finalRMSDs;
    
    wiRMSDrep = np.nanargmin(finalWScoreArray/nArray);
    wRMSDrep = basicMetrics[wiRMSDrep]["wrmsd"];
    wfinalRMSDs = [wRMSDrep*finalWScoreArray[i]/finalWScoreArray[wiRMSDrep] for i in range(len(algorithmOutputList))];
    finalScores["final_wrmsd"] = wfinalRMSDs;
    
    uncendtoend=[]; # this is the combined uncertainty between RMSD and measurement uncertainty.
    for i in range(len(algorithmOutputList)):
         x= np.sqrt( finalScores["wbias"][i]**2 + finalScores["final_wrmsd"][i]**2 )
         uncendtoend.append(x)
    finalScores["unc_end_end"]=uncendtoend
    
    return basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores;
    
    








































    
    
    
    