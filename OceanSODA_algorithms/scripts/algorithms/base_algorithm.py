#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 10:00:46 2020

@author: tom holding
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 08:57:49 2020

@author: tom holding
"""

import pandas as pd;
import utilities;


#Ternon, J.F., Oudot, C., Dessier, A. and Diverres, D., 2000. A seasonal tropical sink for atmospheric CO2 in the Atlantic Ocean: the role of the Amazon River discharge. Marine Chemistry, 68(3), pp.183-201.
#Amazon plume region
class BaseAlgorithm:
    #These should be overwritten with values specific to the implemented algorithm
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = (); #Tuple containing the algorithm coefficients, index 0 should be the intercept.
        self.coefsUncertainty = (); #Tuple containing uncertainties associated with each coefficient. Must be specified in the same order as self.coeffs
        self.rmsd = None; #Output uncertainty from the algorithm's original fit (root mean squared difference)
        self.r = None; #Optional, r from the model fit (correlation coefficient)
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [];
        self.includedRegionsLats = [];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {}; #Dictionary of string column 'common names' mapped to tuples with (min, max) allowed bounds (inclusive)
                                  #    np.nan can be specified to indicate no bound.
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {}; #Dictionary of string column 'common names' mapped to tuples with (min, max) allowed bounds (inclusive)
                              #    np.nan can be specified to indicate no bound.

    #common names input and output variables (see global_settings for definitions of these)
    @staticmethod
    def input_names():
        raise NotImplementedError("BaseAlgorithm.input_names() is not implemented but should be implemented in all derived classes.");
    @staticmethod
    def output_name():
        raise NotImplementedError("BaseAlgorithm.output_name() is not implemented but should be implemented in all derived classes.");
    
    #Subsets the input data to only use values in the range the algorithm was fitted to
    #data provided as a pandas dataframe
    def _subset_for_restrictions(self, data):
        #Apply algorithm specific restrictions
        inRangeRows = pd.Series([True]*len(data), index=data.index);
        
        #For each range restriction, filter rows outside of the range.
        for key, bounds in self.restrictRanges.items():
            if key not in data.keys():
                print("WARNING in %s: %s does not exist in the matchup data provided. Check input name is correct and that you're using a matchup database which includes this variable." % (type(self).__name__, key));
                continue;
            if (bounds[0] != None): #Some algorithms may not have a lower range
                inRangeRows = inRangeRows & (data[key] >= bounds[0]);
            if (bounds[1] != None): #Some algorithms may not have an upper range
                inRangeRows = inRangeRows & (data[key] <= bounds[1]);
                
        #Subset and return the data
        data = data.loc[inRangeRows];
        return data;
    
    #Flags values outside the experimental conditions to the user
    def _report_flagged_values(self, fullData, usedRowIndices):
        #Note: fullData must be first subset by usedRowIndices to give a fullColumn version of the used data.
        #      Then we can follow simialr logic as in self._subset_for_restrictions
        
        print(self.__str__(), "Flagged values not yet implemented.");
        #For each self.flagRange, go to the rows used of the full-column
        #   dataset and report the number of values which fall outside that
        #   range as well as the total rows that are outside the boundary
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        raise NotImplementedError("BaseAlgorithm._kernal() is not implemented but should be implemented in all derived classes.");
    
    #when predict is set to True, only the input data columns will be required
    def __call__(self, fullColumnData, predict=False):
        #Basic filtering to select only usable rows and columns with required data in
        if predict:
            data = fullColumnData[["lon", "lat", "date"]+self.input_names()]; #Subset to select only the columns used
        else: #using output variable, date etc.
            data = fullColumnData[["date", "lon", "lat"]+self.input_names()+[self.output_name()]]; #Subset to select only the columns used
        data = utilities.subset_complete_rows(data); #Select only rows that are complete (no NaNs)
        
        #Filter out any data points outside the regions supported by the algorithm
        if self.settings["usingTestMatchupDataset"] == False: #When using the test matchup dataset there isn't many datapoints, so don't do algorithm specific filtering
            data = utilities.subset_from_inclusive_coord_list(self.includedRegionsLons, self.includedRegionsLats, data);
        
        #Discard data outside the valid range for the algorithm
        if self.settings["usingTestMatchupDataset"] == False: #When using the test matchup dataset there isn't many datapoints, so don't do algorithm specific filtering
            data = self._subset_for_restrictions(data); #Sometimes returns 'None' for some reason...?
        
        #Check that there is some data left
        if (data is None) or (len(data) == 0): #if no complete rows so continue to the next algorithm
            raise ValueError("No data in subsetted data for", self);

        #Report the number of rows with values inside the flag range
        self._report_flagged_values(fullColumnData, data.index);
        
        #Compute output estimate
        modelOutput = self._kernal(data);
        
        #remove any rows which couldn't be estimated for any reason
        modelOutput = modelOutput.dropna();
        data = data.loc[modelOutput.index];
        
        return modelOutput, data; #Return estimates and the data used to calculate them (as this data has been further subsetted)
    
    def __str__(self):
        return "BaseAlgorithm: Serves as the parent of all other algorithm implementations";
