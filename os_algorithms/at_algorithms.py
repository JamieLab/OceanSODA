#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 08:32:17 2020

@author: tom holding
"""

import numpy as np;
import pandas as pd;

from .base_algorithm import BaseAlgorithm;
from .utilities import subset_from_mask, subset_from_inclusive_coord_list;


#Astor, Y.M., Lorenzoni, L., Guzman, L., Fuentes, G., Muller-Karger, F., Varela, R., Scranton, M., Taylor, G.T. and Thunell, R., 2017. Distribution and variability of the dissolved inorganic carbon system in the Cariaco Basin, Venezuela. Marine Chemistry, 195, pp.15-26.
#2004 data, Carribean sea, Eastern sub-basin of the Cariaco basin
class Astor2017c_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Astor2017c_at: A17c(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-774.9, 86.3]; #intersept, salinity slope, fig 7b
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None; #equation 1
        self.r = None; #But see section 3.3 for quoted values?
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #See fig 2, 3 for extents used by authors
        self.includedRegionsLons = [(-65.1, -64.3), #Eastern sub-basin of the Cariaco basin
                                    ];
        self.includedRegionsLats = [(10.0, 11.4),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35.5, 37), #See section 3.3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5d
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Astor, Y.M., Lorenzoni, L., Guzman, L., Fuentes, G., Muller-Karger, F., Varela, R., Scranton, M., Taylor, G.T. and Thunell, R., 2017. Distribution and variability of the dissolved inorganic carbon system in the Cariaco Basin, Venezuela. Marine Chemistry, 195, pp.15-26.
#2008 data, Carribean sea, Easten sub-basin of the Cariaco basin
class Astor2017a_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Astor2017a_at: A17a(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        self.coefs = [831.5, 42.8]; #intersept, SSS, see fig 7b
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None; #equation 1
        self.r = None; #But see section 3.3 for quoted values?
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #See fig 2, 3 for extents used by authors
        self.includedRegionsLons = [(-65.1, -64.3), #Eastern sub-basin of the Cariaco basin
                                    ];
        self.includedRegionsLats = [(10.0, 11.4),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35.5, 37), #See section 3.3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
        
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 7b
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Astor, Y.M., Lorenzoni, L., Guzman, L., Fuentes, G., Muller-Karger, F., Varela, R., Scranton, M., Taylor, G.T. and Thunell, R., 2017. Distribution and variability of the dissolved inorganic carbon system in the Cariaco Basin, Venezuela. Marine Chemistry, 195, pp.15-26.
#2008 data, Carribean sea, Western sub-basin of the Cariaco basin
class Astor2017d_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Astor2017d_at: A17d(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        self.coefs = [581.8, 49.7]; #intersept, SSS, see fig 7a
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None; #equation 1
        self.r = None; #But see section 3.3 for quoted values?
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #See fig 2, 3 for extents used by authors
        self.includedRegionsLons = [(-66.1, -65.1), #Eastern sub-basin of the Cariaco basin
                                    ];
        self.includedRegionsLats = [(10.0, 11.4),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35.5, 37), #See section 3.3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 7a
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Astor, Y.M., Lorenzoni, L., Guzman, L., Fuentes, G., Muller-Karger, F., Varela, R., Scranton, M., Taylor, G.T. and Thunell, R., 2017. Distribution and variability of the dissolved inorganic carbon system in the Cariaco Basin, Venezuela. Marine Chemistry, 195, pp.15-26.
#2009 data, Carribean sea, Easten sub-basin of the Cariaco basin
class Astor2017b_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Astor2017b_at: A17b(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        self.coefs = [-595, 81.7]; #intersept, SSS, see fig 7b
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None; #equation 1
        self.r = None; #But see section 3.3 for quoted values?
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #See fig 2, 3 for extents used by authors
        self.includedRegionsLons = [(-65.1, -64.3), #Eastern sub-basin of the Cariaco basin
                                    ];
        self.includedRegionsLats = [(10.0, 11.4),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35.5, 37), #See section 3.3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 7b
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;
  
#Astor, Y.M., Lorenzoni, L., Guzman, L., Fuentes, G., Muller-Karger, F., Varela, R., Scranton, M., Taylor, G.T. and Thunell, R., 2017. Distribution and variability of the dissolved inorganic carbon system in the Cariaco Basin, Venezuela. Marine Chemistry, 195, pp.15-26.
#2008 data, Carribean sea, Western sub-basin of the Cariaco basin
class Astor2017e_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Astor2017e_at: A17e(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        self.coefs = [1075.8, 94.8]; #intersept, SSS, see fig 7a
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None; #equation 1
        self.r = None; #But see section 3.3 for quoted values?
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #See fig 2, 3 for extents used by authors
        self.includedRegionsLons = [(-66.1, -65.1), #Eastern sub-basin of the Cariaco basin
                                    ];
        self.includedRegionsLats = [(10.0, 11.4),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35.5, 37), #See section 3.3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 7a
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Brewer, P.G., Glover, D.M., Goyet, C. and Shafer, D.K., 1995. The pH of the North Atlantic Ocean: Improvements to the global model for sound absorption in seawater. Journal of Geophysical Research: Oceans, 100(C5), pp.8761-8776.
#Implementation of equation 8 (<250m depth)
class Brewer1995_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Brewer1995_at: B95(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "DO", "SiO4", "PO4", "NO3"];

    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [530.0, 51.132, -0.033, 0.973, 59.3, -4.31]; #intersept, salinity, DO, SiO4, PHO4, NO3 see equation 8
        self.coefsUncertainty = [None, None, None, None, None, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = None;
        self.r = 0.975**0.5; #equation 8
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-100, 30)]; #See fig 1a
        self.includedRegionsLats = [(0, 80)]; #See fig 1a
        #TODO: Mask away pacific corner
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (31, 37.5), #See fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (-2+273.15, 28+273.15), #See fig 5
                           "PO4": (0, 2.2), #See fig 5 (continued, second page)
                           "SiO4": (0, 65), #See fig 5 (continues, 3rd page)
                           "NO3": (0, 34), #See fig 5 (continues, 3rd page)
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 11
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cai, W.J., Hu, X., Huang, W.J., Jiang, L.Q., Wang, Y., Peng, T.H. and Zhang, X., 2010. Alkalinity distribution in the western North Atlantic Ocean margins. Journal of Geophysical Research: Oceans, 115(C8).
#"Tropical North Atlantic Ocean margine surrounding the Amazon River plume" - fig 11 (shown spatially in fig 1)
class Cai2010a_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cai2010a_at: Ca10a(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [320.8, 56.8]; #intersept, salinity slope, see fig 11
        self.coefsUncertainty = [7.1, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = 7.9; #See fig 11
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-60, -40)]; #See "Western Tropical Atlantic" fig 1
        self.includedRegionsLats = [(0, 20)]; #See "Western Tropical Atlantic" fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (24, 37), #See fig 11
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 11
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cai, W.J., Hu, X., Huang, W.J., Jiang, L.Q., Wang, Y., Peng, T.H. and Zhang, X., 2010. Alkalinity distribution in the western North Atlantic Ocean margins. Journal of Geophysical Research: Oceans, 115(C8).
#Northern Gulf of Mexico continental shelf" - fig 9 (shown in fig 1)
class Cai2010b_at(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Cai2010b_at: Ca10b(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [1996.7, 11.7]; #intersept, salinity slope, see fig 9
        self.coefsUncertainty = [31.5, None]; #Uncertainty reported for the coefficients, see fig 9
        self.rmsd = 8.6; #See fig 9
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-100, -82)]; #See "nGMx" fig 1
        self.includedRegionsLats = [(23, 32)]; #See "nGMx" fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (5, 37.5), #See fig 9
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 9
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cai, W.J., Hu, X., Huang, W.J., Jiang, L.Q., Wang, Y., Peng, T.H. and Zhang, X., 2010. Alkalinity distribution in the western North Atlantic Ocean margins. Journal of Geophysical Research: Oceans, 115(C8).
#"Labrador Sea" - fig 3 (shown spatially in fig 1)
class Cai2010c_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cai2010c_at: Ca10c(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [1124.4, 33.0]; #intersept, salinity slope, see fig 3
        self.coefsUncertainty = [103.1, None]; #Uncertainty reported for the coefficients, see fig 3
        self.rmsd = 12.7; #See fig 11
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-60, -45)]; #See "Labrador Sea" fig 1
        self.includedRegionsLats = [(50, 62)]; #See "Labrador Sea" fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (33, 35), #See fig 3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cai, W.J., Hu, X., Huang, W.J., Jiang, L.Q., Wang, Y., Peng, T.H. and Zhang, X., 2010. Alkalinity distribution in the western North Atlantic Ocean margins. Journal of Geophysical Research: Oceans, 115(C8).
#"Gulf of Main: inner shelf (low salinity)" - fig 4 (shown spatially in fig 1)
class Cai2010d_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cai2010d_at: Ca10d(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [75.1, 65.8]; #intersept, salinity slope, see fig 4
        self.coefsUncertainty = [291.2, None]; #Uncertainty reported for the coefficients, see fig 4
        self.rmsd = 9.7; #See fig 11
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-72, -68)]; #See inset of fig 4
        self.includedRegionsLats = [(40, 44)]; #See inset of fig 4
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (30.5, 31.75), #See fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Cai, W.J., Hu, X., Huang, W.J., Jiang, L.Q., Wang, Y., Peng, T.H. and Zhang, X., 2010. Alkalinity distribution in the western North Atlantic Ocean margins. Journal of Geophysical Research: Oceans, 115(C8).
#"Gulf of Main: offshore waters (higher salinity)" - fig 4 (shown spatially in fig 1)
class Cai2010e_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cai2010e_at: Ca10e(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [932.7, 39.1]; #intersept, salinity slope, see fig 4
        self.coefsUncertainty = [16.5, None]; #Uncertainty reported for the coefficients, see fig 4
        self.rmsd = 1.7; #See fig 11
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-72, -68)]; #See inset of fig 4
        self.includedRegionsLats = [(40, 44)]; #See inset of fig 4
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (30.5, 31.75), #See fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cooley, S.R. and Yager, P.L., 2006. Physical and biological contributions to the western tropical North Atlantic Ocean carbon sink formed by the Amazon River plume. Journal of Geophysical Research: Oceans, 111(C8).
class Cooley2006a_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cooley2006a_at: C06(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [369.7, 55.0]; #intersept, salinity slope
        self.coefsUncertainty = [68.1, 2.1]; #Uncertainty reported for the coefficients
        self.rmsd = 9.6; #See section 3.6
        self.r = 983; #See Table 2 "plume"
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-40, -59)]; #See Cooley2006 section 2.1
        self.includedRegionsLats = [(3, 14)]; #See Cooley2006 section 2.1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (25, 35), #Upper bound: see Cooley2006 section 2.3.2 (definition of Amazon plume). Lower bound see Fig 4.
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (20+273.15, 30+273.15), #+273.15 to convert from C to K. See figure 4
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5c
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Copin‐Montégut, C., 1993. Alkalinity and carbon budgets in the Mediterranean Sea. Global Biogeochemical Cycles, 7(4), pp.915-925.
class CopinMontegut1993_at(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "CopinMontegut1993_at: CM93(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-1072.6, 94.85]; #intersept, salinity slope, see first equation in Results
        self.coefsUncertainty = [16.0, 0.4]; #Uncertainty reported for the coefficients, see first equation in Results
        self.rmsd = None;
        self.r = 0.998; #See first paragraph of Results
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-3, 0)]; #See fig 1
        self.includedRegionsLats = [(35, 38)]; #See fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (36.4, 38.4), #See first paragraph of Results and fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {#none reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See first equation in Results
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Copin-Montégut, C. and Bégovic, M., 2002. Distributions of carbonate properties and oxygen along the water column (0–2000 m) in the central part of the NW Mediterranean Sea (Dyfamed site): influence of winter vertical mixing on air–sea CO2 and O2 exchanges. Deep Sea Research Part II: Topical Studies in Oceanography, 49(11), pp.2049-2066.
class CopinMontegut2002_at(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "CopinMontegut2002_at: CM02(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-1304.6, 100.97]; #intersept, salinity slope, see eq 1
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, see eq 1
        self.rmsd = 2.9; #See eq 1
        self.r = 0.977; #See eq 1
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(7, 9)]; #See fig 1, DYFAMED site and surrounding region
        self.includedRegionsLats = [(42.5, 44)]; #See fig 1, DYFAMED site and surrounding region
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (38, 38.6), #See figs 3b, 4, 5b and discussion around eq 1-3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (12+273.15, 25.5+273.15), #See fig 3a
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Corbière, A., Metzl, N., Reverdin, G., Brunet, C. and Takahashi, T., 2007. Interannual and decadal variability of the oceanic carbon sink in the North Atlantic subpolar gyre. Tellus B: Chemical and Physical Meteorology, 59(2), pp.168-178.
class Corbiere2007_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Corbiere2007_at: Co07(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [713.5, 45.808]; #intersept, salinity slope, see eq 1
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, see eq 1
        self.rmsd = 10.3; #See eq 1
        self.r = 0.92**0.5;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-60, -20)]; #See fig 1
        self.includedRegionsLats = [(44, 65)]; #See fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (31, 35.5), #See fig 2
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (1.0+273.15, 16.0+273.15), #See fig 3a
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Gemayel, E., Hassoun, A.E.R., Benallal, M.A., Goyet, C., Rivaro, P., Abboud-Abi Saab, M., Krasakopoulou, E., Touratier, F. and Ziveri, P., 2015. Climatological variations of total alkalinity and total dissolved inorganic carbon in the Mediterranean Sea surface waters. Earth System Dynamics, 6(2), pp.789-800.
class Gemayel2015_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Gemayel2015_at: Ge15(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [2558.4, 49.83, -3.12, -3.89, -1.06]; #intersept, SSS, SSS^2, SST, SST^2, eq 1
        self.coefsUncertainty = [None, None, None, None, None]; #Uncertainty reported for the coefficients
        self.rmsd = 10.34; #See table 2 NOTE: This is validation RMSD, training RMSD is 10.6
        self.r = 0.96; #See eq 2
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #Mediterranean excluding Black Sea, Adriatic Sea and Aegean Sea (fig 1, fig 2)
        self.includedRegionsLons = [(-5, 12),
                                    (12, 22),
                                    (22, 36),
                                    ];
        self.includedRegionsLats = [(30, 44), #Mediterranean up to Italy
                                    (30, 40), #below adriatic sea
                                    (30, 27), #below Aegean sea
                                    ];
        
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SST": (13.0+273.15, None), #See text below eq 1
                               "SSS": (36.3, 39.65), #See text below eq 1
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        SSS = dataToUse["SSS"]-38.2;
        SST = dataToUse["SST"]-18-273.15;
        
        #From equation 1
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*SSS + \
                      self.coefs[2]*(SSS**2) + \
                      self.coefs[3]*SST + \
                      self.coefs[4]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = self.coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = self.coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = self.coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = self.coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );

        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Goyet, C., Adams, R. and Eischeid, G., 1998. Observations of the CO2 system properties in the tropical Atlantic Ocean. Marine Chemistry, 60(1-2), pp.49-61.
class Goyet1998_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Goyet1998_at: G98(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-22.7, 66.141]; #intersept, salinity slope, see near end of section 3
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = 4.9; #See near end of section 3
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-40, -15)]; #See fig 1
        self.includedRegionsLats = [(-32, 8)]; #See near end of section 3 and fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (35, 37.5), #See fig 5
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (20+273.15, 30+273.15), #fig 5
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 11
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Hassoun, A.E.R., Gemayel, E., Krasakopoulou, E., Goyet, C., Saab, M.A.A., Ziveri, P., Touratier, F., Guglielmi, V. and Falco, C., 2015. Modeling of the total alkalinity and the total inorganic carbon in the Mediterranean Sea.
#Whole mediterranean as one region (surface 0-25m) - see table 1 eq2
class Hassoun2015_full_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Hassoun2015_full_at: Ha15_full(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-845.0, 89.0]; #intersept, salinity slope, see table 1 eq2
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = 18.0; #See table 1 eq2
        self.r = 0.96; #See table 1 eq2
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #Mediterranean excluding Black Sea and Aegean Sea (fig 1)
        self.includedRegionsLons = [(-7,  22),
                                    (22, 36),
                                    ];
        self.includedRegionsLats = [(30, 47),
                                    (30, 37),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See fig 11
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Hassoun, A.E.R., Gemayel, E., Krasakopoulou, E., Goyet, C., Saab, M.A.A., Ziveri, P., Touratier, F., Guglielmi, V. and Falco, C., 2015. Modeling of the total alkalinity and the total inorganic carbon in the Mediterranean Sea.
#Mediterranean split into East and West basins (surface 0-25m) - see table 1 eq 6 and 9
class Hassoun2015_basins_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Hassoun2015_basins_at: Ha15b(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        #super().__init__(settings); #Call the parent class's initator
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
    
    
    #Note, separate models for different regions within the single algorithm is ok here, because together they cover the OceanSODA Mediterranean region completely, with no overlap to other regions
    #east Mediterranean basin, iho definition
    def _east_basin(self, data):
        coefs = [-846.0, 89.0]; #intersept, salinity, see table 1 eq 9
        rmsd = 19.0; #See table 1 eq 9
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "east_basin_mask");
        
        ### No condition ranges given
#        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
#                              (dataToUse["SSS"] > 31) &
#                              (dataToUse["SSS"] < 38)
#                              ];
        
        #Equation from table 1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
                      
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #east Mediterranean basin, iho definition
    def _west_basin(self, data):
        coefs = [-956.0, 92.0]; #intersept, salinity, see table 1 eq 6
        rmsd = 12.5; #See table 1 eq 6
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "west_basin_mask");
        
        ### No condition ranges given
#        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
#                              (dataToUse["SSS"] > 31) &
#                              (dataToUse["SSS"] < 38)
#                              ];
        
        #Equation from table 1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each basin
        zoneData, zoneUncertainty, zoneRmsd = self._east_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Hassoun2015_basins_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._west_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Hassoun2015_basins_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;


#Huang, W.J., Cai, W.J., Powell, R.T., Lohrenz, S.E., Wang, Y., Jiang, L.Q. and Hopkinson, C.S., 2012. The stoichiometry of inorganic carbon and nutrient removal in the Mississippi River plume and adjacent continental shelf. Biogeosciences, 9(7), pp.2781-2792.
class Huang2012_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Huang2012_at: Hu12(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [2002.7, 12.113]; #intersept, salinity slope, fig 2c
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = None;
        self.r = 0.956**0.5; #Fig 2c
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-91.5, -89)]; #fig 2d
        self.includedRegionsLats = [(28, 30)]; #fig 2d
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (0.0, 35.0), #See salinity range in fig 2c
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {};

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5d
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Lee, K., Tong, L.T., Millero, F.J., Sabine, C.L., Dickson, A.G., Goyet, C., Park, G.H., Wanninkhof, R., Feely, R.A. and Key, R.M., 2006. Global relationships of total alkalinity with salinity and temperature in surface waters of the world's oceans. Geophysical research letters, 33(19).
#All zones/equations
class Lee2006_at(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Lee2006_at: L06(at)";
    
    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"];
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        #super().__init__(settings); #Call the parent class's initator
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
    
    
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    #Subtropics and tropics (excluding zone 2)
    def _zone1(self, data):
        coefs = [2305, 58.66, 2.32, -1.41, 0.040]; #intersept, salinity, salinity^2, SST, SST^2 see table 1
        rmsd = 8.6; #See table 1 (zone 1)
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SSS"] > 31) &
                              (dataToUse["SSS"] < 38)
                              ];
        
        SSS = dataToUse["SSS"]-35;
        SST = dataToUse["SST"]-20-273.15;
        
        #Equation from table 1 (zone 1)
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SSS + \
                      coefs[2]*(SSS**2) + \
                      coefs[3]*SST + \
                      coefs[4]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #Equatorial upwelling Pacific
    def _zone2(self, data):
        coefs = [2294, 64.88, 0.39, -4.52, -0.232]; #intersept, salinity, salinity^2, SST, SST^2 see table 1
        rmsd = 7.3; #See table 1 (zone 2)
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 18+273.15) &
                              (dataToUse["SSS"] > 31) &
                              (dataToUse["SSS"] < 36.5)
                              ];
        
        #TODO: also allow data points from zone 1 if their SST is (18, 20])
        
        SSS = dataToUse["SSS"]-35;
        SST = dataToUse["SST"]-29-273.15;
        
        #Equation from table 1 (zone 2)
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SSS + \
                      coefs[2]*(SSS**2) + \
                      coefs[3]*SST + \
                      coefs[4]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #North Atlantic
    def _zone3(self, data):
        coefs = [2305, 53.97, 2.74, -1.16, -0.040]; #intersept, salinity, salinity^2, SST, SST^2 see table 1
        rmsd = 6.4; #See table 1 (zone 3)
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 0+273.15) &
                              (dataToUse["SST"] < 20+273.15) &
                              (dataToUse["SSS"] > 31) &
                              (dataToUse["SSS"] < 37)
                              ];
        
        SSS = dataToUse["SSS"]-35;
        SST = dataToUse["SST"]-20-273.15;
        
        #Equation from table 1 (zone 3)
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SSS + \
                      coefs[2]*(SSS**2) + \
                      coefs[3]*SST + \
                      coefs[4]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #North Pacific
    def _zone4(self, data):
        coefs = [2305, 53.23, 1.85, -14.72, -0.158, 0.062]; #intersept, salinity, salinity^2, SST, SST^2, SST*longitude see table 1
        rmsd = 8.7; #See table 1 (zone 4)
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "zone4_mask");
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (dataToUse["SSS"] > 31) &
                              (dataToUse["SSS"] < 35)
                              ];
        
        SSS = dataToUse["SSS"]-35;
        SST = dataToUse["SST"]-20-273.15;
        lon = dataToUse["lon"]+180; #Should longitude be (-180, 180) or (0 360)?
        
        #Equation from table 1 (zone 4)
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SSS + \
                      coefs[2]*(SSS**2) + \
                      coefs[3]*SST + \
                      coefs[4]*(SST**2) + \
                      coefs[5]*SST*lon; #Should longitude be (-180, 180) or (0 360)?
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm5 = coefs[5] * dataToUse["SST_err"]*lon; #B*SST*longitude. Assume 0 uncertainty in longitude.
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #Southern Ocean
    def _zone5(self, data):
        coefs = [2305, 52.48, 2.85, -0.49, 0.086]; #intersept, salinity, salinity^2, SST, SST^2 see table 1
        rmsd = 8.4; #See table 1 (zone 5)
        
        #Subset data to only rows valid for this zone. See Table 1
        dataToUse = subset_from_mask(data, self.regionMasks, "zone5_mask");
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (dataToUse["SSS"] > 33) &
                              (dataToUse["SSS"] < 36)
                              ];
        
        SSS = dataToUse["SSS"]-35;
        SST = dataToUse["SST"]-20-273.15;
        
        #Equation from table 1 (zone 5)
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SSS + \
                      coefs[2]*(SSS**2) + \
                      coefs[3]*SST + \
                      coefs[4]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*SSS; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["SST_err"]; #B*SST
        uterm4 = coefs[4] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    def _kernal(self, dataToUse):
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each zone
        zoneData, zoneUncertainty, zoneRmsd = self._zone1(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone2(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone3(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone4(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone5(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;


#Lefèvre, N., Diverrès, D. and Gallois, F., 2010. Origin of CO2 undersaturation in the western tropical Atlantic. Tellus B: Chemical and Physical Meteorology, 62(5), pp.595-607.
class Lefevre2010_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Lefevre2010_at: Lf10(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [265.0, 58.1]; #intersept, salinity slope, eq 1
        self.coefsUncertainty = [18.0, 0.5]; #Uncertainty reported for the coefficients, see eq 1
        self.rmsd = 11.6; #See eq 1
        self.r = 0.997;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-56, -30)]; #See Lefevre2010 fig 1 and section 4
        self.includedRegionsLats = [(-10, 12)]; #See Lefevre2010 section 2

        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (18.0, 36.5), #See Lefevre2010 fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (25.0+273.15, 30.0+273.15), #See Lefevre2010 fig 4
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See eq 1
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;
    

#Millero, F.J., Lee, K. and Roche, M., 1998. Distribution of alkalinity in the surface waters of the major oceans. Marine Chemistry, 60(1-2), pp.111-130.
class Millero1998_at(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Millero1998_at: M98(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SST", "SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        #import region mask data
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SSS": (33.75, 36), #See fig 7
                           };
    
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    def _zone1(self, data):
        coefs = [2291.0]; #intersept, see table 4
        rmsd = 4.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
        
        #equation from table 4
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone2(self, data):
        coefs = [2291.0, -2.69, -0.046]; #intersept, SST, SST^2 see table 4
        rmsd = 5.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 0+273.15) &
                              (dataToUse["SST"] < 20+273.15)
                              ];
        
        SST = dataToUse["SST"]-(20+273.15);
        
        #equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone3(self, data):
        coefs = [2300, -2.94, -0.058]; #intersept, SST, SST^2 see table 4
        rmsd = 5.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15)
                              ];
        
        SST = dataToUse["SST"]-29-273.15;
        
        #equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone4(self, data):
        coefs = [2300.0]; #intersept, SST, SST^2 see table 4
        rmsd = 6.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone4_mask");
        ### No temperature limits specified
        
        #equation from table 4
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone5(self, data):
        coefs = [2300, -7.0, -0.158]; #intersept, SST, SST^2 see table 4
        rmsd = 5.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone5_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 7+273.15) &
                              (dataToUse["SST"] < 20+273.15)
                              ];
        
        SST = dataToUse["SST"]-20-273.15;
        
        #equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]* + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone6a(self, data):
        coefs = [2291, -2.52, 0.056]; #intersept, SST, SST^2 see table 4
        rmsd = 5.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone6_atlantic_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > -1+273.15) &
                              (dataToUse["SST"] < 20+273.15)
                              ];
        
        SST = dataToUse["SST"]-20-273.15;
        
        #equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _zone6p(self, data):
        coefs = [2300, -2.52, 0.056]; #intersept, SST, SST^2 see table 4
        rmsd = 5.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone6_pacific_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > -1+273.15) &
                              (dataToUse["SST"] < 20+273.15)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;


    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput[:] = rmsds[:] = np.nan;
        
        #Perform calculations for each zone
        zoneData, zoneUncertainty, zoneRmsd = self._zone1(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone2(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone3(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone4(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone5(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone6a(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._zone6p(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Lee06_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        #This algorithm gives AT normalised to 35 PSU salinity. Reverse normalisation:
        
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;



#Schneider, A., Wallace, D.W. and Körtzinger, A., 2007. Alkalinity of the Mediterranean sea. Geophysical Research Letters, 34(15).
class Schneider2007_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Schneider2007_at: Sc07(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-285.7, 73.7]; #intersept, salinity slope, equation 1
        self.coefsUncertainty = [114.94, 3.0]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = 8.2; #equation 1
        self.r = 0.98; #equation 1
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #Mediterranean excluding Black Sea, Adriatic Sea and Aegean Sea (fig 1)
        self.includedRegionsLons = [(-5, 12),
                                    (12, 22),
                                    (22, 36),
                                    ];
        self.includedRegionsLats = [(30, 44), #Mediterranean up to Italy
                                    (30, 40), #below adriatic sea
                                    (30, 27), #below Aegean sea
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (37.25, 39.5), #See fig 3a
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5d
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Tait, V.K., Gershey, R.M. and Jones, E.P., 2000. Inorganic carbon in the Labrador Sea: Estimation of the anthropogenic component. Deep Sea Research Part I: Oceanographic Research Papers, 47(2), pp.295-308.
class Tait2000_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Tait2000_at: Ta00(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [785.3, 43.01, 2.046]; #intersept, SSS, SST (as potential temperature), equation 6
        self.coefsUncertainty = [None, None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = None;
        self.r = None;
        self.standardError = 5.4;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #Estimated from fig 1
        self.includedRegionsLons = [(-58, -42),
                                    ];
        self.includedRegionsLats = [(51, 62),
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (32, 35), #See fig 2
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        SST = dataToUse["SST"]-273.15;
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"] + self.coefs[2]*SST; #From Tait2000 eq. 6
        outputUncertaintyDueToInputUncertainty = np.sqrt( (self.coefs[1]*dataToUse["SSS_err"])**2 + (self.coefs[2]*dataToUse["SST_err"])**2 );
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Rivaro, P., Messa, R., Massolo, S. and Frache, R., 2010. Distributions of carbonate properties along the water column in the Mediterranean Sea: Spatial and temporal variations. Marine Chemistry, 121(1-4), pp.236-245.
#Mediterranean split into East and West basins
class Rivaro2010_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Rivaro2010_at: R10(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        #super().__init__(settings); #Call the parent class's initator
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
    
    
    #east Mediterranean basin, iho definition
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    def _east_basin(self, data):
        coefs = [-499.8, 80.04]; #intersept, salinity, see eq 2
        rmsd = None;
        
        #Subset data to only rows valid for this zone
        dataToUse = subset_from_mask(data, self.regionMasks, "east_basin_mask");
        dataToUse = dataToUse[(dataToUse["SSS"] > 38.0) & #fig 2b
                              (dataToUse["SSS"] < 39.6) #fig 2b
                              ];
        
        #Equation 2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #east Mediterranean basin, iho definition
    def _west_basin(self, data):
        coefs = [-1089.3, 95.25]; #intersept, salinity, eq 1
        rmsd = 3.3; #See eq 1
        
        #Subset data to only rows valid for this zone.
        dataToUse = subset_from_mask(data, self.regionMasks, "west_basin_mask");
        dataToUse = dataToUse[(dataToUse["SSS"] > 36.8) & #fig 2a
                              (dataToUse["SSS"] < 38.8) #fig 2a
                              ];
        
        #Equation 1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        print("*** NOTE: Rivaro2010_at missing east basin RMSD so cannot perform calculation in the Eastern basin...");
        #Perform calculations for each basin
        zoneData, zoneUncertainty, zoneRmsd = self._east_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Rivaro2010_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._west_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Rivaro2010_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
        rmsds[zoneData.index] = zoneRmsd;
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;


#Requires DO
#Sasse, T.P., McNeil, B.I. and Abramowitz, G., 2013. A novel method for diagnosing seasonal to inter-annual surface ocean carbon dynamics from bottle data using neural networks. Biogeosciences, 10(6), pp.4319-4340.
class Sasse2013_at(BaseAlgorithm): 
    #String representation of the algorithm
    def __str__(self):
        return "Sasse2013_at: S13(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SST", "SSS", "DO", "SiO4", "PO4"];
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        #import region mask data
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    #subtropical
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    def _zone1(self, data):
        #For coefficients, see supplemental table T2
        coefs = [2064.66, -0.3, -47.57, 1.54, 0.13, -1.12, 10.1]; #intersept, SST, SSS, SSS^2, DO, Si, PO4
        rmsd = 11.0; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "AT_zone1_subtropical");
        #dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["SSS"]**2) + \
                      coefs[4]*(dataToUse["DO"]) + \
                      coefs[5]*(dataToUse["SiO4"]) + \
                      coefs[6]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3] * 2.0*dataToUse["SSS_err"]*dataToUse["SSS"]; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm4 = coefs[4]*dataToUse["DO_err"]; #B*DO
        uterm5 = coefs[5]*dataToUse["SiO4_err"]; #B*SiO4
        uterm6 = coefs[6]*dataToUse["PO4_err"]; #B*PO4
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #Equatorial Pacific
    def _zone2(self, data):
        #For coefficients, see supplemental table T2
        coefs = [1142.6, -1.39, 0.96, 0.14, -3.51]; #intersept, SST, SSS^2, DO, PO4
        rmsd = 9.4; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "AT_zone2_equatorial_pacific");
        #dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]**2) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*dataToUse["SSS"]; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3]*dataToUse["DO_err"]; #B*DO
        uterm4 = coefs[4]*dataToUse["PO4_err"]; #B*PO4
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2);
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #North Atlantic
    def _zone3(self, data):
        #For coefficients, see supplemental table T2
        coefs = [1543.52, -4.78, 0.64, 0.04, -0.29, -9.04, 0.13]; #intersept, SST, SSS^2, DO, Si, PO4, SSS*SST
        rmsd = 7.9; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "AT_zone3_north_atlantic");
        #dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]**2) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["SiO4"]) + \
                      coefs[5]*(dataToUse["PO4"]) + \
                      coefs[6]*(dataToUse["SSS"]*SST);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SSS_err"]*dataToUse["SSS"]; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3]*dataToUse["DO_err"]; #B*DO
        uterm4 = coefs[4]*dataToUse["SiO4_err"]; #B*SiO4
        uterm5 = coefs[5]*dataToUse["PO4_err"]; #B*PO4
        uterm6 = coefs[6] * ((dataToUse["SSS_err"]/dataToUse["SSS"]) + (dataToUse["SST_err"]/SST))  * (dataToUse["SSS"]*SST); #B*SSS*SST
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
        
    
    #North Pacific
    def _zone4(self, data):
        #For coefficients, see supplemental table T2
        coefs = [721.6, 44.31, 0.09, -7.81, 9.97, 0.24]; #intersept, SSS, DO, Si, PO4, SSS*SiO4
        rmsd = 14.8; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "AT_zone4_north_pacific");
        #dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
    
        
        #See supplemental table T2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]) + \
                      coefs[2]*(dataToUse["DO"]) + \
                      coefs[3]*(dataToUse["SiO4"]) + \
                      coefs[4]*(dataToUse["PO4"]) + \
                      coefs[5]*(dataToUse["DO"]*dataToUse["SiO4"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SSS_err"]; #B*SSS
        uterm2 = coefs[2]*dataToUse["DO_err"]; #B*DO
        uterm3 = coefs[3]*dataToUse["SiO4_err"]; #B*SiO4
        uterm4 = coefs[4]*dataToUse["PO4_err"]; #B*PO4
        uterm5 = coefs[5] * ((dataToUse["DO_err"]/dataToUse["DO"]) + (dataToUse["SiO4_err"]/dataToUse["SiO4"]))  * (dataToUse["DO"]*dataToUse["SiO4"]); #B*DO*SiO4
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #Southern Ocean
    def _zone5(self, data):
        #For coefficients, see supplemental table T2
        coefs = [7661.04, -1.46, -362.53, 5.86, 0.54, -12.17, -6.56, 0.08, 0.44, 0.01]; #intersept, SST, SSS, SSS^2, DO, Si, PO4, SSS*SST, SSS*SiO4, DO*SiO4
        rmsd = 9.4; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "AT_zone5_southern_ocean");
        #dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15)];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T2
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["SSS"]**2) + \
                      coefs[4]*(dataToUse["DO"]) + \
                      coefs[5]*(dataToUse["SiO4"]) + \
                      coefs[6]*(dataToUse["PO4"]) + \
                      coefs[7]*(dataToUse["SSS"]*SST) + \
                      coefs[8]*(dataToUse["SSS"]*dataToUse["SiO4"]) + \
                      coefs[9]*(dataToUse["DO"]*dataToUse["SiO4"]);
        
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3] * 2.0*dataToUse["SSS_err"]*dataToUse["SSS"]; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm4 = coefs[4]*dataToUse["DO_err"]; #B*DO
        uterm5 = coefs[5]*dataToUse["SiO4_err"]; #B*SiO4
        uterm6 = coefs[6]*dataToUse["PO4_err"]; #B*PO4
        uterm7 = coefs[7] * ((dataToUse["SSS_err"]/dataToUse["SSS"]) + (dataToUse["SST_err"]/SST))  * (dataToUse["SSS"]*SST); #B*SSS*SST
        uterm8 = coefs[8] * ((dataToUse["SSS_err"]/dataToUse["SSS"]) + (dataToUse["SiO4_err"]/dataToUse["SiO4"]))  * (dataToUse["SSS"]*dataToUse["SiO4"]); #B*SSS*SiO4
        uterm9 = coefs[9] * ((dataToUse["DO_err"]/dataToUse["DO"]) + (dataToUse["SiO4_err"]/dataToUse["SiO4"]))  * (dataToUse["DO"]*dataToUse["SiO4"]); #B*DO*SiO4
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 + uterm7**2 + uterm8**2 + uterm9**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;


    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Innernal function used to run, check and assign values for each equation/zone
        def run_single_zone(function, data, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds):
            zoneData, zoneUncertainty, zoneRmsd = function(data);
            if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
                raise RuntimeError("Overlapping zones in Sasse2013_at. Something has done wrong!");
            modelOutput[zoneData.index] = zoneData;
            outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
            rmsds[zoneData.index] = zoneRmsd;
        
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each zone
        run_single_zone(self._zone1, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone2, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone3, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone4, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone5, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;


#Requires DO, SiO4, PO4
#Sasse, T.P., McNeil, B.I. and Abramowitz, G., 2013. A novel method for diagnosing seasonal to inter-annual surface ocean carbon dynamics from bottle data using neural networks. Biogeosciences, 10(6), pp.4319-4340.
class Sasse2013_global_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Sasse2013_global_at: S13g(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST", "DO", "SiO4", "PO4"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        #For coefficients see supplemental table T2
        self.coefs = [1972.44, -12.78, -33.44, 1.19, 0.16, 0.39, 6.89, 0.37]; #intersept, SST, SSS, SSS^2, DO, SiO4, PO4, SSS*SST
        self.coefsUncertainty = [None, None, None, None, None, None, None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = 11.1; #RSE from testing dataset (see table 3)
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-180, 180)]; #Global
        self.includedRegionsLats = [(-90, 60)]; #Not higher than 60N, see grey region in Fig 3
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T2
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*SST + \
                      self.coefs[2]*(dataToUse["SSS"]) + \
                      self.coefs[3]*(dataToUse["SSS"]**2) + \
                      self.coefs[4]*(dataToUse["DO"]) + \
                      self.coefs[5]*(dataToUse["SiO4"]) + \
                      self.coefs[6]*(dataToUse["PO4"]) + \
                      self.coefs[7]*(dataToUse["SSS"]*SST);
        
        
        #Propagate uncertainty
        uterm1 = self.coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = self.coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = self.coefs[3] * 2.0*dataToUse["SSS_err"]*dataToUse["SSS"]; #B*SSS^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm4 = self.coefs[4]*dataToUse["DO_err"]; #B*DO
        uterm5 = self.coefs[5]*dataToUse["SiO4_err"]; #B*SiO4
        uterm6 = self.coefs[6]*dataToUse["PO4_err"]; #B*PO4
        uterm7 = self.coefs[7] * ((dataToUse["SSS_err"]/dataToUse["SSS"]) + (dataToUse["SST_err"]/SST))  * (dataToUse["SSS"]*SST); #B*SSS*SST

        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 + uterm7**2 );
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Requires NO3- (as fitted to potential alkalinity)
#Takahashi, T., Sutherland, S.C., Chipman, D.W., Goddard, J.G., Ho, C., Newberger, T., Sweeney, C. and Munro, D.R., 2014. Climatological distributions of pH, pCO2, total CO2, alkalinity, and CaCO3 saturation in the global surface ocean, and temporal changes at selected locations. Marine Chemistry, 164, pp.95-125.
class Takahashi2013_at(BaseAlgorithm): 
    #String representation of the algorithm
    def __str__(self):
        return "Takahashi2013_at: TS13(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "NO3"]; #NO3 used to convert from potential alkalinity to TA
    @staticmethod
    def output_name():
        return "AT";
    
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        #import region mask data
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    #North Atlantic Drift
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions (except osoda_mediterranean).
    #      This algorithm should not be used in the Mediterranean region!
    def _zone7(self, data):
        #For coefficients see table 1
        coefs = [733.0, 45.30]; #intersept, SSS
        rmsd = 6.5; #table 1
        
        #Subset data to only rows valid for this zone. See table 1 and fig 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone7_north_atlantic_drift");
        dataToUse = dataToUse[(dataToUse["SSS"] > 31) & #See fig 4
                              (dataToUse["SSS"] < 36.5)];
        
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #Central Atlantic
    def _zone8(self, data):
        #For coefficients see table 1
        coefs = [270.9, 58.25]; #intersept, SSS
        rmsd = 12.6; #table 1
        
        #Subset data to only rows valid for this zone. See table 1 and fig 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone8_central_atlantic");
        dataToUse = dataToUse[(dataToUse["SSS"] > 31) & #See fig 4
                              (dataToUse["SSS"] < 38)];
        
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;


    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Innernal function used to run, check and assign values for each equation/zone
        def run_single_zone(function, data, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds):
            zoneData, zoneUncertainty, zoneRmsd = function(data);
            if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
                raise RuntimeError("Overlapping zones in Takahashi2013. Something has done wrong!");
            modelOutput[zoneData.index] = zoneData;
            outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
            rmsds[zoneData.index] = zoneRmsd;
        
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each zone
        run_single_zone(self._zone7, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone8, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        
        #Convert from potential alkalinity to TA
        #PALK = TA + NO3-, so calculate TA using:
        modelOutput = modelOutput - dataToUse["NO3"];
        outputUncertaintyDueToInputUncertainty = np.sqrt( outputUncertaintyDueToInputUncertainty**2 + dataToUse["NO3_err"]**2);
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;



#Ternon, J.F., Oudot, C., Dessier, A. and Diverres, D., 2000. A seasonal tropical sink for atmospheric CO2 in the Atlantic Ocean: the role of the Amazon River discharge. Marine Chemistry, 68(3), pp.183-201.
class Ternon2000_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Ternon2000_at: Te00(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [272.3, 58.85]; #intersept, salinity slope, fig 5d
        self.coefsUncertainty = [9.5, 0.29]; #Uncertainty reported for the coefficients, fig5d
        self.rmsd = 20.1; #Ternon2000 Fig 5d, Peter calles this 'rms2' - why?
        self.r = 0.998; #Ternon2000 Fig 5d
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-55.0, -35.0)]; #See Ternon2000 fig 1a
        self.includedRegionsLats = [(-7.5, 7.5)]; #See Ternon2000 fig 1a
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (17.0, 37.0), #Algorithm valid for this range
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (26, 30), #See Ternon2000 fig 2a
                           "SSS": (17, 37), #See Ternon2000 fig 2b, fig 3b, fig 5
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5d
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Touratier, F. and Goyet, C., 2009. Decadal evolution of anthropogenic CO2 in the northwestern Mediterranean Sea from the mid-1990s to the mid-2000s. Deep Sea Research Part I: Oceanographic Research Papers, 56(10), pp.1708-1716.
class Touratier2009_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Touratier2009_at: TG09(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-1238.4, 99.26]; #intersept, SSS equation 4
        self.coefsUncertainty = [4.5, 0.0]; #Uncertainty reported for the coefficients, equation 4
        self.rmsd = None;
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(7, 9)]; #See fig 1, DYFAMED site and surrounding region
        self.includedRegionsLats = [(42.5, 44)]; #See fig 1, DYFAMED site and surrounding region
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (37.8, 38.8), #See fig 4
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From equation 4
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Touratier, F. and Goyet, C., 2011. Impact of the Eastern Mediterranean Transient on the distribution of anthropogenic CO2 and first estimate of acidification for the Mediterranean Sea. Deep Sea Research Part I: Oceanographic Research Papers, 58(1), pp.1-15.
class Touratier2011_at(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Touratier2011_at: TG11(at)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"];
    @staticmethod
    def output_name():
        return "AT";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-0.0000657, 0.0177, 5.93]; #intersept, SSS, SST (as potential temperature), equation 4
        self.coefsUncertainty = [None, None, None]; #Uncertainty reported for the coefficients, equation 4
        self.rmsd = None;
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        #Mediterranean excluding Black Sea, Adriatic Sea and Aegean Sea, see fig 1
        self.includedRegionsLons = [(-5, 12),
                                    (12, 22),
                                    (22, 36),
                                    ];
        self.includedRegionsLats = [(30, 44), #Mediterranean up to Italy
                                    (30, 40), #below adriatic sea
                                    (30, 27), #below Aegean sea
                                    ];
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #equation 4
        modelOutput = 1.0 / (self.coefs[0] + (self.coefs[1])/dataToUse["SSS"] - (self.coefs[2]*np.log(dataToUse["SST"]))/(dataToUse["SST"]**2) );
        
        
        #Propagate uncertainty
        uterm1 = self.coefs[1] / dataToUse["SST_err"]; #B/SST
        
        #Split term2 into numerator and denominator
        ulogSst = dataToUse["SST_err"] / dataToUse["SST"]; #Numerator. d/dx of ln(x) = 1/x
        usstSquared = 2.0*dataToUse["SST_err"]*dataToUse["SST"]; #Denominator. SSS^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm2 = self.coefs[2] * ((ulogSst/np.log(dataToUse["SST"])) + (usstSquared/dataToUse["SST"])) * (np.log(dataToUse["SST"]/(dataToUse["SST"]**2))); #B*ln(SST)/SST^2
        
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

