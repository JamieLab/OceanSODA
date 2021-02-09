#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 08:57:49 2020

@author: tom holding
"""

import pandas as pd;
import numpy as np;

from .base_algorithm import BaseAlgorithm;
from .utilities import subset_from_mask;


#Bakker, D.C., de Baar, H.J. and de Jong, E., 1999. The dependence on temperature and salinity of dissolved inorganic carbon in East Atlantic surface waters. Marine Chemistry, 65(3-4), pp.263-280.
#Low salinity region, North of the Congo outflow, north region 1
class Bakker1999_lcr1_dic(BaseAlgorithm):  
    #String representation of the algorithm
    def __str__(self):
        return "Bakker1999_lcr1_dic: B95_lcr1(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"]; #SST for restricted range
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [-1575, 99.3]; #intersept, SSS, Table 3
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = 8.6; #See table 3
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-30, -8)]; #Approximated from fig 1
        self.includedRegionsLats = [(7.5, 15.0)]; #From table 3
        #TODO: Mask away pacific corner
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (34, 36.5), #See fig 3
                               "SST": (20+273.15, 30+273.15), #See fig 3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };


    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #equation from table 3
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*dataToUse["SSS"];
        
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;
    

#Bakker, D.C., de Baar, H.J. and de Jong, E., 1999. The dependence on temperature and salinity of dissolved inorganic carbon in East Atlantic surface waters. Marine Chemistry, 65(3-4), pp.263-280.
#Low salinity region, North of the Congo outflow, north region 2
class Bakker1999_lcr2_dic(BaseAlgorithm):  
    #String representation of the algorithm
    def __str__(self):
        return "Bakker1999_lcr2_dic: B95_lcr2(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"]; #SST for restricted range
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [-1415, 95.5]; #intersept, SSS, Table 3
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = 5.5; #See table 3
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-30, -0)]; #Approximated from fig 1
        self.includedRegionsLats = [(-1, 4.0)]; #From table 3
        #TODO: Mask away pacific corner
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (34, 36.5), #See fig 3
                               "SST": (20+273.15, 30+273.15), #See fig 3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #equation from table 3
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*dataToUse["SSS"];
        
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Bakker, D.C., de Baar, H.J. and de Jong, E., 1999. The dependence on temperature and salinity of dissolved inorganic carbon in East Atlantic surface waters. Marine Chemistry, 65(3-4), pp.263-280.
#Congo outflow region
class Bakker1999_outflow_dic(BaseAlgorithm):  
    #String representation of the algorithm
    def __str__(self):
        return "Bakker1999_outflow_dic: B95_outflow(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"]; #SST for restricted range
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [-1575, 99.3]; #intersept, SSS, Table 3
        self.coefsUncertainty = [None, None, None, None, None, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = 14.1; #See table 3
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(0.0, 5.0)]; #Approximated from fig 1
        self.includedRegionsLats = [(-10.0, -5.0)]; #From table 3
        #TODO: Mask away pacific corner?
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (33, 36), #See fig 3
                               "SST": (26.2+273.15, 29.9+273.15), #See fig 3
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
    
    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #equation from table 3
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*dataToUse["SSS"];
        
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Brewer, P.G., Glover, D.M., Goyet, C. and Shafer, D.K., 1995. The pH of the North Atlantic Ocean: Improvements to the global model for sound absorption in seawater. Journal of Geophysical Research: Oceans, 100(C5), pp.8761-8776.
#Implementation of equation 10 (<250m depth)
class Brewer1995_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Brewer1995_dic: B95(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST", "DO", "PO4", "NO3"];

    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [944.0, 35.584, -7.099, -0.464, 138.0, -4.62]; #intersept, SSS, SST (as potential temperature), DO, PHO4, NO3 see equation 10
        self.coefsUncertainty = [None, None, None, None, None, None]; #Uncertainty reported for the coefficients, see fig 11
        self.rmsd = None;
        self.r = 0.940**0.5; #equation 10
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-100, 30)]; #See fig 1a
        self.includedRegionsLats = [(0, 80)]; #See fig 1a
        #TODO: Mask away pacific corner
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (31, 37.5), #See fig 4
                               "SST": (-2+273.15, 28+273.15), #See fig 5
                               "PO4": (0, 2.2), #See fig 5 (continued, second page)
                               "NO3": (0, 34), #See fig 5 (continues, 3rd page)
                               "DO": (100, 420), #See fig 5...
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #equation 10
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*dataToUse["SSS"] + \
                      self.coefs[2]*(dataToUse["SST"]-273.15) + \
                      self.coefs[3]*dataToUse["DO"] + \
                      self.coefs[4]*dataToUse["PO4"] + \
                      self.coefs[5]*dataToUse["NO3"];
        
        #equation 10
        outputUncertaintyDueToInputUncertainty = np.sqrt( (self.coefs[1]*dataToUse["SSS_err"])**2 + \
                                          (self.coefs[2]*dataToUse["SST_err"])**2 + \
                                          (self.coefs[3]*dataToUse["DO_err"])**2 + \
                                          (self.coefs[4]*dataToUse["PO4_err"])**2 + \
                                          (self.coefs[5]*dataToUse["NO3_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Cooley, S.R. and Yager, P.L., 2006. Physical and biological contributions to the western tropical North Atlantic Ocean carbon sink formed by the Amazon River plume. Journal of Geophysical Research: Oceans, 111(C8).
class Cooley2006a_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Cooley2006a_dic: C06a(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [378.5, 44.1]; #intersept, salinity slope
        self.coefsUncertainty = [120.4, 3.8]; #Uncertainty reported for the coefficients
        self.rmsd = 11.6; #See section 3.6
        self.r = 912; #See Table2 "plume"
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-59, -40)]; #See Cooley2006 section 2.1
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
class CopinMontegut1993_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "CopinMontegut1993_dic: CM93(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-1628, 102.3]; #intersept, salinity slope, see second equation in Results
        self.coefsUncertainty = [66.0, 1.8]; #Uncertainty reported for the coefficients, see second equation in Results
        self.rmsd = None;
        self.r = 0.985; #See first paragraph of Results
        
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
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See second equation in Results
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Copin-Montégut, C. and Bégovic, M., 2002. Distributions of carbonate properties and oxygen along the water column (0–2000 m) in the central part of the NW Mediterranean Sea (Dyfamed site): influence of winter vertical mixing on air–sea CO2 and O2 exchanges. Deep Sea Research Part II: Topical Studies in Oceanography, 49(11), pp.2049-2066.
#Relationship with salinity above the salinity maximum
class CopinMontegut2002a_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "CopinMontegut2002a_dic: CM02a(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-3662.6, 155.17]; #intersept, salinity slope, see second unnamed equation of seciton 3.4
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = None;
        self.r = None; #See second unnamed equation of seciton 3.4
        
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
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See second unnamed equation of seciton 3.4
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;

#Copin-Montégut, C. and Bégovic, M., 2002. Distributions of carbonate properties and oxygen along the water column (0–2000 m) in the central part of the NW Mediterranean Sea (Dyfamed site): influence of winter vertical mixing on air–sea CO2 and O2 exchanges. Deep Sea Research Part II: Topical Studies in Oceanography, 49(11), pp.2049-2066.
#Relationship with salinity below the salinity maximum
class CopinMontegut2002b_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "CopinMontegut2002b_dic: CM02b(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-555.2, 74.53]; #intersept, salinity slope, see first unnamed equation of seciton 3.4
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = None;
        self.r = 0.81; #See first unnamed equation of seciton 3.4
        
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
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #See first unnamed equation of seciton 3.4
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;
    
    

#Gemayel, E., Hassoun, A.E.R., Benallal, M.A., Goyet, C., Rivaro, P., Abboud-Abi Saab, M., Krasakopoulou, E., Touratier, F. and Ziveri, P., 2015. Climatological variations of total alkalinity and total dissolved inorganic carbon in the Mediterranean Sea surface waters. Earth System Dynamics, 6(2), pp.789-800.
class Gemayel2015_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Gemayel2015_dic: Ge15(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [2234, #intercept
                      38.15, #SSS
                      -4.48, #SSS^2
                      -1.10, #SSS^3
                      -14.38, #SST
                      9.62, #SST^2
                      -4.61, #SST^3
                      -1.43, #SSS*SST
                      3.53, #SST*SSS^2
                      1.47, #SST^2*SSS
                      ]; #See eq 2
        self.coefsUncertainty = [None, None, None, None, None, None, None, None, None, None]; #Uncertainty reported for the coefficients
        self.rmsd = 16.2; #See table 4 NOTE: This is validation RMSD, training RMSD is 10.6
        self.r = np.sqrt(0.90); #See eq 2
        
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
        self.restrictRanges = {"SST": (13.0+273.15, None), #See text below eq 2
                               "SSS": (36.3, 39.65), #See text below eq 2
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #Calculate the uncertainty on the main calculation here, returns the uncertainty on the model output.
    def _uncertainty_kernal(self, dataToUse):
        SSS = dataToUse["SSS"]-38.2; #Simplify things by applying offsets / unit conversion and storing it
        SST = dataToUse["SST"]-273.15-17.7; #Simplify things by applying offsets / unit conversion and storing it
        SSS_err = dataToUse["SSS_err"];
        SST_err = dataToUse["SST_err"];
        
        #calculate each uncertainty term seperately then add in quadrature
        uterm1 = self.coefs[1] * SSS_err; #B*SSS
        uterm2 = self.coefs[2] * 2*SSS_err*SSS; # B*SSS^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = self.coefs[3] * 3*SSS_err*(SSS**2); # B*SSS^3: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm4 = self.coefs[4] * SST_err; #B*SST
        uterm5 = self.coefs[5] * 2*SST_err*SST; # B*SSS^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm6 = self.coefs[6] * 3*SST_err*(SST**2); # B*SSS^3: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm7 = self.coefs[7] * ((SSS_err/SSS) + (SST_err/SST)) * (SSS*SST); #B*SSS*SST
        uterm8 = self.coefs[8] * ((2*SSS/SSS_err) + (SST_err/SST)) * ((SSS**2)*SST); #B*(SSS^2)*SST
        uterm9 = self.coefs[9] * ((SSS/SSS_err) + (2*SST_err/SST)) * (SSS*(SST**2)); #B*SSS*(SST^2)
        
        #Add absolute uncertainty terms in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( (uterm1)**2 + \
                                          (uterm2)**2 + \
                                          (uterm3)**2 + \
                                          (uterm4)**2 + \
                                          (uterm5)**2 + \
                                          (uterm6)**2 + \
                                          (uterm7)**2 + \
                                          (uterm8)**2 + \
                                          (uterm9)**2 );
        return outputUncertaintyDueToInputUncertainty;

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        SSS = dataToUse["SSS"]-38.2; #Simplify things by applying offsets / unit conversion and storing it
        SST = dataToUse["SST"]-273.15-17.7; #Simplify things by applying offsets / unit conversion and storing it
        
        #From equation 2
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*SSS + \
                      self.coefs[2]*SSS**2 + \
                      self.coefs[3]*SSS**3 + \
                      self.coefs[4]*SST + \
                      self.coefs[5]*SST**2 + \
                      self.coefs[6]*SST**3 + \
                      self.coefs[7]*SSS*SST + \
                      self.coefs[8]*(SSS**2)*SST + \
                      self.coefs[9]*SSS*(SST**2);
        
        outputUncertaintyDueToInputUncertainty = self._uncertainty_kernal(dataToUse);
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Hassoun, A.E.R., Gemayel, E., Krasakopoulou, E., Goyet, C., Saab, M.A.A., Ziveri, P., Touratier, F., Guglielmi, V. and Falco, C., 2015. Modeling of the total alkalinity and the total inorganic carbon in the Mediterranean Sea.
#Whole mediterranean as one region (surface 0-25m) - see table 2 eq2
class Hassoun2015_full_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Hassoun2005_full_dic: Ha15_full(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [-198, 63.65]; #intersept, salinity slope, see table 2 eq2
        self.coefsUncertainty = [None, None]; #Uncertainty reported for the coefficients
        self.rmsd = 18.0; #See table 1 eq2
        self.r = 0.93; #See table 1 eq2
        
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
#Mediterranean split into East and West basins (surface 0-25m) - see table 2 eq 5 and 7
class Hassoun2015_basins_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Hassoun2005b_dic: Ha15b(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    def __init__(self, settings):
        #super().__init__(settings); #Call the parent class's initator
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
    
    
    #east Mediterranean basin, iho definition
    #Note: Multiple models for a single algorithm is ok here because together spatial extent OceanSODA Mediterranean region and do not overlap any other regions
    def _east_basin(self, data):
        coefs = [-292.6, 66.0]; #intersept, salinity, see table 2 eq 7
        rmsd = 20.0; #See table 1 eq 9
        #r = 0.68;
        
        #Subset data to only rows valid for this zone. See Table 2
        dataToUse = subset_from_mask(data, self.regionMasks, "east_basin_mask");
        
        #Equation from table 1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]);
        
        outputUncertaintyDueToInputUncertainty = coefs[1]*dataToUse["SSS_err"];
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #east Mediterranean basin, iho definition
    def _west_basin(self, data):
        coefs = [-583.5, 74.0]; #intersept, salinity, see table 1 eq 5
        rmsd = 14.0; #See table 1 eq 5
        #r = 0.94;
        
        #Subset data to only rows valid for this zone. See Table 2
        dataToUse = subset_from_mask(data, self.regionMasks, "west_basin_mask");
        
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
        modelRmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each basin
        zoneData, zoneUncertainty, zoneRmsd = self._east_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Hassoun2015_basins_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneUncertainty.index] = zoneUncertainty;
        modelRmsds[zoneData.index] = zoneRmsd;
        
        zoneData, zoneUncertainty, zoneRmsd = self._west_basin(dataToUse);
        if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
            raise RuntimeError("Overlapping zones in Hassoun2015_basins_at. Something has done wrong!");
        modelOutput[zoneData.index] = zoneData;
        outputUncertaintyDueToInputUncertainty[zoneUncertainty.index] = zoneUncertainty;
        modelRmsds[zoneData.index] = zoneRmsd;
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty; #Update the instance's rmsd to reflect the computation just carried out.
        return modelOutput, outputUncertaintyDueToInputUncertainty, modelRmsds;



#Lee, K., Wanninkhof, R., Feely, R.A., Millero, F.J. and Peng, T.H., 2000. Global relationships of total inorganic carbon with temperature and nitrate in surface seawater. Global Biogeochemical Cycles, 14(3), pp.979-994.
class Lee2000_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Lee2000_dic: L00(dic)";

    #common names of input and output variables (see global_settings for definitions of these)
    @staticmethod
    def input_names():
        return ["SST", "SSS", "NO3"];
    @staticmethod
    def output_name():
        return "DIC";
    
    def __init__(self, settings):
        #super().__init__(settings); #Call the parent class's initator
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
        
        self.northernSummer = np.array([5, 6, 7, 8, 9]); #May-Sept
        self.northernWinter = np.array([10, 11, 12, 1, 2, 3, 4]); #Oct - April
        self.southernSummer = np.array([10, 11, 12, 1, 2, 3, 4]); #Oct - April
        self.southernWinter = np.array([5, 6, 7, 8, 9]); #May-Sept
    
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    def _equation1(self, data):
        coefs = [1940.0, -10.327, -0.451, 7.829]; #intersept, SST, SST^2, NO3. Table 3
        rmsd = 8.1; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation1_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 18+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (dataToUse["NO3"] > 0.5)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 3
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2. Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2);
        
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _equation2(self, data):
        coefs = [1950.0, -7.604, -0.178, 6.883]; #intersept, SST, SST^2, NO3. Table 3 eq 2
        rmsd = 7.0 #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation2_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 18+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (dataToUse["NO3"] > 0.5)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 3
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2);
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 1 equation 3 for SST between 20 and 29 C
    def _equation3a(self, data):
        coefs = [1940.0, -11.976, -0.513]; #intersept, SST, SST^2, Table 3
        rmsd = 6.2; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation3_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        SST = dataToUse["SST"]-273.15-29
        
        #Equation from table 3
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 1 equation 3 for SST > 29 C
    def _equation3b(self, data):
        coefs = [1940.0]; #intersept only, Table 3
        rmsd = 6.5; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation3_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        #Equation from table 3
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Convert to non-normalised DIC
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;

    
    #This is zone 1 equation 4 for SST between 20 and 29 C
    def _equation4a(self, data):
        coefs = [1950.0, -5.549, 0.125]; #intersept, SST, SST^2, Table 3
        rmsd = 6.3; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation4_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 3
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    
    
    #This is zone 1 equation 4 for SST > 29 C
    def _equation4b(self, data):
        coefs = [1950.0]; #intersept only, Table 3
        rmsd = 7.0; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation4_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        #Equation from table 3
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Convert to non-normalised DIC
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 1 equation 5 for SST between 20 and 29 C
    def _equation5a(self, data):
        coefs = [1940.0, -33.385, 2.407]; #intersept, SST, SST^2, Table 3
        rmsd = 7.5; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation5_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 3
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);

        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 1 equation 5 for SST > 29 C
    def _equation5b(self, data):
        coefs = [1940.0]; #intersept only, Table 3
        rmsd = 6.5; #See table 3
        
        #Subset data to only rows valid for this zone. See Table 3
        dataToUse = subset_from_mask(data, self.regionMasks, "zone1_equation5_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15) &
                              (dataToUse["NO3"] < 0.5)
                              ];
        
        #Equation from table 3
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Convert to non-normalised DIC
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 6 for SST between 20 and 29 C (summer)
    def _equation6s(self, data):
        coefs = [1940.0, -3.039, 0.494]; #intersept, SST, SST^2, Table 4
        rmsd = 7.5; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation6_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 6 for SST between 20 and 29 C (winter)
    def _equation6w(self, data):
        coefs = [1940.0, -1.003, 0.372]; #intersept, SST, SST^2, Table 4
        rmsd = 4.9; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation6_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
                      
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 6 for SST > 29 C (all seasons)
    def _equation6b(self, data):
        coefs = [1940.0]; #intersept,Table 4
        rmsd = 6.7; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation6_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15)
                              ];
        
        #Equation from table 4
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 7 for SST between 20 and 29 C (summer)
    def _equation7s(self, data):
        coefs = [1950.0, 0.188, 0.725]; #intersept, SST, SST^2, Table 4
        rmsd = 8.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation7_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
                      
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 7 for SST between 20 and 29 C (winter)
    def _equation7w(self, data):
        coefs = [1950.0, 2.570, 0.640]; #intersept, SST, SST^2, Table 4
        rmsd = 5.2; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation7_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
                      
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 2 equation 7 for SST > 29 C (all seasons)
    def _equation7b(self, data):
        coefs = [1950.0]; #intersept,Table 4
        rmsd = 8.0; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation7_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15)
                              ];
        
        #Equation from table 4
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #This is zone 2 equation 8 for SST between 20 and 29 C
    def _equation8a(self, data):
        coefs = [1940.0, 1.842, 0.468]; #intersept, SST, SST^2, Table 4
        rmsd = 6.2; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation8_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 20+273.15) &
                              (dataToUse["SST"] < 29+273.15)
                              ];
        
        SST = dataToUse["SST"]-273.15-29;
        
        #Equation from table 4
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2);
                      
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #This is zone 2 equation 7 for SST > 29 C
    def _equation8b(self, data):
        coefs = [1940.0]; #intersept,Table 4
        rmsd = 6.1; #See table 4
        
        #Subset data to only rows valid for this zone. See Table 4
        dataToUse = subset_from_mask(data, self.regionMasks, "zone2_equation8_mask");
        dataToUse = dataToUse[(dataToUse["SST"] > 29+273.15)
                              ];
        
        #Equation from table 4
        modelOutput = pd.Series([coefs[0]]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([0.0]*len(dataToUse), index=dataToUse.index);
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 3 equation 9 for SST < 20 C (summer)
    def _equation9s(self, data):
        coefs = [2010.0, -8.633, -0.036, -0.279]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 6.9; #See table 5
        
        #Subset data to only rows valid for this zone. See Table 6
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_equation9_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 3 equation 9 for SST < 20 C (winter)
    def _equation9w(self, data):
        coefs = [1980.0, -14.680, -0.297, -1.152]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 7.5; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_equation9_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 3 equation 10 for SST < 20 C (summer)
    def _equation10s(self, data):
        coefs = [2010, -4.262, -0.013, 5.054]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 5.9; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 6
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_equation10_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 3 equation 10 for SST < 20 C (winter)
    def _equation10w(self, data):
        coefs = [1980.0, -10.864, -0.311, 4.235]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 6.7; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone3_equation10_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 4 equation 11 for SST < 20 C (summer)
    def _equation11s(self, data):
        coefs = [2010, -7.805, 0.069, 3.891]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 7.8; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 6
        dataToUse = subset_from_mask(data, self.regionMasks, "zone4_equation11_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 4 equation 11 for SST < 20 C (winter)
    def _equation11w(self, data):
        coefs = [1980, -13.199, -0.172, 3.983]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 6.7; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone4_equation11_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 4 equation 12 for SST < 20 C (summer)
    def _equation12s(self, data):
        coefs = [2010, -7.415, 0.024, 2.343]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 6.2; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 6
        dataToUse = subset_from_mask(data, self.regionMasks, "zone5_equation12_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #This is zone 4 equation 12 for SST < 20 C (winter)
    def _equation12w(self, data):
        coefs = [1980, -12.884, -0.112, 1.365]; #intersept, SST, SST^2, NO3, Table 5
        rmsd = 6.8; #See table 5
        
        
        #Subset data to only rows valid for this zone. See Table 5
        dataToUse = subset_from_mask(data, self.regionMasks, "zone5_equation12_mask");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        
        dataToUse = dataToUse[(dataToUse["SST"] < 20+273.15) &
                              (seasonalIndices)
                              ];
        
        SST = dataToUse["SST"]-273.15-20;
        
        #Equation from table 5
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(SST**2) + \
                      coefs[3]*dataToUse["NO3"];
        
        #Calculate each uncertainty term serparately, for simplicity
        uterm1 = coefs[1] * dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2] * 2.0*dataToUse["SST_err"]*SST; #B*SST^2: Rearranged form of Taylor eq. 3.10: if x=q^n, then dx = n*u*q^(n-1), where u is uncertainty on q
        uterm3 = coefs[3] * dataToUse["NO3_err"]; #B*NO3
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 );
        
        #Normalised DIC (N_DIC) is N_DIC = DIC*(35/SSS), so DIC = N_DIC/(35/SSS)
        #Convert to non-normalised DIC
        outputUncertaintyDueToInputUncertaintyRatio = (outputUncertaintyDueToInputUncertainty/modelOutput) + (dataToUse["SSS_err"]/dataToUse["SSS"]); #Propagate uncertainty through normalisation
        modelOutput = modelOutput / (35.0/dataToUse["SSS"]);
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertaintyRatio*modelOutput; #two steps to avoid duplicate calculation
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    def _kernal(self, dataToUse):
        #Innernal function used to run, check and assign values for each equation/zone
        def run_single_zone(function, data, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds):
            zoneData, zoneUncertainty, zoneRmsd = function(data);
            if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
                raise RuntimeError("Overlapping zones in Lee00_dic. Something has done wrong!");
            modelOutput[zoneData.index] = zoneData;
            outputUncertaintyDueToInputUncertainty[zoneUncertainty.index] = zoneUncertainty;
            rmsds[zoneData.index] = zoneRmsd;
        
        
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each zone
        run_single_zone(self._equation1, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation2, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation3a, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation3b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation4a, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation4b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation5a, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation5b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation6s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation6w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation6b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation7s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation7w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation7b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation8a, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation8b, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation9s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation9w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation10s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation10w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation11s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation11w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation12s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._equation12w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;




#Lefèvre, N., Diverrès, D. and Gallois, F., 2010. Origin of CO2 undersaturation in the western tropical Atlantic. Tellus B: Chemical and Physical Meteorology, 62(5), pp.595-607.
#Implementation of the Lefevre et al 2010 DIC algorithm, fit in the region of the Amazon plume.
class Lefevre2010_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Lefevre2010: Lf10(dic)";
    
    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [189.0, 50.6]; #intersept, salinity slope. Lefevre2010 Eq. 2
        self.coefsUncertainty = [31.0, 0.9]; #Uncertainty reported for the coefficients. Lefevre2010 Eq. 2
        self.rmsd = 16.2; #Lefevre2010 Eq. 2 (text)
        self.r = 0.999; #Lefevre2010 Eq. 2
        
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
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #Lefevre2010 eq. 2
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Lefèvre, N., Flores Montes, M., Gaspar, F.L., Rocha, C., Jiang, S., De Araújo, M.C. and Ibánhez, J., 2017. Net heterotrophy in the Amazon Continental Shelf changes rapidly to a sink of CO2 in the Outer Amazon Plume. Frontiers in Marine Science, 4, p.278.
class Lefevre2017_dic(BaseAlgorithm):
    #String representation of the algorithm
    def __str__(self):
        return "Lefevre2017: Lf17(dic)";
    
    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [367.75, 45.86]; #intersept, salinity slope. Lefevre2017 Eq. 3
        self.coefsUncertainty = [12.96, 0.85]; #Uncertainty reported for the coefficients. Lefevre2017 Eq. 3
        self.rmsd = 27.1; #Lefevre2017 Eq. 3 (text)
        self.r = 0.998; #Lefevre2017 Eq. 3
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-55, -35)]; #See end of Lefevre2017 Methods/Sampling and figure 1
        self.includedRegionsLats = [(-5.0, 10.0)]; #See end of Lefevre2017 Methods/Sampling and figure 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (1, 37), #See Lefevre2017 fig 4 for upper bound, Section Results/"Variability along the Transect on the Amazon Continental Shelf" for lower bound
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {"SST": (27.8+273.15, None), #All temperatures above 27.8, no upper bound specified but maybe around 30 based on figures?
                           };
    
                           
    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #Lefevre2017 eq. 3
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Ternon, J.F., Oudot, C., Dessier, A. and Diverres, D., 2000. A seasonal tropical sink for atmospheric CO2 in the Atlantic Ocean: the role of the Amazon River discharge. Marine Chemistry, 68(3), pp.183-201.
#Implementation of the Ternon et al 2000 DIC algorithm, fit in the region of the Amazon plume.
class Ternon2000_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Ternon00_dic: Te00(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        self.coefs = [226.8, 49.48]; #intersept, salinity slope
        self.coefsUncertainty = [9.3, 0.28]; #Uncertainty reported for the coefficients
        self.rmsd = 19.7; #Ternon2000 Fig 5c, Peter calles this 'rms2' - why?
        self.r = 0.997; #Ternon2000 Fig 5c
        
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
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #From Ternon2000 fig 5c
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;




#Requires DO, SiO4, PO4
#Sasse, T.P., McNeil, B.I. and Abramowitz, G., 2013. A novel method for diagnosing seasonal to inter-annual surface ocean carbon dynamics from bottle data using neural networks. Biogeosciences, 10(6), pp.4319-4340.
class Sasse2013_dic(BaseAlgorithm): 
    #String representation of the algorithm
    def __str__(self):
        return "Sasse2013_dic: S13(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SST", "SSS", "DO", "NO3", "SiO4", "PO4"];
    @staticmethod
    def output_name():
        return "DIC";
    
    def __init__(self, settings):
        BaseAlgorithm.__init__(self, settings); #Call the parent class's initator
        self.settings = settings;
        
        #import region mask data
        from netCDF4 import Dataset;
        self.regionMasks = Dataset(settings["algorithmSpecificDataPaths"][type(self).__name__], 'r');
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };
        
        self.northernSummer = np.array([5, 6, 7, 8, 9, 10]); #May-Oct
        self.northernWinter = np.array([11, 12, 1, 2, 3, 4]); #Nov-Apr
        self.southernSummer = np.array([11, 12, 1, 2, 3, 4]); #Nov-Apr
        self.southernWinter = np.array([5, 6, 7, 8, 9, 10]); #May-Oct
    
    #Note: Multiple models for a single algorithm is ok here because together spatial extent of each model covers all of the OceanSODA regions
    #north pacific summer
    def _zone1s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1066.9, 24.76, 0.38, 5.07]; #intersept, SSS, DO, NO3
        rmsd = 16.8; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone1_north_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]) + \
                      coefs[2]*(dataToUse["DO"]) + \
                      coefs[3]*(dataToUse["NO3"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SSS_err"])**2 + \
                                          (coefs[2]*dataToUse["DO_err"])**2 + \
                                          (coefs[3]*dataToUse["NO3_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #north pacific winter
    def _zone1w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [868.73, -7.95, 36.54, 4.73, -0.01]; #intersept, SST, SSS, SiO4, SiO4*DO
        rmsd = 16.8; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone1_north_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["SiO4"]) + \
                      coefs[4]*(dataToUse["SiO4"]*dataToUse["DO"]);
        
        #Propagate uncertainty
        utermInteraction = ((dataToUse["SiO4_err"]/dataToUse["SiO4"]) + (dataToUse["DO_err"]/dataToUse["DO"])) * dataToUse["SiO4"]*dataToUse["DO"]; #Absolute uncertainty term for the SiO4*DO interaction
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["SiO4_err"])**2 + \
                                          (utermInteraction)**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #southern ocean summer
    def _zone2s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [698.74, 35.84, 0.25, 0.42, 83.28]; #intersept, SSS, DO, SiO4, PO4
        rmsd = 16.4; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone2_southern_ocean");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]) + \
                      coefs[2]*(dataToUse["DO"]) + \
                      coefs[3]*(dataToUse["SiO4"])+ \
                      coefs[4]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SSS_err"])**2 + \
                                          (coefs[2]*dataToUse["DO_err"])**2 + \
                                          (coefs[3]*dataToUse["SiO4_err"])**2 + \
                                          (coefs[4]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #southern ocean winter
    def _zone2w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1494.14, -48.81, 22.3, -0.31, 0.48, 0.03, 0.92]; #intersept, SST, SSS, DO, SiO4, SST*DO, SST*SSS
        rmsd = 16.4; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone2_southern_ocean");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        SST = dataToUse["SST"]-273.15; #Just simplifies the following equations
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["SiO4"]) + \
                      coefs[5]*(SST*dataToUse["DO"]) + \
                      coefs[6]*(SST*dataToUse["SSS"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3]*dataToUse["DO_err"]; #B*DO
        uterm4 = coefs[4]*dataToUse["SiO4_err"]; #B*SiO4
        uterm5 = coefs[5] * ((dataToUse["SST_err"]/SST) + (dataToUse["DO_err"]/dataToUse["DO"])) * (SST*dataToUse["DO"]); #B*SST*DO
        uterm6 = coefs[6] * ((dataToUse["SST_err"]/SST) + (dataToUse["SSS_err"]/dataToUse["SSS"]))  * (SST*dataToUse["SSS"]); #B*SST*SSS
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #northwest atlantic summer
    def _zone3s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1709.15, 9.13, 9.14]; #intersept, SSS, NO3
        rmsd = 15.5; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone3_northwest_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]) + \
                      coefs[2]*(dataToUse["NO3"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SSS_err"])**2 + \
                                          (coefs[2]*dataToUse["NO3_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #northwest atlantic winter
    def _zone3w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1128.79, -5.93, 28.71, 17.31, -0.05]; #intersept, SST, SSS, NO3, NO3*DO
        rmsd = 15.5; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone3_northwest_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["NO3"]*dataToUse["DO"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3]*dataToUse["NO3_err"]; #B*NO3
        uterm4 = coefs[4] * ((dataToUse["NO3_err"]/dataToUse["NO3"]) + (dataToUse["DO_err"]/dataToUse["DO"])) * (dataToUse["NO3"]*dataToUse["DO"]); #B*NO2*DO
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #northeast atlantic summer
    def _zone4s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1625.17, 12.35, 6.13]; #intersept, SSS, NO3
        rmsd = 15.5; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone4_northeast_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SSS"]) + \
                      coefs[2]*(dataToUse["NO3"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SSS_err"])**2 + \
                                          (coefs[2]*dataToUse["NO3_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #northeast atlantic winter
    def _zone4w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1013.95, 61.48, 32.75, -61.93, -2.05, -12.54, 0.4, -0.3, -1.67, 1.78]; #intersept, SST, SSS, NO3, SiO4, PO4, NO3*SiO4, SST*DO, SST*SSS
        rmsd = 15.5; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone4_northeast_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*SST + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["SiO4"]) + \
                      coefs[5]*(dataToUse["PO4"]) + \
                      coefs[6]*(SST*dataToUse["DO"]) + \
                      coefs[7]*(SST*dataToUse["SSS"]) + \
                      coefs[8]*(dataToUse["NO3"]*dataToUse["SSS"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3]*dataToUse["SiO4_err"]; #B*SiO4
        uterm4 = coefs[4]*dataToUse["PO4_err"]; #B*PO4
        uterm5 = coefs[5] * ((dataToUse["SST_err"]/SST) + (dataToUse["DO_err"]/dataToUse["DO"])) * (SST*dataToUse["DO"]); #B*SST*DO
        uterm6 = coefs[6] * ((dataToUse["SST_err"]/SST) + (dataToUse["SSS_err"]/dataToUse["SSS"])) * (SST*dataToUse["SSS"]); #B*SST*SSS
        uterm7 = coefs[7] * ((dataToUse["NO3_err"]/dataToUse["NO3"]) + (dataToUse["SSS_err"]/dataToUse["SSS"])) * (dataToUse["NO3"]*dataToUse["SSS"]); #B*NO3*SSS
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 + uterm7**2 );
        
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #equatorial pacific
    def _zone5(self, data):
        #For coefficients, see supplemental table T1
        coefs = [467.55, -7.11, 48.4, 2.34, 1.44, 38.85]; #intersept, SST, SSS, NO3, SiO4, PO4
        rmsd = 18.9; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone5_equatorial_pacific");
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["SiO4"]) + \
                      coefs[5]*(dataToUse["PO4"]);
                      
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["NO3_err"])**2 + \
                                          (coefs[4]*dataToUse["SiO4_err"])**2 + \
                                          (coefs[5]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #subtropical north pacific summer
    def _zone6s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [519.77, -9.98, 48.97, 20.25, 1.92]; #intersept, SST, SSS, NO3, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone6_north_subtropical_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["NO3_err"])**2 + \
                                          (coefs[4]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #subtropical north pacific winter
    def _zone6w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [236.82, -5.32, 50.65, 0.38, -1.5, 139.22]; #intersept, SST, SSS, DO, NO3, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone6_north_subtropical_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["NO3"]) + \
                      coefs[5]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["DO_err"])**2 + \
                                          (coefs[4]*dataToUse["NO3_err"])**2 + \
                                          (coefs[5]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #subtropical south pacific summer
    def _zone7s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [147.12, -4.61, 52.63, 0.34, 7.48, 1.67, 72.68]; #intersept, SST, SSS, DO, NO3, SiO4, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone7_south_subtropical_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["NO3"]) + \
                      coefs[5]*(dataToUse["SiO4"]) + \
                      coefs[6]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["DO_err"])**2 + \
                                          (coefs[4]*dataToUse["NO3_err"])**2 + \
                                          (coefs[5]*dataToUse["SiO4_err"])**2 + \
                                          (coefs[6]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #subtropical south pacific winter
    def _zone7w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [643.05, -12.01, 46.68, -1.28, 107.88]; #intersept, SST, SSS, SiO4, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone7_south_subtropical_pacific");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["SiO4"]) + \
                      coefs[4]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["SiO4_err"])**2 + \
                                          (coefs[4]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    
    #indian ocean summer
    def _zone8s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [551.82, -6.59, 45.21, 27.16, 0.14]; #intersept, SST, SSS, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone8_indian_ocean");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["PO4"]) + \
                      coefs[4]*(dataToUse["SSS"]*dataToUse["NO3"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = coefs[3]*dataToUse["PO4_err"]; #B*PO4
        uterm4 = coefs[4] * ((dataToUse["SSS_err"]/dataToUse["SSS"]) + (dataToUse["NO3_err"]/dataToUse["NO3"])) * (dataToUse["SSS"]*dataToUse["NO3"]); #B*SSS*NO3
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #indian ocean winter
    def _zone8w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [1733.55, -1.84, -4.78, 18.13, 2.64, 67.1, 0.17]; #intersept, SST, DO, NO3, SiO4, PO4, DO*SSS
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone8_indian_ocean");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["DO"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["SiO4"]) + \
                      coefs[5]*(dataToUse["PO4"]) + \
                      coefs[6]*(dataToUse["DO"]*dataToUse["SSS"]);
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["DO_err"]; #B*DO
        uterm3 = coefs[3]*dataToUse["NO3_err"]; #B*NO3
        uterm4 = coefs[4]*dataToUse["SiO4_err"]; #B*SiO4
        uterm5 = coefs[5]*dataToUse["PO4_err"]; #B*PO4
        uterm6 = coefs[6] * ((dataToUse["DO_err"]/dataToUse["DO"]) + (dataToUse["SSS_err"]/dataToUse["SSS"])) * (dataToUse["DO"]*dataToUse["SSS"]); #B*DO*SSS
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;

    
    #subtropical north atlantic summer
    def _zone9s(self, data):
        #For coefficients, see supplemental table T1
        coefs = [619.34, -8.18, 46.94, -0.37, -31.07]; #intersept, SST, SSS, DO, NO3
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone9_north_subtropical_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernSummer));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernSummer)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["DO"]) + \
                      coefs[4]*(dataToUse["NO3"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["DO_err"])**2 + \
                                          (coefs[4]*dataToUse["NO3_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #subtropical north atlantic winter
    def _zone9w(self, data):
        #For coefficients, see supplemental table T1
        coefs = [765.0, -7.36, 39.89, -4.88, 109.27]; #intersept, SST, SiO4, PO4
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone9_north_subtropical_atlantic");
        
        #Calculate indices of summer data
        seasonalIndices = (dataToUse["lat"] >= 0) & (dataToUse["date"].dt.month.isin(self.northernWinter));
        seasonalIndices = seasonalIndices | ( (dataToUse["lat"] < 0) & (dataToUse["date"].dt.month.isin(self.southernWinter)) );
        dataToUse = dataToUse[seasonalIndices];
        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["SiO4"]) + \
                      coefs[4]*(dataToUse["PO4"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["SiO4_err"])**2 + \
                                          (coefs[4]*dataToUse["PO4_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #equatorial atlantic (summer and winter)
    def _zone10(self, data):
        #For coefficients, see supplemental table T1
        coefs = [163.5, -8.91, 59.32, -0.17]; #intersept, SST, SSS, DO
        rmsd = 18.9; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone10_equatorial_atlantic");

        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["SSS"]) + \
                      coefs[3]*(dataToUse["DO"]);
        
        #Propagate uncertainty
        outputUncertaintyDueToInputUncertainty = np.sqrt( (coefs[1]*dataToUse["SST_err"])**2 + \
                                          (coefs[2]*dataToUse["SSS_err"])**2 + \
                                          (coefs[3]*dataToUse["DO_err"])**2
                                         );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;
    
    #south subtropical atlantic (summer and winter)
    def _zone11(self, data):
        #For coefficients, see supplemental table T1
        coefs = [2277.89, -6.15, -7.48, -5.08, 74.92, 0.2]; #intersept, SST, DO, NO3, PO4, DO*SSS
        rmsd = 15.2; #RSE from testing dataset (see table 3)
        
        #Subset data to only rows valid for this zone. See fig 3 and supplemental fig F2
        dataToUse = subset_from_mask(data, self.regionMasks, "DIC_zone11_south_subtropical_atlantic");

        
        #See supplemental table T1
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        modelOutput = coefs[0] + \
                      coefs[1]*(dataToUse["SST"]-273.15) + \
                      coefs[2]*(dataToUse["DO"]) + \
                      coefs[3]*(dataToUse["NO3"]) + \
                      coefs[4]*(dataToUse["PO4"]) + \
                      coefs[5]*(dataToUse["DO"]*dataToUse["SSS"]);
        
        
        #Propagate uncertainty
        uterm1 = coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = coefs[2]*dataToUse["DO_err"]; #B*DO
        uterm3 = coefs[3]*dataToUse["NO3_err"]; #B*NO3
        uterm4 = coefs[4]*dataToUse["PO4_err"]; #B*PO4
        uterm5 = coefs[5] * ((dataToUse["DO_err"]/dataToUse["DO"]) + (dataToUse["SSS_err"]/dataToUse["SSS"])) * (dataToUse["DO"]*dataToUse["SSS"]); #B*DO*SSS
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsd;



    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        #Innernal function used to run, check and assign values for each equation/zone
        def run_single_zone(function, data, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds):
            zoneData, zoneUncertainty, zoneRmsd = function(data);
            if np.any(np.isfinite(modelOutput[zoneData.index])==True): #Sanity check for overlaps
                raise RuntimeError("Overlapping zones in Sasse2013_dic. Something has done wrong!");
            modelOutput[zoneData.index] = zoneData;
            outputUncertaintyDueToInputUncertainty[zoneData.index] = zoneUncertainty;
            rmsds[zoneData.index] = zoneRmsd;
        
        #Create empty output array
        modelOutput = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        outputUncertaintyDueToInputUncertainty = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        rmsds = pd.Series([np.nan]*len(dataToUse), index=dataToUse.index);
        
        #Perform calculations for each zone
        run_single_zone(self._zone1s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone1w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone2s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone2w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone3s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone3w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone4s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone4w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone5, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone6s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone6w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone7s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone7w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone8s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone8w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone9s, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone9w, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone10, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        run_single_zone(self._zone11, dataToUse, modelOutput, outputUncertaintyDueToInputUncertainty, rmsds);
        
        
        outputUncertaintyDueToInputUncertainty = outputUncertaintyDueToInputUncertainty;
        return modelOutput, outputUncertaintyDueToInputUncertainty, rmsds;


#Sasse, T.P., McNeil, B.I. and Abramowitz, G., 2013. A novel method for diagnosing seasonal to inter-annual surface ocean carbon dynamics from bottle data using neural networks. Biogeosciences, 10(6), pp.4319-4340.
class Sasse2013_global_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Sasse2013_global_dic: S13g(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS", "SST", "DO", "SiO4", "PO4", "NO3"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;
        #For coefficients see supplemental table T1
        self.coefs = [596.77, -8.21, 45.5, -0.17, 1.12, 0.45, 17.83, 0.01, 1.52]; #intersept, SST, SSS, DO, NO3, SiO4, PO4, SST*DO, SST*PO4
        self.coefsUncertainty = [None, None, None, None, None, None, None, None, None]; #Uncertainty reported for the coefficients, equation 1
        self.rmsd = 15.6; #RSE from testing dataset (see table 3)
        self.r = None;
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-180, 180)]; #Global
        self.includedRegionsLats = [(-90, 70)]; #Not higher than 70N, see grey region in Fig 3
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = { #None reported
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        SST = dataToUse["SST"]-273.15;
        
        #See supplemental table T1
        modelOutput = self.coefs[0] + \
                      self.coefs[1]*SST + \
                      self.coefs[2]*(dataToUse["SSS"]) + \
                      self.coefs[3]*(dataToUse["DO"]) + \
                      self.coefs[4]*(dataToUse["NO3"]) + \
                      self.coefs[5]*(dataToUse["SiO4"]) + \
                      self.coefs[6]*(dataToUse["PO4"]) + \
                      self.coefs[7]*(dataToUse["DO"]*SST) + \
                      self.coefs[8]*(dataToUse["PO4"]*SST);
        
        #Propagate uncertainty
        uterm1 = self.coefs[1]*dataToUse["SST_err"]; #B*SST
        uterm2 = self.coefs[2]*dataToUse["SSS_err"]; #B*SSS
        uterm3 = self.coefs[3]*dataToUse["DO_err"]; #B*DO
        uterm4 = self.coefs[4]*dataToUse["NO3_err"]; #B*NO3
        uterm5 = self.coefs[5]*dataToUse["SiO4_err"]; #B*SiO4
        uterm6 = self.coefs[6]*dataToUse["PO4_err"]; #B*PO4
        uterm7 = self.coefs[7] * ((dataToUse["DO_err"]/dataToUse["DO"]) + (dataToUse["SST_err"]/SST)) * (dataToUse["DO"]*SST); #B*DO*SST
        uterm8 = self.coefs[8] * ((dataToUse["PO4_err"]/dataToUse["PO4"]) + (dataToUse["SST_err"]/SST)) * (dataToUse["PO4"]*SST); #B*PO4*SST
        #Add in quadrature
        outputUncertaintyDueToInputUncertainty = np.sqrt( uterm1**2 + uterm2**2 + uterm3**2 + uterm4**2 + uterm5**2 + uterm6**2 + uterm7**2 + uterm8**2 );
        
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;



#Vangriesheim, A., Pierre, C., Aminot, A., Metzl, N., Baurand, F. and Caprais, J.C., 2009. The influence of Congo River discharges in the surface and deep layers of the Gulf of Guinea. Deep Sea Research Part II: Topical Studies in Oceanography, 56(23), pp.2183-2196.
#Open ocean/non-plume (surface), see second unlabelled equation at the end of section 3.1
class Vangriesheim2009_open_ocean_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Vangriesheim2009_open_ocean_dic: V00(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [221.0, 50.6]; #intercept, salinity slope
        self.coefsUncertainty = [97.0, 3.0]; #Uncertainty reported for the coefficients
        self.rmsd = None; #
        self.r = np.sqrt(0.957); #Ternon2000 Fig 6
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-4.0, 12.0)]; #See fig 1
        self.includedRegionsLats = [(-10, -1.5)]; #See fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (33.0, 36.0), #end of section 3.1
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #see second unlabelled equation at the end of section 3.1
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;


#Vangriesheim, A., Pierre, C., Aminot, A., Metzl, N., Baurand, F. and Caprais, J.C., 2009. The influence of Congo River discharges in the surface and deep layers of the Gulf of Guinea. Deep Sea Research Part II: Topical Studies in Oceanography, 56(23), pp.2183-2196.
#All data (surface), see fig 6
class Vangriesheim2009_all_dic(BaseAlgorithm):    
    #String representation of the algorithm
    def __str__(self):
        return "Vangriesheim2009_all_dic: V00(dic)";

    #common names of input and output variables (see global_settings for definitions of these
    @staticmethod
    def input_names():
        return ["SSS"];
    @staticmethod
    def output_name():
        return "DIC";
    
    #Set algorithm specific variables
    def __init__(self, settings):
        self.settings = settings;

        self.coefs = [355.0, 46.5]; #intersept, salinity slope
        self.coefsUncertainty = [48.0, 1.0]; #Uncertainty reported for the coefficients
        self.rmsd = None; #
        self.r = np.sqrt(0.969); #Ternon2000 Fig 6
        
        #Specify rectangular regions which the algorithm is valid for. Defaults to global when empty.
        self.includedRegionsLons = [(-4.0, 12.0)]; #See fig 1
        self.includedRegionsLats = [(-10, -1.5)]; #See fig 1
        
        #Algorithm will only be applied to values inside these ranges
        self.restrictRanges = {"SSS": (23.0, 36.0), #fig 6
                               };
        
        #If the matchup dataset contains values outside of these ranges they will be flagged to the user
        self.flagRanges = {
                           };

    #The main calculation is performed here, returns the model output
    def _kernal(self, dataToUse):
        modelOutput = self.coefs[0] + self.coefs[1]*dataToUse["SSS"]; #see fig 6
        outputUncertaintyDueToInputUncertainty = self.coefs[1]*dataToUse["SSS_err"];
        return modelOutput, outputUncertaintyDueToInputUncertainty, self.rmsd;









































