#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:46:01 2020

@author: tom holding
"""

from os import path;
from string import Template;

from os_dic_outflow import calc_dic_outflow_from_transects;
from os_dic_outflow import process_dic_outflow_transects;


def main(carbonateParametersTemplate, outputDirectoryRoot, regions, regionMaskPath, gridAreasPath, perimeterRadii=range(1, 25), numSamples=100):
    outputDirectoryTemplate = Template(path.join(outputDirectoryRoot, "${REGION}"));
    calc_dic_outflow_from_transects.calculate_dic_outflow_from_circle_transects(carbonateParametersTemplate, outputDirectoryTemplate, regions, regionMaskPath, gridAreasPath, perimeterRadii, numSamples, verbose=True)
    
    inputDirectoryTemplate = outputDirectoryTemplate; #the input directory is the output directory from the previous step.
    process_dic_outflow_transects.process_dic_outflow_transects(inputDirectoryTemplate, regions, outputDirectoryRoot);
    