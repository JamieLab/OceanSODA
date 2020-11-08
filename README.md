
# Optimal satellite remote sensing of the carbonate system using empirical methods

Methods for evaluating and using empirical approaches for studying the surface marine carbonate system. All work builds on the methods and approaches developed within Land et al., (2019).

# Scripts and directory structure
To run the complete analysis run `osoda_driver.py`. The other scripts in the root directory run each component of the analysis:
* `osoda_global_settings.py` - This contains the global settings used to configure each aspect of the analysis (e.g. data file paths and activating/deactivating different aspects of the analysis).
* `osoda_algorithm_comparison.py` - Uses matchup database to compute algorithm performance metrics for each algorithm and input data combination.
* `osoda_calculate_gridded_predictions.py` - Uses an 'optimal' algorithm table produced by the algorithm comparison, with Earth observation data sets to calculate gridded time series for DIC, AT and other carbonate system parameters.
* `osoda_dic_outflow.py` - Uses guaging station discharge data for the Amazon river and gridded DIC time series prediction data, to estimate DIC outflow from the Amazon.
* `osoda_reef_vulnerability.py` - Uses ReefBase data and gridded DIC time series predictions to estimate reef vulnerability.

All output from each step is written to a separate subdirectory in the `output` directory.


### References

Land PE, Findlay H, Shutler J, Ashton I, Holding T, Grouazel A, GIrard-Ardhuin F, Reul N, Piolle J-F, Chapron B, et al (2019). Optimum satellite remote sensing of the marine carbonate system using empirical algorithms in the Global Ocean, the Greater Caribbean, the Amazon Plume and the Bay of Bengal. Remote Sensing of Environment, doi: 10.1016/j.rse.2019.111469
