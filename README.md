
# Optimal satellite remote sensing of the carbonate system using empirical methods

Methods for evaluating and using empirical approaches for studying the surface marine carbonate system. All work builds on the methods and approaches developed within Land et al., (2019).

# Overview of pipeline and scripts
1) A set of algorithms for predicting DIC or AT are compared, using the 'matchup database' as a test/validation data set (i.e. using the algorithm to predict DIC/AT and comparing to in situ derrived measures in the matchup database).

2) Optimal algorithms are selected, and these are used to calculated gridded time series for AT, DIC and other carbonate system parameters. In this step, the 'gridded prediction data' are used as input to the algorithms. The resulting output is gridded netCDF (.nc) files containing carbonate chemistry fields for each region. A copy of input data sets, and uncertainty data are also included in these files.

3) The gridded time series are used to calculate DIC outflow for Amazon river.

4) The gridded time series are used to assess coral reef vulnerability in each region.

To run the complete analysis run `osoda_driver.py`. This file essentially calls the other scripts in the root directory, each of which run each one of the steps enumerated above:
* `osoda_algorithm_comparison.py` - Uses matchup database to compute algorithm performance metrics for each algorithm and input data combination.
* `osoda_calculate_gridded_predictions.py` - Uses an 'optimal' algorithm table produced by the algorithm comparison, with Earth observation data sets to calculate gridded time series for DIC, AT and other carbonate system parameters.
* `osoda_dic_outflow.py` - Uses guaging station discharge data for the Amazon river and gridded DIC time series prediction data, to estimate DIC outflow from the Amazon.
* `osoda_reef_vulnerability.py` - Uses ReefBase data and gridded DIC time series predictions to estimate reef vulnerability.
* `osoda_global_settings.py` - This contains the global settings used to configure each aspect of the analysis (e.g. data file paths and activating/deactivating different aspects of the analysis).

All output from each step is written to a separate subdirectory in the `output` directory (although this is configurable in the global settings file). Parameters, variable names and filepaths (e.g. for matchup data base, prediction data, region masks etc.) and global options are all defined in osoda_global_settings.py. The driver script makes use of these values, and some parts of the analysis will import the global settings using osoda_global_settings.get_default_settings(). The resulting object is a Python dictionary, so values can be easily overwritten in a custom driver script if required.

# Calculating the algorithm comparison metrics
osoda_algorithm_comparison.py performs the pairwise algorithm comparison using the matchup database and generates a set of weighted and unweighted metrics for each region/input data combination. It uses the methodology of Land et al (2019).

Which algorithms are compared for each region, region boundaries, output directory and other settings are all controlled by the global settings file (osoda_global_settings.py::get_default_settings). 

There is no commandline interface to run this tool and it must be ran using Python. To run, simply import osoda_algorithm_comparison.py and run the main function, passing the global settings object (or your own global settings object) as an argument:

```
import osoda_global_settings
import osoda_algorithm_comparison
settings = osoda_global_settings.get_default_settings()
osoda_algorithm_comparison.main(settings)
```

Or just run the file from the command line to run it with the default global settings
`python osoda_algorithm_comparison.py`

A detailed description of the output files and directory structure is provided in the comments at the top of osoda_algorithm_comparison.py

# Running the gridded time series tool separately
osoda_calculate_gridded_predictions.py can be run as a command line tool. Use -h to see the help information:
`python osoda_calculate_gridded_predictions.py -h`

Required inputs are the path to the selected algorithm table (e.g. generated from osoda_algorithm_comparison) and output path template. The path to the gridded input prediction data can be supplied (as in the example below) but will default to that defined by the global settings file. Start and end years can be supplied also, but will default to those defined in the global settings file if left blank. A list of regions and path to the region mask are optional, if not supplied it will use all the regions and region mask defined in the global settings file. Example usage:
`python osoda_calculate_gridded_predictions.py "output/algo_metrics/overall_best_algos_min_years=8.csv" "output/gridded_predictions_min_year_range/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc" --input_data_root "output/gridded_output_new"`

Note that the output path takes the format of a string.Template definition from the Python standard library. This means that, for example, `${REGION}` is a placeholder which is replaces by the region name. Similarly `${LATRES}`, `${LONRES}` and `${OUTPUTVAR}` are replaced by the resolution and output variable (DIC or AT).

When running the gridded time series calculation for the first time, the tool will try to automatically download and resample the gridded input data, unless it detects it at the path provided (--input_data_root). It takes a long time to download everything (potentially over a day with a fast internet connection), so plan accordingly. When running again it should automatically detect that you have the input data. You may want to delete the raw `downloaded_data` folder after running the tool for the first time to save space, as you only need the resampled files. The specification for the file format, variable names and directory hierarchy are all in the global settings file under the `datasetInfoMap` key.

`datasetInfoMap` in the osoda_global_settings.py defines the mapping between an ocean parameter (e.g. sea surface tempoerature, SST) and the matchup database and prediction data set. To do this it defined a `commonName` (the key used to refer to the ocean parameter), a `datasetName` (a unique name for the specific data set - since multiple data sets for the same ocean parameter can be used). File paths to the netCDF file/s, netCDF variable name, and netCDF variable name for the uncertainty field are also defined here for the matchup database and gridded prediction data used by osoda_calculate_gridded_predictions.py to create the gridded time sereis.

# Running the Amazon DIC outflow calculation
Simply import osoda_dic_outflow.py in python and run main. The main function is commented with descriptions of each argument, all of which can be supplied using values from the global settings file. For example:

```
import osoda_global_settings
import osoda_dic_outflow
settings = osoda_global_settings.get_default_settings()
osodaMasksPath = settings["regionMasksPath"];
precomputedGridAreaPath = settings["gridAreasPath"];
regions = ["oceansoda_amazon_plume"];
carbonateParametersTemplate = settings["longGriddedTimeSeriesPathTemplate"];
outputDir = path.join(settings["outputPathRoot"], "dic_outflow_amazon_best");
osoda_dic_outflow.main(carbonateParametersTemplate, outputDir, regions, osodaMasksPath, precomputedGridAreaPath);
```

# Running the reef vulnerability analysis
Either import osoda_reef_vulnerability and run main with the defualt (or custom) global settings dictionary, or run the python file from command line with no arguments to use the default settings. E.g.
```
import osoda_global_settings
import osoda_reef_vulnerability
settings = osoda_global_settings.get_default_settings()
osoda_reef_vulnerability(settings)
```
is equivailent to running `python osoda_reef_vulnerability.py` in your terminal / command prompt.


# Calculating metrics for custom algorithms
It is possible to add new algorithms which do not conform to the standard 'implemented' algorithm format used in the `os_algorithms` module. The easiest way to do this is to pre-compute model output, RMSD, propagated input data uncertainty and combined (input and model) uncertainty, and add these to the matchup database. Algorithms added to the analysis in this way are referred to as 'custom' algorithms, while algorithms inside the `os_algorithms` module (i.e. with simple Python implementations) are described as 'implemented' algorithms within the code.

The `osoda_algorithm_comparison.py` file contains a function ('custom_algorithm_metrics') for calculating the metrics using a combination of custom and implemented algorithms. The `example_metrics_from_algo_output.py` script (in the root directory) provides a full example of how to use this function and is fully commented. This script first generates four synthesised test data sets and adds them to the matchup database. For each custom algorithm you must provide the model output (calculated from matchup database inputs), model RMSD, propagated input uncertainty and combined (model and input) uncertainty. The script then defines which matchup database variables correspond to each of these fields for each custom algorithm, and runs the `custom_algorith_metrics` function to perform the full algorithm comparison.

`custom_algorith_metrics` requires a list of algorithm info about each custom algorithm (locations of the fields described above in the matchup database), a list of any 'implemented' algorithms you want to also include in the analysis, a region name (which must correspond to one of the regions in the mask file), an SST and SSS input dataset name (which must correspond to one of the data sets defined in the global settings 'datasetInfoMap' dictionary), and an output path. Optionally, the script can generate simple diagnostic plots which plot in situ derrived matchup database output (DIC/AT) against the model output as a quick visual aid. Note that, since a given set of custom algorithm output data must have been generated using a single input data combination, the function can only be ran using one input combination at a time (hence, the SST and SSS data set names are required inputs to the function).


### References

Land PE, Findlay H, Shutler J, Ashton I, Holding T, Grouazel A, GIrard-Ardhuin F, Reul N, Piolle J-F, Chapron B, et al (2019). Optimum satellite remote sensing of the marine carbonate system using empirical algorithms in the Global Ocean, the Greater Caribbean, the Amazon Plume and the Bay of Bengal. Remote Sensing of Environment, doi: 10.1016/j.rse.2019.111469
