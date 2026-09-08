customAlgorithmInfo = [
###########################################################################################
# ETHZ-OHOA TA
###########################################################################################
                       {"name": "ETHZ_OHOA_aragonite", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_aragonite_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '',
                       "variable": 'Aragonite Saturation',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_aragonite_cali", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_aragonite_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '',
                       "variable": 'Aragonite Saturation',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_aragonite_cali_open", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "insitu_aragonite_mean", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '',
                      "variable": 'Aragonite Saturation',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [28,48 ],
                      "lon_limits": [-138,-114],
                      "lon_or": False,
                      "coastal": 'open',
                      "exclude_box": None
                      },
                      {"name": "ETHZ_OHOA_aragonite_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "insitu_aragonite_mean", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '',
                      "variable": 'Aragonite Saturation',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [28,48 ],
                      "lon_limits": [-138,-114],
                      "lon_or": False,
                      "coastal": 'coastal',
                      "exclude_box": None
                      },
                       {"name": "ETHZ_OHOA_aragonite_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_aragonite_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '',
                       "variable": 'Aragonite Saturation',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_aragonite_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "insitu_aragonite_mean", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '',
                      "variable": 'Aragonite Saturation',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [0,30 ],
                      "lon_limits": [[-180,-150],[140,180]],
                      "lon_or": True,
                      "coastal": 'open',
                      "exclude_box": None
                      },
                      {"name": "ETHZ_OHOA_aragonite_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "insitu_aragonite_mean", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '',
                      "variable": 'Aragonite Saturation',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [0,30 ],
                      "lon_limits": [[-180,-150],[140,180]],
                      "lon_or": True,
                      "coastal": 'coastal',
                      "exclude_box": None
                      },
                       {"name": "ETHZ_OHOA_aragonite_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_aragonite_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '',
                       "variable": 'Aragonite Saturation',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_aragonite_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_aragonite_mean",  #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_aragonite_sat", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '',
                       "variable": 'Aragonite Saturation',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       }
]
