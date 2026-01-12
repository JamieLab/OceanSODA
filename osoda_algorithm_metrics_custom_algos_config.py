customAlgorithmInfo = [
###########################################################################################
# ETHZ-OHOA TA
###########################################################################################
                       {"name": "ETHZ_OHOA_at", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "AT", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Total Alkalinity',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_at_cali", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "AT", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Total Alkalinity',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_at_cali_open", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "AT", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '$\mu$mol kg$^{-1}$',
                      "variable": 'Total Alkalinity',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [28,48 ],
                      "lon_limits": [-138,-114],
                      "lon_or": False,
                      "coastal": 'open',
                      "exclude_box": None
                      },
                      {"name": "ETHZ_OHOA_at_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "AT", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '$\mu$mol kg$^{-1}$',
                      "variable": 'Total Alkalinity',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [28,48 ],
                      "lon_limits": [-138,-114],
                      "lon_or": False,
                      "coastal": 'coastal',
                      "exclude_box": None
                      },
                       {"name": "ETHZ_OHOA_at_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "AT", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Total Alkalinity',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_at_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "AT", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '$\mu$mol kg$^{-1}$',
                      "variable": 'Total Alkalinity',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [0,30 ],
                      "lon_limits": [[-180,-150],[140,180]],
                      "lon_or": True,
                      "coastal": 'open',
                      "exclude_box": None
                      },
                      {"name": "ETHZ_OHOA_at_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                      "outputVar": "AT", #DIC or AT
                      "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                      "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                      "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                      "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                      "units": '$\mu$mol kg$^{-1}$',
                      "variable": 'Total Alkalinity',
                      "secondary": None,
                      "thrid": None,
                      "lat_limits": [0,30 ],
                      "lon_limits": [[-180,-150],[140,180]],
                      "lon_or": True,
                      "coastal": 'coastal',
                      "exclude_box": None
                      },
                       {"name": "ETHZ_OHOA_at_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "AT", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Total Alkalinity',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_at_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "AT", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Total Alkalinity',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
###########################################################################################
# ETHZ-OHOA DIC
###########################################################################################
                      {"name": "ETHZ_OHOA_dic", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "DIC", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$mol kg$^{-1}$',
                       "variable": 'Dissolved Inorganic Carbon',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_dic_open", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "DIC", #DIC or AT
                        "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$mol kg$^{-1}$',
                        "variable": 'Dissolved Inorganic Carbon',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'open',
                        "exclude_box": None
                        },
                        {"name": "ETHZ_OHOA_dic_coastal", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": None,
                         "lon_limits": None,
                         "lon_or": False,
                         "coastal": 'coastal',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [0,30 ],
                         "lon_limits": [[-180,-150],[140,180]],
                         "lon_or": True,
                         "coastal": 'all',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [0,30 ],
                         "lon_limits": [[-180,-150],[140,180]],
                         "lon_or": True,
                         "coastal": 'open',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [0,30 ],
                         "lon_limits": [[-180,-150],[140,180]],
                         "lon_or": True,
                         "coastal": 'coastal',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_cali", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [28,48 ],
                         "lon_limits": [-138,-114],
                         "lon_or": False,
                         "coastal": 'all',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_cali_open", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [28,48 ],
                         "lon_limits": [-138,-114],
                         "lon_or": False,
                         "coastal": 'open',
                         "exclude_box": None
                         },
                        {"name": "ETHZ_OHOA_dic_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "DIC", #DIC or AT
                         "matchupVariableName": "ethz_ohoa_median_dic", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
                         "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$mol kg$^{-1}$',
                         "variable": 'Dissolved Inorganic Carbon',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [28,48 ],
                         "lon_limits": [-138,-114],
                         "lon_or": False,
                         "coastal": 'coastal',
                         "exclude_box": None
                         },
###########################################################################################
# ETHZ-OHOA pH
###########################################################################################
                      {"name": "ETHZ_OHOA_ph", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_cali", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_cali_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_ph_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
###########################################################################################
# ETHZ-OHOA pCO2sw
###########################################################################################
                      {"name": "ETHZ_OHOA_pco2", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                       {"name": "ETHZ_OHOA_pco2_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'all',
                        "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                        },
                      {"name": "ETHZ_OHOA_pco2_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_open_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                       },
                      {"name": "ETHZ_OHOA_pco2_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_coastal_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                       },
                      {"name": "ETHZ_OHOA_pco2_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_cali", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_cali_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_pco2_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar": "insitu_pco2w_mean", #DIC or AT
                       "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": '$\mu$atm',
                       "variable":  'pCO$_{2 (sw)}$',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
###########################################################################################
# ETHZ-OHOA resampled to monthly pCO2
###########################################################################################
                       {"name": "ETHZ_OHOA_monthly_pco2", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'all',
                        "exclude_box": None
                        },
                       {"name": "ETHZ_OHOA_monthly_pco2_open", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'open',
                        "exclude_box": None
                        },
                       {"name": "ETHZ_OHOA_monthly_pco2_coastal", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'coastal',
                        "exclude_box": None
                        },
                       {"name": "ETHZ_OHOA_monthly_pco2_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": [0,30 ],
                        "lon_limits": [[-180,-150],[140,180]],
                        "lon_or": True,
                        "coastal": 'all',
                        "exclude_box": None
                        },
                       {"name": "ETHZ_OHOA_monthly_pco2_cali", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": [28,48 ],
                        "lon_limits": [-138,-114],
                        "lon_or": False,
                        "coastal": 'all',
                        "exclude_box": None
                        },
###########################################################################################
# ETHZ-OHOA resampled to monthly pH
###########################################################################################
                      {"name": "ETHZ_OHOA_monthly_ph", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_monthly_ph_open", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'open',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_monthly_ph_coastal", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": None,
                       "lon_limits": None,
                       "lon_or": False,
                       "coastal": 'coastal',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_monthly_ph_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [0,30 ],
                       "lon_limits": [[-180,-150],[140,180]],
                       "lon_or": True,
                       "coastal": 'all',
                       "exclude_box": None
                       },
                      {"name": "ETHZ_OHOA_monthly_ph_cali", #Human readable name, this can be set to anything and is only used as a label
                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                       "units": 'Total Scale',
                       "variable": 'pH',
                       "secondary": None,
                       "thrid": None,
                       "lat_limits": [28,48 ],
                       "lon_limits": [-138,-114],
                       "lon_or": False,
                       "coastal": 'all',
                       "exclude_box": None
                       },
###########################################################################################
# CMEMS pH
###########################################################################################
                       {"name": "CMEMS_ph", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #propagated input data uncertainty
                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": 'Total Scale',
                        "variable": 'pH',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'all',
                        "exclude_box": None
                        },
                        {"name": "CMEMS_ph_open", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": 'Total Scale',
                         "variable": 'pH',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": None,
                         "lon_limits": None,
                         "lon_or": False,
                         "coastal": 'open',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_ph_coastal", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": 'Total Scale',
                         "variable": 'pH',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": None,
                         "lon_limits": None,
                         "lon_or": False,
                         "coastal": 'coastal',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_ph_cali", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": 'Total Scale',
                         "variable": 'pH',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [28,48 ],
                         "lon_limits": [-138,-114],
                         "lon_or": False,
                         "coastal": 'all',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_ph_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": 'Total Scale',
                         "variable": 'pH',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [0,30 ],
                         "lon_limits": [[-180,-150],[140,180]],
                         "lon_or": True,
                         "coastal": 'all',
                         "exclude_box": None
                         },
###########################################################################################
# CMEMS pCO2sw
###########################################################################################
                       {"name": "CMEMS_pco2", #Human readable name, this can be set to anything and is only used as a label
                        "outputVar": "insitu_pco2w_mean", #DIC or AT
                        "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                        "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                        "inputUncertaintyName": None, #propagated input data uncertainty
                        "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                        "units": '$\mu$atm',
                        "variable":  'pCO$_{2 (sw)}$',
                        "secondary": None,
                        "thrid": None,
                        "lat_limits": None,
                        "lon_limits": None,
                        "lon_or": False,
                        "coastal": 'all',
                        "exclude_box": None
                        },
                        {"name": "CMEMS_pco2_open", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "insitu_pco2w_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$atm',
                         "variable":  'pCO$_{2 (sw)}$',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": None,
                         "lon_limits": None,
                         "lon_or": False,
                         "coastal": 'open',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_pco2_coastal", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "insitu_pco2w_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$atm',
                         "variable":  'pCO$_{2 (sw)}$',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": None,
                         "lon_limits": None,
                         "lon_or": False,
                         "coastal": 'coastal',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_pco2_cali", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "insitu_pco2w_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$atm',
                         "variable":  'pCO$_{2 (sw)}$',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [28,48 ],
                         "lon_limits": [-138,-114],
                         "lon_or": False,
                         "coastal": 'all',
                         "exclude_box": None
                         },
                        {"name": "CMEMS_pco2_pacific_islands", #Human readable name, this can be set to anything and is only used as a label
                         "outputVar": "insitu_pco2w_mean", #DIC or AT
                         "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                         "inputUncertaintyName": None, #propagated input data uncertainty
                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                         "units": '$\mu$atm',
                         "variable":  'pCO$_{2 (sw)}$',
                         "secondary": None,
                         "thrid": None,
                         "lat_limits": [0,30 ],
                         "lon_limits": [[-180,-150],[140,180]],
                         "lon_or": True,
                         "coastal": 'all',
                         "exclude_box": None
                         },
###########################################################################################
# PML-OHOA TA
###########################################################################################
                            {"name": "PML_OHOA_at", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'all',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_cali", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [28,48 ],
                            "lon_limits": [-138,-114],
                            "lon_or": False,
                            "coastal": 'all',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_cali_open", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [28,48 ],
                            "lon_limits": [-138,-114],
                            "lon_or": False,
                            "coastal": 'open',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [28,48 ],
                            "lon_limits": [-138,-114],
                            "lon_or": False,
                            "coastal": 'coastal',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [0,30 ],
                            "lon_limits": [[-180,-150],[140,180]],
                            "lon_or": True,
                            "coastal": 'all',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [0,30 ],
                            "lon_limits": [[-180,-150],[140,180]],
                            "lon_or": True,
                            "coastal": 'open',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": [0,30 ],
                            "lon_limits": [[-180,-150],[140,180]],
                            "lon_or": True,
                            "coastal": 'coastal',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_coastal", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'coastal',
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_at_open", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "AT", #DIC or AT
                            "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$mol kg$^{-1}$',
                            "variable": 'Total Alkalinity',
                            "secondary": None,
                            "thrid": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'open',
                            "exclude_box": None
                            },
###########################################################################################
# PML-OHOA TA compared to ETHZ-OHOA TA
###########################################################################################
                           {"name": "PML_OHOA_at_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'all',
                           "secondary": "ethz_ohoa_mean_alk",
                           "exclude_box": None
                           },
                           {"name": "PML_OHOA_at_cali_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": [28,48 ],
                           "lon_limits": [-138,-114],
                           "lon_or": False,
                           "coastal": 'all',
                           "secondary": "ethz_ohoa_mean_alk",
                           "exclude_box": None
                           },
                           {"name": "PML_OHOA_at_pacificislands_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": [0,30 ],
                           "lon_limits": [[-180,-150],[140,180]],
                           "lon_or": True,
                           "coastal": 'all',
                           "secondary": "ethz_ohoa_mean_alk",
                           "exclude_box": None
                           },
                           {"name": "PML_OHOA_at_coastal_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'coastal',
                           "secondary": "ethz_ohoa_mean_alk",
                           "exclude_box": None
                           },
                           {"name": "PML_OHOA_at_open_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "PML_OHOA_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'open',
                           "secondary": "ethz_ohoa_mean_alk",
                           "exclude_box": None
                           },
#######################################################################################################
                           {"name": "ETHZ_OHOA_at_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'all',
                           "secondary": "PML_OHOA_alk",
                           "exclude_box": None
                           },
                           {"name": "ETHZ_OHOA_at_cali_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": [28,48 ],
                           "lon_limits": [-138,-114],
                           "lon_or": False,
                           "coastal": 'all',
                           "secondary": "PML_OHOA_alk",
                           "exclude_box": None
                           },
                           {"name": "ETHZ_OHOA_at_pacificislands_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": [0,30 ],
                           "lon_limits": [[-180,-150],[140,180]],
                           "lon_or": True,
                           "coastal": 'all',
                           "secondary": "PML_OHOA_alk",
                           "exclude_box": None
                           },
                           {"name": "ETHZ_OHOA_at_coastal_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'coastal',
                           "secondary": "PML_OHOA_alk",
                           "exclude_box": None
                           },
                           {"name": "ETHZ_OHOA_at_open_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                           "outputVar": "AT", #DIC or AT
                           "matchupVariableName": "ethz_ohoa_mean_alk", #netCDF variable name of the model output (algorithm prediction)
                           "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                           "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                           "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                           "units": '$\mu$mol kg$^{-1}$',
                           "variable": 'Total Alkalinity',
                           "secondary": None,
                           "thrid": None,
                           "lat_limits": None,
                           "lon_limits": None,
                           "lon_or": False,
                           "coastal": 'open',
                           "secondary": "PML_OHOA_alk",
                           "exclude_box": None
                           },
###########################################################################################
# CMEMS comparisions to ETHZ-OHOA with the monthly CMEMS version
###########################################################################################
                           {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS+OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_pco2w_mean", #DIC or AT
                            "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": '$\mu$atm',
                            "variable":  'pCO$_{2 (sw)}$',
                            "secondary": "cmems_biocarbon_mean_spco2",
                            "thrid":"ethz_ohoa_mean_spco2",
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'all',
                            "exclude_box": None
                            },
                            {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS+OHOA_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "insitu_pco2w_mean", #DIC or AT
                             "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$atm',
                             "variable":  'pCO$_{2 (sw)}$',
                             "secondary": "cmems_biocarbon_mean_spco2",
                             "thrid":"ethz_ohoa_mean_spco2",
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                             },
                            {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS+monthly_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "insitu_pco2w_mean", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$atm',
                             "variable":  'pCO$_{2 (sw)}$',
                             "secondary": "cmems_biocarbon_mean_spco2",
                             "thrid":"ETHZ_OHOA_month_spco2",
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "exclude_box": None
                             },
                            {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS+monthly_OHOA_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "insitu_pco2w_mean", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$atm',
                             "variable":  'pCO$_{2 (sw)}$',
                             "secondary": "cmems_biocarbon_mean_spco2",
                             "thrid":"ETHZ_OHOA_month_spco2",
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                             },
                            {"name": "CMEMS_pco2_consistent_to_OHOA_monthly+OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "insitu_pco2w_mean", #DIC or AT
                             "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #propagated input data uncertainty
                             "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$atm',
                             "variable":  'pCO$_{2 (sw)}$',
                             "secondary": "ETHZ_OHOA_month_spco2",
                             "thrid": "ethz_ohoa_mean_spco2",
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "exclude_box": None
                             },
                             {"name": "CMEMS_pco2_consistent_to_OHOA_monthly+OHOA_no_baltic", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "insitu_pco2w_mean", #DIC or AT
                              "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #propagated input data uncertainty
                              "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$atm',
                              "variable":  'pCO$_{2 (sw)}$',
                              "secondary": "ETHZ_OHOA_month_spco2",
                              "thrid": "ethz_ohoa_mean_spco2",
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": [ [10,30,53,62], [16,30,54,67] ]
                              },
                             {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS+OHOA", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total scale',
                              "variable":  'pH',
                              "secondary": "cmems_biocarbon_mean_ph",
                              "thrid":"ethz_ohoa_mean_ph_total",
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                              {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS+monthly_OHOA", #Human readable name, this can be set to anything and is only used as a label
                               "outputVar": "insitu_ph_mean", #DIC or AT
                               "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                               "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                               "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
                               "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                               "units": 'Total scale',
                               "variable":  'pH',
                               "secondary": "cmems_biocarbon_mean_ph",
                               "thrid":"ETHZ_OHOA_month_ph_total",
                               "lat_limits": None,
                               "lon_limits": None,
                               "lon_or": False,
                               "coastal": 'all',
                               "exclude_box": None
                               },
                              {"name": "CMEMS_ph_consistent_to_OHOA_monthly+OHOA", #Human readable name, this can be set to anything and is only used as a label
                               "outputVar": "insitu_ph_mean", #DIC or AT
                               "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
                               "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
                               "inputUncertaintyName": None, #propagated input data uncertainty
                               "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                               "units": 'Total scale',
                               "variable":  'pH',
                               "secondary": "ETHZ_OHOA_month_ph_total",
                               "thrid": "ethz_ohoa_mean_ph_total",
                               "lat_limits": None,
                               "lon_limits": None,
                               "lon_or": False,
                               "coastal": 'all',
                               "exclude_box": None
                               },
###########################################################################################
# PML-OHOA DIC
###########################################################################################
                              {"name": "PML_OHOA_ct", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_cali", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_cali_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'open',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'open',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
                              {"name": "PML_OHOA_ct_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar": "DIC", #DIC or AT
                              "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": '$\mu$mol kg$^{-1}$',
                              "variable": 'Dissolved Inorganic Carbon',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'open',
                              "exclude_box": None
                              },
###########################################################################################
# PML-OHOA pH
###########################################################################################
                             {"name": "PML_OHOA_ph", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName":"PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'open',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": None,
                              "lon_limits": None,
                              "lon_or": False,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_pacificislands", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_pacificislands_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'open',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_pacificislands_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [0,30 ],
                              "lon_limits": [[-180,-150],[140,180]],
                              "lon_or": True,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_cali", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'all',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_cali_open", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'open',
                              "exclude_box": None
                              },
                             {"name": "PML_OHOA_ph_cali_coastal", #Human readable name, this can be set to anything and is only used as a label
                              "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
                              "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                              "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
                              "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
                              "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                              "units": 'Total Scale',
                              "variable": 'pH',
                              "secondary": None,
                              "thrid":None,
                              "lat_limits": [28,48 ],
                              "lon_limits": [-138,-114],
                              "lon_or": False,
                              "coastal": 'coastal',
                              "exclude_box": None
                              },
###########################################################################################
# PML-OHOA DIC consistent to ETHZ-OHOA
###########################################################################################
                             {"name": "PML_OHOA_ct_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "secondary": "ethz_ohoa_mean_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "PML_OHOA_ct_cali_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": [28,48 ],
                             "lon_limits": [-138,-114],
                             "lon_or": False,
                             "coastal": 'all',
                             "secondary": "ethz_ohoa_mean_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "PML_OHOA_ct_pacificislands_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": [0,30 ],
                             "lon_limits": [[-180,-150],[140,180]],
                             "lon_or": True,
                             "coastal": 'all',
                             "secondary": "ethz_ohoa_mean_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "PML_OHOA_ct_coastal_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'coastal',
                             "secondary": "ethz_ohoa_mean_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "PML_OHOA_ct_open_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "PML_OHOA_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'open',
                             "secondary": "ethz_ohoa_mean_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "ETHZ_OHOA_ct_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'all',
                             "secondary": "PML_OHOA_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "ETHZ_OHOA_dic_cali_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": [28,48 ],
                             "lon_limits": [-138,-114],
                             "lon_or": False,
                             "coastal": 'all',
                             "secondary": "PML_OHOA_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "ETHZ_OHOA_dic_pacificislands_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": [0,30 ],
                             "lon_limits": [[-180,-150],[140,180]],
                             "lon_or": True,
                             "coastal": 'all',
                             "secondary": "PML_OHOA_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "ETHZ_OHOA_dic_coastal_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'coastal',
                             "secondary": "PML_OHOA_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
                             {"name": "ETHZ_OHOA_ct_open_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                             "outputVar": "DIC", #DIC or AT
                             "matchupVariableName": "ethz_ohoa_mean_dic", #netCDF variable name of the model output (algorithm prediction)
                             "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                             "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                             "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                             "units": '$\mu$mol kg$^{-1}$',
                             "variable": 'Dissolved Inorganic Carbon',
                             "secondary": None,
                             "lat_limits": None,
                             "lon_limits": None,
                             "lon_or": False,
                             "coastal": 'open',
                             "secondary": "PML_OHOA_dic",
                             "thrid":None,
                             "exclude_box": None
                             },
###########################################################################################
# PML-OHOA pH consistent to ETHZ-OHOA
###########################################################################################
                            {"name": "PML_OHOA_ph_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'all',
                            "secondary": "ethz_ohoa_mean_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_ph_cali_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": [28,48 ],
                            "lon_limits": [-138,-114],
                            "lon_or": False,
                            "coastal": 'all',
                            "secondary": "ethz_ohoa_mean_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_ph_pacificislands_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": [0,30 ],
                            "lon_limits": [[-180,-150],[140,180]],
                            "lon_or": True,
                            "coastal": 'all',
                            "secondary": "ethz_ohoa_mean_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_ph_coastal_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'coastal',
                            "secondary": "ethz_ohoa_mean_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "PML_OHOA_ph_open_consistent_to_ETHZ_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "PML_OHOA_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'open',
                            "secondary": "ethz_ohoa_mean_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            #######################################################################################################
                            {"name": "ETHZ_OHOA_ph_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'all',
                            "secondary": "PML_OHOA_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "ETHZ_OHOA_ph_cali_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": [28,48 ],
                            "lon_limits": [-138,-114],
                            "lon_or": False,
                            "coastal": 'all',
                            "secondary": "PML_OHOA_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "ETHZ_OHOA_ph_pacificislands_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": [0,30 ],
                            "lon_limits": [[-180,-150],[140,180]],
                            "lon_or": True,
                            "coastal": 'all',
                            "secondary": "PML_OHOA_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "ETHZ_OHOA_ph_coastal_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'coastal',
                            "secondary": "PML_OHOA_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                            {"name": "ETHZ_OHOA_ph_open_consistent_to_PML_OHOA", #Human readable name, this can be set to anything and is only used as a label
                            "outputVar": "insitu_ph_mean", #DIC or AT
                            "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
                            "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
                            "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
                            "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
                            "units": 'Total Scale',
                            "variable": 'pH',
                            "secondary": None,
                            "lat_limits": None,
                            "lon_limits": None,
                            "lon_or": False,
                            "coastal": 'open',
                            "secondary": "PML_OHOA_ph_total",
                            "thrid":None,
                            "exclude_box": None
                            },
                    ];






#
#
# #######################################################################
#                        {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ethz_ohoa_mean_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'all'
#                        },
# #######################################################################
#                        {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_monthly", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ETHZ_OHOA_month_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'all'
#                        },
# ########################################
#                        {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS_coastal", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'coastal'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_coastal", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ethz_ohoa_mean_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'coastal'
#                        },
# ########################################
#                        {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS_coastal", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'coastal'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_monthly_coastal", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ETHZ_OHOA_month_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'coastal'
#                        },
# ########################################
#                        {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS_open", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'open'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_open", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ethz_ohoa_mean_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'open'
#                        },
# ########################################
#                        {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS_open", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": None,
#                       "lon_limits": None,
#                       "lon_or": False,
#                       "coastal": 'open'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_monthly_open", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ETHZ_OHOA_month_ph_total",
#                        "thrid": None,
#                        "lat_limits": None,
#                        "lon_limits": None,
#                        "lon_or": False,
#                        "coastal": 'open'
#                        },
# ###########################################
#                        {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS_cali", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": [28,48 ],
#                       "lon_limits": [-138,-114],
#                       "lon_or": False,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_cali", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ethz_ohoa_mean_ph_total",
#                        "thrid": None,
#                        "lat_limits": [28,48 ],
#                        "lon_limits": [-138,-114],
#                        "lon_or": False,
#                        "coastal": 'all'
#                        },
# ###########################################
#                        {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS_cali", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": [28,48 ],
#                       "lon_limits": [-138,-114],
#                       "lon_or": False,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_monthly_cali", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ETHZ_OHOA_month_ph_total",
#                        "thrid": None,
#                        "lat_limits": [28,48 ],
#                        "lon_limits": [-138,-114],
#                        "lon_or": False,
#                        "coastal": 'all'
#                        },
# ##############################################
#                        {"name": "ETHZ_OHOA_ph_consistent_to_CMEMS_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ethz_ohoa_mean_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": [0,30 ],
#                       "lon_limits": [[-180,-150],[140,180]],
#                       "lon_or": True,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ethz_ohoa_mean_ph_total",
#                        "thrid": None,
#                        "lat_limits": [0,30 ],
#                        "lon_limits": [[-180,-150],[140,180]],
#                        "lon_or": True,
#                        "coastal": 'all'
#                        },
# ##############################################
#                        {"name": "ETHZ_OHOA_monthly_ph_consistent_to_CMEMS_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                       "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                       "matchupVariableName": "ETHZ_OHOA_month_ph_total", #netCDF variable name of the model output (algorithm prediction)
#                       "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                       "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       "units": 'Total Scale',
#                       "variable": 'pH',
#                       "secondary": 'cmems_biocarbon_mean_ph',
#                       "thrid": None,
#                       "lat_limits": [0,30 ],
#                       "lon_limits": [[-180,-150],[140,180]],
#                       "lon_or": True,
#                       "coastal": 'all'
#                       },
#                        {"name": "CMEMS_ph_consistent_to_OHOA_monthly_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "cmems_biocarbon_mean_ph", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #propagated input data uncertainty
#                        "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        "units": 'Total Scale',
#                        "variable": 'pH',
#                        "secondary":"ETHZ_OHOA_month_ph_total",
#                        "thrid": None,
#                        "lat_limits": [0,30 ],
#                        "lon_limits": [[-180,-150],[140,180]],
#                        "lon_or": True,
#                        "coastal": 'all'
#                        },
# ##########################################
#                        {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar": "insitu_pco2w_mean", #DIC or AT
#                         "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                         "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         "units": '$\mu$atm',
#                         "variable":  'pCO$_{2 (sw)}$',
#                         "secondary": "cmems_biocarbon_mean_spco2",
#                         "thrid": None,
#                         "lat_limits": None,
#                         "lon_limits": None,
#                         "lon_or": False,
#                         "coastal": 'all'
#                         },
#                         {"name": "CMEMS_pco2_consistent_to_OHOA", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #propagated input data uncertainty
#                          "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "ethz_ohoa_mean_spco2",
#                          "thrid": None,
#                          "lat_limits": None,
#                          "lon_limits": None,
#                          "lon_or": False,
#                          "coastal": 'all'
#                          },
# ##########################################
#                        {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar": "insitu_pco2w_mean", #DIC or AT
#                         "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                         "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         "units": '$\mu$atm',
#                         "variable":  'pCO$_{2 (sw)}$',
#                         "secondary": "cmems_biocarbon_mean_spco2",
#                         "thrid": None,
#                         "lat_limits": None,
#                         "lon_limits": None,
#                         "lon_or": False,
#                         "coastal": 'all'
#                         },
#                         {"name": "CMEMS_pco2_consistent_to_OHOA_monthly", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #propagated input data uncertainty
#                          "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "ETHZ_OHOA_month_spco2",
#                          "thrid": None,
#                          "lat_limits": None,
#                          "lon_limits": None,
#                          "lon_or": False,
#                          "coastal": 'all'
#                          },
# ######################################
#                         {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS_coastal", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                          "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "cmems_biocarbon_mean_spco2",
#                          "thrid": None,
#                          "lat_limits": None,
#                          "lon_limits": None,
#                          "lon_or": False,
#                          "coastal": 'coastal'
#                          },
#                          {"name": "CMEMS_pco2_consistent_to_OHOA_coastal", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #propagated input data uncertainty
#                           "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "ethz_ohoa_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": None,
#                           "lon_limits": None,
#                           "lon_or": False,
#                           "coastal": 'coastal'
#                           },
# ######################################
#                         {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS_coastal", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                          "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "cmems_biocarbon_mean_spco2",
#                          "thrid": None,
#                          "lat_limits": None,
#                          "lon_limits": None,
#                          "lon_or": False,
#                          "coastal": 'coastal'
#                          },
#                          {"name": "CMEMS_pco2_consistent_to_OHOA_monthly_coastal", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #propagated input data uncertainty
#                           "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "ETHZ_OHOA_month_spco2",
#                           "thrid": None,
#                           "lat_limits": None,
#                           "lon_limits": None,
#                           "lon_or": False,
#                           "coastal": 'coastal'
#                           },
# #############################
#                          {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS_open", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                           "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "cmems_biocarbon_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": None,
#                           "lon_limits": None,
#                           "lon_or": False,
#                           "coastal": 'open'
#                           },
#                           {"name": "CMEMS_pco2_consistent_to_OHOA_open", #Human readable name, this can be set to anything and is only used as a label
#                            "outputVar": "insitu_pco2w_mean", #DIC or AT
#                            "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                            "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                            "inputUncertaintyName": None, #propagated input data uncertainty
#                            "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                            "units": '$\mu$atm',
#                            "variable":  'pCO$_{2 (sw)}$',
#                            "secondary": "ethz_ohoa_mean_spco2",
#                            "thrid": None,
#                            "lat_limits": None,
#                            "lon_limits": None,
#                            "lon_or": False,
#                            "coastal": 'open'
#                            },
# #############################
#                          {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS_open", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                           "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "cmems_biocarbon_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": None,
#                           "lon_limits": None,
#                           "lon_or": False,
#                           "coastal": 'open'
#                           },
#                           {"name": "CMEMS_pco2_consistent_to_OHOA_monthly_open", #Human readable name, this can be set to anything and is only used as a label
#                            "outputVar": "insitu_pco2w_mean", #DIC or AT
#                            "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                            "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                            "inputUncertaintyName": None, #propagated input data uncertainty
#                            "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                            "units": '$\mu$atm',
#                            "variable":  'pCO$_{2 (sw)}$',
#                            "secondary": "ETHZ_OHOA_month_spco2",
#                            "thrid": None,
#                            "lat_limits": None,
#                            "lon_limits": None,
#                            "lon_or": False,
#                            "coastal": 'open'
#                            },
# #############################
#                         {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS_cali", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                          "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "cmems_biocarbon_mean_spco2",
#                          "thrid": None,
#                          "lat_limits": [28,48 ],
#                          "lon_limits": [-138,-114],
#                          "lon_or": False,
#                          "coastal": 'all'
#                          },
#                          {"name": "CMEMS_pco2_consistent_to_OHOA_cali", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #propagated input data uncertainty
#                           "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "ethz_ohoa_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": [28,48 ],
#                           "lon_limits": [-138,-114],
#                           "lon_or": False,
#                           "coastal": 'all'
#                           },
# #############################
#                         {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS_cali", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar": "insitu_pco2w_mean", #DIC or AT
#                          "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                          "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          "units": '$\mu$atm',
#                          "variable":  'pCO$_{2 (sw)}$',
#                          "secondary": "cmems_biocarbon_mean_spco2",
#                          "thrid": None,
#                          "lat_limits": [28,48 ],
#                          "lon_limits": [-138,-114],
#                          "lon_or": False,
#                          "coastal": 'all'
#                          },
#                          {"name": "CMEMS_pco2_consistent_to_OHOA_monthly_cali", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #propagated input data uncertainty
#                           "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "ETHZ_OHOA_month_spco2",
#                           "thrid": None,
#                           "lat_limits": [28,48 ],
#                           "lon_limits": [-138,-114],
#                           "lon_or": False,
#                           "coastal": 'all'
#                           },
# #################################################
#                          {"name": "ETHZ_OHOA_pco2_consistent_to_CMEMS_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "ethz_ohoa_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                           "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "cmems_biocarbon_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": [0,30 ],
#                           "lon_limits": [[-180,-150],[140,180]],
#                           "lon_or": True,
#                           "coastal": 'all'
#                           },
#                           {"name": "CMEMS_pco2_consistent_to_OHOA_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                            "outputVar": "insitu_pco2w_mean", #DIC or AT
#                            "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                            "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                            "inputUncertaintyName": None, #propagated input data uncertainty
#                            "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                            "units": '$\mu$atm',
#                            "variable":  'pCO$_{2 (sw)}$',
#                            "secondary": "ethz_ohoa_mean_spco2",
#                            "thrid": None,
#                            "lat_limits": [0,30 ],
#                            "lon_limits": [[-180,-150],[140,180]],
#                            "lon_or": True,
#                            "coastal": 'all'
#                            },
# #################################################
#                          {"name": "ETHZ_OHOA_monthly_pco2_consistent_to_CMEMS_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                           "outputVar": "insitu_pco2w_mean", #DIC or AT
#                           "matchupVariableName": "ETHZ_OHOA_month_spco2", #netCDF variable name of the model output (algorithm prediction)
#                           "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                           "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                           "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                           "units": '$\mu$atm',
#                           "variable":  'pCO$_{2 (sw)}$',
#                           "secondary": "cmems_biocarbon_mean_spco2",
#                           "thrid": None,
#                           "lat_limits": [0,30 ],
#                           "lon_limits": [[-180,-150],[140,180]],
#                           "lon_or": True,
#                           "coastal": 'all'
#                           },
#                           {"name": "CMEMS_pco2_consistent_to_OHOA_monthly_pacificisland", #Human readable name, this can be set to anything and is only used as a label
#                            "outputVar": "insitu_pco2w_mean", #DIC or AT
#                            "matchupVariableName": "cmems_biocarbon_mean_spco2", #netCDF variable name of the model output (algorithm prediction)
#                            "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                            "inputUncertaintyName": None, #propagated input data uncertainty
#                            "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                            "units": '$\mu$atm',
#                            "variable":  'pCO$_{2 (sw)}$',
#                            "secondary": "ETHZ_OHOA_month_spco2",
#                            "thrid": None,
#                            "lat_limits": [0,30 ],
#                            "lon_limits": [[-180,-150],[140,180]],
#                            "lon_or": True,
#                            "coastal": 'all'
#                            },
# #############################################
