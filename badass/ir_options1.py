################################## Fit Options #################################
# Fitting Parameters
fit_options={
"fit_reg"    : (48000,280000), # Fitting region; Note: Indo-US Library=(3460,9464)
"good_thresh": 0.0, # percentage of "good" pixels required in fig_reg for fit.
"mask_bad_pix": False, # mask pixels SDSS flagged as 'bad' (careful!)
"mask_emline" : False, # automatically mask lines for continuum fitting.
"mask_metal": False, # mask metal absorption lines for high-z spectra
"fit_stat": "RCHI2", # fit statistic; ML = Max. Like. , OLS = Ordinary Least Squares, RCHI2 = reduced chi2
"n_basinhop": 10, # Number of consecutive basinhopping thresholds before solution achieved
"test_outflows": False, # only test for outflows; stops after test
"test_line": {"bool":False,
              "line":"na_FeII_5"},
"max_like_niter": 0, # number of maximum likelihood iterations
"output_pars": False, # only output free parameters of fit and stop code (diagnostic)
"cosmology": {"H0":70.0, "Om0": 0.30}, # Flat Lam-CDM Cosmology
}
################################################################################

########################### MCMC algorithm parameters ##########################
mcmc_options={
"mcmc_fit"    : False, # Perform robust fitting using emcee
"nwalkers"    : 100,  # Number of emcee walkers; min = 2 x N_parameters
"auto_stop"   : False, # Automatic stop using autocorrelation analysis
"conv_type"   : "all", # "median", "mean", "all", or (tuple) of parameters
"min_samp"    : 1000,  # min number of iterations for sampling post-convergence
"ncor_times"  : 1.0,  # number of autocorrelation times for convergence
"autocorr_tol": 10.0,  # percent tolerance between checking autocorr. times
"write_iter"  : 100,   # write/check autocorrelation times interval
"write_thresh": 100,   # iteration to start writing/checking parameters
"burn_in"     : 1500, # burn-in if max_iter is reached
"min_iter"    : 2500, # min number of iterations before stopping
"max_iter"    : 2500, # max number of MCMC iterations
}
################################################################################

############################ Fit component op dtions #############################
comp_options={
"fit_opt_feii"     : False, # optical FeII
"fit_uv_iron"      : False, # UV Iron 
"fit_balmer"       : False, # Balmer continuum (<4000 A)
"fit_losvd"        : False, # stellar LOSVD
"fit_host"         : False, # host template
"fit_power"        : True, # AGN power-law
"fit_poly"         : True, # Add polynomial continuum component
"fit_narrow"       : True, # narrow lines
"fit_broad"        : False, # broad lines
"fit_outflow"      : False, # outflow lines
"fit_absorp"       : False, # absorption lines
"tie_line_disp"    : False, # tie line widths
"tie_line_voff"    : False, # tie line velocity offsets
"na_line_profile"  : "gaussian",     # narrow line profile
"br_line_profile"  : "gaussian",     # broad line profile
"out_line_profile" : "gaussian",     # outflow line profile
"abs_line_profile" : "gaussian",     # absorption line profile
"n_moments"        : 4, # number of Gauss-Hermite moments for Gauss-Hermite line profiles
                        # must be >2 and <10 for higher-order moments (default = 4)
}
################################################################################

########################### Emission Lines & Options ###########################
# If not specified, defaults to SDSS-QSO Emission Lines (http://classic.sdss.org/dr6/algorithms/linestable.html)
################################################################################
user_lines = {
    # MIRI
    # CH 1
    'na_FeII_4': {'center':48890., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S8_5': {'center':50530., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeII_5': {'center':53400., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeVIII_5': {'center':54470., 'line_type':'user', 'line_profile':'gaussian'},
    'na_MgVII_5': {'center':55030., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S7_5': {'center':55110., 'line_type':'user', 'line_profile':'gaussian'},
    'na_MgV_5': {'center':56100., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S6_6': {'center':61090., 'line_type':'user', 'line_profile':'gaussian'},
    'na_PAH62_6': {'center':62000., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArIII_6': {'center':63670., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NiII_6': {'center':66360., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeII_6': {'center':67210., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S5_6': {'center':69090., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArII_6': {'center':69850., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NaIII_7': {'center':73180., 'line_type':'user', 'line_profile':'gaussian'},
    # CH 2
    'na_HI63_7': {'center':74600., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NeVI_7': {'center':76520., 'line_type':'user', 'line_profile':'gaussian'},
    'na_PAH77_7': {'center':77000., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeVII_7': {'center':78140., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArV_7': {'center':79020., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S4_8': {'center':80260., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArIII_8': {'center':89910., 'line_type':'user', 'line_profile':'gaussian'},
    'na_MgVII_9': {'center':90090., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeVII_9': {'center':95270., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S3_9': {'center':96650., 'line_type':'user', 'line_profile':'gaussian'},
    'na_SIV_10': {'center':105100., 'line_type':'user', 'line_profile':'gaussian'},
    'na_PAH113_11': {'center':113000., 'line_type':'user', 'line_profile':'gaussian'},
    # CH 3
    'na_ClIV_11': {'center':117630., 'line_type':'user', 'line_profile':'gaussian'},
    'na_SIII_12': {'center':120000., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S2_12': {'center':122800., 'line_type':'user', 'line_profile':'gaussian'},
    'na_HI76_12': {'center':123700., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NeII_12': {'center':128100., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArV_13': {'center':131000., 'line_type':'user', 'line_profile':'gaussian'},
    'na_MgV_13': {'center':135200., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NeV_14': {'center':143200., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ClII_14': {'center':143700., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NeIII_15': {'center':155600., 'line_type':'user', 'line_profile':'gaussian'},
    'na_H200S1_17': {'center':170300., 'line_type':'user', 'line_profile':'gaussian'},
    # CH 4
    'na_FeII_17': {'center':179400., 'line_type':'user', 'line_profile':'gaussian'},
    'na_SIII_18': {'center':187100., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeVI_19': {'center':195527., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ClIV_20': {'center':203200., 'line_type':'user', 'line_profile':'gaussian'},
    'na_ArIII_21': {'center':218250., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeIII_22': {'center':229250., 'line_type':'user', 'line_profile':'gaussian'},
    'na_NeV_24': {'center':243200., 'line_type':'user', 'line_profile':'gaussian'},
    'na_OIV_25': {'center':258900., 'line_type':'user', 'line_profile':'gaussian'},
    'na_FeII_25': {'center':259900., 'line_type':'user', 'line_profile':'gaussian'},
}
user_constraints = [
#     ("NA_OIII_5007_AMP","OUT_OIII_5007_AMP"),
]
# User defined masked regions (list of tuples)
user_mask = [
# (6250,6525),
]

combined_lines = {

}
########################## LOSVD Fitting & Options ##################
# For direct fitting of the stellar kinematics (stellar LOSVD), one can 
# specify a stellar template library (Indo-US, Vazdekis 2010, or eMILES).
# One can also hold velocity or dispersion constant to avoid template
# convolution during the fitting process.
################################################################################

losvd_options = {
"library"   : "IndoUS", # Options: IndoUS, Vazdekis2010, eMILES
"vel_const" : {"bool":False, "val":0.0},
"disp_const": {"bool":False, "val":250.0},
"losvd_apoly": {"bool":True , "order":3},
}


########################## SSP Host Galaxy Template & Options ##################
# The default is zero velocity, 100 km/s dispersion 10 Gyr template from 
# the eMILES stellar library. 
################################################################################

host_options = {
"age"       : [1.0,5.0,10.0], # Gyr; [0.09 Gyr - 14 Gyr] 
"vel_const" : {"bool":False, "val":0.0},
"disp_const": {"bool":False, "val":250.0}
}

########################### AGN power-law continuum & Options ##################
# The default is a simple power law.
################################################################################

power_options = {
"type" : "simple" # alternatively, "broken" for smoothly-broken power-law
}

########################### Polynomial Continuum Options #######################
# Disabled by default.  Options for a power series polynomial continuum, 
# additive legendre polynomial, or multiplicative polynomial to be included in 
# the fit.
################################################################################

poly_options = {
"ppoly" : {"bool": False, "order": 3}, # positive definite additive polynomial 
"apoly" : {"bool": True , "order": 3}, # Legendre additive polynomial 
"mpoly" : {"bool": False, "order": 3}, # Legendre multiplicative polynomial 
}

############################### Optical FeII options ###############################
# Below are options for fitting FeII.  For most objects, you don't need to 
# perform detailed fitting on FeII (only fit for amplitudes) use the 
# Veron-Cetty 2004 template ('VC04') (2-6 free parameters)
# However in NLS1 objects, FeII is much stronger, and sometimes more detailed 
# fitting is necessary, use the Kovacevic 2010 template 
# ('K10'; 7 free parameters).

# The options are:
# template   : VC04 (Veron-Cetty 2004) or K10 (Kovacevic 2010)
# amp_const  : constant amplitude (default False)
# disp_const : constant disp (default True)
# voff_const : constant velocity offset (default True)
# temp_const : constant temp ('K10' only)

opt_feii_options={
"opt_template"  :{"type":"VC04"}, 
"opt_amp_const" :{"bool":False, "br_opt_feii_val":1.0   , "na_opt_feii_val":1.0},
"opt_disp_const":{"bool":False, "br_opt_feii_val":3000.0, "na_opt_feii_val":500.0},
"opt_voff_const":{"bool":False, "br_opt_feii_val":0.0   , "na_opt_feii_val":0.0},
}
# or
# opt_feii_options={
# "opt_template"  :{"type":"K10"},
# "opt_amp_const" :{"bool":False,"f_feii_val":1.0,"s_feii_val":1.0,"g_feii_val":1.0,"z_feii_val":1.0},
# "opt_disp_const":{"bool":False,"opt_feii_val":1500.0},
# "opt_voff_const":{"bool":False,"opt_feii_val":0.0},
# "opt_temp_const":{"bool":True,"opt_feii_val":10000.0},
# }
################################################################################

############################### UV Iron options ################################
uv_iron_options={
"uv_amp_const"  :{"bool":False, "uv_iron_val":1.0},
"uv_disp_const" :{"bool":False, "uv_iron_val":3000.0},
"uv_voff_const" :{"bool":False, "uv_iron_val":0.0},
"uv_legendre_p" :{"bool":False , "uv_iron_val":3},
}
################################################################################

########################### Balmer Continuum options ###########################
# For most purposes, only the ratio R, and the overall amplitude are free paramters
# but if you want to go crazy, you can fit everything.
balmer_options = {
"R_const"          :{"bool":False, "R_val":1.0}, # ratio between balmer continuum and higher-order balmer lines
"balmer_amp_const" :{"bool":False, "balmer_amp_val":1.0}, # amplitude of overall balmer model (continuum + higher-order lines)
"balmer_disp_const":{"bool":False,  "balmer_disp_val":5000.0}, # broadening of higher-order Balmer lines
"balmer_voff_const":{"bool":False,  "balmer_voff_val":0.0}, # velocity offset of higher-order Balmer lines
"Teff_const"       :{"bool":True,  "Teff_val":15000.0}, # effective temperature
"tau_const"        :{"bool":True,  "tau_val":1.0}, # optical depth
}

################################################################################

############################### Plotting options ###############################
plot_options={
"plot_param_hist"    : False,# Plot MCMC histograms and chains for each parameter
"plot_flux_hist"     : False,# Plot MCMC hist. and chains for component fluxes
"plot_lum_hist"      : False,# Plot MCMC hist. and chains for component luminosities
"plot_eqwidth_hist"  : False, # Plot MCMC hist. and chains for equivalent widths 
"plot_HTML"          : True,# make interactive plotly HTML best-fit plot
}
################################################################################

################################ Output options ################################
output_options={
"write_chain" : False, # Write MCMC chains for all paramters, fluxes, and
                         # luminosities to a FITS table We set this to false 
                         # because MCMC_chains.FITS file can become very large, 
                         # especially  if you are running multiple objects.  
                         # You only need this if you want to reconstruct chains 
                         # and histograms. 
"verbose"     : True,  # prints steps of fitting process in Notebook
}
################################################################################
