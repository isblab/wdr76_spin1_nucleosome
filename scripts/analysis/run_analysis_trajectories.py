import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change this to the location of your PMI_analysis folder
sys.path.append("/home/shreyas/imp-clean/pmi_analysis/pyext/src")

from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 120
top_dir = sys.argv[1]
dir_head = sys.argv[2]  # How are the trajectories dir names
out_dirs = glob.glob(top_dir + "/" + dir_head + "*/")

analys_dir = os.getcwd() + "/model_analysis/"

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)


################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {"dsso1": 35.0, "dsso2": 35.0}

# Load module
AT = AnalysisTrajectories(
    out_dirs,
    dir_name=dir_head,
    analysis_dir=analys_dir,
    burn_in_fraction=0.2,
    nskip=200,
    nproc=nproc,
)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs=XLs_cutoffs)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
# AT.set_analyze_score_only_restraint(handle = 'MinimumPairDistanceBindingRestraint', short_name = 'MPDBR', do_sum = True)
# AT.set_analyze_EM_restraint()
# AT.restraint_names["GaussianEMRestraint"] = "GaussianEMRestraint"
# ConnectivityRestraint_None

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.get_psi_stats()

# What scores do we cluster on?
AT.hdbscan_clustering(["EV_sum", "XLs_sum"])  # , 'GaussianEMRestraint_None'
AT.summarize_XLs_info()
# exit()
