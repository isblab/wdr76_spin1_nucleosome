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
top_dir = sys.argv[1]
# How are the trajectories dir names
dir_head = sys.argv[2]
c = sys.argv[3]  # Cluster we wish to extract
filter_applied = sys.argv[4]
nproc = 30
analys_dir = os.getcwd() + "/model_analysis/"

out_dirs = glob.glob(top_dir + "/" + dir_head + "*/")
print(out_dirs)

################################
# Extract frames
################################

# Load module
AT = AnalysisTrajectories(
    out_dirs, dir_name=dir_head, analysis_dir=analys_dir, nproc=nproc
)

# Point to the selected_models file
if filter_applied == "True":
    HA = AT.get_models_to_extract(
        analys_dir + "good_scoring_models_A_cluster" + str(c) + "_detailed.csv"
    )
    HB = AT.get_models_to_extract(
        analys_dir + "good_scoring_models_B_cluster" + str(c) + "_detailed.csv"
    )
else:
    HA = AT.get_models_to_extract(
        analys_dir + "selected_models_A_cluster" + str(c) + "_detailed.csv"
    )
    HB = AT.get_models_to_extract(
        analys_dir + "selected_models_B_cluster" + str(c) + "_detailed.csv"
    )

rmf_file_out_A = "A_models_clust" + str(c) + ".rmf3"
rmf_file_out_B = "B_models_clust" + str(c) + ".rmf3"

# Arguments for do_extract_models_single_rmf:
# HA :: Dataframe object from AT.get_models_to_extract()
# file_out :: The output rmf file
AT.do_extract_models_single_rmf(
    HA,
    rmf_file_out_A,  # RMF file outputted for Sample A
    top_dir,  # Top directory containing the PMI output folders
    analys_dir,  # Analysis directory to write RMF and score files
    scores_prefix="A_models_clust" + str(c),
)  # Prefix for the scores file

AT.do_extract_models_single_rmf(
    HB, rmf_file_out_B, top_dir, analys_dir, scores_prefix="B_models_clust" + str(c)
)
