import os
import sys

import IMP
import RMF
import IMP.core
import IMP.atom
import IMP.rmf
import argparse
import numpy as np
import concurrent.futures
import dmaps_functions
from matplotlib import pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Get contact maps for a given protein pair"
    )

    parser.add_argument(
        "--inputA",
        "-ia",
        dest="ia",
        help="cluster list of sample A RMFs. (eg. cluster.0.sample_A.txt)",
        required=True,
    )
    parser.add_argument(
        "--inputB",
        "-ib",
        dest="ib",
        help="cluster list of sample B RMFs. (eg. cluster.0.sample_B.txt)",
        required=True,
    )
    parser.add_argument(
        "--rmfA",
        "-ra",
        dest="ra",
        help="rmfA file. (eg. A_gsm_clust0.rmf3)",
        required=True,
    )
    parser.add_argument(
        "--rmfB",
        "-rb",
        dest="rb",
        help="rmfB file. (eg. B_gsm_clust0.rmf3)",
        required=True,
    )
    parser.add_argument(
        "--textA",
        "-ta",
        dest="ta",
        help="Text file associated with rmfA. (eg. A_gsm_clust0.txt)",
        required=True,
    )
    parser.add_argument(
        "--nprocs",
        "-p",
        dest="nprocs",
        help="Number of processes to be launched",
        required=True,
    )

    return parser.parse_args()


args = parse_args()


###################################################################################
########################### Get protein names and sizes ###########################
###################################################################################
mdl0 = IMP.Model()
rmf_fh0 = RMF.open_rmf_file_read_only(args.ra)
hier0 = IMP.rmf.create_hierarchies(rmf_fh0, mdl0)

all_proteins = dmaps_functions.get_protein_names(hier0)
sizes_dict = dmaps_functions.get_protein_sizes(hier0, all_proteins)


# Create list of model indices for sampleA, sample_B
nA = dmaps_functions.get_nmodels_in_A(args.ta)
RESOLUTION = 1

sample_A_models = []
sample_B_models = []

with open(args.ia, "r") as iaf:
    for ln in iaf.readlines():
        sample_A_models.append(int(ln.strip()))
with open(args.ib, "r") as ibf:
    for ln in ibf.readlines():
        sample_B_models.append(int(ln.strip()) - nA)

sample_A_models.sort()
sample_B_models.sort()
print(f"Total number of models:\t {len(sample_A_models) + len(sample_B_models)}")

sample_A_models = np.array_split(sample_A_models, int(args.nprocs))
print(f"Maximum number of models handled by a process:\t {len(sample_A_models[0])}")


for p1 in all_proteins:
    for p2 in all_proteins:
        if p1 == p2:
            continue

        s1, s2 = sizes_dict[p1], sizes_dict[p2]

        (
            distances,
            mean_distances,
        ) = dmaps_functions.get_mean_bead_distances_for_two_proteins(
            p1,
            p2,
            s1,
            s2,
            rmf_fh0,
            list(sample_A_models[0]),  # [:5],
            mdl0,
            hier0,
            resolution=RESOLUTION,
        )

        print(mean_distances)
        plt.figure(dpi=400)
        plt.imshow(mean_distances, cmap="hot")
        plt.xlabel(p2)
        plt.ylabel(p1)
        plt.colorbar()
        plt.savefig("tmp.png")
        break
    break
