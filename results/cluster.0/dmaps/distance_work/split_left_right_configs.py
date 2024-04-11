import os
import sys
import shutil
from sklearn.cluster import KMeans
from jenkspy import jenks_breaks
import IMP
import RMF
import IMP.rmf
import glob
from tqdm import tqdm
import numpy as np
import concurrent.futures
import matplotlib.pyplot as plt


sys.path.append("/home/shreyas/Projects/distance_maps_v3_temp/")
import dmaps_functions


def measure_beadwise_distances(
    p1: str,
    p2: str,
    s1: int,
    s2: int,
    rmf_file: str,
    resolution: int,
    selection1: range,
    selection2: range,
) -> np.ndarray:
    """Reads the concatednated RMF3 file and measures p1-p2 distances for all models"""

    rmf_file_handle = RMF.open_rmf_file_read_only(rmf_file)
    mdl = IMP.Model()
    hier = IMP.rmf.create_hierarchies(rmf_file_handle, mdl)

    distances = []
    for frame_id in tqdm(range(rmf_file_handle.get_number_of_frames())):
        IMP.rmf.load_frame(rmf_file_handle, frame_id)
        mdl.update()

        sel1 = IMP.atom.Selection(
            hier,
            residue_indexes=selection1,
            resolution=resolution,
            molecule=p1.split(".")[0],
            copy_index=int(p1.split(".")[1]),
        )
        sel2 = IMP.atom.Selection(
            hier,
            residue_indexes=selection2,
            resolution=resolution,
            molecule=p2.split(".")[0],
            copy_index=int(p2.split(".")[1]),
        )

        temp1, temp2 = 0, 0
        for _ in selection1:
            temp1 += 1
        for _ in selection2:
            temp2 += 1

        distances_in_frame = np.ones((temp1 + 2, temp2 + 2, 1)) * 1000
        for bead1 in sel1.get_selected_particles():
            for bead2 in sel2.get_selected_particles():
                dist = IMP.core.get_distance(IMP.core.XYZR(bead1), IMP.core.XYZR(bead2))
                if dist < 0:
                    dist = 0

                _, start1, end1 = dmaps_functions.get_bead_name(bead1)
                _, start2, end2 = dmaps_functions.get_bead_name(bead2)

                for r1 in range(start1, end1 + 1):
                    for r2 in range(start2, end2 + 1):
                        distances_in_frame[r1, r2] = dist

        distances.append(distances_in_frame)

    del rmf_file_handle
    out_distance_matrix = np.concatenate(distances, axis=2)
    out_distance_matrix = np.delete(
        np.delete(out_distance_matrix, 0, axis=0), 0, axis=1
    )

    return out_distance_matrix


def main():
    rmf_fname = "sampcon_curr_cluster_models.rmf3"
    nprocs = 120
    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)

    all_proteins = dmaps_functions.get_protein_names(hier)
    sizes_dict = dmaps_functions.get_protein_sizes(hier, all_proteins)
    nmodels = rmf_fh.get_number_of_frames()
    mdl_ids = [i for i in range(nmodels)]
    del rmf_fh

    dmaps_functions.split_models_into_subsets(rmf_fname, mdl_ids, nprocs)

    poi1 = "H3.0"
    poi2 = "H3.1"
    target_prot = "WDR76.0"
    selection_target = range(1, 276)
    selection_h3 = range(1, 36)

    target_prot_sz = sizes_dict[target_prot]

    distances = {}
    concatenated_rmfs = glob.glob("concatenated_models/*rmf3")
    for prot in (poi1, poi2):
        print(f"Processing {prot}:{target_prot} pair...")
        all_distances = []

        with concurrent.futures.ProcessPoolExecutor(nprocs) as executor:
            for dist in executor.map(
                measure_beadwise_distances,
                [target_prot for _ in range(nprocs)],
                [prot for _ in range(nprocs)],
                [target_prot_sz for _ in range(nprocs)],
                [sizes_dict[prot] for _ in range(nprocs)],
                concatenated_rmfs,
                [1 for _ in range(nprocs)],
                [selection_target for _ in range(nprocs)],
                [selection_h3 for _ in range(nprocs)],
            ):
                all_distances.append(dist)

        all_distances = np.concatenate(all_distances, axis=2)
        distances[f"dist_{prot}_{target_prot}"] = all_distances.min(axis=(0, 1))

    i = 1
    for k, val in distances.items():
        plt.hist(x=val, bins=100, histtype="step", color=f"C{i}", label=k)
        i += 1

    plt.legend()
    plt.xlabel("Distances")
    plt.savefig("wdrN_h3N_min_histogram.png")
    plt.close()

    print("\nInitiating clustering...")
    coords = []
    keys = tuple(distances.keys())
    for i in range(len(distances[keys[0]])):
        coords.append((distances[keys[0]][i], distances[keys[1]][i]))

    matrix = np.array(coords)
    print(f"Number of points to cluster: {len(coords)}")
    kmeans = KMeans(n_clusters=3)
    clusters = kmeans.fit(matrix)

    print("Clustering done...")
    x, y, c = [], [], []
    for i, mat in enumerate(coords):
        x.append(mat[0])
        y.append(mat[1])
        c.append(f"C{clusters.labels_[i]+1}")

    plt.scatter(x, y, c=c, alpha=0.6)

    plt.xlabel(f"Distance from {poi1}")
    plt.ylabel(f"Distance from {poi2}")
    plt.savefig("wdrN_h3N_min_scatter.png")


main()
shutil.rmtree(os.path.join(os.getcwd(), "concatenated_models"))
