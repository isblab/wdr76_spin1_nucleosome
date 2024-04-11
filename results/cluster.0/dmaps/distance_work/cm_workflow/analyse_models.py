import sys
import IMP
import RMF
import math
import pickle
import hdbscan
import IMP.rmf
import IMP.atom
import IMP.core
import numpy as np
import IMP.algebra
from tqdm import tqdm
import matplotlib.pyplot as plt


def get_all_particle_xyzm(selection):
    xyzm = []
    for leaf in selection:
        p = IMP.core.XYZR(leaf)
        coords = list(p.get_coordinates())
        mass = IMP.atom.Mass(leaf)
        mass = float(str(mass).split(" ")[-1])
        coords.append(mass)
        xyzm.append(coords)

    return xyzm


def get_center_of_mass(all_xyzm):
    mass_x, mass_y, mass_z, total_mass = 0, 0, 0, 0

    for i in all_xyzm:
        mass_x += i[0] * i[-1]
        mass_y += i[1] * i[-1]
        mass_z += i[2] * i[-1]
        total_mass += i[-1]

    cm = (mass_x / total_mass, mass_y / total_mass, mass_z / total_mass)
    return cm


def measure_distance(coords1, coords2):
    x = (coords1[0] - coords2[0]) ** 2
    y = (coords1[1] - coords2[1]) ** 2
    z = (coords1[2] - coords2[2]) ** 2
    return math.sqrt(x + y + z)


rmf_fname = sys.argv[1]

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
mdl.update()

target_proteins = {"WDR76": 1, "H3": 2}
target_proteins_selection = {"WDR76": range(286), "H3": range(37)}

pairs = []
distances = {}

for frame_id in tqdm(
    range(rmf_fh.get_number_of_frames()),
    desc="Computing center of mass for the target proteins",
):
    c_masses_inframe = {}
    for prot, n_copies in target_proteins.items():
        IMP.rmf.load_frame(rmf_fh, frame_id)

        if n_copies == 1:
            # sel0 = IMP.atom.Selection(
            #     hier, molecule=prot, copy_index=0, resolution=1
            # ).get_selected_particles()
            sel0 = IMP.atom.Selection(
                hierarchy=hier,
                molecule=prot,
                residue_indexes=target_proteins_selection[prot],
                copy_index=0,
                resolution=1,
            ).get_selected_particles()
            xyzm = get_all_particle_xyzm(sel0)
            center_of_mass = get_center_of_mass(xyzm)
            c_masses_inframe[f"{prot}.{0}"] = center_of_mass

        elif n_copies > 1:
            for copy in range(n_copies):
                # sel0 = IMP.atom.Selection(
                #     hier, molecule=prot, copy_index=copy, resolution=1
                # ).get_selected_particles()
                sel0 = IMP.atom.Selection(
                    hierarchy=hier,
                    molecule=prot,
                    copy_index=copy,
                    residue_indexes=target_proteins_selection[prot],
                    resolution=1,
                ).get_selected_particles()
                xyzm = get_all_particle_xyzm(sel0)
                center_of_mass = get_center_of_mass(xyzm)
                c_masses_inframe[f"{prot}.{copy}"] = center_of_mass

    if len(pairs) == 0:
        keys = list(c_masses_inframe.keys())
        for k1 in keys:
            for k2 in keys:
                if k1.split(".")[0] != k2.split(".")[0]:
                    if (k2, k1) not in pairs:
                        pairs.append((k1, k2))

    for pair in pairs:
        if pair not in distances:
            distances[pair] = []
        dist = measure_distance(c_masses_inframe[pair[0]], c_masses_inframe[pair[1]])
        distances[pair].append(dist)

with open("distances", "wb") as pklf:
    pickle.dump(distances, pklf)

for pair in distances:
    plt.hist(distances[pair], bins=100, histtype="step", label=f"{pair[0]}-{pair[1]}")
plt.xlabel("Distances")
plt.legend()
plt.savefig("center of mass histogram.png", dpi=400)
plt.close()


### Do HDBScan clustering
print("Initializing HDBScan clustering...")
keys = list(distances.keys())
dmat = np.zeros((len(distances[keys[0]]), 2))
dmat[:, 0] = distances[keys[0]]
dmat[:, 1] = distances[keys[1]]

clusterer = hdbscan.HDBSCAN(min_cluster_size=50)
cluster_labels = clusterer.fit_predict(dmat)

colors = []
for idx, coords in enumerate(dmat):
    if cluster_labels[idx] != -1:
        color = f"C{cluster_labels[idx]}"
    else:
        color = "#aeaeae"
    colors.append(color)


plt.scatter(
    x=dmat[:, 0],
    y=dmat[:, 1],
    c=colors,
    alpha=0.8,
)
plt.xlabel(f"Distance from {keys[0]}")
plt.ylabel(f"Distance from {keys[1]}")
plt.savefig("center of mass clustered scatter.png")
plt.close()

plt.plot(
    [i for i in range(len(distances[keys[0]]))],
    dmat[:, 0],
    c="C1",
    alpha=0.8,
    label=keys[0],
)
plt.plot(
    [i for i in range(len(distances[keys[0]]))],
    dmat[:, 1],
    c="C2",
    alpha=0.6,
    label=keys[1],
)
plt.xlabel("Model ID")
plt.ylabel("Distance")
plt.legend()
plt.savefig("center of mass plot.png")
