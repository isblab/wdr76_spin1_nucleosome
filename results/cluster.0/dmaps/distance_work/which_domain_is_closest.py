import sys
import IMP
import RMF
import pickle
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
import pandas as pd
from tqdm import tqdm
import plotly.express as px


def get_winner(minimum_distances: list[dict[str, float]]):
    winners = []
    win_counts = {"Both": 0}
    for frame in minimum_distances:
        values = tuple(frame.items())
        p1, d1 = values[0][0], values[0][1]
        p2, d2 = values[1][0], values[1][1]

        if p1 not in win_counts.keys():
            win_counts[p1] = 0
        if p2 not in win_counts.keys():
            win_counts[p2] = 0

        if d1 == d2:
            winners.append("Both")
            win_counts["Both"] += 1
        elif d1 < d2:
            winners.append(p1)
            win_counts[p1] += 1
        elif d2 < d1:
            winners.append(p2)
            win_counts[p2] += 1
        else:
            print("Something really bad has happened...")
    return winners, win_counts


rmf_fname = "/home/shreyas/Projects/washburn/wdr_work/v3/analysis/sampcon/cluster.0/dmaps/test/sampcon_curr_cluster_models.rmf3"

ref_prot = {"H3": 2}
target_proteins_copies = {"SPIN1": 1, "WDR76": 1}

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
mdl.update()

nmodels = rmf_fh.get_number_of_frames()
mdl_ids = [i for i in range(nmodels)]


all_minimum_distances = []
for frame_id in tqdm(range(rmf_fh.get_number_of_frames())):
    IMP.rmf.load_frame(rmf_fh, frame_id)
    sel_ref = IMP.atom.Selection(
        hierarchy=hier, molecule=tuple(ref_prot.keys())[0]
    ).get_selected_particles()

    sel1 = IMP.atom.Selection(
        hierarchy=hier, molecule=tuple(target_proteins_copies.keys())[0]
    ).get_selected_particles()
    sel2 = IMP.atom.Selection(
        hierarchy=hier, molecule=tuple(target_proteins_copies.keys())[1]
    ).get_selected_particles()

    proteins = {
        tuple(target_proteins_copies.keys())[0]: sel1,
        tuple(target_proteins_copies.keys())[1]: sel2,
    }

    distances = {}
    interactor = None
    for prot, selection in proteins.items():
        temp_dists = []
        for ref_leaf in sel_ref:
            for leaf in selection:
                dist = IMP.core.get_distance(
                    IMP.core.XYZR(ref_leaf), IMP.core.XYZR(leaf)
                )
                dist = 0 if dist < 0 else dist
                temp_dists.append(dist)
        distances[prot] = min(temp_dists)

    all_minimum_distances.append(distances)

with open("min_distances.pkl", 'wb') as pklf:
    pickle.dump(all_minimum_distances, pklf)

interactors, interaction_counts = get_winner(all_minimum_distances)
print(interaction_counts)

fig = px.pie(
    values=list(interaction_counts.values()), names=list(interaction_counts.keys())
)


fig.show()
fig.write_html("Minimum distance comparision pie.html")
