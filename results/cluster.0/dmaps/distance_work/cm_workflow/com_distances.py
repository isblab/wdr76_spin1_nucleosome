import IMP
import RMF
import math
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
from tqdm import tqdm
import matplotlib.pyplot as plt
import plotly.express as px

###############################################################################
################################# User Inputs #################################
###############################################################################

rmf_fname = "/home/shreyas/Projects/washburn/wdr_work/v3/analysis/sampcon/cluster.0/dmaps/test/sampcon_curr_cluster_models.rmf3"
target_proteins = {"WDR76": 1, "H3": 2}
target_proteins_selection = {
    "WDR76": {
        "WDR_N": range(0, 285),
        "WD1": range(286, 327),
        "WD2": range(328, 392),
        "WD3": range(393, 438),
        "WD4": range(439, 485),
        "WD5": range(486, 532),
        "WD6": range(533, 581),
        "WD7": range(582, 621),
        "WDR_C": range(622, 626),
    },
    "SPIN1": {
        "SPIN1_N": range(1, 53),
        "SPIN1_Tudor1": range(54, 103),
        "SPIN1_Td12": range(104, 132),
        "SPIN1_Tudor2": range(133, 182),
        "SPIN1_Td23": range(183, 213),
        "SPIN1_Tudor3": range(214, 262),
    },
}
nucleosome_proteins = ["H2A", "H2B", "H3", "H4"]

###############################################################################
################################## Functions ##################################
###############################################################################


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


def get_center_of_mass_wrapper(
    rmf_fh,
    frame_id,
    target_selection,
    nucleosome_com,
):
    IMP.rmf.load_frame(rmf_fh, frame_id)

    min_dist = 1_000_000
    interacting_domain = None

    for prot_name, val in target_selection.items():
        for k, res_range in val.items():
            sel0 = IMP.atom.Selection(
                hierarchy=hier,
                molecule=prot_name,
                residue_indexes=res_range,
                copy_index=0,
                resolution=1,
            ).get_selected_particles()

            xyzm = get_all_particle_xyzm(sel0)
            sel_com = get_center_of_mass(xyzm)
            dist = measure_distance(nucleosome_com, sel_com)
            if dist < min_dist:
                min_dist = dist
                interacting_domain = k

    return (interacting_domain, min_dist)


def measure_distance(coords1, coords2):
    x = (coords1[0] - coords2[0]) ** 2
    y = (coords1[1] - coords2[1]) ** 2
    z = (coords1[2] - coords2[2]) ** 2
    return math.sqrt(x + y + z)


###############################################################################
################################## Main work ##################################
###############################################################################

mdl = IMP.Model()
rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
mdl.update()

# Get COM of fixed nucleosome
IMP.rmf.load_frame(rmf_fh, 0)
sel_nucleosome = IMP.atom.Selection(
    hierarchy=hier,
    molecules=nucleosome_proteins,
    resolution=1,
).get_selected_particles()
xyzm = get_all_particle_xyzm(sel_nucleosome)
nucleosome_com = get_center_of_mass(xyzm)

# Iterate over all frames to get the interactions
interacting_domains = []
for frame_id in tqdm(
    range(rmf_fh.get_number_of_frames()),
    desc="Computing center of mass for the target proteins",
):
    interaction = get_center_of_mass_wrapper(
        rmf_fh, frame_id, target_proteins_selection, nucleosome_com
    )
    interacting_domains.append(interaction)


actual_interacting_domains = {}
for i in interacting_domains:
    if i[0] not in actual_interacting_domains.keys():
        actual_interacting_domains[i[0]] = 1
    else:
        actual_interacting_domains[i[0]] += 1

print(actual_interacting_domains)

# plt.pie(
#     list(actual_interacting_domains.values()),
#     labels=list(actual_interacting_domains.keys()),
# )
# plt.show()

fig = px.pie(
    values=list(actual_interacting_domains.values()),
    names=list(actual_interacting_domains.keys()),
)
# fig.show()
fig.write_html("Minimum distance domains pie.html")
