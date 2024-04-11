##########################################################################################
################# IMP Modeling Script for WDR76/SPIN1/Nucleosome Complex #################
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
import os
import sys

# Imports needed to use ProtocolOutput
import IMP.pmi.mmcif
import ihm
import ihm.location
import ihm.model
import ihm.cross_linkers


runType = sys.argv[1]  # Specify test or prod
runID = sys.argv[2]  # Specify the number of runs
run_output_dir = "run_" + str(runID)

num_frames = 0
if runType == "test":
    num_frames = 2000
elif runType == "prod":
    num_frames = 15_000
elif runType == "mmcif":
    num_frames = 5

rex_max_temp = 2.6

# Identify data files
data_dir = "../../data/"

topology_file = os.path.join(data_dir, "topology.txt")
xl_data_set1 = os.path.join(data_dir, "xlinks", "modeling_xlfile_sheetA.dat")
xl_data_set2 = os.path.join(data_dir, "xlinks", "modeling_xlfile_sheetD.dat")

xl_weight = 10

# Actual modeling begins here
mdl = IMP.Model()
t = IMP.pmi.topology.TopologyReader(topology_file)
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# Add deposition information
po = IMP.pmi.mmcif.ProtocolOutput()
bs.system.add_protocol_output(po)
po.system.title = "Integrative structure of the human NuDe complex"
# po.system.citations.append(ihm.Citation.from_pubmed_id(000000)) #TODO

root_hier, dof = bs.execute_macro(
    max_rb_trans=5,
    max_rb_rot=0.1,
    max_bead_trans=4,
    max_srb_trans=0.01,
    max_srb_rot=0.04,
)

molecules = t.get_components()


################################################################################
########################## Fixing Particles ####################################
################################################################################

fixed_particles = []
for prot in ["H2A", "H2B", "H3", "H4"]:
    for cp in [0, 1]:
        fixed_particles += IMP.atom.Selection(
            root_hier, molecule=prot, copy_index=cp
        ).get_selected_particles()

fixed_beads, fixed_rbs = dof.disable_movers(
    fixed_particles, [IMP.core.RigidBodyMover, IMP.pmi.TransformMover]
)


##### Uncomment the following lines to get test.rmf file to visualise the system representation
# IMP.atom.show_with_representations(root_hier)
# fname = "test.rmf"
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)

# exit()

#####################################################
##################### RESTRAINTS ####################
#####################################################

output_objects = []

# CONNECTIVITY RESTRAINT
for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# EXCLUDED VOLUME RESTRAINT
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    included_objects=[root_hier], resolution=1000
)
output_objects.append(evr)

print("Excluded volume restraint applied")


# CROSSLINKING RESTRAINT
# -- Set1
xldbkc1 = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc1.set_standard_keys()

xldb1 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb1.create_set_from_file(file_name=xl_data_set1, converter=xldbkc1)
xlr1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,
    database=xldb1,
    length=25,
    resolution=1,
    slope=0.0001,
    label="dsso1",
    weight=xl_weight,
    linker=ihm.cross_linkers.dsso,
)

# -- Set2
xldbkc2 = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc2.set_standard_keys()

xldb2 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb2.create_set_from_file(file_name=xl_data_set2, converter=xldbkc2)
xlr2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,
    database=xldb2,
    length=25,
    resolution=1,
    slope=0.0001,
    label="dsso2",
    weight=xl_weight,
    linker=ihm.cross_linkers.dsso,
)

output_objects.append(xlr1)
output_objects.append(xlr2)

print("Cross-linking restraint applied")


#####################################################
###################### SAMPLING #####################
#####################################################
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))

IMP.pmi.tools.shuffle_configuration(
    root_hier,
    max_translation=50,
    excluded_rigid_bodies=fixed_rbs,
    hierarchies_included_in_collision=fixed_particles,
)

dof.optimize_flexible_beads(500)

evr.add_to_model()
xlr1.add_to_model()
xlr2.add_to_model()

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

rex = IMP.pmi.macros.ReplicaExchange0(
    mdl,
    root_hier=root_hier,
    monte_carlo_temperature=1.0,
    replica_exchange_minimum_temperature=1.0,
    replica_exchange_maximum_temperature=rex_max_temp,
    monte_carlo_sample_objects=dof.get_movers(),
    global_output_directory=run_output_dir,
    output_objects=output_objects,
    monte_carlo_steps=10,
    number_of_best_scoring_models=0,
    number_of_frames=num_frames,
)

rex.execute_macro()

# -----------------------------
# Finalize the protocol output
po.finalize()
s = po.system
import ihm.dumper

with open("model_init_wdr76-spin1-nucleosome.cif", "w") as fh:
    ihm.dumper.write(fh, [s])

# Datasets for XL-MS restraint
for r in s.restraints:
    if isinstance(r, ihm.restraint.CrossLinkRestraint):
        print("XL-MS dataset at:", r.dataset.location.path)
        print("Details:", r.dataset.location.details)


last_step = s.orphan_protocols[-1].steps[-1]
last_step.num_models_end = (
    750_000  # 15,000 models per run and 50 independent runs (8 cores per run)
)

protocol = po.system.orphan_protocols[-1]
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis_nmodels = 27314  # nModels in cluster0
analysis.steps.append(
    ihm.analysis.ClusterStep(
        feature="RMSD",
        num_models_begin=last_step.num_models_end,
        num_models_end=analysis_nmodels,
    )
)

mg = ihm.model.ModelGroup(name="Cluster 0")
po.system.state_groups[-1][-1].append(mg)
e = ihm.model.Ensemble(
    model_group=mg,
    num_models=analysis_nmodels,
    post_process=analysis.steps[-1],
    name="Cluster 0",
    clustering_method="Density based threshold-clustering",
    clustering_feature="RMSD",
    precision="24",
)
po.system.ensembles.append(e)

Uniprot = {
    "WDR76.0": "Q9H967",
    "SPIN1.0": "Q9Y657",
    "H2A.0": "P0C0S8",
    "H2A.1": "P0C0S8",
    "H2B.0": "P62807",
    "H2B.1": "P62807",
    "H3.0": "P68431",
    "H3.1": "P68431",
    "H4.0": "P62805",
    "H4.1": "P62805",
}
lpep = ihm.LPeptideAlphabet()


for prot, entry in Uniprot.items():
    ref = ihm.reference.UniProtSequence.from_accession(entry)
    ref.alignments.append(ihm.reference.Alignment())
    po.asym_units[prot].entity.references.append(ref)

m = IMP.Model()
inf1 = RMF.open_rmf_file_read_only("../../results/cluster.0/cluster_center_model.rmf3")
h = IMP.rmf.create_hierarchies(inf1, m)[0]
IMP.rmf.link_hierarchies(inf1, [h])
IMP.rmf.load_frame(inf1, RMF.FrameID(0))
m.update()

model = po.add_model(e.model_group)

#! TODO?
# repo = ihm.location.Repository(
#     doi="10.5281/zenodo.6674232",
#     root="../..",
#     top_directory="nurd_zenodo",
#     url="https://zenodo.org/record/6674232/files/nurd_zenodo.zip",
# )

loc_density_list = {
    "WDR76.0": [
        "LPD_wdr76_N",
        "LPD_wdr76_wd1",
        "LPD_wdr76_wd2",
        "LPD_wdr76_wd3",
        "LPD_wdr76_wd4",
        "LPD_wdr76_wd5",
        "LPD_wdr76_wd6",
        "LPD_wdr76_wd7",
    ],
    "SPIN1.0": [
        "LPD_spin1_N",
        "LPD_spin1_td12",
        "LPD_spin1_td23",
        "LPD_spin1_tudor1",
        "LPD_spin1_tudor2",
        "LPD_spin1_tudor3",
    ],
    "H2A.0": ["LPD_H2A0", "LPD_H2A0_n", "LPD_H2A0_c"],
    "H2A.1": ["LPD_H2A1", "LPD_H2A1_n", "LPD_H2A1_c"],
    "H2B.0": ["LPD_H2B0", "LPD_H2B0_n"],
    "H2B.1": ["LPD_H2B1", "LPD_H2B1_n"],
    "H3.0": ["LPD_H30", "LPD_H30_n"],
    "H3.1": ["LPD_H31", "LPD_H31_n"],
    "H4.0": ["LPD_H40", "LPD_H40_n"],
    "H4.1": ["LPD_H41", "LPD_H41_n"],
}

for prot, density in loc_density_list.items():
    asym = po.asym_units[prot]
    for domain_density in density:
        loc = ihm.location.OutputFileLocation(
            "../../results/cluster.0/" + domain_density + ".mrc"
        )
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
        e.densities.append(den)

# po.system.update_locations_in_repositories([repo])
po.finalize()

with open("model_wdr76-spin1-nucleosome.cif", "w") as fh:
    ihm.dumper.write(fh, [po.system])

import ihm.reader

with open("model_wdr76-spin1-nucleosome.cif") as fh:
    (s,) = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)
