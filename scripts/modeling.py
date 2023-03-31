
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
import ihm
import ihm.cross_linkers
import os
import sys

runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = 'run_' + str(runID)

num_frames = 0
if runType == "test":
    num_frames = 2500
elif runType == "prod":
    num_frames = 10_000

rex_max_temp = 2.24

# Identify data files
data_dir = "/home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/data/"

topology_file = os.path.join(data_dir,"topology.txt")
xl_data = os.path.join(data_dir,"xlinks",'modeling_xlfile_sheetA.dat')

xl_weight = 100

# Actual modeling begins here 
mdl = IMP.Model()
t = IMP.pmi.topology.TopologyReader(topology_file)
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

root_hier, dof = bs.execute_macro(max_rb_trans= 5,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 4,
                                  max_srb_trans= 0.01,
                                  max_srb_rot=0.04)

molecules = t.get_components()


################################################################################
########################## Fixing Particles ####################################
################################################################################

fixed_particles=[]
for prot in ["H2A","H2B","H3","H4"]:
    for cp in [0,1]:
        fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot,copy_index=cp).get_selected_particles()

fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover, #type: ignore
                                          IMP.pmi.TransformMover])



##### Uncomment the following lines to get test.rmf file to visualise the system representation
# IMP.atom.show_with_representations(root_hier)
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)



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
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)

print("Excluded volume restraint applied")


# CROSSLINKING RESTRAINT
xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name=xl_data,
                              converter=xldbkc)
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,
                database=xldb,
                length=25,
                resolution=1,
                slope=0.0001,
                label="dsso",
                weight=xl_weight,
                linker=ihm.cross_linkers.dsso)

output_objects.append(xlr)

print("Cross-linking restraint applied")


#####################################################
###################### SAMPLING #####################
#####################################################
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))

IMP.pmi.tools.shuffle_configuration(root_hier,          #type: ignore
                                    max_translation=50,
                                    excluded_rigid_bodies=fixed_rbs,
                                    hierarchies_included_in_collision=fixed_particles)

dof.optimize_flexible_beads(500)

evr.add_to_model()
xlr.add_to_model()

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=run_output_dir,
        output_objects=output_objects,
        monte_carlo_steps=10,
        number_of_best_scoring_models=0,
        number_of_frames=num_frames)

rex.execute_macro()

