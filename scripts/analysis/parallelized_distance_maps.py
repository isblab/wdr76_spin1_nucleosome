import numpy as np
import os


print("\n<-----------Creating Contact maps----------->")
cmap_path = os.getcwd()

protein1 = [
    "WDR76.0",
    "SPIN1.0",
    "H2A.0",
    "H2A.1",
    "H2B.0",
    "H2B.1",
    "H3.0",
    "H3.1",
    "H4.0",
    "H4.1",
]
protein2 = [
    "WDR76.0",
    "SPIN1.0",
    "H2A.0",
    "H2A.1",
    "H2B.0",
    "H2B.1",
    "H3.0",
    "H3.1",
    "H4.0",
    "H4.1",
]

prots = []
for index1, prot1 in enumerate(protein1):
    for index2, prot2 in enumerate(protein2):
        if prot1 == prot2:
            continue
        elif f"{prot2},{prot1}" in prots:
            continue
        else:
            if prot1.startswith("H") and prot2.startswith("H"):
                continue
            else:
                prots.append(f"{prot1},{prot2}")


print(f"Contact maps will be created for the following protein pairs -- \n{prots}")
print(f"No. of Contact maps created = {len( prots )}")


for prot in prots:
    os.system(
        f"~/imp-clean/build/setup_environment.sh \
		python /home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/scripts/analysis/contact_maps_all_pairs_surface.py \
		-ia ../../cluster.0.sample_A.txt \
		-ib ../../cluster.0.sample_B.txt \
		-ra ../../../model_analysis/A_models_clust1.rmf3 \
		-rb ../../../model_analysis/B_models_clust1.rmf3 \
		-c ../../cluster.0/cluster_center_model.rmf3 \
		-ta ../../../model_analysis/A_models_clust1.txt \
		-p {prot} &"
    )
