import IMP
import RMF
import math
import hdbscan
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


class Protein:
    def __init__(
        self,
        prot_name: str,
        cp_id: int,
        resrange: dict[str, range],
        hierarchy: IMP.atom.Hierarchy,
    ):
        self.prot = prot_name
        self.prot_copy = cp_id
        self.res_range = resrange
        self.hierarchy = hierarchy
        self.selection = {}
        self.do_selection()
        self.com = {}
        self.get_all_particle_xyzm()

    def do_selection(self):
        for domain, res_range in self.res_range.items():
            sel0 = IMP.atom.Selection(
                hierarchy=self.hierarchy,
                molecule=self.prot,
                copy_index=self.prot_copy,
                residue_indexes=res_range,
            ).get_selected_particles()
            self.selection[domain] = sel0

    def get_all_particle_xyzm(self):
        def get_center_of_mass(all_xyzm):
            mass_x, mass_y, mass_z, total_mass = 0, 0, 0, 0

            for i in all_xyzm:
                mass_x += i[0] * i[-1]
                mass_y += i[1] * i[-1]
                mass_z += i[2] * i[-1]
                total_mass += i[-1]

            cm = (mass_x / total_mass, mass_y / total_mass, mass_z / total_mass)
            return cm

        for k, leaf in self.selection.items():
            xyzm = []
            for particle in leaf:
                p = IMP.core.XYZR(particle)
                coords = list(p.get_coordinates())
                mass = IMP.atom.get_mass(particle)
                coords.append(mass)
                xyzm.append(coords)

            self.com[k] = get_center_of_mass(xyzm)

    def __str__(self) -> str:
        return f"{self.prot}, {self.prot_copy}"


def do_protein_selections(
    rmf_fh0: RMF.FileConstHandle,
    frame_id: int,
    hier0: IMP.atom.Hierarchy,
    ref_prot_copies: dict,
    ref_prot_resrange: dict,
    target_prot_copies: dict,
    target_prot_resrange: dict,
) -> tuple[list, list]:
    IMP.rmf.load_frame(rmf_fh0, frame_id)

    ref_proteins = []
    for pname, copies in ref_prot_copies.items():
        for cp in range(copies):
            ref_proteins.append(
                Protein(
                    prot_name=pname,
                    cp_id=cp,
                    resrange=ref_prot_resrange[pname],
                    hierarchy=hier0,
                )
            )

    target_proteins = []
    for pname, copies in target_prot_copies.items():
        for cp in range(copies):
            target_proteins.append(
                Protein(
                    prot_name=pname,
                    cp_id=cp,
                    resrange=target_prot_resrange[pname],
                    hierarchy=hier0,
                )
            )

    return ref_proteins, target_proteins


if __name__ == "__main__":
    rmf_fname = "/home/shreyas/Projects/washburn/wdr_work/v3/analysis/sampcon_v2/cluster.0/dmaps/distance_work/sampcon_curr_cluster_models.rmf3"

    mdl = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
    hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)
    mdl.update()

    ref_prot_cps = {"H3": 2}
    target_prot_cps = {"SPIN1": 1, "WDR76": 1}
    ref_prot_res_range = {"H3": {"full": range(137)}}
    target_prot_res_range = {
        "SPIN1": {
            "N": range(54),
            "Tudor1": range(54, 104),
            "Td12": range(104, 133),
            "Tudor2": range(133, 183),
            "Td23": range(183, 214),
            "Tudor3": range(214, 263),
        },
        "WDR76": {
            "N": range(286),
            "WD1": range(286, 328),
            "WD2": range(328, 393),
            "WD3": range(393, 439),
            "WD4": range(439, 486),
            "WD5": range(486, 533),
            "WD6": range(533, 582),
            "WD7": range(582, 627),
        },
    }

    ref_proteins_all_mdls, target_proteins_all_mdls = [], []
    for frame_id in tqdm(range(rmf_fh.get_number_of_frames())):
        ref, target = do_protein_selections(
            rmf_fh0=rmf_fh,
            frame_id=frame_id,
            hier0=hier,
            ref_prot_copies=ref_prot_cps,
            ref_prot_resrange=ref_prot_res_range,
            target_prot_copies=target_prot_cps,
            target_prot_resrange=target_prot_res_range,
        )
        ref_proteins_all_mdls.append(tuple(ref))
        target_proteins_all_mdls.append(tuple(target))
        # if frame_id == 100:
        #     break

    all_distances = {}
    temp_ref_prots, temp_target_prots = (
        ref_proteins_all_mdls[0],
        target_proteins_all_mdls[0],
    )
    for ref_prot in temp_ref_prots:
        for target_prot in temp_target_prots:
            if (
                f"{target_prot.prot}.{target_prot.prot_copy}-{ref_prot.prot}.{ref_prot.prot_copy}"
                not in all_distances
            ):
                all_distances[
                    f"{target_prot.prot}.{target_prot.prot_copy}-{ref_prot.prot}.{ref_prot.prot_copy}"
                ] = []

    wins = {}
    min_distances = {}
    for ref_in_mdl, target_in_mdl in zip(
        ref_proteins_all_mdls, target_proteins_all_mdls
    ):
        target_min = ("Dummy", 1_000_000)
        for target_protein in target_in_mdl:
            min_dist = ("Dummy", 1_000_000)
            for ref_protein in ref_in_mdl:
                ref_prot_mindist = 1_000_000
                for domain1, com1 in ref_protein.com.items():
                    for domain2, com2 in target_protein.com.items():
                        dist = math.dist(com1, com2)
                        dist = 0 if dist < 0 else dist

                        if dist < min_dist[1]:
                            min_dist = (f"{domain1}-{domain2}", dist)
                        if dist < ref_prot_mindist:
                            ref_prot_mindist = dist

                all_distances[
                    f"{target_protein.prot}.{target_protein.prot_copy}-{ref_protein.prot}.{ref_protein.prot_copy}"
                ].append(ref_prot_mindist)

            if target_protein.prot not in min_distances:
                min_distances[target_protein.prot] = []

            min_distances[target_protein.prot].append(min_dist[1])

            if min_dist[1] == target_min[1]:
                target_min = ("Both", min_dist[1])
            elif min_dist[1] < target_min[1]:
                target_min = (target_protein.prot, min_dist[1])

        if target_min[0] not in wins:
            wins[target_min[0]] = 1
        else:
            wins[target_min[0]] += 1

    print(wins)
    plt.pie(x=list(wins.values()), labels=list(wins.keys()))
    plt.savefig("wins_pie.png", bbox_inches="tight", dpi=1200)
    plt.close()

    for i, prot in enumerate(min_distances):
        plt.plot(min_distances[prot], label=prot, alpha=0.8)
        mean = np.mean(min_distances[prot])
        plt.plot(
            [mean for _ in range(len(min_distances[prot]))],
            linestyle="dashed",
            c=f"C{i}",
            label=f"Average {prot}-H3 distance",
        )

    plt.xlabel("Model index")
    plt.ylabel("Distances")
    plt.legend(bbox_to_anchor=(1.15, 1.0), loc="upper left")
    plt.savefig("mdlike_distances.png", bbox_inches="tight", dpi=1200)
    plt.close()

    ## Do HDBScan clustering
    print("Initiating HDBScan clustering")
    clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
    dmat = list(zip(min_distances["SPIN1"], min_distances["WDR76"]))
    dmat = np.array(dmat)

    clusters = clusterer.fit_predict(dmat)
    plt.scatter(
        x=min_distances["SPIN1"], y=min_distances["WDR76"], c=clusters, alpha=0.8
    )
    plt.xlabel("Distance from SPIN1")
    plt.ylabel("Distance from WDR76")
    plt.savefig("hdbscan_scatter.png", bbox_inches="tight", dpi=1200)
    plt.close()

    for k, val in all_distances.items():
        if "SPIN" in k:
            cp, cpid = None, int(k.split("-")[1].split(".")[1])
            if cpid > 1:
                print(f"Err: {cpid}")
                raise ValueError("More than two copies found")
            else:
                if cpid == 0:
                    cp = "first copy"
                elif cpid == 1:
                    cp = "second copy"

            k = f"{k.split('.')[0]}-{k.split('-')[1].split('.')[0]} {cp}"
            plt.plot(val, label=k, alpha=0.8)

    plt.legend(bbox_to_anchor=(1.15, 1.0), loc="upper left")
    plt.xlabel("Model index")
    plt.ylabel("Distances")
    plt.savefig("spin-h3_distances_mdlike.png", bbox_inches="tight", dpi=1200)
    plt.close()

    for k, val in all_distances.items():
        if "WDR" in k:
            cp, cpid = None, int(k.split("-")[1].split(".")[1])
            if cpid > 1:
                print(f"Err: {cpid}")
                raise ValueError("More than two copies found")
            else:
                if cpid == 0:
                    cp = "first copy"
                elif cpid == 1:
                    cp = "second copy"

            k = f"{k.split('.')[0]}-{k.split('-')[1].split('.')[0]} {cp}"
            plt.plot(val, label=k, alpha=0.8)

    plt.legend(bbox_to_anchor=(1.15, 1.0), loc="upper left")
    plt.xlabel("Model index")
    plt.ylabel("Distances")
    plt.savefig("wdr-h3_distances_mdlike.png", bbox_inches="tight", dpi=1200)
    plt.close()
