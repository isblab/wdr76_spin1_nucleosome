import IMP
import IMP.rmf
import IMP.core
import IMP.atom
import numpy as np
from tqdm import tqdm
import concurrent.futures


def get_nmodels_in_A(ta_file):
    with open(ta_file, "r") as taf:
        ln_count = 0
        for ln in taf.readlines():
            ln_count += 1
    return ln_count


def get_bead_name(particle):
    """
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    """

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number = IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    if IMP.atom.Fragment.get_is_setup(
        particle
    ):  # CG bead                      ###################### Did not understand get_is_setup
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = (
            mol_name
            + ":"
            + str(copy_number)
            + ":"
            + str(min(residues_in_bead))
            + ":"
            + str(max(residues_in_bead))
        )
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = (
            mol_name
            + ":"
            + str(copy_number)
            + ":"
            + residue_in_bead
            + ":"
            + residue_in_bead
        )

    return bead_name


def get_protein_names(hier: IMP.atom.Hierarchy) -> list[str]:
    proteins = []
    sel0 = IMP.atom.Selection(hierarchy=hier).get_selected_particles()

    for particle in sel0:
        name = get_bead_name(particle)

        protein_name = ".".join(name.split(":")[0:2])
        if protein_name not in proteins:
            proteins.append(protein_name)
    return proteins


def get_protein_sizes(hier: IMP.atom.Hierarchy, all_proteins: list[str]) -> dict:
    sizes_dict = {}
    for protein in all_proteins:
        if protein not in sizes_dict:
            sel0 = IMP.atom.Selection(
                hier, molecule=protein.split(".")[0]
            ).get_selected_particles()

            if "bead" in str(sel0[-1]):
                last_res = str(sel0[-1]).split("_")[0].split("-")[1]
            else:
                last_res = str(sel0[-1]).strip()[1:-1]

            sizes_dict[protein] = int(last_res)
    return sizes_dict


def get_all_res_from_cg_bead(bead: IMP.Particle) -> list[int]:
    res_range = []
    bead_str = str(bead)[1:-1]
    start_res = bead_str.split("_")[0].split("-")[0]
    end_res = bead_str.split("_")[0].split("-")[1]
    for res in range(int(start_res), int(end_res) + 1):
        res_range.append(res)
    return res_range


def measure_beadwise_distances(
    p1: str,
    p2: str,
    s1: int,
    s2: int,
    rmf_file_handle,
    frame_id: int,
    mdl: IMP.Model,
    hier: IMP.atom.Hierarchy,
    resolution: int,
) -> tuple[int, np.ndarray]:
    """Given a frame and a pair of proteins, return the all vs all distances for all beads in the two proteins in that frame"""

    distances = np.ones(shape=(s1 + 1, s2 + 1))
    IMP.rmf.load_frame(rmf_file_handle, frame_id)
    mdl.update()

    sel1 = IMP.atom.Selection(
        hier,
        resolution=resolution,
        molecule=p1.split(".")[0],
        copy_index=int(p1.split(".")[1]),
    )
    sel2 = IMP.atom.Selection(
        hier,
        resolution=resolution,
        molecule=p2.split(".")[0],
        copy_index=int(p2.split(".")[1]),
    )

    for bead1 in sel1.get_selected_particles():
        for bead2 in sel2.get_selected_particles():
            dist = IMP.core.get_distance(IMP.core.XYZR(bead1), IMP.core.XYZR(bead2))
            if dist < 0:
                dist = 0

            if "bead" not in str(bead1) and "bead" not in str(bead2):
                distances[int(str(bead1)[1:-1]), int(str(bead2)[1:-1])] = dist

            elif ("bead" in str(bead1)) and ("bead" not in str(bead2)):
                res_range = get_all_res_from_cg_bead(bead1)
                for res in res_range:
                    distances[res, int(str(bead2)[1:-1])] = dist

            elif ("bead" not in str(bead1)) and ("bead" in str(bead2)):
                res_range = get_all_res_from_cg_bead(bead2)
                for res in res_range:
                    distances[int(str(bead1)[1:-1]), res] = dist

            else:
                res_range1 = get_all_res_from_cg_bead(bead1)
                res_range2 = get_all_res_from_cg_bead(bead2)

                for res1 in res_range1:
                    for res2 in res_range2:
                        distances[res1, res2] = dist

    return frame_id, distances


def get_mean_bead_distances_for_two_proteins(
    p1: str,
    p2: str,
    s1: int,
    s2: int,
    rmf_file_handle,
    sample_models: list,
    mdl: IMP.Model,
    hier: IMP.atom.Hierarchy,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Given a pair of proteins a list of models, get average pairwise distances between them"""

    distances = np.ones(shape=(s1 + 1, s2 + 1, len(sample_models)))
    mean_distances = np.ones(shape=(s1 + 1, s2 + 1))
    distances_in_frame = np.zeros(shape=(s1 + 1, s2 + 1))

    mdl_id = 0
    for mdl_id, frame in enumerate(tqdm(sample_models)):
        distances_in_frame = measure_beadwise_distances(
            p1, p2, s1, s2, rmf_file_handle, frame, mdl, hier, resolution
        )

    distances[:, :, mdl_id] = distances_in_frame
    mean_distances: np.ndarray = np.mean(distances, axis=2)
    return distances, mean_distances
