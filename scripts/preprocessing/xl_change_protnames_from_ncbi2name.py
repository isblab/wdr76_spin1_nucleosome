import os
import yaml
import pandas as pd

xl_file = (
    "/home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/data/xlinks/xls_sheet4.csv"
)

out_xlfile = os.path.join(
    "/".join(xl_file.split("/")[:-1]), "modeling_xlfile_sheetD.dat"
)

xl_df = pd.read_csv(xl_file)

#! Pull the necessary columns out
imp_columns = ["Protein1", "PepPos1", "LinkPos1", "Protein2", "PepPos2", "LinkPos2"]
xl_df = xl_df[imp_columns]
xl_df["Residue1"] = xl_df["PepPos1"] + xl_df["LinkPos1"] - 1
xl_df["Residue2"] = xl_df["PepPos2"] + xl_df["LinkPos2"] - 1

xl_df = xl_df[["Protein1", "Residue1", "Protein2", "Residue2"]]

xl_df.to_csv(
    "/home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/data/xlinks/xls_p1r1p2r2.csv",
    index=False,
)


xlinks = []
with open(
    "/home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/data/xlinks/S4_ncbi_map.yaml",
    "r",
) as mapf:
    mapping = yaml.safe_load(mapf)["mapping_from_sheet5"]

with open(
    "/home/shreyas/Projects/washburn/wdr76_spin1_nucleosome/data/xlinks/xls_p1r1p2r2.csv",
    "r",
) as xlf:
    for ln in xlf.readlines():
        ln1 = ln.split(",")
        if ln1[0] in list(mapping.keys()):
            ln1[0] = mapping[ln1[0]]
        if ln1[2] in list(mapping.keys()):
            ln1[2] = mapping[ln1[2]]

        xlinks.append(",".join(ln1))


def xl_filter(xls: list, mapping: dict) -> list[str]:
    outl = []
    for link in xls:
        ln = link.split(",")
        if not ln[0] == "Protein1":
            if ln[0] in list(mapping.values()) and ln[2] in list(mapping.values()):
                outl.append(link)
        else:
            outl.append(link)
    return outl


xls = xl_filter(xlinks, mapping)

with open(out_xlfile, "w") as outxlf:
    for ln in xls:
        outxlf.write(ln)
