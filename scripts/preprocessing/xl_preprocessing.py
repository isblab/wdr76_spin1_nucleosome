import yaml
import os,sys


xlfile = sys.argv[1]
uniprot_mapping_file = sys.argv[2]
outxlf = os.path.join(os.getcwd(),"modeling_xlfile.dat")

with open(uniprot_mapping_file,'r') as mapf:
    mapping = yaml.safe_load(mapf)


def xl_filter(xl_file_name: str, uniprot_mapping: dict) -> list[str]:
    xlinks = []
    with open(xl_file_name,'r') as inxlf:
        for ln in inxlf.readlines():
            ln1 = ln.split(',')
            p1, p2 = ln1[0], ln1[4]
            if p1 in uniprot_mapping["mapping_from_xlfile"].keys() and p2 in uniprot_mapping["mapping_from_xlfile"].keys():
                p1name = uniprot_mapping["protein_names_to_use"][p1]
                p2name = uniprot_mapping["protein_names_to_use"][p2]
                ln1[0], ln1[4] = p1name, p2name
                ln1 = ','.join(ln1)
                xlinks.append(ln1)
    return xlinks


def write_modeling_xls_file(outxlf: str, filtered_xlinks: list[str]) -> None:
    with open(outxlf, 'w') as outf:
        outf.write("Protein1,Residue1,Protein2,Residue2\n")
        for ln in filtered_xlinks:
            ln1 = ln.split(',')
            p1, p2 = ln1[0], ln1[4]
            res1 = int(ln1[1]) + int(ln1[3]) - 1
            res2 = int(ln1[5]) + int(ln1[7]) - 1
            outf.write(f"{p1},{res1},{p2},{res2}\n")


filtered_xls = xl_filter(xlfile,mapping)
write_modeling_xls_file(outxlf, filtered_xls)