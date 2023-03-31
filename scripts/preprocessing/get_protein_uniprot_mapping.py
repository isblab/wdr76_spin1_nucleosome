"""
Given an input XL file and a list of proteins 
"""
import yaml
import argparse
import uniprot_helper
import pandas as pd
import concurrent.futures
    

############################################################
####################### User Inputs ########################
parser = argparse.ArgumentParser()
parser.add_argument("-x", "--xl_file", type=str, required=True, help="Path to the XL file")     
# xl_file = "/home/shreyas/Dropbox/washburn_wdr_spin/xls_sheet1.csv"
parser.add_argument("-m","--mode", type=str, default="prot_names", help="prot_names or mapping")
parser.add_argument("-p","--poi", type=str, default=None, help="File containing the names of proteins of interest")
args = parser.parse_args()

df = pd.read_csv(args.xl_file)
todo = args.mode
poi_file = args.poi


def get_all_proteins_from_xl_and_fetch_pnames(xls_df) -> list[tuple[str, str]]:
    all_proteins = uniprot_helper.get_unique_protein_ids(xls_df)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:      #! Do not increase the max_workers parameter value, as higher values return HTTP error 429 on Entrez
        prot_details = []
        failed_searches = []
        for pdeet in executor.map(uniprot_helper.get_protein_name, all_proteins):
        # for i in all_proteins:
        #     pdeet,f_search = uniprot_helper.get_protein_name(i)
            if pdeet[0]!=pdeet[1]:
                prot_details.append(pdeet)
            else:
                failed_searches.append(pdeet[0])
    return prot_details, failed_searches


def get_sequences(prot_ids):
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        sequences = []
        for prot_seq in executor.map(uniprot_helper.get_sequence, list(prot_ids.keys())):
            sequences.append(prot_seq)

    with open("sequences.fasta", "w") as seqf:
        for ln in sequences:
            seqf.write(f"{ln}\n")



if __name__=="__main__":
    protein_details, failed_searches = get_all_proteins_from_xl_and_fetch_pnames(df)
    print(protein_details)
    print(failed_searches)

    if todo == "prot_details":
        print("From the above output, add the names of proteins of interest to a file and rerun this script with mode=mapping")
        exit()

    elif todo == "mapping":
        poi = []
        with open(poi_file,'r') as poif:
            for ln in poif.readlines():
                if ln.strip() not in poi:
                    poi.append(ln.strip())

        necessary_uniprot_ids = {"mapping_from_xlfile": {}}
        with open('uniprot_mapping.yaml','w') as umf:
            for uni_id, pname in protein_details:
                if pname in poi:
                    if uni_id not in necessary_uniprot_ids:
                        necessary_uniprot_ids["mapping_from_xlfile"][uni_id] = pname
            
            yaml.dump(necessary_uniprot_ids,umf)
            
        get_sequences(necessary_uniprot_ids["mapping_from_xlfile"])