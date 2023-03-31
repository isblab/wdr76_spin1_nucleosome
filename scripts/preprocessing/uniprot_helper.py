def get_unique_protein_ids(xls_df) -> list:
    proteins = []
    for col in ["Protein1", "Protein2"]:
        for prot in xls_df[col]:
            proteins.append(prot.strip())
    return list(set(proteins))


def get_protein_name(prot_id) -> tuple[str, str]:
    # input is a tuple with protein_name, failed_searches list
    
    from bioservices import UniProt
    from Bio import Entrez

    uniprot_service = UniProt()
    Entrez.email = 'shreyasarvindekar@ncbs.res.in'
    db = "protein"
    
    failed = None
    prot_id = prot_id.strip()
    if  "_" in prot_id:
        try:
            eSearch = Entrez.esummary(db=db, id=prot_id)
            entry_name = Entrez.read(eSearch)[0]["Title"]
        except Exception as err:
            # print(prot_id,'\t',err,'\tin Entrez search')
            entry_name = prot_id    
    else:
        try:
            result = uniprot_service.quick_search(prot_id)
            entry_name = result[prot_id]['Entry name']
        except Exception as err:
            # print(prot_id,'\t',err,'\tin UniProt search')
            entry_name = prot_id
    
    return (prot_id, entry_name)
    


def get_sequence(uniprot_id: str) -> str:
    from bioservices import UniProt
    uniprot_service = UniProt()
    sequence: str = uniprot_service.search(uniprot_id, frmt="fasta") # type: ignore
    return sequence