def get_unique_uniprot_ids(xls_df) -> list:
    proteins = []
    for col in ["Protein1", "Protein2"]:
        for prot in xls_df[col]:
            proteins.append(prot.strip())
    return list(set(proteins))


def get_protein_name(uniprot_id) -> tuple[str, str]:
    from bioservices import UniProt
    uniprot_service = UniProt()
    
    result = uniprot_service.quick_search(uniprot_id)
    entry_name = result[uniprot_id]['Entry name'] # type: ignore
    return (uniprot_id, entry_name)


def get_sequence(uniprot_id: str) -> str:
    from bioservices import UniProt
    uniprot_service = UniProt()
    sequence: str = uniprot_service.search(uniprot_id, frmt="fasta") # type: ignore
    return sequence