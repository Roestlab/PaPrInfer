
import sqlite3

# noinspection PyUnresolvedReferences
from typing import Any, Dict, List, Set, Union

from pyopenms import IDFilter, IdXMLFile


def main(epifany_file: str, idpicker_file: str, pyprophet_file: str,
         threshold: str, remove_decoy: bool, return_qvalue: bool,
         return_distinct_protein: bool, return_mapping: bool):

    num_epifany_distinct_protein = len(get_epifany_result(epifany_file,
                                                          threshold,
                                                          remove_decoy,
                                                          return_qvalue,
                                                          return_distinct_protein,
                                                          return_mapping))

    num_idpicker_distinct_protein = get_idpicker_numbers(idpicker_file,
                                                         threshold)

    num_pyprophet_distinct_protein = get_pyprophet_numbers(pyprophet_file,
                                                           threshold)

    return num_epifany_distinct_protein, \
        num_idpicker_distinct_protein, \
        num_pyprophet_distinct_protein


def get_epifany_result(epifany_file: str, threshold: str, remove_decoy: bool,
                       return_qvalue: bool, return_distinct_protein: bool,
                       return_mapping: bool) -> \
        Union[Set[str], List[int], Dict[str, int]]:
    # getting all protein from epifany
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    all_protein_groups = []

    protein_group_probability = {}

    if remove_decoy:
        remove_decoy_protein_groups(prot_ids)

    # read the protein groups
    for protein_id in prot_ids:

        # because I put all the target groups in getProteinGroup
        if remove_decoy:
            idxml_protein_group_list = protein_id.getProteinGroups()
        else:
            idxml_protein_group_list = protein_id.getIndistinguishableProteins()

        for group in idxml_protein_group_list:

            # list of bytes representing the accessions
            accession_list_bytes = group.accessions

            # a list of str representing the accessions
            accession_list_str = [accession_bytes.decode("utf-8") for
                                  accession_bytes in accession_list_bytes]

            # sort it and strip it
            sorted_accession_list_str = sorted(accession_list_str)
            stripped_sorted_accession_list_str = [s.strip() for s in sorted_accession_list_str]

            # this means use this string ' ' to join the list
            accessions = ' '.join(stripped_sorted_accession_list_str)

            all_protein_groups.append(accessions)
            protein_group_probability.setdefault(accessions, []).append(
                group.probability)

    # a list of string, each string is a protein group
    # separate by a whitespace
    # all_protein_groups

    # only protein group below the threshold
    epifany_proteins_groups = []

    qvalue_list = []

    # protein_group_probability is a dict
    for groups in protein_group_probability:
        probability = min(protein_group_probability[groups])

        qvalue_list.append(probability)

        if probability <= float(threshold):
            epifany_proteins_groups.append(groups)

    epifany_distinct_proteins_groups = set(epifany_proteins_groups)
    print("epifany groups", len(epifany_distinct_proteins_groups))


    if return_qvalue:
        return qvalue_list
    elif return_distinct_protein:
        return epifany_distinct_proteins_groups
    elif return_mapping:
        return protein_group_probability
    else:
        print("you are returning nothing")


def remove_decoy_protein_groups(prot_ids) -> None:
    # first remove all the decoy protein hits, it remove from getHits()
    IDFilter().removeDecoyHits(prot_ids)
    for protein_id in prot_ids:
        # then get the protein groups
        groups = protein_id.getIndistinguishableProteins()
        # then get all the (now is updated to only target) protein hits
        hits = protein_id.getHits()
        # based on the target protein hits, remove all decoy protein group
        IDFilter().updateProteinGroups(groups, hits)

        # add the protein group back in

        # in epifany output protein groups is empty anyway
        # check that
        assert protein_id.getProteinGroups() == []

        # append only target protein groups into "protein groups"
        for group in groups:
            protein_id.insertProteinGroup(group)


def get_idpicker_numbers(idpicker_file: str, threshold: str) -> int:
    """now for idpicker """
    con = sqlite3.connect(idpicker_file)

    c = con.cursor()

    # PEP
    # c.execute(
    #     """SELECT PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID, SCORE_PROTEIN_GROUP.PEP
    #     FROM PROTEIN_GROUP
    #     INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
    #     WHERE PROTEIN_GROUP.DECOY = 0 AND PEP <= :threshold""", {'threshold': float(threshold)}
    # )

    # q_value
    c.execute(
        """SELECT DISTINCT PROTEIN_GROUP.PROTEIN_GROUP_ID, SCORE_PROTEIN_GROUP.QVALUE
        FROM PROTEIN_GROUP
        INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
        WHERE PROTEIN_GROUP.DECOY = 0 AND SCORE_PROTEIN_GROUP.QVALUE <= :threshold""",
        {'threshold': float(threshold)}
    )

    idpicker_q_values = []

    for row in c.fetchall():
        q_value = row[1]

        idpicker_q_values.append(q_value)

    con.commit()
    con.close()

    print("IDpicker", len(idpicker_q_values))

    return len(idpicker_q_values)


def get_idpicker_accessions(idpicker_file: str, threshold: str) -> Set[str]:

    con = sqlite3.connect(idpicker_file)

    c = con.cursor()

    # I cannot remember why I need distinct, was it like inner join
    # produce multiple identical row or something?
    c.execute(
        """SELECT DISTINCT PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID
        FROM PROTEIN_GROUP
        INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
        WHERE PROTEIN_GROUP.DECOY = 0 AND SCORE_PROTEIN_GROUP.QVALUE <= :threshold""",
        {'threshold': float(threshold)}
    )

    protein_group_to_accessions = {}

    for row in c.fetchall():
        protein_group_id = row[0]
        protein_accession = row[1]

        # if the protein_group_id already exist in dict, then get its value
        # if not give a empty list and append the protein accession
        protein_group_to_accessions.setdefault(protein_group_id, []).append(
            protein_accession)

    con.commit()
    con.close()


    all_accessions = []
    for protein_group, protein_accession_list in protein_group_to_accessions.items():

        # sort and strip them
        sorted_protein_accession_list = sorted(protein_accession_list)
        stripped_sorted_protein_accession_list = [s.strip() for s in sorted_protein_accession_list]

        # this means use this string ' ' to join the list
        accessions = ' '.join(stripped_sorted_protein_accession_list)

        accessions.strip()
        all_accessions.append(accessions)

    return set(all_accessions)



def get_pyprophet_numbers(pyprophet_file: str, threshold: str) -> int:
    # also need the file with swissprot instead of uniprot

    con = sqlite3.connect(pyprophet_file)

    c = con.cursor()

    c.execute(
        """ 
        SELECT SCORE_PROTEIN.PROTEIN_ID
        FROM SCORE_PROTEIN
        INNER JOIN PROTEIN ON PROTEIN.ID = SCORE_PROTEIN.PROTEIN_ID
        WHERE SCORE_PROTEIN.CONTEXT = 'global'
        AND QVALUE <= :threshold
        AND PROTEIN.DECOY = 0
        """, {'threshold': float(threshold)}
    )

    pyprophet_proteins = []

    for row in c.fetchall():
        protein = row[0]

        pyprophet_proteins.append(protein)

    con.commit()
    con.close()

    pyprophet_distinct_protein = set(pyprophet_proteins)
    print("pyprophet", len(pyprophet_distinct_protein))

    return len(pyprophet_distinct_protein)


def get_pyprophet_accessions(pyprophet_file: str, threshold: str) -> Set[str]:
    con = sqlite3.connect(pyprophet_file)

    c = con.cursor()

    c.execute(
        """select PROTEIN_ACCESSION
        from PROTEIN 
        inner join SCORE_PROTEIN on PROTEIN.ID = SCORE_PROTEIN.PROTEIN_ID
        where DECOY = 0 and QVALUE <=  :threshold""",
        {'threshold': float(threshold)}
    )

    pyprophet_proteins = []

    for row in c.fetchall():
        protein_accession = row[0]

        protein_accession.strip()

        pyprophet_proteins.append(protein_accession)

    con.commit()
    con.close()

    pyprophet_distinct_protein = set(pyprophet_proteins)

    return pyprophet_distinct_protein



# if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # epifany output file, idpicker output file, pyprophet file q-value threshold
    # whether to remove decoys
    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
