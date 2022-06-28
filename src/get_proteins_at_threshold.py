
import sqlite3

# noinspection PyUnresolvedReferences
from typing import Any, List, Set, Union

from pyopenms import IDFilter, IdXMLFile


def main(epifany_file: str, idpicker_file: str, pyprophet_file: str,
         threshold: str, remove_decoy: bool, return_qvalue: bool):
    num_epifany_distinct_protein = len(get_epifany_result(epifany_file,
                                                          threshold,
                                                          remove_decoy,
                                                          return_qvalue))

    num_idpicker_distinct_protein = get_idpicker_numbers(idpicker_file,
                                                         threshold)

    num_pyprophet_distinct_protein = get_pyprophet_result(pyprophet_file,
                                                          threshold)

    return num_epifany_distinct_protein, \
        num_idpicker_distinct_protein, \
        num_pyprophet_distinct_protein


def get_epifany_result(epifany_file: str, threshold: str, remove_decoy: bool, return_qvalue: bool) -> \
        Union[List[Any], Set[Any]]:
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

            # list of bytes
            accession_list_bytes = group.accessions

            accession_list_str = [accession_bytes.decode("utf-8") for
                                  accession_bytes in accession_list_bytes]

            accessions = ' '.join(accession_list_str)

            all_protein_groups.append(accessions)
            protein_group_probability.setdefault(accessions, []).append(
                group.probability)

    # _ = plt.hist(epifany_q_values, bins='auto')
    # plt.title("epifany qvalue")
    # plt.xlabel("Q value")
    # plt.ylabel("Number of proteins")
    # plt.show()

    # _ = plt.hist(epifany_pep, bins='auto')
    # plt.title("epifany pep")
    # plt.xlabel("Posterior Error Probability")
    # plt.ylabel("Number of proteins")
    # plt.show()


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
    else:
        return epifany_distinct_proteins_groups


def remove_decoy_protein_groups(prot_ids):
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

    # TODO maybe I should just use set default to group protein into groups
    #   and then venn diagram compare that to pyprophet protein

    # _ = plt.hist(idpicker_pep, bins='auto')
    # plt.title("Idpicker pep")
    # plt.xlabel("Posterior Error Probability")
    # plt.ylabel("Number of Proteins")
    # plt.show()
    #
    # _ = plt.hist(idpicker_q_values, bins='auto')
    # plt.title("Idpicker qvalue")
    # plt.xlabel("Q Value")
    # plt.ylabel("Number of Proteins")
    # plt.show()

    print("IDpicker", len(idpicker_q_values))

    # in_both = idpicker_distinct_protein.intersection(epifany_distinct_proteins)
    #
    # print("in both", len(in_both))
    #
    # only_in_epifany = epifany_distinct_proteins.difference(idpicker_distinct_protein)
    # only_in_idpicker = idpicker_distinct_protein.difference(epifany_distinct_proteins)
    #
    # print("num only in epifany", len(only_in_epifany))
    # print("num only in idpicker", len(only_in_idpicker))

    return len(idpicker_q_values)


def get_idpicker_accessions(idpicker_file: str, threshold: str):

    con = sqlite3.connect(idpicker_file)

    c = con.cursor()

    c.execute(
        """SELECT DISTINCT PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID
        FROM PROTEIN_GROUP
        INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
        WHERE PROTEIN_GROUP.DECOY = 0 AND SCORE_PROTEIN_GROUP.QVALUE <= :threshold""",
        {'threshold': float(threshold)}
    )

    protein_group_to_accession = {}

    for row in c.fetchall():
        protein_group_id = row[0]
        protein_accession = row[1]

        protein_group_to_accession.setdefault(protein_group_id, []).append(
            protein_accession)

    con.commit()
    con.close()

    all_accessions = []
    for protein_group, protein_accession_list in protein_group_to_accession.items():
        accessions = ' '.join(protein_accession_list)
        all_accessions.append(accessions)

    return all_accessions



def get_pyprophet_result(pyprophet_file: str, threshold: str) -> int:
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


# if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # epifany output file, idpicker output file, pyprophet file q-value threshold
    # whether to remove decoys
    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
