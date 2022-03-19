# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 10:07:46 2021

@author: kren
"""

import collections
import sqlite3
import sys

from pyopenms import IDFilter, IdXMLFile


def main(epifany_file: str, idpicker_file: str, pyprophet_file: str,
         threshold: str):
    num_epifany_distinct_protein = get_epifany_result(epifany_file, threshold)

    num_idpicker_distinct_protein = get_idpicker_result(idpicker_file,
                                                        threshold)

    num_pyprophet_distinct_protein = get_pyprophet_result(pyprophet_file,
                                                          threshold)

    return num_epifany_distinct_protein, \
           num_idpicker_distinct_protein, \
           num_pyprophet_distinct_protein


def get_epifany_result(epifany_file: str, threshold: str) -> int:
    # getting all protein from epifany
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    all_protein_hits = collections.deque()

    all_protein_groups = []

    protein_hit_scores = {}
    protein_group_probability = {}

    # there is 1 protein identification
    # ~200000 hits
    # each hit is one protein accession
    # why do some hits have the same protein accession, but different q_value
    # that is because when I convert from OSW i forgot to remove duplciates

    # proteins, should be protein identification list load from idXML
    # protrun, is then just the a protein identification
    # I need getIndistinguishableProteins() instead of getIndistinguishableProteinGroups()

    # IDFilter().removeDecoyHits(proteins)
    # protrun = proteins[0]
    # grps = protrun.getIndistinguishableProteinGroups()
    # IDFilter().updateIndistinguishableProteinGroups(grps,
    #                                                 protrun.getHits())
    # protrun.setIndistinguishableProteinGroups(grps)
    # proteins[0] = protrun

    score_type = "None"

    epifany_q_values = []

    # need to check if this filter in place (i.e. mutates, and not return a copy)
    IDFilter().removeDecoyHits(prot_ids)

    for protein_id in prot_ids:
        groups = protein_id.getIndistinguishableProteins()
        hits = protein_id.getHits()
        IDFilter().updateProteinGroups(groups, hits)

        # in epifany output protein groups is empty anyway
        # check that
        assert protein_id.getProteinGroups() == []

        # append only target protein groups into "protein groups"
        for group in groups:
            protein_id.insertProteinGroup(group)

    for protein_id in prot_ids:

        for group in protein_id.getProteinGroups():

            # list of bytes
            accession_list_bytes = group.accessions

            accession_list_str = [accession_bytes.decode("utf-8") for
                                  accession_bytes in accession_list_bytes]

            accessions = ''.join(accession_list_str)

            all_protein_groups.append(accessions)
            protein_group_probability.setdefault(accessions, []).append(
                group.probability)

        # for hit in protein_id.getHits():
        #
        #     if hit.getMetaValue("target_decoy") == "target":
        #         accession = hit.getAccession()
        #         all_protein_hits.append(accession)
        #         epifany_q_values.append(hit.getScore())
        #         # epifany_pep.append(1 - hit.getScore())
        #         protein_hit_scores.setdefault(accession, []).append(
        #             hit.getScore())

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

    epifany_proteins = []

    if score_type == "Posterior Probability":
        for accession in protein_hit_scores:
            pp = min(protein_hit_scores[accession])
            pep = 1 - pp
            if pep <= float(threshold):
                epifany_proteins.append(accession)

    elif score_type == "q-value":
        for accession in protein_hit_scores:
            q_value = min(protein_hit_scores[accession])

            # only if it is less than or equal to the threshold, we reject the null
            # hypothesis i.e. accept it as inferred protein
            if q_value <= float(threshold):
                epifany_proteins.append(accession)
    else:
        print("Unrecognized Score Type")

    # epifany_distinct_proteins = set(epifany_proteins)
    # print("epifany", len(epifany_distinct_proteins))

    # return len(epifany_distinct_proteins)

    """protein groups"""
    epifany_proteins_groups = []

    for groups in protein_group_probability:
        probability = min(protein_group_probability[groups])

        if probability <= float(threshold):
            epifany_proteins_groups.append(groups)

    epifany_distinct_proteins_groups = set(epifany_proteins_groups)
    print("epifany groups", len(epifany_distinct_proteins_groups))

    return len(epifany_distinct_proteins_groups)


def get_idpicker_result(idpicker_file: str, threshold: str) -> int:
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
        group_id = row[0]
        # protein_id = row[1]
        q_value = row[1]

        idpicker_q_values.append(q_value)

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

    pyprophet_distinct_protein = set(pyprophet_proteins)
    print("pyprophet", len(pyprophet_distinct_protein))

    return len(pyprophet_distinct_protein)


if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # print("""usage: osw_idXML_converter.py <sql file path> <q-value threshold for peptides>""")
    # main(sys.argv[1], sys.argv[2])
    # epifany output file, idpicker output file, pyprophet file q-value threshold
    # pyprophet output file is the same as the idpicker
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
