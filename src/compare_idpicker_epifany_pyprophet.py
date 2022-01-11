# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 10:07:46 2021

@author: kren
"""

import collections
import sqlite3
import sys

from pyopenms import IdXMLFile


def main(epifany_file: str, idpicker_file: str, threshold: str):
    # getting all protein from epifany
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    all_protein_hits = collections.deque()

    protein_hit_scores = {}

    # there is 1 protein identification
    # ~200000 hits
    # each hit is one protein accession
    # why do some hits have the same protein accession, but different q_value
    # that is because when I convert from OSW i forgot to remove duplciates

    score_type = "None"

    epifany_q_values = []
    epifany_pep = []

    for protein_id in prot_ids:

        score_type = protein_id.getScoreType()

        # print(protein_id.getScoreType())

        for hit in protein_id.getHits():

            if hit.getMetaValue("target_decoy") == "target":
                accession = hit.getAccession()
                all_protein_hits.append(accession)
                epifany_q_values.append(hit.getScore())
                # epifany_pep.append(1 - hit.getScore())
                protein_hit_scores.setdefault(accession, []).append(hit.getScore())

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
            if q_value <= float(threshold):
                epifany_proteins.append(accession)
    else:
        print("Unrecognized Score Type")

    epifany_distinct_proteins = set(epifany_proteins)
    print("epifany", len(epifany_distinct_proteins))




    """now for idpicker """
    con = sqlite3.connect(idpicker_file)

    c = con.cursor()

    #PEP
    # c.execute(
    #     """SELECT PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID, SCORE_PROTEIN_GROUP.PEP
    #     FROM PROTEIN_GROUP
    #     INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
    #     WHERE PROTEIN_GROUP.DECOY = 0 AND PEP <= :threshold""", {'threshold': float(threshold)}
    # )

    #q_value
    c.execute(
        """SELECT PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID, SCORE_PROTEIN_GROUP.QVALUE, SCORE_PROTEIN_GROUP.PEP
        FROM PROTEIN_GROUP
        INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
        WHERE PROTEIN_GROUP.DECOY = 0 AND QVALUE <= :threshold""", {'threshold': float(threshold)}
    )

    idpicker_protein = []
    idpicker_q_values = []
    idpicker_pep = []

    for row in c.fetchall():
        group_id = row[0]
        protein_id = row[1]
        q_value = row[2]
        pep = row[3]

        idpicker_protein.append(protein_id)
        idpicker_q_values.append(q_value)
        idpicker_pep.append(pep)

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



    idpicker_distinct_protein = set(idpicker_protein)
    print("IDpicker", len(idpicker_distinct_protein))

    in_both = idpicker_distinct_protein.intersection(epifany_distinct_proteins)

    # print("in both", len(in_both))

    only_in_epifany = epifany_distinct_proteins.difference(idpicker_distinct_protein)
    only_in_idpicker = idpicker_distinct_protein.difference(epifany_distinct_proteins)

    # print("num only in epifany", len(only_in_epifany))
    # print("num only in idpicker", len(only_in_idpicker))



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


    return len(idpicker_distinct_protein), len(epifany_distinct_proteins),\
           len(pyprophet_distinct_protein)


if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # print("""usage: osw_idXML_converter.py <sql file path> <q-value threshold for peptides>""")
    # main(sys.argv[1], sys.argv[2])
    # epifany output file, idpicker output file, q-value threshold
    # pyprophet output file is the same as the idpicker
    main(sys.argv[1], sys.argv[2], sys.argv[3])





