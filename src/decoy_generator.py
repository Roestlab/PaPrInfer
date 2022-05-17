import sqlite3
import sys
from typing import Any, Dict, Tuple


def main(pqp_file: str) -> None:
    """
    make decoy proteins
    :param pqp_file: the transition pqp file
    :return: none
    """

    # a.	For every target protein, make a decoy protein entry,
    # b.	Find the decoy counterpart of peptides of this target protein
    # c.	That is the decoy peptides for this decoy protein
    # d.	And then it should work, since epifany only want peptide PEP
    # i.	    And protein name
    con = sqlite3.connect(pqp_file)
    c = con.cursor()
    delete_old_decoys(c)
    sequence_id, decoy_protein_peptide = get_mapping(c)
    add_decoy_protein_peptide(c, sequence_id, decoy_protein_peptide)
    con.commit()
    c.close()


def delete_old_decoys(c) -> None:

    # first delete all score protein entries that are decoy proteins
    c.execute(
        """DELETE FROM SCORE_PROTEIN
        WHERE SCORE_PROTEIN.PROTEIN_ID IN
        (SELECT PROTEIN_ID FROM SCORE_PROTEIN 
        INNER JOIN PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.ID
        WHERE PROTEIN.DECOY = 1)"""
    )

    # second delete all decoy proteins
    c.execute(
        """DELETE FROM PROTEIN
        WHERE PROTEIN.DECOY = 1"""
    )


def get_mapping(c) -> Tuple[Dict[Any, Any], Dict[Any, Any]]:
    # get a dictionary that maps peptide unmod sequence to peptide id
    # this is 1 to 1
    peptide_seq_id_dict = {}

    c.execute(
        """SELECT PEPTIDE.UNMODIFIED_SEQUENCE, PEPTIDE.ID
        FROM PEPTIDE
        WHERE PEPTIDE.DECOY = 1"""
    )

    for row in c.fetchall():
        peptide_seq = row[0]
        peptide_id = row[1]
        peptide_seq_id_dict[peptide_seq] = peptide_id

    # read all target proteins map to its peptides, as a dictionary
    # this is 1 to n
    protein_accession_peptide_seq_dict = {}

    c.execute(
        """SELECT DISTINCT PROTEIN.PROTEIN_ACCESSION, PEPTIDE.UNMODIFIED_SEQUENCE
        FROM PROTEIN 
        INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID
        INNER JOIN PEPTIDE ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID
        WHERE PROTEIN.DECOY = 0"""
    )

    for row in c.fetchall():
        protein_accession = row[0]
        peptide_seq = row[1]
        protein_accession_peptide_seq_dict.setdefault(protein_accession, []).append(peptide_seq)

    return peptide_seq_id_dict, protein_accession_peptide_seq_dict



def add_decoy_protein_peptide(c, peptide_seq_id_dict, protein_accession_peptide_seq_dict):

    # for loop over the dictionary
    # based on the dictionary keys, insert all decoy protein entries
    # I just need protein accession (just use decoy + target protein accession)
    # a new ID, and set decoy = 1

    # get the id that is numerically greatest
    c.execute(
        """SELECT max(id)
        FROM PROTEIN"""
    )

    max_id = int(c.fetchone()[0])

    current_id = 1
    for protein_accession, peptide_sequence_list in protein_accession_peptide_seq_dict.items():

        protein_id = max_id + current_id
        decoy_protein_accession = 'DECOY ' + protein_accession
        decoy = 1

        c.execute(
            """INSERT INTO PROTEIN(ID, PROTEIN_ACCESSION, DECOY) 
            VALUES(:id, :protein_accession, :decoy)""",
            {'id': protein_id, 'protein_accession': decoy_protein_accession, 'decoy': decoy}
        )

        # then add in score protein
        # do global context, null run id, protein id, 0 score, 0 pvalue, 0 qvalue
        # 0 pep
        context = 'global'
        run_id = 0
        score = 0
        pvalue = 0
        qvalue = 0
        pep = 0

        c.execute(
            """INSERT INTO SCORE_PROTEIN(CONTEXT, RUN_ID, PROTEIN_ID, SCORE,
            PVALUE, QVALUE, PEP) VALUES(:context, :run_id, :protein_id, :score,
            :pvalue, :qvalue, :pep)""",
            {'context': context, 'run_id': run_id, "protein_id": protein_id,
             'score': score, 'pvalue': pvalue, 'qvalue': qvalue, 'pep': pep}
        )

        # at the same time, take the peptide sequence, and pseudo reverse it
        for peptide_sequence in peptide_sequence_list:

            pseudo_reverse(c, peptide_seq_id_dict, peptide_sequence, protein_id)

        current_id += 1

        # Also, I think I remove the issue of multiple peptide entries have
        # thew same seqeucen but different id

    # I dont understand the notation of unimod
    # but I think that if we have M(UniMod:35)S
    # the modification is on M
    # and pseudo reverse does not change which amino acid the modification
    # also my protein inference does not consider modification, so
    # I think we can make decoy peptides without them for now


    pass


def pseudo_reverse(c, peptide_seq_id_dict, peptide_sequence, protein_id):

    # ignore the last letter, reverse the rest, convert k/r to r/k

    # with the pseudo reserve sequence
    # first check if it is in the dictionary of sequence to id
    # get its peptide id
    # write mapping entires

    last_aa = peptide_sequence[len(peptide_sequence) - 1]
    except_last = peptide_sequence[0:-1]
    reversed_except_last = except_last[::-1]
    if last_aa == 'K':
        pseudo_reverse_seq = reversed_except_last + 'R'
    elif last_aa == 'R':
        pseudo_reverse_seq = reversed_except_last + 'K'
    else:
        print("Not sure if this is right")
        pseudo_reverse_seq = reversed_except_last
    if pseudo_reverse_seq in peptide_seq_id_dict:
        decoy_peptide_id = peptide_seq_id_dict[pseudo_reverse_seq]

        c.execute(
            """INSERT INTO PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID, PROTEIN_ID) 
            VALUES(:peptide_id, :protein_id)""",
            {'peptide_id': decoy_peptide_id, 'protein_id': protein_id}
        )


if __name__ == "__main__":
    main(sys.argv[1])
