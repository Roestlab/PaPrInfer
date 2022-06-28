import sqlite3
from collections import Counter
import sys


"""READ ME VERY IMPORTANT
this was a complete waste of time, since I wanted to compare the pep of PROTEIN
not peptides, between the more decoy and less decoy, not the PEPTIDES, since
they are the same"""


def main(more_decoy_file, less_decoy_file):
    less_decoy_pep_list = get_pep_from_sqlite_db(less_decoy_file)
    more_decoy_pep_list = get_pep_from_sqlite_db(more_decoy_file)

    print(len(less_decoy_pep_list))
    print(len(more_decoy_pep_list))

    # less = [i for i in less_decoy_pep_list if i <= 0.001]
    # more = [i for i in more_decoy_pep_list if i <= 0.001]
    #
    # print('less decoy')
    # print(Counter(less))
    # print('more decoy')
    # print(Counter(more))


def get_pep_from_sqlite_db(idpicker_file):

    con = sqlite3.connect(idpicker_file)

    c = con.cursor()

    c.execute(
        """select PEP
        from SCORE_PEPTIDE
        inner join PEPTIDE on PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
        where DECOY = 1"""
    )

    pep_list = []
    for row in c.fetchall():
        pep = row[0]
        pep_list.append(pep)

    con.commit()
    con.close()

    return pep_list

if __name__ == "__main__":
    # the epifany output file and a threshold
    main(sys.argv[1], sys.argv[2])

