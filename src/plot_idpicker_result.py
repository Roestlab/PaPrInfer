import sqlite3
import sys

import matplotlib.pyplot as plt

if __name__ == "__main__":

    decoy_pep_list = []

    db_name = sys.argv[1]

    con = sqlite3.connect(db_name)

    c = con.cursor()

    c.execute(
        """SELECT SCORE_PROTEIN_GROUP.QVALUE
        FROM SCORE_PROTEIN_GROUP""")

    for row in c.fetchall():
        posterior_error_probability = row[0]
        decoy_pep_list.append(posterior_error_probability)

    print(decoy_pep_list)

    _ = plt.hist(decoy_pep_list, bins='auto', label="target")
    plt.title("AFTER idpicker qvalue distribution in histogram")
    plt.xlabel("qvalue")
    plt.ylabel("Number of Protein")
    plt.legend()
    plt.show()

