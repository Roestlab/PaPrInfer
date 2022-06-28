import sqlite3
import sys

import matplotlib.pyplot as plt

if __name__ == "__main__":

    decoy_qvalue_list = []

    db_name = sys.argv[1]

    con = sqlite3.connect(db_name)

    c = con.cursor()

    c.execute(
        """SELECT SCORE_PROTEIN_GROUP.QVALUE
        FROM SCORE_PROTEIN_GROUP""")

    for row in c.fetchall():
        qvalue = row[0]
        decoy_qvalue_list.append(qvalue)

    print(decoy_qvalue_list)

    _ = plt.hist(decoy_qvalue_list, bins=100, label="target")
    plt.title("AFTER idpicker qvalue distribution in histogram")
    plt.xlabel("qvalue")
    plt.ylabel("Number of Protein")
    plt.legend()
    plt.show()

