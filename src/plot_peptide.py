
import sqlite3, sys
import matplotlib.pyplot as plt

if __name__ == "__main__":

    decoy_pep_list = []
    target_pep_list = []

    db_name = sys.argv[1]

    con = sqlite3.connect(db_name)

    c = con.cursor()

    # note that there is no global context in the 20180911 file
    c.execute(
        """SELECT Score_Peptide.PEP
        FROM Peptide
        INNER JOIN Score_Peptide ON Score_Peptide.peptide_id = Peptide.id
        WHERE DECOY = 1
        AND SCORE_PEPTIDE.CONTEXT = 'global'""")

    for row in c.fetchall():
        posterior_error_probability = row[0]
        decoy_pep_list.append(posterior_error_probability)

    c.execute(
        """SELECT Score_Peptide.PEP
        FROM Peptide
        INNER JOIN Score_Peptide ON Score_Peptide.peptide_id = Peptide.id
        WHERE DECOY = 0
        AND SCORE_PEPTIDE.CONTEXT = 'global'""")

    for row in c.fetchall():
        posterior_error_probability = row[0]
        target_pep_list.append(posterior_error_probability)

    _ = plt.hist(decoy_pep_list, bins='auto', label="decoy", alpha=0.5)
    _ = plt.hist(target_pep_list, bins='auto', label="target", alpha=0.5)
    plt.title("Before epifany PEP distribution in histogram")
    plt.xlabel("PEP")
    plt.ylabel("Number of Peptides")
    plt.legend()
    plt.show()

