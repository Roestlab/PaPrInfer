import sqlite3
import sys

from matplotlib import pyplot as plt
from matplotlib_venn import venn2


def main(swissprot_library, uniprot_library):

    swissprot_peptides = set()
    uniprot_peptides = set()

    get_all_peptide_sequence(swissprot_library, swissprot_peptides)
    get_all_peptide_sequence(uniprot_library, uniprot_peptides)

    venn2([swissprot_peptides, uniprot_peptides], ('swissprot peptides', 'uniprot peptides'))
    plt.savefig('peptide venn.pdf')
    plt.show()


def get_all_peptide_sequence(library, peptide_set):
    con = sqlite3.connect(library)
    c = con.cursor()
    c.execute(
        """select UNMODIFIED_SEQUENCE 
        from PEPTIDE""")
    for row in c.fetchall():
        peptide_set.add(row[0])
    con.commit()
    con.close()


if __name__ == "__main__":
    # swissprot and uniprot osw file
    main(sys.argv[1], sys.argv[2])
