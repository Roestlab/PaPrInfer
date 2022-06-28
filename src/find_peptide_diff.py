import sqlite3
import sys

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

from collections import Counter


def main(swissprot_library, uniprot_library):

    swissprot_peptides = set()
    uniprot_peptides = set()

    get_all_peptide_sequence(swissprot_library, swissprot_peptides)
    get_all_peptide_sequence(uniprot_library, uniprot_peptides)

    # find the differences
    swissprot_only = swissprot_peptides.difference(uniprot_peptides)
    uniprot_only = uniprot_peptides.difference(swissprot_peptides)

    # analyse the differences
    library_unique_list = list(swissprot_only)

    library_unique_list.extend(list(uniprot_only))

    library_unique_list.sort()

    first_aa_list_swissprot = []

    first_aa_list_uniprot = []

    for pep in swissprot_only:
        first_aa_list_swissprot.append(pep[0])

    print('swissprot', len(first_aa_list_swissprot))
    print(Counter(first_aa_list_swissprot))

    for pep in uniprot_only:
        first_aa_list_uniprot.append(pep[0])

    print('uniprot', len(first_aa_list_uniprot))
    print(Counter(first_aa_list_uniprot))

    with open('figures and files/library unique peptides.txt', 'w') as f:
        for line in library_unique_list:
            f.write(line)
            f.write('\n')

    venn2([swissprot_peptides, uniprot_peptides], ('swissprot peptides', 'uniprot peptides'))
    plt.savefig('figures and files/peptide venn.pdf')


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
