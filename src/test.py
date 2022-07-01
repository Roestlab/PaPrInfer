import sys

import get_proteins_at_threshold


if __name__ == "__main__":
    uniprot_epifany_file = sys.argv[1]
    uniprot_idpicker_file = sys.argv[2]
    threshold = sys.argv[3]

    uniprot_epifany_proteins_accession = \
        get_proteins_at_threshold.get_epifany_result(uniprot_epifany_file,
                                                     remove_decoy=True,
                                                     return_qvalue=False,
                                                     threshold=threshold)

    with open('figures and files/epifany_accessions.txt', 'w') as f:
        for line in list(uniprot_epifany_proteins_accession):
            f.write(line)
            f.write('\n')

    uniprot_idpicker_proteins_accession = \
        get_proteins_at_threshold.get_idpicker_accessions(uniprot_idpicker_file,
                                                          threshold=threshold)

    with open('figures and files/idpicker_accessions.txt', 'w') as f:
        for line in list(uniprot_idpicker_proteins_accession):
            f.write(line)
            f.write('\n')

