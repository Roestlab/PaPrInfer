import sys

from Bio import SeqIO


def main(swissprot_fasta_file, uniprot_fasta_file):
    swissprot_all_accession_list = list()
    uniprot_all_accession_list = list()

    get_all_accession(swissprot_all_accession_list, swissprot_fasta_file)
    get_all_accession(uniprot_all_accession_list, uniprot_fasta_file)

    swissprot_all_accession_set = set(swissprot_all_accession_list)
    uniprot_all_accession_set = set(uniprot_all_accession_list)

    only_swissprot = swissprot_all_accession_set.difference(uniprot_all_accession_set)

    print(only_swissprot)

    with open('figures and files/only_in_swissprot_protein_accessions.txt', 'w') as f:
        for line in list(only_swissprot):
            f.write(line)
            f.write('\n')

    return only_swissprot


def get_all_accession(all_accession, fasta_file):
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        all_accession.append(name)

    return all_accession


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
