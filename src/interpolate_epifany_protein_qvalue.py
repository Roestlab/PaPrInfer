from scipy.interpolate import interp1d
import sys

from src.get_proteins_at_threshold import get_epifany_result


def main(epifany_file_pep, epifany_file_qvalue):
    protein_group_pep_dict = get_epifany_result(epifany_file_pep, threshold='1',
                                                remove_decoy=False, return_qvalue=False,
                                                return_distinct_protein=False, return_mapping=True)

    protein_group_qvalue_dict = get_epifany_result(epifany_file_qvalue, threshold='1',
                                                   remove_decoy=False, return_qvalue=False,
                                                   return_distinct_protein=False, return_mapping=True)

    # this is a list of list of list, each middle sublist has 3 elements,
    # a list of protein group's accessions
    # a list of posterior error probabilities
    # a list of qvalue
    protein_group_pep_qvalue_list = []
    for protein_group, pep in protein_group_pep_dict.items():
        qvalue = protein_group_qvalue_dict[protein_group]

        protein_group_pep_qvalue_list.append([protein_group, pep, qvalue])

    with open('figures and files/protein_group_pep_qvalue.txt', 'w') as f:
        for protein_group in protein_group_pep_qvalue_list:
            accessions = protein_group[0]
            pep = protein_group[1]
            qvalue = protein_group[2]

            f.write(accessions)
            f.write(pep)
            f.write(qvalue)
            f.write('\n')


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
