from scipy.interpolate import interp1d
import sys

from get_proteins_at_threshold import get_epifany_result





def main(epifany_file_pep, epifany_file_qvalue):

    # first get the accession_pep_qvalue_list of all, target only, decoy only
    all_accession_pp_qvalue_list = get_accession_qvalue_pep(
        epifany_file_pep, epifany_file_qvalue, remove_decoy=False)

    target_accession_pp_qvalue_list = get_accession_qvalue_pep(
        epifany_file_pep, epifany_file_qvalue, remove_decoy=True)

    # then use the 2nd element of each protein group, which is pp, to sort
    all_accession_pp_qvalue_list.sort(key=lambda x: x[1], reverse=True)
    target_accession_pp_qvalue_list.sort(key=lambda x: x[1], reverse=True)

    write_into_txt_file(all_accession_pp_qvalue_list)

    sorted_target_accession_list, sorted_target_pp_list, sorted_target_qvalue_list = split_into_three_list(
        target_accession_pp_qvalue_list)

    # find the lowest pep for each qvalue, that is what need to use here for interpolate
    pep_lowest_pep, qvalue_lowest_pep = find_lowest_pep_each_qvalue(target_accession_pp_qvalue_list)

    interpolation_function = interp1d(pep_lowest_pep, qvalue_lowest_pep)

    sorted_target_pep_list = convert_to_pp(sorted_target_pp_list)

    new_qvalue = interpolation_function(sorted_target_pep_list)

    # what to do with the newqvalue

    write_into_txt_file_2(sorted_target_accession_list, sorted_target_pep_list, new_qvalue)


def find_lowest_pep_each_qvalue(target_accession_pp_qvalue_list):
    lowest_pep_list = []

    # make a mapping from each qvalue to a pep

    qvalue_pep = {}
    for protein_group in target_accession_pp_qvalue_list:
        pp = protein_group[1]
        qvalue = protein_group[2]
        pep = 1 - pp

        # only replace when the pep is lower
        # i think this could be written with set_default?
        if qvalue not in qvalue_pep:
            qvalue_pep[qvalue] = pep
        elif qvalue in qvalue_pep and pep < qvalue_pep[qvalue]:
            qvalue_pep[qvalue] = pep

    for qvalue, pep in qvalue_pep.items():
        lowest_pep_list.append([pep, qvalue])

    lowest_pep_list.sort()

    # then split into 2 list
    pep_lowest_pep, qvalue_lowest_pep = [], []

    for pep_qvalue in lowest_pep_list:
        pep_lowest_pep.append(pep_qvalue[0])
        qvalue_lowest_pep.append(pep_qvalue[1])

    # # this is quick and bad solution
    # lowest_pp_protein_group = target_accession_pp_qvalue_list[-1]
    # qvalue_highest_pep = lowest_pp_protein_group[1]
    # highest_pep = 1 - lowest_pp_protein_group[2]

    # pep_lowest_pep.append(highest_pep)
    # qvalue_lowest_pep.append(qvalue_highest_pep)

    return pep_lowest_pep, qvalue_lowest_pep


def split_into_three_list(accession_pep_qvalue_list):
    sorted_accession_list = []
    sorted_pep_list = []
    sorted_qvalue_list = []
    for protein_group in accession_pep_qvalue_list:
        accessions_list = protein_group[0]
        pep = protein_group[1]
        qvalue = protein_group[2]
        sorted_accession_list.append(accessions_list)
        sorted_pep_list.append(pep)
        sorted_qvalue_list.append(qvalue)
    return sorted_accession_list, sorted_pep_list, sorted_qvalue_list,


def get_accession_qvalue_pep(epifany_file_pep, epifany_file_qvalue,
                             remove_decoy):

    protein_group_pep_dict = get_epifany_result(epifany_file_pep, threshold='1',
                                                remove_decoy=remove_decoy,
                                                return_qvalue=False,
                                                return_distinct_protein=False,
                                                return_mapping=True)

    protein_group_qvalue_dict = get_epifany_result(epifany_file_qvalue,
                                                   threshold='1',
                                                   remove_decoy=remove_decoy,
                                                   return_qvalue=False,
                                                   return_distinct_protein=False,
                                                   return_mapping=True)

    # this is a list of list of list, each middle sublist has 3 elements,
    # a list of protein group's accessions
    # a list of posterior error probabilities
    # a list of qvalue (but both pep and qvalue really is a single element list)
    accession_pep_qvalue_list = []
    # a protein group only has 1 pep and 1 qvalue anyway,
    # so convert list into float
    for protein_group_accessions, pep_list in protein_group_pep_dict.items():

        qvalue = protein_group_qvalue_dict[protein_group_accessions][0]
        pep = pep_list[0]
        accession_pep_qvalue_list.append([protein_group_accessions, pep, qvalue])

    return accession_pep_qvalue_list


def write_into_txt_file(accession_pep_qvalue_list):
    with open('figures and files/protein_group_pep_qvalue.txt', 'w') as f:
        for protein_group in accession_pep_qvalue_list:
            accessions_list = protein_group[0]
            pep = protein_group[1]
            qvalue = protein_group[2]

            for accessions in accessions_list:
                f.write(str(accessions))
            f.write('\n')

            f.write('pp: ')
            f.write(str(pep))
            f.write('\n')

            f.write('qvalue: ')
            f.write(str(qvalue))
            f.write('\n')
            f.write('\n')

def convert_to_pp(pp_list):
    pep_list = []
    for pp in pp_list:
        pep_list.append(1 - pp)

    return pep_list


def write_into_txt_file_2(sorted_target_accession_list, sorted_target_pp_list, new_qvalue):
    with open('figures and files/protein_group_pep_new_qvalue.txt', 'w') as f:
        for i in range(len(sorted_target_accession_list)):
            accessions_list = sorted_target_accession_list[i]
            pep = sorted_target_pp_list[i]
            qvalue = new_qvalue[i]

            for accessions in accessions_list:
                f.write(str(accessions))
            f.write('\n')

            f.write('pep: ')
            f.write(str(pep))
            f.write('\n')

            f.write('qvalue: ')
            f.write(str(qvalue))
            f.write('\n')
            f.write('\n')

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
