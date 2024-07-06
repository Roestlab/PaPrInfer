

# write one list into text file
def write_into_txt_file(accession_pep_qvalue_list, text_file_name):
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


# write three list into text file
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
