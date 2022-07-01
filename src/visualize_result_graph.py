import sys
from typing import Any, Dict, List, Tuple

from pyopenms import IdXMLFile
from get_proteins_at_threshold import remove_decoy_protein_groups
import sqlite3


"""this does 3 things, 
1. write all result of idpicker and epifany into a text file
so it is easy to text search it
2. (Not finished) compute the number of inferred protein per component for both epifany and idpicker
3. (not started) plot the entire graph into cytoscape"""


# import py2cytoscape
# from py2cytoscape.data.cyrest_client import CyRestClient


def get_graph_epifany(epifany_file) -> Tuple[
        Dict[str, List[Any]], Dict[Any, Any]]:
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    epifany_all_protein_peptide_dict = {}

    # use peptide evidence to get a mapping of all protein to its peptides
    for peptide_id in pep_ids:
        # PeptideHits
        for hit in peptide_id.getHits():
            if hit.getMetaValue("target_decoy") == "target":

                peptide_sequence = str(hit.getSequence())

                protein_list = [ev.getProteinAccession() for ev in
                                hit.getPeptideEvidences()]

                for protein in protein_list:
                    epifany_all_protein_peptide_dict.setdefault(protein,
                                                                []).append(
                        peptide_sequence)

    # the peptide_sequences are unique (in each list) since the list of
    # peptide_hits is unique

    remove_decoy_protein_groups(prot_ids)

    # select the subset where the protein is inferred
    protein_peptide_dict = {}

    # read the protein groups, after remove all decoy protein groups
    for protein_id in prot_ids:

        for group in protein_id.getProteinGroups():

            # list of bytes
            accession_list_bytes = group.accessions

            accession_list_str = [accession_bytes.decode("utf-8") for
                                  accession_bytes in accession_list_bytes]

            # all the peptide that map to this protein group
            all_mapped_peptide = []
            for accession in accession_list_str:
                all_mapped_peptide.extend(
                    epifany_all_protein_peptide_dict[accession])

            all_mapped_peptide = list(set(all_mapped_peptide))

            accessions = ' '.join(accession_list_str)

            aaccessions_fdr = accessions + ' ' + str(group.probability)

            protein_peptide_dict[aaccessions_fdr] = all_mapped_peptide

    peptide_protein_dict = {}
    for protein_accessions, peptide_list in protein_peptide_dict.items():
        protein_list = protein_accessions.split()

        for peptide in peptide_list:
            peptide_protein_dict.setdefault(peptide, []).extend(
                protein_list
            )

    for a, b in protein_peptide_dict.items():
        if len(b) == 0:
            print(a)

    return protein_peptide_dict, peptide_protein_dict


def get_graph_idpicker(idpicker_file) -> Tuple[Dict[str, Any], Dict[Any, Any]]:
    con = sqlite3.connect(idpicker_file)
    c = con.cursor()

    c.execute("""select PROTEIN_GROUP.PROTEIN_GROUP_ID, PROTEIN_GROUP.PROTEIN_ID
     from PROTEIN_GROUP
     where PROTEIN_GROUP.DECOY = 0""")

    protein_group_protein_dict = {}
    for row in c.fetchall():
        protein_group_id = row[0]
        protein_id = row[1]
        protein_group_protein_dict.setdefault(protein_group_id, []).append(
            protein_id)

    c.execute(
        """select distinct PROTEIN_GROUP_PEPTIDE_MAPPING.PROTEIN_GROUP_ID, PEPTIDE.UNMODIFIED_SEQUENCE
        from PROTEIN_GROUP_PEPTIDE_MAPPING
        inner join PEPTIDE on PROTEIN_GROUP_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID""")

    protein_group_peptide_dict = {}
    for row in c.fetchall():
        protein_group_id = row[0]
        peptide_sequence = row[1]
        protein_group_peptide_dict.setdefault(protein_group_id, []).append(
            peptide_sequence)

    # make a mapping from protein group to peptide
    protein_peptide_dict = {}
    for protein_group_id, proteins in protein_group_protein_dict.items():
        protein_accessions = ' '.join(proteins)

        protein_peptide_dict[protein_accessions] = \
            protein_group_peptide_dict[protein_group_id]

    # make a mapping from peptide to protein group
    peptide_protein_dict = {}
    for protein_accessions, peptide_list in protein_peptide_dict.items():

        for peptide in peptide_list:
            peptide_protein_dict.setdefault(peptide, []).append(
                protein_accessions
            )

    return protein_peptide_dict, peptide_protein_dict


def write_into_text_file(epifany_protein_peptide_dict,
                         idpicker_protein_peptide_dict):
    lines = []
    for protein, peptide in epifany_protein_peptide_dict.items():
        # print("epifany", protein)
        # print("epifany", peptide)
        peptide_sequence = ' '.join(peptide)

        lines.append("epifany")
        lines.append(protein)
        lines.append(peptide_sequence)

    for protein, peptide in idpicker_protein_peptide_dict.items():
        # print("idpicker", protein)
        # print("idpicker", peptide)
        peptide_sequence = ' '.join(peptide)

        lines.append("idpicker")
        lines.append(protein)
        lines.append(peptide_sequence)

    with open('figures and files/all_result.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')


def is_white(protein, discovered_protein, explored_protein):
    return (protein not in discovered_protein) \
           and (protein not in explored_protein)


def count_proteins_both(
        epifany_protein_peptide_dict, epifany_peptide_protein_dict,
        idpicker_peptide_protein_dict
):
    # each element is the number of protein inferred by a method in a component
    num_epifany_protein_component = []
    num_idpicker_protein_component = []
    discovered = {}
    explored = {}

    all_protein_list = []
    for protein in epifany_protein_peptide_dict:
        all_protein_list.append(protein)

    all_peptide_list = []
    for peptide in epifany_peptide_protein_dict:
        all_peptide_list.append(peptide)

    for protein, peptide_list in epifany_protein_peptide_dict.items():
        if is_white(protein, discovered, explored):
            # find all peptides and all protein inferred by epifany for
            # this component
            component_protein_list_epifany = []
            component_peptide_list = []
            for peptide in peptide_list:
                search_in_protein(peptide, component_protein_list_epifany,
                                  component_peptide_list,
                                  epifany_protein_peptide_dict,
                                  epifany_peptide_protein_dict, discovered,
                                  explored)

            # now find what is the idpicker protein of this component
            component_protein_list_idpicker = []
            for peptide in component_peptide_list:
                proteins = idpicker_peptide_protein_dict[peptide]
                component_protein_list_idpicker.append(proteins)

            num_epifany_protein_component.append(
                len(component_protein_list_epifany))
            num_idpicker_protein_component.append(
                len(component_protein_list_idpicker))

    return num_epifany_protein_component, num_idpicker_protein_component


def search_in_protein(a_node, component_protein_list, component_peptide_list,
                      protein_peptide_dict, peptide_protein_dict,
                      discovered, explored):
    discovered[a_node] = ''

    neighbour_list = []
    if a_node in protein_peptide_dict:  # is protein
        neighbour_list = protein_peptide_dict[a_node]
        component_protein_list.append(a_node)
    elif a_node in peptide_protein_dict:  # is peptide
        neighbour_list = peptide_protein_dict[a_node]
        component_peptide_list.append(a_node)

    # there are 1000 epifany proteins that has no peptides
    # is that right?
    if len(neighbour_list) == 0:
        print(a_node)

    # base case is that no neighbour are white, so then this loop
    # does not happen
    for neighbour in neighbour_list:
        if is_white(neighbour, discovered, explored):
            search_in_protein(neighbour, component_protein_list,
                              component_peptide_list,
                              protein_peptide_dict, peptide_protein_dict,
                              discovered, explored)

    explored[a_node] = ''


def analyse_comparison(num_epifany_protein_component,
                       num_idpicker_protein_component):
    assert len(num_epifany_protein_component) == len(
        num_idpicker_protein_component)

    how_many_more = []
    for i in range(len(num_epifany_protein_component)):
        how_many_more.append(
            num_epifany_protein_component[i] / num_idpicker_protein_component[i]
        )

    print("max", max(how_many_more))
    print("min", min(how_many_more))
    print("average", sum(how_many_more) / len(how_many_more))


if __name__ == "__main__":
    epifany_protein_peptide_dict, epifany_peptide_protein_dict = get_graph_epifany(
        sys.argv[1])
    idpicker_protein_peptide_dict, idpicker_peptide_protein_dict = get_graph_idpicker(
        sys.argv[2])
    write_into_text_file(epifany_protein_peptide_dict,
                         idpicker_protein_peptide_dict)

    # num_epifany_protein_per_component_list, num_idpicker_protein_per_component_list = \
    #     count_proteins_both(epifany_protein_peptide_dict,
    #                         epifany_peptide_protein_dict,
    #                         idpicker_peptide_protein_dict)

    # analyse_comparison(num_epifany_protein_per_component_list,
    #                    num_idpicker_protein_per_component_list)
