from __future__ import annotations

from typing import Any, Dict, List, Tuple

import sqlite3
import sys
# import igraph
from graph import Graph
from components import Component
from node import Peptide, Protein

"""main"""


def main(input_file: str, q_limit_pep: str, context: str, run_id: str) -> None:
    assert context == 'global' or context == 'run-specific'

    con = begin_connection(input_file)

    run_id = int(run_id)

    q_limit_pep = int(q_limit_pep)

    protein_peptide_graph = Graph()

    initialize(protein_peptide_graph, con, q_limit_pep, context, run_id)

    collapse(protein_peptide_graph)

    # visualize(protein_peptide_graph)

    component_list = separate(protein_peptide_graph)

    reduce(component_list, con)

    end_connection(con)


"""the helper function for main"""


def initialize(protein_peptide_graph: Graph, con,
               q_limit: int, context: str, run_id: int) -> None:
    peptide_id_list = get_all_peptide(con, q_limit, context, run_id)
    print("got all peptide")

    protein_id_list = get_all_protein_id(con, context, run_id)
    print("got all protein")

    linked_peptide_dict = get_all_link_for_protein(con, context, q_limit,
                                                   run_id)
    print("got all link from protein to peptide")

    protein_accession_dict = get_all_protein_accession(con)
    print("got all accession per protein sqlite id")

    linked_protein_dict = get_all_link_for_peptide(con)
    print("got all link from peptide to protein")

    protein_peptide_graph.add_peptide(peptide_id_list)
    print("added all peptide")

    protein_peptide_graph.add_protein(protein_id_list, linked_peptide_dict,
                                      protein_accession_dict)
    print("added all protein with edges")

    protein_peptide_graph.make_edges_from_peptide(linked_protein_dict,
                                                  protein_accession_dict)
    print("added edges from peptide")

    protein_peptide_graph.get_sort_keys()
    print("keys sorted by their id")

    # TODO test
    print("after just initialization")
    for key, value_list in protein_peptide_graph.node_dict.items():
        if key.get_first_id() == 'sp|Q96Q89-3|KI20B_HUMAN':
            value_accession_list = []
            for value in value_list:
                accession = value.get_first_id()
                value_accession_list.append(accession)
            print(value_accession_list)

    print("initialized")

    num_t_protein = 0
    num_d_peptide = 0
    num_d_protein = 0
    num_t_peptide = 0
    num_edges = 0
    for node in protein_peptide_graph.node_dict:
        num_edges += len(protein_peptide_graph.node_dict[node])
        if isinstance(node, Protein) and node.get_target_decoy() == '0':
            num_t_protein += 1
        elif isinstance(node, Peptide) and node.get_target_decoy() == '0':
            num_t_peptide += 1
        elif isinstance(node, Protein) and node.get_target_decoy() == '1':
            num_d_protein += 1
        elif isinstance(node, Peptide) and node.get_target_decoy() == '1':
            num_d_peptide += 1

    print("total nodes", len(protein_peptide_graph.get_node_dict()))
    print("total edges", num_edges)
    print("target protein", num_t_protein)
    print("target peptide", num_t_peptide)
    print("decoy protein", num_d_protein)
    print("decoy peptide", num_d_peptide)


def collapse(protein_peptide_graph: Graph) -> None:
    # protein_peptide_graph.collapse_graph_old()
    protein_peptide_graph.collapse_graph()


# def visualize(protein_peptide_graph):
#     pass
#     my_graph = protein_peptide_graph
#
#     their_graph = igraph.Graph()
#
#     # my_graph is a object with 1 dict
#     # one for protein, one for peptide
#     # the protein dict has protein as key
#     # and peptide as values
#
#     my_node_dict = my_graph.get_node_dict()
#
#     for vertex in my_node_dict:
#         if vertex.get_first_id() not in my_graph.node_to_delete:
#             my_string = ', '.join(map(str, vertex.get_id()))
#             their_graph.add_vertex(my_string)
#
#     print("vertex done")
#
#     protein_number = 1
#     peptide_number = 1
#     for key, value in my_node_dict.items():
#         if (key.get_first_id() not in my_graph.node_to_delete) \
#                 and isinstance(key, Protein):
#             vertex_string_1 = ', '.join(map(str, key.get_id()))
#             print("protein number", protein_number, "got protein vertex")
#             for element in value:
#                 vertex_string_2 = ', '.join(map(str, element.get_id()))
#                 their_graph.add_edge(vertex_string_1, vertex_string_2)
#                 peptide_number += 1
#             protein_number += 1
#
#     print("edges done")
#
#     layout = their_graph.layout_auto()
#     print("layout done")
#     plot = igraph.plot(their_graph, vertex_label=their_graph.vs["name"],
#                        layout=layout)
#     print("plot done")
#     plot.show()


def separate(protein_peptide_graph: Graph) -> List[Component]:
    # for all white peptide nodes, explore them
    protein_peptide_graph.make_accession_object_dict()

    if 'sp|Q9UI08-5|EVL_HUMAN' + '0' in protein_peptide_graph.node_to_delete:
        print('-5 it is deleted')
    else:
        print('-5 is not deleted')

    component_list = []

    for node in protein_peptide_graph.node_dict.keys():

        # if it is white and not deleted
        if protein_peptide_graph.is_white(node) and \
                (node.get_first_id() + node.get_target_decoy()) \
                not in protein_peptide_graph.node_to_delete:

            a_component = Component()
            protein_peptide_graph.dfs(node, a_component)
            component_list.append(a_component)

    print("separated")
    return component_list


# def visualize_component(component):
#     """
#     all nodes in a component are assumed to be not deleted
#     """
#     pass
#
#     their_graph = igraph.Graph()
#
#     # my_graph is a object with 1 dict
#     # one for protein, one for peptide
#     # the protein dict has protein as key
#     # and peptide as values
#
#     for vertex in component._protein_dict:
#         my_string = ', '.join(map(str, vertex.get_id()))
#         their_graph.add_vertex(my_string, type=False)
#     print("protein done")
#
#     for vertex in component._peptide_dict:
#         my_string = ', '.join(map(str, vertex.get_id()))
#         their_graph.add_vertex(my_string, type=True)
#     print("protein done")
#
#     protein_number = 1
#     peptide_number = 1
#     for key, value in component._protein_dict.items():
#         vertex_string_1 = ', '.join(map(str, key.get_id()))
#
#         for element in value:
#
#             vertex_string_2 = ', '.join(map(str, element.get_id()))
#             their_graph.add_edge(vertex_string_1, vertex_string_2)
#             peptide_number += 1
#         print("protein number", protein_number, "got all edges")
#         protein_number += 1
#
#     print("edges done")
#
#     layout = their_graph.layout_bipartite()
#     print("layout done")
#     plot = igraph.plot(their_graph, vertex_label=their_graph.vs["name"],
#                        layout=layout)
#     print("plot done")


def reduce(component_list: List[Component], con) -> None:
    min_pro_list = []

    for component in component_list:

        # this is a list of list of protein accession
        # sublist is one meta-protein vertex, contain multiple protein accession
        component_accession_list = component.make_protein_list()

        if len(component_accession_list) == 0:
            continue

        # this just a list of list of list of protein accession
        # sublist of it all belongs to the same component
        # sublist of the sublist belong to the same meta-protein vertex
        min_pro_list.append(component_accession_list)

    create_table_protein_group(con)
    create_table_protein_group_peptide_mapping(con)
    protein_group_data_entry(con, min_pro_list)

    print("reduced")


"""function that uses the Sqlite file directly"""


def begin_connection(db_name: str):
    con = sqlite3.connect(db_name)
    return con


def get_all_protein_id(con, context: str, run_id: int) -> List[Tuple[str, int]]:
    """
    choose all protein from the context
    :param run_id: the run id if the context is 'run-specific'
    :param context: the context, either global or 'run-specific'
    :param con: the connection to the sqlite database
    :return: a list of all proteins
    """
    c = con.cursor()
    # TODO: there must be a better solution
    if context == 'global':
        c.execute(
            """SELECT PROTEIN_ID, DECOY
            FROM SCORE_PROTEIN 
            INNER JOIN PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.ID
            WHERE CONTEXT='global'""")
    elif context == 'run-specific':
        c.execute(
            """SELECT PROTEIN_ID, DECOY
            FROM SCORE_PROTEIN 
            INNER JOIN PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.ID
            WHERE CONTEXT='run-specific' AND RUN_ID=:run_id""",
            {"run_id": run_id})
    all_protein_id_list = []
    for row in c.fetchall():
        all_protein_id_list.append((str(row[0]), row[1]))
    c.close()
    return all_protein_id_list


def get_all_peptide(con, q_limit: int, context: str, run_id: int) \
        -> List[Tuple[str, float, int]]:
    c = con.cursor()
    # choose all peptide from the global context and less than the threshold
    # TODO: there should be a more clever solution
    if context == 'global':
        c.execute(
            """SELECT PEPTIDE_ID, SCORE, DECOY
            FROM SCORE_PEPTIDE
            INNER JOIN PEPTIDE ON PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
            WHERE QVALUE<:q_limit AND SCORE_PEPTIDE.CONTEXT = 'global'""",
            {'q_limit': q_limit}
        )
    elif context == 'run-specific':
        c.execute(
            """SELECT PEPTIDE_ID, SCORE, DECOY
            FROM SCORE_PEPTIDE
            INNER JOIN PEPTIDE ON PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
            WHERE QVALUE<:q_limit AND SCORE_PEPTIDE.CONTEXT = 'run-specific' 
            AND RUN_ID=:run_id""",
            {'q_limit': q_limit, 'run_id': run_id}
        )
    all_peptide_id_list = []
    for row in c.fetchall():
        # each row is a tuple, since c.fetchall() returns a list of tuples
        all_peptide_id_list.append((str(row[0]), row[1], row[2]))

    c.close()
    return all_peptide_id_list

    # if I ever want to select only peptide or protein that fit a certain
    #  properties, use c.execute("SELECT * WHERE value=3 AND keyword=1")
    #  something like that, note this data is in SCORE_PROTEIN or SCORE_PEPTIDE


def get_all_link_for_protein(con, context, q_limit, run_id) -> Dict[
    str, List[Tuple[str, float, int]]]:
    c = con.cursor()
    if context == 'global':
        c.execute(
            """SELECT PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID,
            PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID, SCORE_PEPTIDE.SCORE, 
            PEPTIDE.DECOY
            FROM PEPTIDE_PROTEIN_MAPPING
            INNER JOIN SCORE_PEPTIDE 
            ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
            INNER JOIN PEPTIDE ON PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
            WHERE QVALUE<=:q_limit AND SCORE_PEPTIDE.CONTEXT = 'global'""",
            {'q_limit': q_limit})
    elif context == 'run-specific':
        c.execute(
            """SELECT PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID,
            PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID, SCORE_PEPTIDE.SCORE, 
            PEPTIDE.DECOY
            FROM PEPTIDE_PROTEIN_MAPPING
            INNER JOIN SCORE_PEPTIDE 
            ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
            INNER JOIN PEPTIDE ON PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
            WHERE QVALUE<=:q_limit AND SCORE_PEPTIDE.CONTEXT = 'run-specific' 
            AND RUN_ID=:run_id""",
            {'q_limit': q_limit, 'run_id': run_id})
    linked_peptide_dict = {}
    for row in c.fetchall():
        protein_id = str(row[0])
        peptide_id = str(row[1])
        peptide_score = row[2]
        peptide_decoy = row[3]
        peptide_info = (peptide_id, peptide_score, peptide_decoy)
        # if key was not in the dict, setdefault return the default value,
        # empty list here. if it was then it returns the value
        linked_peptide_dict.setdefault(protein_id, []).append(peptide_info)
    c.close()
    return linked_peptide_dict


def get_all_protein_accession(con) -> Dict[str, List[str]]:
    """
    returns a dictionary that pairs id with protein accession (as a list)
    """
    c = con.cursor()
    protein_accession_dict = {}
    c.execute(
        """SELECT ID, PROTEIN_ACCESSION 
            FROM PROTEIN"""
    )
    # if there are any row, split the row (which are text) into List[str]
    for row in c.fetchall():
        protein_sqlite_id = str(row[0])
        accession_sublist = row[1].split(",")

        stripped_accession_sublist = [s.strip() for s in accession_sublist]
        # if key was not in the dict, setdefault return the default value,
        # empty list here. if it was then it returns the value
        protein_accession_dict. \
            setdefault(protein_sqlite_id, []).extend(stripped_accession_sublist)

    c.close()

    for accession_list in protein_accession_dict.keys():
        protein_accession_dict[accession_list] \
            = list(set(protein_accession_dict[accession_list]))
    return protein_accession_dict


def get_all_link_for_peptide(con) -> Dict[str, List[Tuple[str, int]]]:
    """
    return a dict that serves as a mapping from peptide to protein
    """
    c = con.cursor()
    c.execute(
        """SELECT PEPTIDE_ID, PROTEIN_ID, PROTEIN.DECOY
        FROM PEPTIDE_PROTEIN_MAPPING
        INNER JOIN PROTEIN ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID"""
    )
    linked_protein_dict = {}
    for row in c.fetchall():
        peptide_id = str(row[0])
        protein_id = str(row[1])
        protein_decoy = row[2]
        protein_info = (protein_id, protein_decoy)
        # if key was not in the dict, setdefault return the default value,
        # empty list here. if it was then it returns the value
        linked_protein_dict.setdefault(peptide_id, []).append(protein_info)
    c.close()
    return linked_protein_dict


def create_table_protein_group(con) -> None:
    """mpl = minimal protein list"""
    c = con.cursor()
    c.execute("""DROP TABLE IF EXISTS PROTEIN_GROUP""")
    c.execute("""CREATE TABLE PROTEIN_GROUP(
              COMPONENT_ID INTEGER,
              PROTEIN_GROUP_ID INTEGER, 
              PROTEIN_ID INTEGER,
              SCORE REAL,
              DECOY INT);""")
    con.commit()
    c.close()


def create_table_protein_group_peptide_mapping(con) -> None:
    c = con.cursor()
    c.execute("""DROP TABLE IF EXISTS PROTEIN_GROUP_PEPTIDE_MAPPING""")
    c.execute("""CREATE TABLE PROTEIN_GROUP_PEPTIDE_MAPPING(
              PROTEIN_GROUP_ID INTEGER, 
              PEPTIDE_ID INTEGER
              );""")
    con.commit()
    c.close()


def protein_group_data_entry(con, min_pro_list: List[List[Tuple[
    List[str], float, str, Any]]]) -> None:
    c = con.cursor()
    component_id = 0
    protein_group_id = 0
    for component in min_pro_list:
        for protein_group in component:
            protein_accession_list = protein_group[0]
            score = protein_group[1]
            decoy = protein_group[2]
            peptide_id_list = protein_group[3]
            for protein_accession in protein_accession_list:
                c.execute("""INSERT INTO PROTEIN_GROUP(COMPONENT_ID, 
                PROTEIN_GROUP_ID, PROTEIN_ID, SCORE, DECOY) VALUES(:component_id
                , :protein_group_id, :protein_accession, :score, :decoy)""",
                          {'component_id': component_id,
                           'protein_group_id': protein_group_id,
                           'protein_accession': protein_accession,
                           'score': score,
                           'decoy': decoy}
                          )
                con.commit()
            unique_peptide_id_list = list(set(peptide_id_list))
            for peptide_id in unique_peptide_id_list:
                c.execute("""INSERT INTO PROTEIN_GROUP_PEPTIDE_MAPPING(
                PROTEIN_GROUP_ID, PEPTIDE_ID) VALUES(:protein_group_id, :peptide_id)
                """, {'protein_group_id': protein_group_id,
                      'peptide_id': peptide_id})
                con.commit()
            protein_group_id += 1
        component_id += 1
    c.close()


def end_connection(con) -> None:
    con.close()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("""usage: ParsimoniousProteinInferencer.py <sql file path>
              <q-value threshold for peptides> <context> <run id>""")
    # input_file path, q_limit_pep, context, run_id
    # if context is global, run_id can be anything
    # don't run global on 20180911_TIMS2_12-2_AnBr_SA_diaPASEF_Test10_42eV_1_A1_01_2927_test_lib_overridegroupid.osw
    # it does not have it
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
