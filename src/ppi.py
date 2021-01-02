from __future__ import annotations

from typing import Dict, List, Tuple
from abc import abstractmethod

import itertools
import random
import sqlite3
import sys
import igraph

"""main"""


def main(input_file: str, q_limit_pep: str) -> None:
    con = begin_connection(input_file)

    q_limit = (int(q_limit_pep))

    protein_peptide_graph = Graph()

    initialize(protein_peptide_graph, con, q_limit)

    collapse(protein_peptide_graph)

    # visualize(protein_peptide_graph)

    component_list = separate(protein_peptide_graph)

    # for component in component_list:
    #     visualize_component(component)

    reduce(component_list, con)

    end_connection(con)


"""the helper function for main"""


def initialize(protein_peptide_graph: Graph, con,
               q_limit: int) -> None:
    peptide_id_list = get_all_peptide(con, q_limit)
    print("got all peptide")

    protein_id_list = get_all_protein_id(con)
    print("got all protein")

    protein_peptide_graph.add_peptide(peptide_id_list)
    print("added all peptide")

    protein_peptide_graph.add_protein(con, protein_id_list)
    print("added all protein with edges")

    protein_peptide_graph.make_edges_from_peptide(con)
    print("added edges from peptide")

    print("initialized")

    protein_num = 0
    peptide_num = 0
    for node in protein_peptide_graph.node_dict:
        if isinstance(node, Protein):
            protein_num += 1
        elif isinstance(node, Peptide):
            peptide_num += 1
        else:
            print("what")

    print(protein_num)
    print(peptide_num)


def collapse(protein_peptide_graph: Graph) -> None:
    protein_peptide_graph.collpase_graph()


def visualize(protein_peptide_graph):
    my_graph = protein_peptide_graph

    their_graph = igraph.Graph()

    # my_graph is a object with 1 dict
    # one for protein, one for peptide
    # the protein dict has protein as key
    # and peptide as values

    my_node_dict = my_graph.get_node_dict()

    for vertex in my_node_dict:
        if vertex.get_first_id() not in my_graph.node_to_delete:
            my_string = ', '.join(map(str, vertex.get_id()))
            their_graph.add_vertex(my_string)

    print("vertex done")

    protein_number = 1
    peptide_number = 1
    for key, value in my_node_dict.items():
        if key.get_first_id() not in my_graph.node_to_delete and isinstance(key,
                                                                            Protein):
            vertex_string_1 = ', '.join(map(str, key.get_id()))
            print("protein number", protein_number, "got protein vertex")
            for element in value:
                vertex_string_2 = ', '.join(map(str, element.get_id()))
                their_graph.add_edge(vertex_string_1, vertex_string_2)
                peptide_number += 1
            protein_number += 1

    print("edges done")

    layout = their_graph.layout_auto()
    print("layout done")
    plot = igraph.plot(their_graph, vertex_label=their_graph.vs["name"],
                       layout=layout)
    print("plot done")
    plot.show()


def separate(protein_peptide_graph: Graph) -> List[Component]:
    # for all white peptide nodes, explore them
    component_list = []
    component_counter = 0

    print(protein_peptide_graph.discovered_nodes)

    for node in protein_peptide_graph.node_dict.keys():
        # if it is white and not deleted)
        if protein_peptide_graph.is_white(node) and \
                node.get_first_id() not in protein_peptide_graph.node_to_delete:

            a_component = Component()
            protein_peptide_graph.dfs(node, a_component)
            component_list.append(a_component)

            print("component number", component_counter, "is seperated")
            component_counter += 1

    print("separated")
    return component_list


def visualize_component(component):
    """
    all nodes in a component are assumed to be not deleted
    :param component_list:
    :return:
    """

    their_graph = igraph.Graph()

    # my_graph is a object with 1 dict
    # one for protein, one for peptide
    # the protein dict has protein as key
    # and peptide as values

    for vertex in component._protein_dict:
        my_string = ', '.join(map(str, vertex.get_id()))
        their_graph.add_vertex(my_string, type = False)
    print("protein done")

    for vertex in component._peptide_dict:
        my_string = ', '.join(map(str, vertex.get_id()))
        their_graph.add_vertex(my_string, type = True)
    print("protein done")



    protein_number = 1
    peptide_number = 1
    for key, value in component._protein_dict.items():
        vertex_string_1 = ', '.join(map(str, key.get_id()))

        for element in value:

            vertex_string_2 = ', '.join(map(str, element.get_id()))
            their_graph.add_edge(vertex_string_1, vertex_string_2)
            peptide_number += 1
        print("protein number", protein_number, "got all edges")
        protein_number += 1

    print("edges done")

    layout = their_graph.layout_bipartite()
    print("layout done")
    plot = igraph.plot(their_graph, vertex_label=their_graph.vs["name"],
                       layout=layout)
    print("plot done")


def reduce(component_list: List[Component], con) -> None:
    min_pro_list = []
    component_accession_list_counter = 0
    for component in component_list:

        # this is a list of list of protein accession
        # sublist is one meta-protein vertex, contain multiple protein accession
        component_accession_list = component.make_protein_list()

        # this just a list of list of list of protein accession
        # sublist of it all belongs to the same component
        # sublist of the sublist belong to the same meta-protein vertex
        min_pro_list.append(component_accession_list)

        print("component number", component_accession_list_counter,
              "is reduced")
        component_accession_list_counter += 1

    create_table_protein_group(con)
    protein_data_entry(con, min_pro_list)

    print("reduced")


"""function that uses the Sqlite file directly"""


def begin_connection(db_name: str):
    con = sqlite3.connect(db_name)
    return con


def get_all_protein_id(con) -> List[int]:
    """
    choose all peptide from the global context
    :param con: the connection to the sqlite database
    :return: a list of all proteins
    """
    c = con.cursor()
    c.execute(
        """SELECT PROTEIN_ID 
        FROM SCORE_PROTEIN 
        WHERE CONTEXT='global'""")
    all_protein_id_list = []
    for row in c.fetchall():
        all_protein_id_list.append(row[0])
    c.close()
    return all_protein_id_list


def get_protein_accession(con, protein_id: int) -> List[str]:
    """
    returns the protein accession that has this id, as a list
    :param con:
    :param protein_id:
    :return:
    """
    c = con.cursor()
    protein_accession_list = []
    c.execute(
        """SELECT PROTEIN_ACCESSION 
            FROM PROTEIN 
            WHERE ID=:protein_sqlite_id AND DECOY=0""",
        {"protein_sqlite_id": protein_id}
    )
    # if there are any row, split the row (which are text) into List[str]
    for row in c.fetchall():
        accession_sublist = row[0].split(";")
        protein_accession_list.extend(accession_sublist)

    c.close()
    return protein_accession_list


def get_all_peptide(con, q_limit: int) -> List[int]:
    c = con.cursor()
    # choose all peptide from the global context and less than the threshold
    c.execute(
        """SELECT PEPTIDE_ID FROM SCORE_PEPTIDE
        INNER JOIN PEPTIDE ON PEPTIDE.ID = SCORE_PEPTIDE.PEPTIDE_ID
        WHERE QVALUE<:q_limit AND DECOY=0""",
        {'q_limit': q_limit}
    )
    all_peptide_id_list = []

    for row in c.fetchall():
        # each row is a tuple, since c.fetchall() returns a list of tuples
        all_peptide_id_list.append(row[0])

    c.close()
    return all_peptide_id_list

    # TODO if i ever want to select only peptide or protein that fit a certain
    #  properties, use c.execute("SELECT * WHERE value=3 AND keyword=1")
    #  smth like that, note this data and result is in SCORE_PROTEIN or SCORE_PEPTIDE


def get_link_for_protein(con, protein_sqlite_id: int) -> List[int]:
    c = con.cursor()
    c.execute(
        """SELECT PEPTIDE_ID 
        FROM PEPTIDE_PROTEIN_MAPPING 
        WHERE PROTEIN_ID=:protein_sqlite_id""",
        {'protein_sqlite_id': protein_sqlite_id})
    peptide_id_list = []
    for row in c.fetchall():
        peptide_id_list.append(row[0])
    c.close()
    return peptide_id_list


def get_link_for_peptide(con, peptide_id: int) -> List[int]:
    c = con.cursor()
    c.execute(
        """SELECT PROTEIN_ID 
        FROM PEPTIDE_PROTEIN_MAPPING 
        where PEPTIDE_ID=:peptide_id""",
        {"peptide_id": peptide_id})
    protein_id_list = []
    for row in c.fetchall():
        protein_id_list.append(row[0])
    c.close()
    return protein_id_list


def create_table_protein_group(con) -> None:
    """mpl = minimal protein list"""
    c = con.cursor()
    c.execute("""DROP TABLE IF EXISTS PROTEIN_GROUPS""")
    c.execute("""CREATE TABLE PROTEIN_GROUPS(
              COMPONENT_ID INTEGER,
              PROTEIN_GROUP_ID INTEGER, 
              PROTEIN_ID INTEGER);""")
    c.close()


def protein_data_entry(con, min_pro_list: List[List[List[str]]]) -> None:
    c = con.cursor()
    component_id = 0
    protein_group_id = 0
    for component in min_pro_list:
        for protein_group in component:
            for protein_accession in protein_group:
                c.execute("""INSERT INTO PROTEIN_GROUPS(COMPONENT_ID, 
                PROTEIN_GROUP_ID, PROTEIN_ID) VALUES(:component_id, 
                :protein_group_id, :protein_accession)""",
                          {'component_id': component_id,
                           'protein_group_id': protein_group_id,
                           'protein_accession': protein_accession}
                          )
                con.commit()
            print("protein_group", protein_group_id, "is entered")
            protein_group_id += 1
        print("component", component_id, "is entered")
        component_id += 1
    c.close()


def end_connection(con) -> None:
    con.close()


class Node:
    """=== Private Attributes ===

    _merge_number: same number means that these nodes are to be merged into
    a meta node, initially set as 0
    """

    def __init__(self) -> None:
        pass

    @abstractmethod
    def get_id(self):
        pass

    @abstractmethod
    def add_id(self, accession_list):
        pass

    @abstractmethod
    def get_first_id(self):
        pass


class Protein(Node):
    _selected: bool
    _id_list: List[str]
    _first_id: str
    _sqlite_ids: List[int]
    _first_sqlite_id: int

    """initially, the protein group only has 1 accession number in it
    first accession serves as the immutable hash value for the dict key
    the protein_ids is just for making edges"""

    def __init_subclass__(cls, **kwargs) -> None:
        pass

    def __init__(self, id_list: List[str], protein_sqlite_id: int) -> None:
        super().__init__()
        self._selected = False
        self._id_list = id_list
        self._first_id = id_list[0]
        self._first_sqlite_id = protein_sqlite_id
        self._sqlite_ids = [protein_sqlite_id]

    def __hash__(self) -> int:
        return hash(self.get_first_id())

    def get_both_first_accession(self, other: Protein) -> Tuple[str, str]:
        this_acc = self.get_first_id()

        that_acc = other.get_first_id()

        return this_acc, that_acc

    def __eq__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc == that_acc

    def __ne__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc != that_acc

    def __lt__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc < that_acc

    def __le__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc <= that_acc

    def __gt__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc > that_acc

    def __ge__(self, other: Protein) -> bool:
        this_acc, that_acc = self.get_both_first_accession(other)
        return this_acc >= that_acc

    def is_selected(self) -> bool:
        return self._selected

    def set_selected(self) -> None:
        self._selected = True

    def get_first_id(self) -> str:
        """
        this is for building the graph, because initially, there is only one id
        :return:
        """
        return self._first_id

    def add_id(self, accession_list: List[str]):
        self._id_list.extend(accession_list)

    def get_id(self) -> List[str]:
        return self._id_list

    def get_sqlite_id(self):
        return self._sqlite_ids

    def add_sqlite_id(self, protein_sqlite_id: int) -> None:
        self._sqlite_ids.append(protein_sqlite_id)

    def get_first_sqlite_id(self) -> int:
        return self._first_sqlite_id


class Peptide(Node):
    _covered: bool
    _first_id: int
    _ids: List[int]

    """initially, the protein group only has 1 id in it"""

    def __init__(self, peptide_id_list: List[int]) -> None:
        super().__init__()
        self._covered = False
        self._ids = peptide_id_list
        self._first_id = peptide_id_list[0]

    def __hash__(self) -> int:
        return hash(self.get_first_id())

    def get_both_ids(self, other: Peptide) -> Tuple[int, int]:
        this_id = self.get_first_id()

        that_id = other.get_first_id()

        return this_id, that_id

    def __eq__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id == that_id

    def __ne__(self, other) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id != that_id

    def __lt__(self, other) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id < that_id

    def __le__(self, other) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id <= that_id

    def __gt__(self, other) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id > that_id

    def __ge__(self, other) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id >= that_id

    def is_covered(self) -> bool:
        return self._covered

    def set_covered(self) -> None:
        self._covered = True

    def get_first_id(self) -> int:
        return self._first_id

    def get_id(self) -> List[int]:
        return self._ids

    def add_id(self, id_list: List[int]):
        self._ids.extend(id_list)


class Component:
    _protein_dict: Dict[Protein, List[Peptide]]
    _peptide_dict: Dict[Peptide, List[Protein]]
    _covered_peptide: dict

    def __init__(self):
        self._protein_dict = {}
        self._peptide_dict = {}
        self._covered_peptide = {}

    def add_peptide(self, current_peptide: Peptide,
                    neighbours) -> None:

        self._peptide_dict[current_peptide] = neighbours

    def add_protein(self, current_protein: Protein,
                    neighbours) -> None:

        self._protein_dict[current_protein] = neighbours

    def make_protein_list(self) -> List[List[str]]:

        component_min_pro_list = []

        while not self.all_component_peptides_covered():
            # select the protein with the most edges to uncovered peptide
            currently_selected_protein = self._find_most_uncovered_protein()
            component_min_pro_list.append(currently_selected_protein)
            currently_selected_protein.set_selected()

            # set the peptide that has edges to this selected protein as covered
            # that is add it to the covered peptide dict
            will_cover_peptide = self._protein_dict[currently_selected_protein]
            for current_peptide in will_cover_peptide:
                self._covered_peptide[current_peptide.get_first_id()] = ''

        component_accession_list = []
        # component_min_pro_list is a list of protein that are selected
        # each protein contain their own accession list (list of str)
        # collecting them into the component_accession_list
        # would a list of list of str
        # where each sublist is the belong to the same protein object
        for current_protein in component_min_pro_list:
            component_accession_list.append(current_protein.get_id())

        return component_accession_list

    def _find_most_uncovered_protein(self) -> Protein:

        # pick any protein from the a_component
        most_edges_protein = random.choice(list(self._protein_dict.keys()))

        # record its number of edges
        most_uncovered = self.find_num_uncovered_peptides(most_edges_protein)

        # loop over the proteins of a a_component
        for current_protein in self._protein_dict.keys():

            # of all the protein that are not selected, find the one with
            # the most edges (to peptides)
            if current_protein.is_selected():
                continue
            else:
                current_protein_num_uncovered = \
                    self.find_num_uncovered_peptides(current_protein)
                if current_protein_num_uncovered > most_uncovered:
                    most_uncovered = current_protein_num_uncovered
                    most_edges_protein = current_protein

        most_edges_protein.set_selected()
        return most_edges_protein

    def find_num_uncovered_peptides(self, current_protein: Protein) -> int:
        num_uncovered_peptides = 0
        peptide_neighbours = self._protein_dict[current_protein]
        for current_peptide in peptide_neighbours:
            if current_peptide.get_first_id() in self._covered_peptide:
                continue
            else:
                num_uncovered_peptides += 1
        return num_uncovered_peptides

    def all_component_peptides_covered(self) -> bool:
        all_covered = True

        # find if any peptide is not yet covered
        for peptide in self._peptide_dict:
            if peptide.get_first_id() not in self._covered_peptide:
                all_covered = False

        return all_covered


class Graph:
    """
    The node_to_delete, discovered_nodes, explored_nodes with key as the first id
    of the node, and an empty string as value
    """
    node_dict: Dict[Node, List[Node]]
    node_to_delete: Dict[int, str]
    # node_to_add_id: dict
    discovered_nodes: Dict[int, str]
    explored_nodes: Dict[int, str]

    def __init__(self):
        self.node_dict = {}
        self.node_to_delete = {}
        # self.node_to_add_id = {}
        self.discovered_nodes = {}
        self.explored_nodes = {}

    def get_node_dict_keys(self):
        return self.node_dict.keys()

    def get_node_dict(self):
        return self.node_dict

    """methods for step 1: initialize"""

    def add_peptide(self, peptide_id_list: List[int]) -> None:

        # for every peptide that is a hit
        for peptide_id in peptide_id_list:
            # make a Peptide node
            current_peptide = Peptide([peptide_id])
            # store that peptide into the graph
            self.node_dict[current_peptide] = []

    def add_protein(self, con, protein_id_list: List[int]) -> None:
        """
        the protein with this accession may already be in the graph
        because in the sqlite database, 1 protein sqlite id correspond to
        multiple protein accessions, and they may overlap (a protein accession
        can exist in multiple protein sqlite ids)
        :param con:
        :param protein_id_list:
        :return:
        """
        for protein_id in protein_id_list:
            protein_accession_list = get_protein_accession(con, protein_id)
            for protein_accession in protein_accession_list:
                current_protein = Protein([protein_accession], protein_id)
                if self.node_in_graph(current_protein):
                    self.make_edge_from_protein_id(con, current_protein)
                else:
                    self.node_dict[current_protein] = []
                    self.make_edge_from_protein_id(con, current_protein)

    def make_edge_from_protein_id(self, con, the_protein: Protein) -> None:
        """
        use the protein id given, add peptides as edges to the current protein
        :param con:
        :param the_protein:
        :return:
        """

        protein_id = the_protein.get_first_sqlite_id()
        neighbour_list_pep = get_link_for_protein(con, protein_id)
        # neighbour list is just a list of peptide ids
        for neighbour_id in neighbour_list_pep:
            a_peptide = Peptide([neighbour_id])
            if self.node_in_graph(a_peptide):
                self.node_dict[the_protein].append(a_peptide)

    def make_edges_from_peptide(self, con) -> None:
        """
        iterate through the peptide keys, and for each, get their link (most
        times there is only one link), for each link, get all the accessions.
        for each accession, make a protein object, check if is a protein key,
        then add it to this peptide keys's value
        :param con: connection to the sqlite database
        :return:
        """
        for node in self.node_dict:
            if isinstance(node, Peptide):
                neighbour_list_pro = get_link_for_peptide(con,
                                                          node.get_first_id())
                for neighbour_id in neighbour_list_pro:
                    accession_list = get_protein_accession(con, neighbour_id)
                    for accession in accession_list:
                        a_protein = Protein([accession], neighbour_id)
                        self.node_dict[node].append(a_protein)

    def node_in_graph(self, a_node: Node):
        """
        check if this node is in the graph
        :param a_node: this node to be check
        :return: whether or not the node is in the graph (boolean)
        """
        if self.get_node_dict().get(a_node) is None:
            in_graph = False
        else:
            in_graph = True

        return in_graph

    """methods for step 2: collapse"""

    def collpase_graph(self) -> None:

        # for each node
        progress_count = 1
        for key, value in self.node_dict.items():

            # if the node has been delete, skip
            if key.get_first_id() in self.node_to_delete:
                progress_count += 1
                continue

            # current neighbours: neighbours of the current node
            current_neighbours = value

            # so uh, actually some nodes are deleted right
            # I could just skip checking its neighbour?

            # if there are more than 2 neighbours
            if len(current_neighbours) >= 2:

                # acually reorganize the neighbours
                neighbours_reorganized = self.reorder_neighbours(
                    current_neighbours)

                # then use that to check mergeability
                self.check_for_merging(neighbours_reorganized)

            progress_count += 1

            print(progress_count, len(self.node_dict))

    """methods for step 2a: reordering"""

    def find_node_most_edges(self, node_list: List[Node]) -> int:
        max_edges = 0

        for current_node in node_list:
            neighbour_list = self.node_dict[current_node]
            if len(neighbour_list) > max_edges:
                max_edges = len(neighbour_list)

        return max_edges

    def num_neighbour(self, node: Node) -> int:
        return len(self.node_dict[node])

    def reorder_neighbours(self, current_neighbours: List[Node]) \
            -> List[List[Node]]:
        """
        take the list of protein, current_neighbours and reformat it into a
        list of list of protein where every protein in each sublist has the same
        number of neighbours, which is returned
        :param current_neighbours:
        :return: neighbours_reorganized
        """

        neighbours_reorganized = []

        num_most_edges = self.find_node_most_edges(current_neighbours)

        # append as many list as num_most_edges
        # where num_most_edges is the greatest number of
        # neighbours a protein has
        for i in range(0, num_most_edges + 1):
            neighbours_reorganized.append([])

        # for all neighbour
        for current_protein in current_neighbours:
            # find num of neighbours for current node
            num_peptide_neighbours = self.num_neighbour(
                current_protein)
            # append them to the list of list appropriately
            neighbours_reorganized[
                num_peptide_neighbours].append(current_protein)

        return neighbours_reorganized

    """methods for step 2b: merging"""

    def check_for_merging(self, neighbours_reorganized: List[List[Node]]) \
            -> None:
        # all the protein neighbours in the same protein_neighbour_list
        # should have the same length
        for neighbour_list in neighbours_reorganized:
            for neighbour_pair in itertools.combinations(neighbour_list, 2):
                # if both nodes are not deleted, then compare them.
                if (
                        neighbour_pair[
                            0].get_first_id() not in self.node_to_delete
                        and
                        neighbour_pair[
                            1].get_first_id() not in self.node_to_delete
                ) and self.compare_neighbours(neighbour_pair):

                    self.delete_node(neighbour_pair[1])

                    id_list = neighbour_pair[1].get_id()

                    neighbour_pair[0].add_id(id_list)

                    self.key_add_id(neighbour_pair[0], id_list)

    def compare_neighbours(self, node_pair: Tuple[Node, ...]) \
            -> bool:
        """
        if either one is None, return false
        :param node_pair: the pair of node to be compared
        :return: True if the neighbour of the node pair are the same
        """
        this_node = node_pair[0]
        that_node = node_pair[1]
        these_neighbours = self.node_dict.get(this_node)
        those_neighbours = self.node_dict.get(that_node)
        if these_neighbours is None or those_neighbours is None:
            return False
        else:
            these_neighbours.sort()
            those_neighbours.sort()

            return these_neighbours == those_neighbours

    def delete_node(self, current_node: Node) -> None:

        # store the id of the node into the "node_to_delete" dict
        # with an empty string
        self.node_to_delete[current_node.get_first_id()] = ''

    def key_add_id(self, same_protein: Node, accession_list):
        for key in self.node_dict.keys():
            if key == same_protein:
                key.add_id(accession_list)

    """methods for step 3: separate"""

    def set_discovered(self, start_node: Node) -> None:
        self.discovered_nodes[start_node.get_first_id()] = ''

    def set_explored(self, start_node: Node) -> None:
        self.explored_nodes[start_node.get_first_id()] = ''

    def is_white(self, node: Node) -> bool:
        # print("discovered", node.get_first_id() in self.discovered_nodes)
        # print("explored", node.get_first_id in self.explored_nodes)

        return (node.get_first_id() not in self.discovered_nodes) \
               and (node.get_first_id not in self.explored_nodes)

    def is_discoverd(self, node) -> bool:
        return node.get_first_id() in self.discovered_nodes

    def is_explored(self, node) -> bool:
        return node.get_first_id() in self.explored_nodes

    def dfs(self, start_node: Node, a_component: Component) -> None:

        self.set_discovered(start_node)

        # if protein, get a peptide list, if peptide get a protein list
        neighbour_list = self.node_dict.get(start_node)

        # for all neighbouring white node, explore them
        for current_neighbour in neighbour_list:
            if self.is_white(current_neighbour) \
                    and current_neighbour.get_first_id() not in self.node_to_delete:
                self.dfs(current_neighbour, a_component)

        self.set_explored(start_node)

        # after the nodes are blacken,
        # add the protein or peptide into the a_component
        # with its neighbours
        if isinstance(start_node, Protein):
            a_component.add_protein(start_node, neighbour_list)
        elif isinstance(start_node, Peptide):
            a_component.add_peptide(start_node, neighbour_list)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("""usage: ParsimoniousProteinInferencer.py <sql file path>
              <q-value threshold for peptides>""")
    main(sys.argv[1], sys.argv[2])
