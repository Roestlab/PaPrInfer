from __future__ import annotations

from typing import Dict, List, Tuple

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
    print("added all protein")
    print("made edges from protein to peptide")

    protein_peptide_graph.make_edges_from_peptide(con)
    print("made edges from peptide to protein")

    print("initialized")


def collapse(protein_peptide_graph: Graph) -> None:
    protein_peptide_graph.reorganize_protein_neighbours()
    print("protein_reorganized")
    protein_peptide_graph.reorganize_peptide_neighbours()
    print("peptide reorganized")
    print("collapsed")


def visualize(protein_peptide_graph):
    my_graph = protein_peptide_graph

    their_graph = igraph.Graph()

    # my_graph is a object with 2 dict
    # one for protein, one for peptide
    # the protein dict has protein as key
    # and peptide as values

    my_protein_dict = my_graph.get_protein_dict()
    my_peptide_dict = my_graph.get_peptide_dict()

    for protein_vertex in my_protein_dict:
        accession_string = ', '.join(map(str, protein_vertex.get_accession()))
        their_graph.add_vertex(accession_string)
        print(accession_string)

    for peptide_vertex in my_peptide_dict:
        id_string = ', '.join(map(str, peptide_vertex.get_id()))
        their_graph.add_vertex(id_string)
        print(id_string)

    print(their_graph.vs["name"])

    print("vertex done")

    protein_number = 1
    peptide_number = 1
    for protein_key in my_protein_dict.keys():
        accession_string = ', '.join(map(str, protein_key.get_accession()))
        print("protein number", protein_number)
        for peptide_value in my_protein_dict[protein_key]:
            id_string = ', '.join(map(str, peptide_value.get_id()))
            print("peptide number", peptide_number)
            their_graph.add_edge(accession_string, id_string)
            peptide_number += 1
        protein_number += 1


    print("edges done")

    layout = my_graph.layout_bipartite()
    print("layout done")
    plot = igraph.plot(their_graph, vertex_label=their_graph.vs["name"],
                       layout=layout)
    print("plot done")
    plot.show()


def separate(protein_peptide_graph: Graph) -> List[Component]:
    # for all white peptide nodes, explore them
    component_list = []
    component_counter = 0
    for current_peptide in protein_peptide_graph.get_peptide_dict_keys():
        if current_peptide.get_color() == 0:  # if it is white
            a_component = Component()
            protein_peptide_graph.dfs(current_peptide, a_component)
            component_list.append(a_component)
        print("component number", component_counter, "is done")
        component_counter += 1
    print("separated")
    return component_list


def reduce(component_list: List[Component], con) -> None:
    min_pro_list = []
    for component in component_list:

        # this is a list of list of protein accession
        # sublist is one meta-protein vertex, contain multiple protein accession
        component_accession_list = component.make_protein_list()
        for vertex in component_accession_list:
            if "Q8NE01" in vertex:
                print(vertex)
                print(component_accession_list)

        # this just a list of list of list of protein accession
        # sublist of it all belongs to the same component
        # sublist of the sublist belong to the same meta-protein vertex
        min_pro_list.append(component_accession_list)

    create_table_protein_group(con)
    protein_data_entry(con, min_pro_list)

    print("reduced")


# TODO i am not sure about this, but since I read from the Sqlite table directly
#  I do not need pyopenms's parser for idXML, so I dont even need pyopenms, the
#  module itself, I mean that is low coupling, so that is good right???

"""function that uses the Sqlite file directly"""


def begin_connection(db_name: str):
    con = sqlite3.connect(db_name)
    return con


def get_all_protein_id(con) -> List[int]:
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
            WHERE ID=:protein_id AND DECOY=0""",
            {"protein_id": protein_id}
    )
    # if there are any row, split the row (which are text) into List[str]
    for row in c.fetchall():
        accession_sublist = row[0].split(";")
        protein_accession_list.extend(accession_sublist)

    c.close()
    return protein_accession_list


def get_all_peptide(con, q_limit: int) -> List[int]:
    c = con.cursor()
    # choose all protein from the global context and less than the threshold
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
    #  smth like that, note this data is in SCORE_PROTEIN or SCORE_PEPTIDE


def get_link_for_protein(con, protein_id: int) -> List[int]:
    c = con.cursor()
    c.execute(
        """SELECT PEPTIDE_ID 
        FROM PEPTIDE_PROTEIN_MAPPING 
        WHERE PROTEIN_ID=:protein_id""",
        {'protein_id': protein_id})
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
            protein_group_id += 1
        component_id += 1
    c.close()


def end_connection(con) -> None:
    con.close()


class Node:
    """=== Private Attributes ===

    _merge_number: same number means that these nodes are to be merged into
    a meta node, initially set as 0
    """
    _color: int

    def __init__(self) -> None:
        self._color = 0

    def get_color(self) -> int:
        return self._color

    def is_white(self) -> bool:
        if self.get_color() == 0:
            return True
        else:
            return False

    def set_color(self, color: int) -> None:
        if not (color == 1 or color == 2):
            self._color = 0
        else:
            self._color = color

    def set_discovered(self) -> None:
        self._color = 1

    def set_explored(self) -> None:
        self._color = 2


class Protein(Node):
    _selected: bool
    _accession_list: List[str]
    _first_accession: str
    _ids: List[int]
    _first_id: int

    """initially, the protein group only has 1 accession number in it
    first accession serves as the immutable hash value for the dict key
    the protein_ids is just for making edges"""

    def __init_subclass__(cls, **kwargs) -> None:
        pass

    def __init__(self, accession_list: List[str], protein_id: int) -> None:
        super().__init__()
        self._selected = False
        self._accession_list = accession_list
        self._first_accession = accession_list[0]
        self._first_id = protein_id
        self._ids = [protein_id]

    def __hash__(self) -> int:

        return hash(self.get_first_accession())

    def get_both_first_accession(self, other: Protein) -> Tuple[str, str]:

        this_acc = self.get_first_accession()

        that_acc = other.get_first_accession()

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

    def get_first_accession(self) -> str:
        """
        this is for building the graph, because initially, there is only one id
        :return:
        """
        return self._first_accession

    def get_id(self):
        return self._ids

    def add_accession(self, accession_list: List[str]):
        self._accession_list.extend(accession_list)

    def get_accession(self) -> List[str]:
        return self._accession_list

    def add_id(self, protein_id: int) -> None:
        self._ids.append(protein_id)

    def get_first_id(self) -> int:
        return self._first_id


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

    def __init__(self):
        self._protein_dict = {}
        self._peptide_dict = {}

    def add_peptide(self, current_peptide: Peptide,
                    neighbours: List[Protein]) -> None:

        self._peptide_dict[current_peptide] = neighbours

    def add_protein(self, current_protein: Protein,
                    neighbours: List[Peptide]) -> None:

        self._protein_dict[current_protein] = neighbours

    def make_protein_list(self) -> List[List[str]]:

        component_min_pro_list = []

        while not self.all_component_peptides_covered():
            # select the protein with the most edges to uncovered peptide
            currently_selected_protein = self._find_most_uncovered_protein()
            component_min_pro_list.append(currently_selected_protein)

            # set the peptide that has edges to this selected protein as covered
            will_cover_peptide = self._protein_dict[currently_selected_protein]
            for current_peptide in will_cover_peptide:
                current_peptide.set_covered()

        component_accession_list = []
        # component_min_pro_list is a list of protein that are selected
        # each protein contain their own accession list (list of str)
        # collecting them into the component_accession_list
        # would a list of list of str
        # where each sublist is the belong to the same protein object
        for current_protein in component_min_pro_list:
            component_accession_list.append(current_protein.get_accession())

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
            if current_peptide.is_covered():
                continue
            else:
                num_uncovered_peptides += 1
        return num_uncovered_peptides

    def all_component_peptides_covered(self) -> bool:
        all_covered = True

        # find if any peptide is not yet covered
        for current_peptide in self._peptide_dict.keys():
            if not current_peptide.is_covered():
                all_covered = False

        return all_covered


class Graph:
    _protein_dict: Dict[Protein, List[Peptide]]
    _peptide_dict: Dict[Peptide, List[Protein]]

    def __init__(self):
        self._protein_dict = {}
        self._peptide_dict = {}

    def get_peptide_dict_keys(self):
        return self._peptide_dict.keys()

    def get_protein_dict(self):
        return self._protein_dict

    def get_peptide_dict(self):
        return self._peptide_dict

    """methods for step 1: initialize"""

    def add_protein(self, con, protein_id_list: List[int]) -> None:
        for protein_id in protein_id_list:
            protein_accession_list = get_protein_accession(con, protein_id)
            for protein_accession in protein_accession_list:
                current_protein = Protein([protein_accession], protein_id)
                if self.is_protein_key(current_protein):
                    self.make_edge_for_this_protein_id(con, current_protein)
                else:
                    self._protein_dict[current_protein] = []
                    self.make_edge_for_this_protein_id(con, current_protein)

    def add_peptide(self, peptide_id_list: List[int]) -> None:

        # for every peptide that is a hit
        for peptide_id in peptide_id_list:
            # make a Peptide node
            current_peptide = Peptide([peptide_id])

            # store that peptide into the graph with no edges yet
            self._peptide_dict[current_peptide] = []

    def make_edge_for_this_protein_id(self, con, the_protein: Protein) -> None:
        protein_id = the_protein.get_first_id()
        neighbour_list_pep = get_link_for_protein(con, protein_id)
        # neighbour list is just a list of peptide ids
        for neighbour_id in neighbour_list_pep:
            a_peptide = Peptide([neighbour_id])
            if self.is_peptide_key(a_peptide):
                self._protein_dict[the_protein].append(a_peptide)

    def make_edges_from_peptide(self, con) -> None:
        """
        iterate through the peptide keys, and for each, get their link (most
        times there is only one link), for each link, get all the accessions.
        for each accession, make a protein object, check if is a protein key,
        then add it to this peptide keys's value
        :param con:
        :return:
        """
        for peptide in self._peptide_dict.keys():
            neighbour_list_pro = get_link_for_peptide(con,
                                                      peptide.get_first_id())
            # by the SQLite format, len(neighbour_list_pro) = 1, (most times)
            for neighbour_id in neighbour_list_pro:
                accession_list = get_protein_accession(con, neighbour_id)
                for accession in accession_list:
                    a_protein = Protein([accession], neighbour_id)
                    if self.is_protein_key(a_protein):
                        self._peptide_dict[peptide].append(a_protein)

    def is_protein_key(self, a_protein: Protein) -> bool:
        """
        check if this id is in a key in the protein dict
        :return:
        """
        if self._protein_dict.get(a_protein) is None:
            truth = False
        else:
            truth = True

        return truth

    def is_peptide_key(self, a_peptide: Peptide) -> bool:
        """
        checks if this id is in a key in the peptide dict
        :param a_peptide:
        :return:
        """
        if self._peptide_dict.get(a_peptide) is None:
            truth = False
        else:
            truth = True

        return truth

    """methods for step 2: collapse"""

    def reorganize_protein_neighbours(self) -> None:

        # for each peptide node
        for current_peptide in self._peptide_dict.keys():

            # current neighbours: neighbours of the current peptide node
            # which are proteins
            current_neighbours = self._peptide_dict[current_peptide]

            # if there are more than 2 protein neighbours
            if len(current_neighbours) >= 2:

                # actually reorganize the neighbours
                neighbours_reorganized = self.reformat_pro(current_neighbours)

                # then use that to check mergeability
                self.check_for_mergeable_protein(neighbours_reorganized)

            print("peptide done", current_peptide.get_first_id())

    def compare_protein_neighbours(self, protein_pair: Tuple[Protein, ...])\
            -> bool:
        """
        if either one is None, return false
        :param protein_pair:
        :return:
        """
        this_node = protein_pair[0]
        that_node = protein_pair[1]
        these_neighbours = self._protein_dict.get(this_node)
        those_neighbours = self._protein_dict.get(that_node)
        if these_neighbours is None or those_neighbours is None:
            return False
        else:
            these_neighbours.sort()
            those_neighbours.sort()
            return these_neighbours == those_neighbours

    def num_protein_neighbour(self, peptide: Peptide) -> int:
        return len(self._peptide_dict[peptide])

    # right now, compare_protein_neighbours assumes that protein neighbour pair
    # is a tuples of variable length, though it really does not,
    # because 2 is specified?
    def check_for_mergeable_protein(self, protein_neighbours_reorganized:
                                    List[List[Protein]]) -> None:
        # all the protein neighbours in the same protein_neighbour_list
        # should have the same length
        for protein_neighbour_list in protein_neighbours_reorganized:
            for protein_neighbour_pair in itertools.combinations(
                    protein_neighbour_list, 2):
                if self.compare_protein_neighbours(protein_neighbour_pair):
                    self.delete_protein(protein_neighbour_pair[1])
                    accession_list = protein_neighbour_pair[1].get_accession()
                    protein_neighbour_pair[0].add_accession(accession_list)
                    self.add_accession_pro_key(protein_neighbour_pair[0], accession_list)

    def add_accession_pro_key(self, same_protein: Protein, accession_list):
        for protein_key in self._protein_dict:
            if same_protein.get_first_accession() == protein_key.get_first_accession():
                protein_key.add_accession(accession_list)


    def delete_protein(self, current_protein: Protein) -> None:

        # remove it as any peptide's neighbour
        peptide_neighbour_list = self._protein_dict[current_protein]
        for peptide_neighbour in peptide_neighbour_list:
            protein_neighbour_list = self._peptide_dict[peptide_neighbour]
            protein_neighbour_list.remove(current_protein)

        # remove protein node
        self._protein_dict.pop(current_protein)

    def reformat_pro(self, current_neighbours: List[Protein]) \
            -> List[List[Protein]]:
        """
        take the list of protein, current_neighbours and reformat it into a
        list of list of protein where every protein in each sublist has the same
        number of neighbours, which is returned
        :param current_neighbours:
        :return: neighbours_reorganized
        """

        neighbours_reorganized = []

        num_most_edges = self.find_protein_most_edges(current_neighbours)

        # append as many list as num_most_edges
        # where num_most_edges is the greatest number of
        # neighbours a protein has
        for i in range(0, num_most_edges + 1):
            neighbours_reorganized.append([])

        # for all protein neighbour
        for current_protein in current_neighbours:
            # find num of peptide neighbours for current protein
            num_peptide_neighbours = self.num_peptide_neighbour(
                current_protein)
            # append them to the list of list appropriately
            neighbours_reorganized[
                num_peptide_neighbours].append(current_protein)

        return neighbours_reorganized

    def reorganize_peptide_neighbours(self) -> None:
        """ peptide_neighbours: neighbouring node that are peptide
        note that any peptide only neighbour proteins
        that is a protein's neighbour is a peptide, and a peptide's neighbour is
        a protein
        """

        # for each protein node
        for current_protein in self._protein_dict.keys():
            current_neighbours = self._protein_dict[current_protein]

            # if there are more than 2 peptide neighbours
            if len(current_neighbours) >= 2:
                neighbour_reformatted = self.reformat_pep(current_neighbours)

                self.check_for_mergeable_peptide(neighbour_reformatted)

            print("protein done", current_protein.get_first_accession())

    def compare_peptide_neighbours(self, peptide_pair: Tuple[Peptide, ...])\
            -> bool:
        """
        if either one is None, return false
        :param peptide_pair:
        :return:
        """
        this_node = peptide_pair[0]
        that_node = peptide_pair[1]

        these_neighbours = self._peptide_dict.get(this_node)
        those_neighbours = self._peptide_dict.get(that_node)

        if these_neighbours is None or those_neighbours is None:
            return False
        else:
            these_neighbours.sort()
            those_neighbours.sort()
            return these_neighbours == those_neighbours

    def num_peptide_neighbour(self, protein: Protein) -> int:
        return len(self._protein_dict[protein])

    def find_protein_most_edges(self, protein_list: List[Protein]) -> int:
        max_edges = 0

        for current_protein in protein_list:
            peptide_neighbour_list = self._protein_dict[current_protein]
            if len(peptide_neighbour_list) > max_edges:
                max_edges = len(peptide_neighbour_list)

        return max_edges

    def find_peptide_most_edges(self, peptide_list: List[Peptide]) -> int:
        max_edges = 0

        for current_peptide in peptide_list:
            peptide_neighbour_list = self._peptide_dict[current_peptide]
            if len(peptide_neighbour_list) > max_edges:
                max_edges = len(peptide_neighbour_list)

        return max_edges

    def reformat_pep(self, current_neighbours: List[Peptide]) \
            -> List[List[Peptide]]:
        """
        take all the peptides in current_neighbours and store then in a
        different format, namely a list of list where every peptide in each
        sublist, has the same number of neighbours
        :param current_neighbours:
        :return:
        """

        neighbours_reformatted = []

        num_most_edges = self.find_peptide_most_edges(current_neighbours)

        # append as many list as num_most_edges
        # where num_most_edges is the greatest number of
        # neighbours a protein has
        for i in range(0, num_most_edges + 1):
            neighbours_reformatted.append([])

        # for all neighbour (which are peptides)
        for current_peptide in current_neighbours:
            # find num of peptide neighbours for current protein
            num_peptide_neighbours = self.num_protein_neighbour(current_peptide)
            # append them to the list of list appropriately
            neighbours_reformatted[
                num_peptide_neighbours].append(current_peptide)

        return neighbours_reformatted

    def check_for_mergeable_peptide(self, peptide_neighbour_reorganized: List[
                                    List[Peptide]]) -> None:
        """
        check every pairwise combination to see if there are indistinguishable
        proteins, if it does then call collapse, and combine id lists
        :param peptide_neighbour_reorganized:
        :return: return nothing
        """
        for peptide_neighbour_list in peptide_neighbour_reorganized:
            if len(peptide_neighbour_list) >= 2:
                for peptide_neighbour_pair in itertools.combinations(
                        peptide_neighbour_list, 2):
                    if self.compare_peptide_neighbours(peptide_neighbour_pair):
                        self.delete_peptide(peptide_neighbour_pair[1])
                        id_list = peptide_neighbour_pair[1].get_id()
                        peptide_neighbour_pair[0].add_id(id_list)
                        self.add_id_pep_key(peptide_neighbour_pair[0], id_list)

    def add_id_pep_key(self, same_peptide: Peptide, id_list):
        """
        it should work, because the one peptide that has matching first id
        with the deleted peptide, is deleted
        :param same_peptide:
        :param id_list:
        :return:
        """
        for peptide_key in self._peptide_dict:
            if same_peptide.get_first_id() == peptide_key.get_first_id():
                peptide_key.add_id(id_list)

    def delete_peptide(self, current_peptide: Peptide) -> None:
        """
        actually remove the protein, from the graph
        :param current_peptide:
        :return:
        """
        # get the protein neighbours of the peptide to be collapsed
        protein_neighbour_list = self._peptide_dict.get(current_peptide)
        # then for every protein neighbor
        for protein_neighbour in protein_neighbour_list:
            # remove the to-be-collapsed-peptide as their neighbour
            peptide_neighbour_list = self._protein_dict[protein_neighbour]
            peptide_neighbour_list.remove(current_peptide)
            self._protein_dict[protein_neighbour] = peptide_neighbour_list

        # pop the to-be-collapsed-peptide from the graph
        self._peptide_dict.pop(current_peptide)

        # assert current_peptide not in self._peptide_dict.keys()
        # assert current_peptide in self._peptide_dict.keys()

    """method for step 3: separate"""

    def set_discovered(self, start_node: Node) -> None:
        self.set_status(start_node, 1)

    def set_explored(self, start_node: Node) -> None:
        self.set_status(start_node, 2)

    def set_status(self, start_node: Node, status: int) -> None:
        start_node.set_color(status)

        peptide_neighbour_list = []
        protein_neighbour_list = []
        if isinstance(start_node, Protein):
            peptide_neighbour_list = self._protein_dict.get(start_node)
            # find it in the dictionary
            for protein in self._protein_dict:
                if protein == start_node:
                    protein.set_color(status)

        elif isinstance(start_node, Peptide):
            protein_neighbour_list = self._peptide_dict.get(start_node)
            # find it in the dictionary
            for peptide in self._peptide_dict:
                if peptide == start_node:
                    peptide.set_color(status)

        # find it as a neighbour of other node
        # sibling: share a node
        for neighbour in peptide_neighbour_list:
            neighbour_neighbour_list = self._peptide_dict[neighbour]
            for sibling in neighbour_neighbour_list:
                if sibling == start_node:
                    sibling.set_color(status)

        for neighbour in protein_neighbour_list:
            neighbour_neighbour_list = self._protein_dict[neighbour]
            for sibling in neighbour_neighbour_list:
                if sibling == start_node:
                    sibling.set_color(status)

    def dfs(self, start_node: Node, a_component: Component) -> None:
        neighbour_list = []

        self.set_discovered(start_node)
        # need to set other start node as discovered too
        #  because it has a different memory address
        #  it can exist as a neighbour of other node
        #  or a key node

        # if protein, get a peptide list, if peptide get a protein list
        if isinstance(start_node, Protein):
            neighbour_list = self._protein_dict.get(start_node)
        elif isinstance(start_node, Peptide):
            neighbour_list = self._peptide_dict.get(start_node)

        # for all neighbouring white node, explore them
        for current_neighbouring_node in neighbour_list:
            if current_neighbouring_node.is_white():
                self.dfs(current_neighbouring_node, a_component)

        self.set_explored(start_node)

        # after the nodes are blacken, add the protein into the a_component
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
