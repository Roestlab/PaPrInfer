import bisect
import itertools
from typing import Dict, List, Tuple, Any

from components import Component
from node import Node, Protein, Peptide


class Graph:
    """
    The node_to_delete, discovered_nodes, explored_nodes are dict
    with key as the first id of the node, and an empty string as value
    """
    node_dict: Dict[Node, List[Node]]
    node_to_delete: Dict[str, str]
    # node_to_add_id: dict
    discovered_nodes: Dict[str, str]
    explored_nodes: Dict[str, str]
    sorted_keys = List[Node]
    sorted_keys_id = List[str]

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

    def add_peptide(self,
                    peptide_id_list: List[Tuple[str, float, int]]) -> None:

        # for every peptide and score pair that is a hit
        for peptide_id, peptide_score, decoy in peptide_id_list:
            # make a Peptide node
            current_peptide = Peptide([peptide_id], peptide_score, decoy)
            # store that peptide into the graph
            self.node_dict[current_peptide] = []

    def add_protein(self, protein_id_list: List[Tuple[str, int]],
                    linked_peptide_dict: Dict[str, List[Tuple[str, float, int]]]
                    , protein_accession_dict: Dict[str, List[str]]) \
            -> None:
        """
        the protein with this accession may already be in the graph
        because in the sqlite database, 1 protein sqlite id correspond to
        multiple protein accessions, and they may overlap (a protein accession
        can exist in multiple protein sqlite ids)
        """
        for protein_id, decoy in protein_id_list:

            protein_accession_list = protein_accession_dict[protein_id]

            for protein_accession in protein_accession_list:
                current_protein = Protein([protein_accession], protein_id,
                                          decoy)
                if not self.node_in_graph(current_protein):
                    self.node_dict[current_protein] = []
                self.make_edge_from_protein_id(current_protein,
                                               linked_peptide_dict)

    def make_edge_from_protein_id(self,
                                  the_protein: Protein,
                                  linked_peptide_dict:
                                  Dict[str, List[
                                      Tuple[str, float, int]]]) \
            -> None:
        """
        based the protein given, add peptides as edges to the current protein
        """

        protein_id = the_protein.get_first_sqlite_id()
        neighbour_list_pep = linked_peptide_dict[protein_id]
        # neighbour list is just a list of peptide ids, with score and decoy
        at_least_one = False
        for neighbour_id, neighbour_score, decoy in neighbour_list_pep:
            a_peptide = Peptide([neighbour_id], neighbour_score, decoy)
            if self.node_in_graph(a_peptide):
                self.node_dict[the_protein].append(a_peptide)
                at_least_one = True

        if not at_least_one:
            self.delete_node(the_protein)

    def make_edges_from_peptide(self,
                                linked_protein_dict: Dict[str,
                                                          List[Tuple[str, int]]]
                                , protein_accession_dict: Dict[str, List[str]]) \
            -> None:
        """
        iterate through the peptide keys, and for each, get their link (most
        times there is only one link), for each link, get all the accessions.
        for each accession, make a protein object, check if is a protein key,
        then add it to this peptide keys's value
        """
        for node in self.node_dict:
            if isinstance(node, Peptide):
                peptide_id = node.get_first_id()
                # a list of protein sqlite id that links to this peptide
                neighbour_list_pro = linked_protein_dict[peptide_id]
                for protein_id, decoy in neighbour_list_pro:
                    accession_list = protein_accession_dict[protein_id]
                    for accession in accession_list:
                        a_protein = Protein([accession], protein_id, decoy)
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

    def store(self) -> None:
        pass
        # with open('initialized_graph.csv', 'w', newline='\n') as csvfile:
        #     fieldnames = []
        #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        #     writer.writeheader()
        #     writer.writerows(self.node_dict)

    def get_sort_keys(self) -> None:

        """
        :param
        """

        # convert dict key to list (since dict key cannot be sorted)
        # sort, then store
        keys_list = list(self.node_dict.keys())
        keys_list.sort()
        sorted_keys_list = keys_list
        self.sorted_keys = sorted_keys_list

        # get all id (first id + target_decoy) of those keys (in sorted order)
        sorted_id_list = []
        for keys in sorted_keys_list:
            sorted_id_list.append(keys.get_first_id() + keys.get_target_decoy())

        # sorted_id_list is then in the same order as sorted keys
        self.sorted_keys_id = sorted_id_list

    """methods for step 2: collapse"""

    def collapse_graph(self) -> None:

        # for each node
        progress_count = 1
        for key, value in self.node_dict.items():

            # skipping peptide for now
            if isinstance(key, Peptide):
                progress_count += 1
                continue

            # if the node has been delete, skip
            if key.get_first_id() in self.node_to_delete:
                progress_count += 1
                continue

            # current neighbours: neighbours of the current node
            current_neighbours = value

            # if there are more than 2 neighbours
            if len(current_neighbours) >= 2:

                # actually reorganize the neighbours
                neighbours_reorganized = self.reorder_neighbours(
                    current_neighbours)

                # for neighbours_list in neighbours_reorganized:
                #     for neighbour in neighbours_list:
                #         print("peptide id", neighbour.get_first_id())

                print(list(map(len, neighbours_reorganized)))

                # further group the neighbours
                # as long as this also is not quadratic, its better
                for index in range(len(neighbours_reorganized)):
                    # skip the first one, since that neighbour list only maps
                    # to 0 protein and the second one since that map to 1 protein
                    if index == 0 or index == 1:
                        continue

                    neighbours_list = neighbours_reorganized[index]
                    print("index", index)

                    further_group = self.group_again(neighbours_list, key.get_first_id())
                    print(list(map(len, further_group)))

                    # then use that to check mergeability
                    self.check_for_merging(further_group)

            progress_count += 1
            print("nodes collapsed", progress_count, len(self.node_dict),
                  "node is peptide", isinstance(key, Peptide))

    """methods for step 2a: reordering"""

    def find_node_most_edges(self, node_list: List[Node]) -> int:
        max_edges = 0

        for current_node in node_list:
            neighbour_list = self.node_dict[current_node]
            if len(neighbour_list) > max_edges:
                max_edges = len(neighbour_list)

        return max_edges

    def num_neighbour(self, node: Node) -> int:
        return len(self.node_dict.get(node))

    def reorder_neighbours(self, current_neighbours: List[Node]) \
            -> List[List[Node]]:
        """
        take the list of protein, current_neighbours and reformat it into a
        list of list of protein where every protein in each sublist has the same
        number of neighbours, which is returned
        :param current_neighbours:
        :return: neighbours_reorganized
        [0, 0, 28] means of the current neighbours of this protein (which are peptide).
        0 map to exactly 0 protein, 0 map to exactly 1 protein,
        28 map to exactly 2 proteins
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
            num_peptide_neighbours = self.num_neighbour(current_protein)
            # append them to the list of list appropriately
            neighbours_reorganized[
                num_peptide_neighbours].append(current_protein)

        return neighbours_reorganized

    def group_again(self, all_peptides: List[Node], first_protein_id) \
            -> List[List[Peptide]]:
        """
        terminology assumes we are trying to merge peptides
        :param all_peptides is all the peptides that have the same number
        of protein that they map to
        :param first_protein_id is the protein that all these peptides definitely
        maps to
        """

        peptide_same_2nd_protein_dict = {}
        for peptide in all_peptides:
            # is a neighbour and a protein
            neighbour_and_protein = self.node_dict.get(peptide)

            # print("num neighbours 2", len(neighbour_and_protein),
            #       "peptide id", peptide.get_first_id(),
            #       "peptide target decoy", peptide.get_target_decoy())
            # needs to remove the protein object has this first_protein_id
            for protein in neighbour_and_protein:
                if protein.get_first_id() == first_protein_id:
                    neighbour_and_protein.remove(protein)

            neighbour_and_protein.sort()

            # first neighbours is the one that was removed
            second_neighbour = neighbour_and_protein[0]

            # if key was not in the dict, setdefault return the default value,
            # empty list here. if key is in, then function returns the value
            peptide_same_2nd_protein_dict\
                .setdefault(second_neighbour, []).append(peptide)


        peptide_same_2nd_protein_list = []
        # key is second neighbour (protein object)
        # value is peptide that map to this second neighbour
        for key, value in peptide_same_2nd_protein_dict.items():
            peptide_same_2nd_protein_list.append(value)

        return peptide_same_2nd_protein_list

    """methods for step 2b: merging"""

    def check_for_merging(self, neighbours_reorganized: List[List[Node]]) \
            -> None:
        # all the protein neighbours in the same protein_neighbour_list
        # should have the same length
        for neighbour_list in neighbours_reorganized:
            for neighbour_pair in itertools.combinations(neighbour_list, 2):
                if (
                # the 2 both neighbour have not been delete yet
                    neighbour_pair[0].get_first_id() not in self.node_to_delete
                    and
                    neighbour_pair[1].get_first_id() not in self.node_to_delete

                # then compare their neighbours
                ) and self.compare_neighbours(neighbour_pair):

                    # either both target and both decoy
                    assert neighbour_pair[0].get_target_decoy() == \
                           neighbour_pair[1].get_target_decoy()

                    self.delete_node(neighbour_pair[1])
                    id_list = neighbour_pair[1].get_id()
                    score = neighbour_pair[1].get_score()
                    neighbour_pair[0].add_id(id_list)
                    # this only keeps the max score of them two
                    neighbour_pair[0].add_score(score)
                    # TODO: bottleneck
                    self.key_add_id(neighbour_pair[0], id_list)

    def compare_neighbours(self, node_pair: Tuple[Node, ...]) \
            -> bool:
        """
        if either one is None, return false
        :param node_pair: the pair of node to be compared
        :return: True if the neighbour of the node pair are the same
        """
        # if their first id and target+decoy is the same, then it is the same
        this_node = node_pair[0]
        that_node = node_pair[1]
        these_neighbours = self.node_dict.get(this_node)
        those_neighbours = self.node_dict.get(that_node)
        if these_neighbours is None or those_neighbours is None:
            return False
        else:
            these_neighbours.sort()
            those_neighbours.sort()

        assert len(these_neighbours) == len(those_neighbours)

        # I dont know if python ordered list comparison for equality has
        # early stop, such that when it reaches an index that is different
        # it exits the equality comparison

        is_the_same = True
        for node_id in range(len(these_neighbours)):
            if these_neighbours[node_id] != those_neighbours[node_id]:
                is_the_same = False
                break

        return is_the_same

    def delete_node(self, current_node: Node) -> None:

        # store the id of the node into the "node_to_delete" dict
        # with an empty string
        self.node_to_delete[current_node.get_first_id()
                            + current_node.get_target_decoy()] = ''

    def key_add_id(self, same_protein: Node, accession_list) -> None:

        # find this protein in the id_list, with binary search
        protein_id = same_protein.get_first_id() + same_protein.get_target_decoy()
        if protein_id in self.sorted_keys_id:
            # find the index of the node in sorted keys, based on its
            # (first id + target/decoy)
            index = bisect.bisect_left(self.sorted_keys_id, protein_id)

            # it is the same index for the protein object in the sorted_key_list
            # so then access it from sorted keys using list indexing
            # this gets me the protein
            equivalent_protein = self.sorted_keys[index]
            equivalent_protein.add_id(accession_list)
            equivalent_protein.add_score(same_protein.get_score())

    """methods for step 3: separate"""

    def set_discovered(self, start_node: Node) -> None:
        self.discovered_nodes[start_node.get_first_id()
                              + start_node.get_target_decoy()] = ''

    def set_explored(self, start_node: Node) -> None:
        self.explored_nodes[start_node.get_first_id()
                            + start_node.get_target_decoy()] = ''

    def is_white(self, node: Node) -> bool:
        return (node.get_first_id() + node.get_target_decoy()
                not in self.discovered_nodes) \
               and (node.get_first_id() + node.get_target_decoy()
                    not in self.explored_nodes)

    def is_discovered(self, node) -> bool:
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
                    and (current_neighbour.get_first_id()
                         not in self.node_to_delete):
                self.dfs(current_neighbour, a_component)

        self.set_explored(start_node)

        # after the nodes are blacken,
        # add the protein or peptide into the a_component
        # with its neighbours
        if isinstance(start_node, Protein):
            a_component.add_protein(start_node, neighbour_list)
        elif isinstance(start_node, Peptide):
            a_component.add_peptide(start_node, neighbour_list)
