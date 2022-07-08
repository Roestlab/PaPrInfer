import bisect
import itertools
from typing import Dict, List, Tuple

from components import Component
from node import Node, Protein, Peptide


class Graph:
    """
    The node_to_delete, discovered_nodes, explored_nodes are dict
    with key as the first id of the node, and an empty string as value
    """
    node_dict: Dict[Node, List[Node]]
    node_to_delete: Dict[str, str]
    accession_object: Dict[str, Node]
    discovered_nodes: Dict[str, str]
    explored_nodes: Dict[str, str]
    sorted_keys = List[Node]
    sorted_keys_id = List[str]

    def __init__(self):
        self.node_dict = {}
        self.node_to_delete = {}
        self.accession_object = {}
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
                    linked_peptide_dict: Dict[str, List[Tuple[str, float, int]]],
                    protein_accession_dict: Dict[str, List[str]]) \
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

                # if node is not in graph, make empty list first
                if not self.node_in_graph(current_protein):
                    self.node_dict[current_protein] = []
                self.make_edge_from_protein_id(current_protein,
                                               linked_peptide_dict)

            # TODO testing
            if 'sp|Q9UI08-5|EVL_HUMAN' in protein_accession_list:
                print('what', protein_id, decoy)


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

        if protein_id not in linked_peptide_dict:
            return

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
                                linked_protein_dict: Dict[str, List[Tuple[str, int]]],
                                protein_accession_dict: Dict[str, List[str]]) \
            -> None:
        """
        iterate through the peptide keys, and for each, get their link (most
        times there is only one link), for each link, get all the accessions.
        for each accession, make a protein object, check if is a protein key,
        then add it to this peptide keys' value
        """
        for node in self.node_dict:
            if isinstance(node, Peptide):

                peptide_id = node.get_first_id()

                if peptide_id not in linked_protein_dict:
                    continue

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
        :param a_node: this node to be checked
        :return: whether the node is in the graph (boolean)
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

        # take all nodes (keys)
        key_list = list(self.node_dict)

        # then group them based on their degree as a list of lists
        all_nodes_reorganized = self.reorder_neighbours(key_list)

        # TODO: Better name
        max_degree = len(all_nodes_reorganized)

        # the value is n when the nodes from degree n to 1 are finished
        progress_count_in_degrees = 1

        # i.e. for every degree
        for num_edges in range(len(all_nodes_reorganized)):

            # skip the group of protein that only map 0 proteins
            # i.e. skip the degree of zero
            if num_edges == 0:
                continue

            # these are all the nodes with the same degree
            same_edges_node_list = all_nodes_reorganized[num_edges]

            # if there are no node that have exactly num_edges degree
            if len(same_edges_node_list) == 0:
                print(progress_count_in_degrees, max_degree, "no node")
                progress_count_in_degrees += 1
                continue

            # start this by first subset (and clone) the node dict
            specific_dict = {}
            self.build_specific_node_dict(same_edges_node_list, specific_dict)
            # then pass the neighbour list and the dict into recursion
            # first protein to remove is key, since all neighbours in
            # neighbour_list share this
            same_edges_node_list.sort()
            self.grouping_recursion(same_edges_node_list, specific_dict)

            print(progress_count_in_degrees, max_degree)
            progress_count_in_degrees += 1

    def collapse_graph_old(self) -> None:
        # for each node
        progress_count = 1
        for key, value in self.node_dict.items():

            # if the node has been delete, skip
            if key.get_first_id() in self.node_to_delete:
                progress_count += 1
                continue

            # code for skipping proteins
            # if isinstance(key, Protein):
            #     progress_count += 1
            #     continue

            # current neighbours: neighbours of the current node
            current_neighbours = value

            # if there are more than 2 neighbours
            if len(current_neighbours) >= 2:

                # actually reorganize the neighbours
                neighbours_reorganized = self.reorder_neighbours(
                    current_neighbours)

                # print(list(map(len, neighbours_reorganized)))

                # the other method of merging
                # not using anymore since it is a lot slower
                # it was used to check mergeability directly
                # if isinstance(key, Peptide):
                #     self.check_for_merging(neighbours_reorganized)
                #     continue

                # use recursion to group protein and then merge them
                # TODO, change to range (1, len(...))?
                for num_protein_mapped_to in range(len(neighbours_reorganized)):
                    # skip the first one, since that peptide in that neighbour
                    # list maps to exactly 0 protein
                    if num_protein_mapped_to == 0:
                        continue

                    neighbours_list = neighbours_reorganized[
                        num_protein_mapped_to]

                    # there are no peptide that map to num_protein_mapped_to protein
                    if len(neighbours_list) == 0:
                        continue

                    # specific_dict is just a subset (and clone) of node dict
                    specific_dict = {}
                    self.build_specific_node_dict(neighbours_list,
                                                  specific_dict)
                    # then pass the neighbour list and the dict into recursion
                    # first protein to remove is key, since all neighbours in
                    # neighbour_list share this
                    self.grouping_recursion(neighbours_list, specific_dict)

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
        :return: neighbours_reorganized [0, 0, 28] means of the current
        neighbours of this protein (which are peptide). 0 map to exactly 0
        protein, 0 map to exactly 1 protein, 28 map to exactly 2 proteins
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

    def build_specific_node_dict(self,
                                 same_degree_protein: List[Node],
                                 protein_peptide_dict: Dict[Node, List[Node]]) \
            -> None:
        """
        terminology assumes we are trying to merge peptides
        :param protein_peptide_dict:
        :param same_degree_protein is all the peptides that have the same number
        of protein that they map to
        :return type should be peptide map to a list of protein since I am only
        using this for collapsing peptides
        """

        for peptide in same_degree_protein:
            # is a neighbour and a protein
            neighbour_and_protein = self.node_dict.get(peptide, [])
            if len(neighbour_and_protein) == 0:
                print('something is wrong')
                exit('something is wrong')
            protein_peptide_dict[peptide] = sorted(neighbour_and_protein)

    def grouping_recursion(self, same_degree_protein: List[Node],
                           protein_peptide_dict: Dict[Node, List[Node]]):
        """terminology assume we are merging protein
        first build a dict based on the protein given, but without the
        protein that is these peptides have in common for sure (this is done by
        build_specific_node_dict method)

        second check how many protein do each of the peptides have left
        this should all be the same number since the initial peptide inputted
        map to the same number of protein

        Third if there are no protein left, then merge them
        else group them by their current first protein
        and pass each group with the current first protein
        into recursion
        """

        num_peptide = []

        # for each protein with the same degree
        for protein in same_degree_protein:
            # get its neighbouring peptides
            protein_list = protein_peptide_dict.get(protein)
            # and count
            num_peptide.append(len(protein_list))

        # the length of the set is the number of unique element in num_peptide (
        # i.e. how many peptides do these same_degree_protein connect to
        # which should be 1, meaning that they all have the same degree
        assert len(set(num_peptide)) < 2

        # this is the number of protein that not yet been used for grouping
        # it is initially the degree
        num_protein_all = list(set(num_peptide))[0]


        # TODO: test
        all_id = []
        for protein in same_degree_protein:
            accession = protein.get_first_id()
            if not 'DECOY' in accession and not accession.isnumeric():
                all_id.append(protein.get_first_id())

        if 'sp|Q9UI08-5|EVL_HUMAN' in all_id:
            print('before merging', print(num_protein_all))
            print('all id', all_id)

        # the protein in same_degree_protein all were part of a list that
        # all map to the same set of peptide

        # base case (we have gotten to final protein)
        # if there are no other peptide that it map to, then this means
        # this protein only map to this set
        # this means that all the protein with the same degree has been grouped
        # into groups where within each group all the protein's peptides are the same
        # each function call that end up in the base case only is given a set of
        # protein that share all its peptide
        if num_protein_all == 0:

            # TODO when it was sort() there was nothing wrong (problem)
            sorted_same_degree_protein = sorted(same_degree_protein)
            # merge the protein that share the same set of peptides
            # TODO: refactor and extract this as a method
            # choose the first object in the list as the one to merge to
            final_protein_group = sorted_same_degree_protein[0]

            # TODO: remove test
            if 'sp|Q9UI08-2|EVL_HUMAN' in final_protein_group.get_id():
                print('protein of interest', final_protein_group.get_id())
                print('interest first id', final_protein_group.get_first_id())
                for protein in sorted_same_degree_protein:
                    print('other protein', protein.get_id())
                    print('other protein first id', protein.get_first_id())

            # for all the protein in this set, skipping first one
            for index in range(1, len(sorted_same_degree_protein)):
                protein = sorted_same_degree_protein[index]
                self.delete_node(protein)

                id_list = protein.get_id()
                score = protein.get_score()

                final_protein_group.add_id(id_list)
                # this only keeps the max score of them two
                final_protein_group.add_score(score)

            if 'sp|Q9UI08-2|EVL_HUMAN' in final_protein_group.get_id():
                print('final protein group of interest', final_protein_group.get_id())

            # TODO supposed bottleneck
            # final_id_list = final_protein_group.get_id()
            # self.key_add_id(final_protein_group, final_id_list)

        # recursive case, not 0
        else:
            protein_same_next_peptide_dict = {}

            # for every entry in this group of protein
            for protein in same_degree_protein:

                # get that protein
                peptide = protein_peptide_dict.get(protein, [])

                assert peptide != []

                # sort their peptides
                sorted_peptide = sorted(peptide)
                next_same_peptide = sorted_peptide[0]

                # group the peptide by their next same protein
                protein_same_next_peptide_dict \
                    .setdefault(next_same_peptide, []).append(protein)

            # now in protein_same_next_peptide_dict, protein in each sublist
            # have the same first peptide

            # recursion
            # each protein list share the first peptide on the list
            for next_same_peptide, protein_list in protein_same_next_peptide_dict.items():
                for protein in protein_list:
                    # get each protein
                    protein = protein_peptide_dict.get(protein, [])
                    assert protein != []
                    # remove that protein
                    protein.remove(next_same_peptide)

                protein_list.sort()
                self.grouping_recursion(protein_list, protein_peptide_dict)

    """methods for step 2b: merging"""

    def check_for_merging_old(self, neighbours_reorganized: List[List[Node]]) \
            -> None:
        # all the protein neighbours in the same protein_neighbour_list
        # should have the same length
        for neighbour_list in neighbours_reorganized:
            for neighbour_pair in itertools.combinations(neighbour_list, 2):
                if (
                        # only if the 2 both neighbour have not been deleted yet
                        neighbour_pair[0].get_first_id() not in self.node_to_delete
                        and
                        neighbour_pair[1].get_first_id() not in self.node_to_delete
                        # then compare their neighbours
                ) and self.compare_neighbours_old(neighbour_pair):

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

    def compare_neighbours_old(self, node_pair: Tuple[Node, ...]) \
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

        # I don't know if python ordered list comparison for equality has
        # early-stop, such that when it reaches an index that is different
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

    def make_accession_object_dict(self) -> None:

        # for all node
        for node in self.node_dict:
            # if not deleted
            if (node.get_first_id() + node.get_target_decoy()) not in self.node_to_delete:
                # map its first id to the node object itself
                self.accession_object[node.get_first_id()] = node

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

        # for all neighbouring white node that is not deleted, explore them
        for current_neighbour in neighbour_list:
            if self.is_white(current_neighbour) \
                    and (current_neighbour.get_first_id() + current_neighbour.get_target_decoy()) \
                    not in self.node_to_delete:
                self.dfs(current_neighbour, a_component)

        self.set_explored(start_node)

        # after the nodes are blackened,
        # add the protein or peptide into the a_component
        # with its neighbours
        # but get the object that is key, instead of the value
        # through the dict accession_object

        component_node = self.accession_object[start_node.get_first_id()]

        # TODO test
        # it is in in
        if start_node.get_first_id() == 'sp|Q9UI08-2|EVL_HUMAN':
            print(component_node.get_id())

        if start_node.get_first_id() == 'sp|Q9UI08-5|EVL_HUMAN':
            print(component_node.get_id())

        if isinstance(component_node, Protein):
            a_component.add_protein(component_node, neighbour_list)
        elif isinstance(component_node, Peptide):
            a_component.add_peptide(component_node, neighbour_list)
