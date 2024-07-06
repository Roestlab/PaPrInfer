from typing import Dict, List, Optional, Tuple, Any
from node import Protein, Peptide
import random


class Component:
    _protein_dict: Dict[Protein, List[Peptide]]
    _peptide_dict: Dict[Peptide, List[Protein]]
    _covered_peptide: Dict[str, str]
    _target_dict: Dict[str, str]

    def __init__(self):
        self._protein_dict = {}
        self._peptide_dict = {}
        self._covered_peptide = {}
        self._target_dict = {}

    def add_peptide(self, current_peptide: Peptide,
                    neighbours) -> None:

        self._peptide_dict[current_peptide] = neighbours

    def add_protein(self, current_protein: Protein,
                    neighbours) -> None:

        self._protein_dict[current_protein] = neighbours

    def make_protein_list(self) -> \
            List[Tuple[List[str], float, str, Any]]:

        component_min_pro_list = []
        claimed_peptide = {}

        if len(self._protein_dict) == 0:
            return []

        while not self.all_component_peptides_covered():

            # find unselected protein with the most edges to uncovered peptide
            currently_selected_protein = self.find_most_uncovered_protein()

            # this should not happen, but sometimes it returns none
            if currently_selected_protein is None:
                break

            # set as selected, so subset is unselected ones
            currently_selected_protein.set_selected()

            # do 3 things, and check 2 things
            # 1. do set uncovered peptide as covered
            # 2. do find max score among uncovered peptide
            # 3. do map those uncovered peptide to the protein
            # 1. check that protein and peptide is same type
            # 2. check that max score after is not negative infinity
            max_score = float('-inf')
            could_cover_peptide = self._protein_dict[currently_selected_protein]

            for current_peptide in could_cover_peptide:
                # decoy should only match up with decoy, same for target
                same_type = currently_selected_protein.get_target_decoy() == current_peptide.get_target_decoy()
                assert same_type

                not_covered = current_peptide.get_first_id() not in self._covered_peptide
                if not_covered:
                    current_score = current_peptide.get_score()
                    if current_score > max_score:
                        max_score = current_score
                    self._covered_peptide[current_peptide.get_first_id()] = ''
                    # TODO: I do have extend now, putting all peptide to the
                    #  same list, but if we need to have peptide group in osw, then
                    #  need to change to append
                    claimed_peptide.setdefault(
                        currently_selected_protein.get_first_id(), []).extend(
                            current_peptide.get_id()
                        )

            assert max_score != float('-inf')

            currently_selected_protein.set_score(max_score)
            component_min_pro_list.append(currently_selected_protein)


        component_accession_list = []
        # component_min_pro_list is a list of protein that are selected
        # each protein contain their own accession list (list of str)
        # and score
        # claimed peptide is a dict of protein id to list of peptide id
        for current_protein in component_min_pro_list:
            component_accession_list.append(
                (
                    current_protein.get_id(), current_protein.get_score(),
                    current_protein.get_target_decoy(),
                    claimed_peptide[current_protein.get_first_id()]
                )
            )

        return component_accession_list

    def find_most_uncovered_protein(self) -> Optional[Protein]:
        """
        find an unselected protein that links to the most uncovered peptides
        :return:
        """

        most_edges_protein = random.choice(list(self._protein_dict.keys()))

        # record its number of edges
        most_uncovered = self.find_num_uncovered_peptides(most_edges_protein)

        # loop over the proteins of a_component
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

        if most_uncovered == 0:
            return None


        return most_edges_protein

    def find_num_uncovered_peptides(self, current_protein: Protein) -> int:
        num_uncovered_peptides = 0
        peptide_neighbours = self._protein_dict[current_protein]
        for current_peptide in peptide_neighbours:
            # if it is covered or if it not same type, skip
            not_covered = current_peptide.get_first_id() not in self._covered_peptide
            same_type = current_peptide.get_target_decoy() == current_protein.get_target_decoy()
            if not_covered and same_type:
                num_uncovered_peptides += 1
            else:
                continue
        return num_uncovered_peptides

    def all_component_peptides_covered(self) -> bool:
        all_covered = True

        # find if any peptide is not yet covered
        for peptide in self._peptide_dict:
            if peptide.get_first_id() not in self._covered_peptide:
                all_covered = False
                break

        return all_covered
