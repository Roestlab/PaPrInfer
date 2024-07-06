from __future__ import annotations

from abc import abstractmethod
from typing import List, Tuple


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

    @abstractmethod
    def add_score(self, other_score: float):
        pass

    @abstractmethod
    def get_score(self):
        pass

    @abstractmethod
    def get_target_decoy(self):
        pass


class Protein(Node):
    _selected: bool
    _id_list: List[str]
    _first_id: str
    _sqlite_ids: List[str]
    _first_sqlite_id: str
    _score: float
    _decoy: int

    """initially, the protein group only has 1 accession number in it
    first accession serves as the immutable hash value for the dict key
    the protein_ids is just for making edges
    this is called protein, but it is actually referring to protein groups
    so initially it is a protein group with 1 protein"""

    def __init_subclass__(cls, **kwargs) -> None:
        pass

    def __init__(self, id_list: List[str], protein_sqlite_id: str, decoy: int) -> None:
        super().__init__()
        self._selected = False
        self._id_list = id_list
        self._first_id = id_list[0]
        self._first_sqlite_id = protein_sqlite_id
        self._sqlite_ids = [protein_sqlite_id]
        self._score = float('-inf')
        self._decoy = decoy

    def __hash__(self) -> int:
        return hash(self.get_first_id() + str(self._decoy))

    # bottleneck
    def get_both_first_accession(self, other: Protein) -> Tuple[str, str]:
        this_acc = self.get_first_id() + self.get_target_decoy()

        that_acc = other.get_first_id() + self.get_target_decoy()

        return this_acc, that_acc

    # bottleneck
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

    # TODO: bottleneck
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

    def add_sqlite_id(self, protein_sqlite_id: str) -> None:
        self._sqlite_ids.append(protein_sqlite_id)

    def get_first_sqlite_id(self) -> str:
        return self._first_sqlite_id

    def set_score(self, score: float):
        self._score = score

    def add_score(self, other_score: float):
        """because protein group's score depends on the peptide it claims,
        so it is not changed during the merging of proteins """
        pass

    def get_score(self) -> float:
        return self._score

    def get_target_decoy(self) -> str:
        return str(self._decoy)


class Peptide(Node):
    _covered: bool
    _first_id: str
    _ids: List[str]
    _score: float
    _decoy: int

    """initially, the protein group only has 1 id in it"""

    def __init__(self, peptide_id_list: List[str], score: float,
                 decoy: int) -> None:
        """
        When initializing, the peptide_id_list only has 1 peptide id
        """
        super().__init__()
        self._covered = False
        self._ids = peptide_id_list
        self._first_id = peptide_id_list[0]
        self._score = score
        self._decoy = decoy

    def __hash__(self) -> int:
        return hash(self.get_first_id() + str(self._decoy))

    # bottleneck
    def get_both_ids(self, other: Peptide) -> Tuple[str, str]:
        this_id = self.get_first_id() + self.get_target_decoy()

        that_id = other.get_first_id() + self.get_target_decoy()

        return this_id, that_id

    # bottleneck
    def __eq__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id == that_id

    def __ne__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id != that_id

    def __lt__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id < that_id

    def __le__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id <= that_id

    def __gt__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id > that_id

    def __ge__(self, other: Peptide) -> bool:
        this_id, that_id = self.get_both_ids(other)
        return this_id >= that_id

    def is_covered(self) -> bool:
        return self._covered

    def set_covered(self) -> None:
        self._covered = True

    # bottleneck
    def get_first_id(self) -> str:
        return self._first_id

    def get_id(self) -> List[str]:
        return self._ids

    def add_id(self, id_list: List[str]):
        self._ids.extend(id_list)

    def add_score(self, other_score: float):
        """
        only keep the bigger score
        """
        self._score = max(self._score, other_score)

    def get_score(self) -> float:
        return self._score

    def get_target_decoy(self) -> str:
        return str(self._decoy)
