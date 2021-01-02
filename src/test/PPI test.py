import unittest
from src import ppi


class TestInitialization(unittest.TestCase):
    def test_something(self):

        example_graph = ppi.Graph()
        accession_List = [
            "1", "2", "3", "4"
        ]

        self.assertEqual(True, False)


class TestCollpase(unittest.TestCase):
    pass


class TestSeperate(unittest.TestCase):
    pass

class TestPPI(unittest.TestCase):

    def test_collapse(self):
        graph_1 = ppi.Graph()
        pro1 = ppi.Protein([1])
        pro2 = ppi.Protein([2])
        pro3 = ppi.Protein([3])
        pro4 = ppi.Protein([4])
        pro5 = ppi.Protein([5])
        pro6 = ppi.Protein([6])
        pro7 = ppi.Protein([7])
        pro8 = ppi.Protein([8])
        pro9 = ppi.Protein([9])

        pep1 = ppi.Peptide([1])
        pep2 = ppi.Peptide([2])
        pep3 = ppi.Peptide([3])
        pep4 = ppi.Peptide([4])
        pep5 = ppi.Peptide([5])
        pep6 = ppi.Peptide([6])
        pep7 = ppi.Peptide([7])
        pep8 = ppi.Peptide([8])
        pep9 = ppi.Peptide([9])
        pep10 = ppi.Peptide([10])

        graph_1._protein_dict[pro1] = [pep4, pep3, pep7, pep9, pep8]
        graph_1._protein_dict[pro2] = [pep8]
        graph_1._protein_dict[pro3] = [pep6]
        graph_1._protein_dict[pro4] = [pep2, pep10]
        graph_1._protein_dict[pro5] = [pep4, pep8]
        graph_1._protein_dict[pro6] = [pep2, pep6]
        graph_1._protein_dict[pro7] = [pep1, pep5]
        graph_1._protein_dict[pro8] = [pep8]
        graph_1._protein_dict[pro9] = [pep2, pep10]

        graph_1._peptide_dict[pep1] = [pro7]
        graph_1._peptide_dict[pep2] = [pro4, pro6, pro9]
        graph_1._peptide_dict[pep3] = [pro1]
        graph_1._peptide_dict[pep4] = [pro1, pro5]
        graph_1._peptide_dict[pep5] = [pro7]
        graph_1._peptide_dict[pep6] = [pro3, pro6]
        graph_1._peptide_dict[pep7] = [pro1]
        graph_1._peptide_dict[pep8] = [pro1, pro5, pro2, pro8]
        graph_1._peptide_dict[pep9] = [pro1]
        graph_1._peptide_dict[pep10] = [pro4, pro9]

        """testing reorganized_protein_neighbours and peptides
        important thing to check about reorganize neighbours: that it does
        reorganize by the number of neighbours of the neighbours"""

        # neighbours of pep8
        protein_list = [pro1, pro5, pro2, pro8]
        # neighbours of pro1
        peptide_list = [pep4, pep3, pep7, pep9, pep8]

        reformatted_protein_neighbours = graph_1.reorder_neighbours(protein_list)
        reformatted_peptide_neighbours = graph_1.reformat_pep(peptide_list)

        for i in range(1, len(reformatted_protein_neighbours)):
            protein_neighbour_list = reformatted_protein_neighbours[i]
            for protein_neighbour in protein_neighbour_list:
                num_neighbours = len(graph_1._protein_dict[protein_neighbour])
                self.assertEqual(num_neighbours, i)

        for j in range(1, len(reformatted_peptide_neighbours)):
            peptide_neighbour_list = reformatted_peptide_neighbours[j]
            for peptide_neighbour in peptide_neighbour_list:
                num_neighbours = len(graph_1._peptide_dict[peptide_neighbour])
                self.assertEqual(num_neighbours, j)

        """test compare_peptide_neighbours"""

        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep1, pep5))
        )
        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep3, pep7))
        )
        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep7, pep9))
        )
        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep7, pep3))
        )
        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep3, pep9))
        )
        self.assertFalse(
            graph_1.compare_peptide_neighbours((pep2, pep4))
        )
        self.assertFalse(  # two neighbours that is the same, one extra for 8
            graph_1.compare_peptide_neighbours((pep8, pep4))
        )
        self.assertTrue(
            graph_1.compare_protein_neighbours((pro2, pro8))
        )
        self.assertTrue(
            graph_1.compare_protein_neighbours((pro4, pro9))
        )
        self.assertFalse(
            graph_1.compare_protein_neighbours((pro3, pro6))
        )
        self.assertTrue(
            graph_1.compare_peptide_neighbours((pep1, pep5))
        )

        """test collapse_peptide"""

        graph_1.delete_protein(pro8)
        # it is not in the protein dict
        self.assertNotIn(pro8, graph_1._protein_dict.keys())
        # nor is it the list of neighbours of pep8 (which used to be one of
        # the neighbours of pro8
        self.assertNotIn(pro8, graph_1._peptide_dict[pep8])

        graph_1.delete_protein(pro9)
        self.assertNotIn(pro9, graph_1._protein_dict.keys())
        self.assertNotIn(pro8, graph_1._peptide_dict[pep2])
        self.assertNotIn(pro8, graph_1._peptide_dict[pep10])

        graph_1.delete_peptide(pep7)
        self.assertNotIn(pep7, graph_1._peptide_dict.keys())
        self.assertNotIn(pep7, graph_1._protein_dict[pro1])

        graph_1.delete_peptide(pep9)
        self.assertNotIn(pep9, graph_1._protein_dict.keys())
        self.assertNotIn(pep7, graph_1._protein_dict[pro1])

    def test_separate(self):
        graph_1 = ppi.Graph()
        pro_1 = ppi.Protein([1])
        pro_2_8 = ppi.Protein([2, 8])
        pro_5 = ppi.Protein([5])
        pro_3 = ppi.Protein([3])
        pro_4_9 = ppi.Protein([49])
        pro_6 = ppi.Protein([6])
        pro_7 = ppi.Protein([7])

        pep_3_7_9 = ppi.Peptide([3, 7, 9])
        pep_4 = ppi.Peptide([4])
        pep_8 = ppi.Peptide([8])
        pep_2 = ppi.Peptide([2])
        pep_6 = ppi.Peptide([6])
        pep_10 = ppi.Peptide([10])
        pep_1_5 = ppi.Peptide([1, 5])

        graph_1._protein_dict[pro_1] = [pep_3_7_9, pep_4, pep_8]
        graph_1._protein_dict[pro_2_8] = [pep_8]
        graph_1._protein_dict[pro_3] = [pep_6]
        graph_1._protein_dict[pro_4_9] = [pep_2, pep_10]
        graph_1._protein_dict[pro_5] = [pep_4, pep_8]
        graph_1._protein_dict[pro_6] = [pep_2, pep_6]
        graph_1._protein_dict[pro_7] = [pep_1_5]

        graph_1._peptide_dict[pep_3_7_9] = [pro_1]
        graph_1._peptide_dict[pep_4] = [pro_1, pro_5]
        graph_1._peptide_dict[pep_8] = [pro_1, pro_5, pro_2_8]
        graph_1._peptide_dict[pep_2] = [pro_6, pro_4_9]
        graph_1._peptide_dict[pep_6] = [pro_3, pro_6]
        graph_1._peptide_dict[pep_10] = [pro_4_9]
        graph_1._peptide_dict[pep_1_5] = [pro_7]

        com_1 = graph_1.dfs(pro_1)
        com_2 = graph_1.dfs(pep_2, )
        com_3 = graph_1.dfs(pro_7, )

        """then just check if the correct proteins and peptide are in the 
        respective a_component"""
        com_1_pro = com_1._protein_dict.keys()
        com_1_pep = com_1._peptide_dict.keys()
        com_2_pro = com_2._protein_dict.keys()
        com_2_pep = com_2._peptide_dict.keys()
        com_3_pro = com_3._protein_dict.keys()
        com_3_pep = com_3._peptide_dict.keys()

        self.assertEqual(len(com_1_pro), 3)
        self.assertIn(pro_1, com_1_pro)
        self.assertIn(pro_2_8, com_1_pro) # failed
        self.assertIn(pro_5, com_1_pro)

        self.assertEqual(len(com_1_pep), 3)
        self.assertIn(pep_3_7_9, com_1_pep)
        self.assertIn(pep_4, com_1_pep)
        self.assertIn(pep_8, com_1_pep)

        self.assertEqual(len(com_2_pro), 3)
        self.assertIn(pro_3, com_2_pro)
        self.assertIn(pro_4_9, com_2_pro) # failed, because forogt between pro49 and pep2
        self.assertIn(pro_6, com_2_pro)


        self.assertEqual(len(com_2_pep), 3)
        self.assertIn(pep_2, com_2_pep)
        self.assertIn(pep_6, com_2_pep)
        self.assertIn(pep_10, com_2_pep) # failed, because forogt between pro49 and pep2

        self.assertEqual(len(com_3_pro), 1)
        self.assertIn(pro_7, com_3_pro)

        self.assertEqual(len(com_3_pep), 1)
        self.assertIn(pep_1_5, com_3_pro)

        # component_list = []
        #
        # ppi.separate(graph_1, component_list)
        #
        # self.assertIn(component_1 in compo)

    def test_reduce(self):

        component_1 = ppi.Component()
        component_2 = ppi.Component()
        component_3 = ppi.Component()

        pro_1 = ppi.Protein([1])
        pro_2_8 = ppi.Protein([2, 8])
        pro_5 = ppi.Protein([5])
        pro_3 = ppi.Protein([3])
        pro_4_9 = ppi.Protein([4, 9])
        pro_6 = ppi.Protein([6])
        pro_7 = ppi.Protein([7])

        pep_3_7_9 = ppi.Peptide([3, 7, 9])
        pep_4 = ppi.Peptide([4])
        pep_8 = ppi.Peptide([8])
        pep_2 = ppi.Peptide([2])
        pep_6 = ppi.Peptide([6])
        pep_10 = ppi.Peptide([10])
        pep_1_5 = ppi.Peptide([15])

        component_1.add_protein(pro_1, [pep_3_7_9, pep_4, pep_8])
        component_1.add_protein(pro_2_8, [pep_8])
        component_1.add_protein(pro_5, [pep_8, pep_4])
        component_1.add_peptide(pep_3_7_9, [pro_1])
        component_1.add_peptide(pep_4, [pro_1, pro_5])
        component_1.add_peptide(pep_8, [pro_1, pro_5, pro_2_8])

        component_2.add_protein(pro_3, [pep_6])
        component_2.add_protein(pro_4_9, [pep_2, pep_10])
        component_2.add_protein(pro_6, [pep_2, pep_6])
        component_2.add_peptide(pep_2, [pro_4_9, pro_6])
        component_2.add_peptide(pep_6, [pro_3, pro_6])
        component_2.add_peptide(pep_10, [pro_4_9])

        component_3.add_protein(pro_7, [pep_1_5])
        component_3.add_peptide(pep_1_5, [pro_7])

        """a_component 1"""
        self.assertEqual(component_1.find_num_uncovered_peptides(pro_1), 3)
        self.assertEqual(component_1._find_most_uncovered_protein(), pro_1)

        self.assertEqual(component_1.make_protein_list(), [pro_1.get_sqlite_id()])

        # now the peptide neighbours of pro1 should be set as covered
        self.assertTrue(pep_3_7_9.is_covered())
        self.assertTrue(pep_4.is_covered())
        self.assertTrue(pep_8.is_covered())

        """a_component 2"""
        self.assertEqual(component_2.find_num_uncovered_peptides(pro_4_9), 2)

        self.assertEqual(component_2.find_num_uncovered_peptides(pro_6), 2)

        self.assertTrue(
            component_2._find_most_uncovered_protein() == pro_4_9 or
            component_2._find_most_uncovered_protein() == pro_6
        )
        # TODO protein 6 returns

        # since pro3 only has 1 peptide neighbor, it should have the same
        # number of uncovered peptide before and after making protein list
        self.assertEqual(component_2.find_num_uncovered_peptides(pro_3), 1)

        self.assertEqual(
            component_2.make_protein_list().sort(),
            [pro_4_9.get_sqlite_id(), pro_6.get_sqlite_id()].sort()
        )

        self.assertTrue(pep_2.is_covered())
        self.assertTrue(pep_6.is_covered())
        self.assertTrue(pep_10.is_covered())

        """a_component 3"""
        self.assertEqual(component_3.find_num_uncovered_peptides(pro_7), 1)
        self.assertEqual(component_3._find_most_uncovered_protein(), pro_7)

        self.assertEqual(component_3.make_protein_list(), [pro_7.get_sqlite_id()])

        self.assertTrue(pep_1_5.is_covered())


if __name__ == '__main__':
    unittest.main()
