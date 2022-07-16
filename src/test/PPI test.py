import unittest


# class TestInitialization(unittest.TestCase):
#     def test_(self):
#
#         example_graph = ppi.Graph()
#         accession_List = [
#             "1", "2", "3", "4"
#         ]
#
#         self.assertEqual(True, False)
#
#


# noinspection PyArgumentList
class TestCollapse(unittest.TestCase):
    """this class contains test for method used in collapse
    """

    def setUp(self) -> None:
        self.graph_1 = ppi.Graph()
        self.pro1 = ppi.Protein(['1'], 1)
        self.pro2 = ppi.Protein(['2'], 1)
        self.pro3 = ppi.Protein(['3'], 1)
        self.pro4 = ppi.Protein(['4'], 1)
        self.pro5 = ppi.Protein(['5'], 1)
        self.pro6 = ppi.Protein(['6'], 1)
        self.pro7 = ppi.Protein(['7'], 1)
        self.pro8 = ppi.Protein(['8'], 1)
        self.pro9 = ppi.Protein(['9'], 1)

        self.pep1 = ppi.Peptide([1])
        self.pep2 = ppi.Peptide([2])
        self.pep3 = ppi.Peptide([3])
        self.pep4 = ppi.Peptide([4])
        self.pep5 = ppi.Peptide([5])
        self.pep6 = ppi.Peptide([6])
        self.pep7 = ppi.Peptide([7])
        self.pep8 = ppi.Peptide([8])
        self.pep9 = ppi.Peptide([9])
        self.pep10 = ppi.Peptide([10])

        self.graph_1.node_dict[self.pro1] = [self.pep4, self.pep3, self.pep7,
                                             self.pep9, self.pep8]
        self.graph_1.node_dict[self.pro2] = [self.pep8]
        self.graph_1.node_dict[self.pro3] = [self.pep6]
        self.graph_1.node_dict[self.pro4] = [self.pep2, self.pep1]
        self.graph_1.node_dict[self.pro5] = [self.pep4, self.pep8]
        self.graph_1.node_dict[self.pro6] = [self.pep2, self.pep6]
        self.graph_1.node_dict[self.pro7] = [self.pep1, self.pep5]
        self.graph_1.node_dict[self.pro8] = [self.pep8]
        self.graph_1.node_dict[self.pro9] = [self.pep2, self.pep1]

        self.graph_1.node_dict[self.pep1] = [self.pro7]
        self.graph_1.node_dict[self.pep2] = [self.pro4, self.pro6, self.pro9]
        self.graph_1.node_dict[self.pep3] = [self.pro1]
        self.graph_1.node_dict[self.pep4] = [self.pro1, self.pro5]
        self.graph_1.node_dict[self.pep5] = [self.pro7]
        self.graph_1.node_dict[self.pep6] = [self.pro3, self.pro6]
        self.graph_1.node_dict[self.pep7] = [self.pro1]
        self.graph_1.node_dict[self.pep8] = [self.pro1, self.pro5, self.pro2,
                                             self.pro8]
        self.graph_1.node_dict[self.pep9] = [self.pro1]
        self.graph_1.node_dict[self.pep10] = [self.pro4, self.pro9]

        # neighbours of pep8
        self.protein_list = [self.pro1, self.pro5, self.pro2, self.pro8]
        # neighbours of pro1
        self.peptide_list = [self.pep4, self.pep3, self.pep7, self.pep9,
                             self.pep8]

    def test_reorder_neighbours_1(self):

        """testing reorganized_protein_neighbours and peptides
        important thing to check about reorganize neighbours: that it does
        reorganize by the number of neighbours of the neighbours"""

        reformatted_protein_neighbours = self.graph_1.reorder_neighbours(
            self.protein_list)

        for i in range(1, len(reformatted_protein_neighbours)):
            protein_neighbour_list = reformatted_protein_neighbours[i]
            for protein_neighbour in protein_neighbour_list:
                num_neighbours = len(self.graph_1.node_dict[protein_neighbour])
                self.assertEqual(num_neighbours, i)

    def test_reorder_neighbours_2(self):

        reformatted_peptide_neighbours = self.graph_1.reorder_neighbours(
            self.peptide_list)

        for j in range(1, len(reformatted_peptide_neighbours)):
            peptide_neighbour_list = reformatted_peptide_neighbours[j]
            for peptide_neighbour in peptide_neighbour_list:
                num_neighbours = len(self.graph_1.node_dict[peptide_neighbour])
                self.assertEqual(num_neighbours, j)

    def test_compare_neighbours_1(self):
        self.assertTrue(
            self.graph_1.compare_neighbours_old((self.pep1, self.pep5))
        )

    def test_compare_neighbours_2(self):
        self.assertTrue(
            self.graph_1.compare_neighbours_old((self.pep3, self.pep9))
        )

    def test_compare_neighbours_3(self):
        self.assertFalse(
            self.graph_1.compare_neighbours_old((self.pep2, self.pep4))
        )

    def test_compare_neighbours_4(self):
        self.assertFalse(
            # two neighbours that is the same, one extra for 8
            self.graph_1.compare_neighbours_old((self.pep8, self.pep4))
        )

    def test_compare_neighbours_5(self):
        self.assertTrue(
            self.graph_1.compare_neighbours_old((self.pro2, self.pro8))
        )

    def test_compare_neighbours_6(self):
        self.assertTrue(
            self.graph_1.compare_neighbours_old((self.pro4, self.pro9))
        )

    def test_compare_neighbours_7(self):
        self.assertFalse(
            self.graph_1.compare_neighbours_old((self.pro3, self.pro6))
        )

    def test_compare_neighbours_8(self):
        self.assertTrue(
            self.graph_1.compare_neighbours_old((self.pep1, self.pep5))
        )

    def test_delete_node_1(self):

        self.graph_1.delete_node(self.pro8)
        self.assertTrue(
            self.pro8.get_first_id() in self.graph_1.node_to_delete
        )

    def test_delete_node_2(self):

        self.graph_1.delete_node(self.pro9)
        self.assertTrue(
            self.pro9.get_first_id() in self.graph_1.node_to_delete
        )

    def test_delete_node_3(self):

        self.graph_1.delete_node(self.pep7)
        self.assertTrue(
            self.pep7.get_first_id() in self.graph_1.node_to_delete
        )

    def test_delete_node_4(self):

        self.graph_1.delete_node(self.pep9)
        self.assertTrue(
            self.pep9.get_first_id() in self.graph_1.node_to_delete
        )

        # TODO: I need to write a test so that whenever any method go uses
        #  the delete node, it does not work


# noinspection PyArgumentList,PyUnresolvedReferences
class TestSeparate(unittest.TestCase):
    """this class contains test for method used in Separate
    """

    def setUp(self) -> None:
        self.graph_1 = ppi.Graph()
        self.pro_1 = ppi.Protein(['1'], 1)
        self.pro_28 = ppi.Protein(['2', '8'], 1)
        self.pro_3 = ppi.Protein(['3'], 1)
        self.pro_49 = ppi.Protein(['4', '9'], 1)
        self.pro_5 = ppi.Protein(['5'], 1)
        self.pro_6 = ppi.Protein(['6'], 1)
        self.pro_7 = ppi.Protein(['7'], 1)

        self.pep_3_7_9 = ppi.Peptide([3, 7, 9])
        self.pep_4 = ppi.Peptide([4])
        self.pep_8 = ppi.Peptide([8])
        self.pep_2 = ppi.Peptide([2])
        self.pep_6 = ppi.Peptide([6])
        self.pep_10 = ppi.Peptide([10])
        self.pep_1_5 = ppi.Peptide([1, 5])

        self.graph_1.node_dict[self.pro_1] = [self.pep_3_7_9, self.pep_4,
                                              self.pep_8]
        self.graph_1.node_dict[self.pro_28] = [self.pep_8]
        self.graph_1.node_dict[self.pro_3] = [self.pep_6]
        self.graph_1.node_dict[self.pro_49] = [self.pep_2, self.pep_10]
        self.graph_1.node_dict[self.pro_5] = [self.pep_4, self.pep_8]
        self.graph_1.node_dict[self.pro_6] = [self.pep_2, self.pep_6]
        self.graph_1.node_dict[self.pro_7] = [self.pep_1_5]

        self.graph_1.node_dict[self.pep_3_7_9] = [self.pro_1]
        self.graph_1.node_dict[self.pep_4] = [self.pro_1, self.pro_5]
        self.graph_1.node_dict[self.pep_8] = [self.pro_1, self.pro_5,
                                              self.pro_28]
        self.graph_1.node_dict[self.pep_2] = [self.pro_6, self.pro_49]
        self.graph_1.node_dict[self.pep_6] = [self.pro_3, self.pro_6]
        self.graph_1.node_dict[self.pep_10] = [self.pro_49]
        self.graph_1.node_dict[self.pep_1_5] = [self.pro_7]

        self.com_1 = ppi.Component()
        self.com_2 = ppi.Component()
        self.com_3 = ppi.Component()

        self.graph_1.dfs(self.pro_1, self.com_1)
        self.graph_1.dfs(self.pep_2, self.com_2)
        self.graph_1.dfs(self.pro_7, self.com_3)

        """then just check if the correct proteins and peptide are in the 
        respective a_component"""
        self.com_1_pro = self.com_1.protein_dict.keys()
        self.com_1_pep = self.com_1.peptide_dict.keys()
        self.com_2_pro = self.com_2.protein_dict.keys()
        self.com_2_pep = self.com_2.peptide_dict.keys()
        self.com_3_pro = self.com_3.protein_dict.keys()
        self.com_3_pep = self.com_3.peptide_dict.keys()

    def test_DFS(self):
        """test whether or not dfs put the correct protein and peptide into the
        correct connected components
        """
        self.assertEqual(len(self.com_1_pro), 3)
        self.assertIn(self.pro_1, self.com_1_pro)
        self.assertIn(self.pro_28, self.com_1_pro)  # failed
        self.assertIn(self.pro_5, self.com_1_pro)

        self.assertEqual(len(self.com_1_pep), 3)
        self.assertIn(self.pep_3_7_9, self.com_1_pep)
        self.assertIn(self.pep_4, self.com_1_pep)
        self.assertIn(self.pep_8, self.com_1_pep)

        self.assertEqual(len(self.com_2_pro), 3)
        self.assertIn(self.pro_3, self.com_2_pro)
        self.assertIn(self.pro_49, self.com_2_pro)
        # failed, because forgot between pro49 and pep2
        self.assertIn(self.pro_6, self.com_2_pro)

        self.assertEqual(len(self.com_2_pep), 3)
        self.assertIn(self.pep_2, self.com_2_pep)
        self.assertIn(self.pep_6, self.com_2_pep)
        self.assertIn(self.pep_10, self.com_2_pep)
        # failed, because forgot between pro49 and pep2

        self.assertEqual(len(self.com_3_pro), 1)
        self.assertIn(self.pro_7, self.com_3_pro)

        self.assertEqual(len(self.com_3_pep), 1)
        self.assertIn(self.pep_1_5, self.com_3_pep)

        # component_list = []
        #
        # ppi.separate(graph_1, component_list)
        #
        # self.assertIn(component_1 in compo)


# noinspection PyArgumentList,PyUnresolvedReferences
class TestPPI(unittest.TestCase):

    def test_reduce(self):
        component_1 = ppi.Component()
        component_2 = ppi.Component()
        component_3 = ppi.Component()

        pro_1 = ppi.Protein(['1'], 1)
        pro_2_8 = ppi.Protein(['2', '8'], 1)
        pro_5 = ppi.Protein(['5'], 1)
        pro_3 = ppi.Protein(['3'], 1)
        pro_4_9 = ppi.Protein(['4', '9'], 1)
        pro_6 = ppi.Protein(['6'], 1)
        pro_7 = ppi.Protein(['7'], 1)

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

        self.assertEqual(component_1.make_protein_list(),
                         [pro_1.get_sqlite_id()])

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

        self.assertEqual(component_3.make_protein_list(),
                         [pro_7.get_sqlite_id()])

        self.assertTrue(pep_1_5.is_covered())


if __name__ == '__main__':
    unittest.main()
