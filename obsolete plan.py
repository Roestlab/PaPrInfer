# THOUGHTS ONLY
    # do I want adjacency list or adjacency matrix?
    # I think the graph is dense, that is |E| is close to |V|^2
    # so there is |V| nodes, half is protein, half is peptides
    # for each protein node, there can be at most 'num_pepties' node
    # that means |E| = |V|/2 * |V|/2 = (|V|/2)^2 = (|V|^2)/4

    # but that is the worse case, I dont think it what is normally though?

    # I think I am going to need 2 pieces of knowledge
    # Typically the list of protein is all possible proteins right, so what fraction
    # of that is present in the sample
    # and also on average, how many peptide does each protein fragment into
    # though, i have a feeling, it is kind of hard to know these without this algorithm

    # I am thinking that, it is at most a quarter of V sqaured, maybe typically it would
    # be much less?
    # and plus, with adjacency list, I can do a quick check for length in the collapse step
    # and also, adjacency matrix allow fast look up of edges for any 2 given vertices\
    # but I need to to check if these do that like every possible edge of 1 vertices
    # then compare it with another vertex

    # EXECUTION
    # first take the list of proteins, and list of peptides
    # for each protein, make a protein node (node's __init__)
    # for each peptides, make a peptide node (node's __init__)
    # add each protein into a dict with node as key and empty list (initially)
    # also do that with each peptide (graph's __init__)

    # Edges (i am not sure)
    # call graph.make_edges
    # maybe loop through the protein dict, for each protein node, loop through the
    # peptide dict, check if it contains parts of the peptide sequence (compare)
    # maybe there is a function in OpenMS that?
    # so if it does, then for the corresponding list, add the peptide node to it
    # also for that peptide node's corresponding list, add the protein node to it.


# compare (needs 1 argument: the other node, a peptide. return bool)
# compare 2 nodes's sequence
# check if the peptide's seq have a match in the protein's se






# # Iterate over its PeptideEvidences objects
            # for evidence in peptide_hit.getPeptideEvidences():
            #     # and find the protein corresponding to evidence, by accession
            #     for current_protein in self.protein_dict.keys():
            #         if current_protein.get_accession() == evidence.getProteinAccession():
            #             # then append it to the protein edges of this peptide
            #             self.peptide_dict[current_peptide].append(current_protein)
            #             # also append this peptide to the peptide edges of this protein
            #             self.protein_dict[current_protein].append(current_peptide)
            #             break

            # actually what if there is a protein in peptideEvidence, but not
            # in the protein files given?
            # i.e. the peptide has a edge to a protein that is "not part of the
            # graph"
            # I mean I think I just add that protein to the protein_dict
            # since the purpose of the algorithm is trying to generate a
            # minimal protein list, that is trying to know what are the most likely
            # protein that was in the sample based on the peptide identified

            # because it is obvious if it the reverse/opposite situation: that is
            # a protein in the protein file does not match any peptide
            # i.e. a protein in the graph has no edges to any peptide
            # obviously you just not touch that protein in the subsequent steps
                # collapse wont affect it
                # neither will seperate
                # nor reduce


# def dfs_protein(self, start_node: Protein, a_component: Component) -> None:
#     start_node.set_discovered()
#
#     a_component.add_protein()
#
#     # for all white neighbouring peptides, explore them
#     for current_neighbouring_peptide in self.protein_dict[start_node]:
#         if current_neighbouring_peptide.get_color() == 0:
#             self.dfs_peptide(current_neighbouring_peptide, a_component)
#
#     start_node.set_explored()



# take the osw file (which is SQlite) and pass it through PercolatorAdapter
# and select output as idXML
# PPI can then read this idXML file
# after I am done, write into a SQLite file

# it seems like that is not what I am meant to do??
