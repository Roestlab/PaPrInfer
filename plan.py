

# first read the protein id that is lower than the threshold
# for each protein id
    # read the protein accessions one by one from the PROTEIN table
    # make protein objects out of them
    # protein will take 2 arguments
    # one is the first protein accession (which serves as the unique immutable dict key)
    # and have a list of protein accession
    # and have protein id (which will not be unique)

# for peptides, read the peptide id that is lower than threshold

# as for the edges, read from protein peptide mapping
# go through the protein keys (and using their ids) and find all corresponding peptides
# check if they are indeed peptide keys, then add them into the values

# just like before go through the peptide keys and find all corresponding protein id
# but now for each, protein id, i need to parse the string into a list of substrings




# note that change the *, will let me select specific columns
# like c.execute("SELECT ID FROM peptide") vs
# c.execute("SELECT * FROM peptide")

# class node
    # === attributes ===
        # name
        # color
            # 0 = white = undiscovered + unexplored,
            # 1 = gray = discovered + unexplored,
            # 2 = black = discovered + explored
            # initially, should be white
        # same
        # component_number: the connected a_component that this node belong to.
            # Before doing step 3 (and after the node is created) this attribute
            # does not mean anything

    # === methods ===
        # __init__ (needs input of name, seq)
            # set up name and seq
            # color = 0
            # a_component = -1
            # same = false

# class protein extends nodes
    # === attributes ===
        # selected(a)
            # this is for step 4: reduce: it says whether or not the protein
            # in question has been selected for the minimal protein list)
        # accession number
    # === methods ===
        # compare_protein_connections (2 arguments needed: first_protein, second
        # protein
            # find first_protein (a key)'s corresponding value
            # in the protein dict
            # also with second_protein
            # they should both be a list, sort that list
            # check if they are equal
            # return the result (which is bool)

# I AM NOT SURE THIS
# but I remember a list of objects are equal, if the class objects in
# same index are equal


# class peptides extends nodes
    # === attributes ===
        # PSM data
        # covered
            # this is for 4: reduce, it say whether or not
            # the peptide is covered by a protein that is selected)
    # === methods ===
        # compare_peptide connections

# class graph
# it is two dict, each with node as key, and a list of nodes that this nodes has
# edges to, as value
# one for protein, one for peptides

    # === attributes ===
        # protein_dict
            # dict consisting of protein node as key
            # and a list of peptide nodes as values
        # peptide_dict
            # dict consisting of peptide node as key
            # and a list of protein nodes as values

    # === methods ===
        # make edges (0 arguments)

        # collpase_protein (0 arguments)
            # loop through peptide dict's keys
                # for each peptide node, if it has more than or equal to 2 edges
                # (to a protein or meta-protein) loop through the protein
                    # for each protein
                    # call protein.compare_protein_connections
                    # if made up of same pepties
                        # exit the loop
                # merge the 2 proteins
                # if it still has more than 2 edges
                    # do that again

        # collpase_peptide (0 arguments)
            # similar logic to collpase_protein, but peptide and protein switches

        # DFS: input a node, output a connected a_component
            # loop over the protein dict
            # if node is not black
                # call explore (which is a recursive function)
            # set node to black

        # explore
            # set node to gray
            # loop over the nodes that it has edges to
            # call discover on the node

        # make_protein_list

        # find_most_edges_protein

        # all_peptides_covered

# class a_component extends graph
    # === methods ===
        # make_protein_list: input a connected a_component, output a minimal protein list


# if __name__ == '__main__':

# 1. initialize

    # I am not sure about the ProteinIdentification object's data structure
    # but it seems like it contain mutiple protein hits objects
    # same thing applies to PeptideIdentification and peptide hits

    # Protein inference algorithm is not meant to search the match between
    # peptide and protein through their sequence, right?
    # I can just use peptide evidence object in the peptide hits object which
    # is in peptide identification objects?

    # each peptide evidence object correspond to one protein

    # make new object of graph

    # so retrive all protein hits object
    # read their accession number
    # and make a protein node
    # append it to the graph's protein dict
    # do that for all protein hits objects

    # for all peptide hits
    # read their PSM data (Peptide-spectrum match)
    # and make a peptide node
    # append it to the graph's peptide dict
    # for each peptide hits, loop over its peptide evidence objects
    # then use get protein accession method to get a string
    # then search for that string in protein dict
    # add that protein node into the value of this peptide



# 2. collapse
    # call collpase_protein of class graph. which:

    # and then call collpase_peptide (which uses similar logic)

    # give a meta node number (at initialization)
    # if they have the same number
    # then they belong to the same meta node


# hold on, what if one peptide neighbour was remove from
# the graph, by collapse_peptide, but that peptide was
# compared in another iteration of the inner for loop
# by compare_peptide_neighbours?
# ok, so I will do a fix
# currently compare_peptide_neighbours tries to access the
# dictionary, it may give an keyerror exception
# but i want it to return none
# so use the get method in dict

# nvm, this fix does not work, because i still need to compare that removed
# node

# but it return none, it is possible (I THINK IT IS POSSIBLE??)
# that both gets return none
# which will return true

# also add sort in the comparing of the compare peptide neighbors




# 3. separate
    # call DFS from peptide

    # I guess i can for each a_component number
    # make a a_component class object
    # and only if there is more than 1 node in the a_component, then call make
    # protein list

    # all nodes are initially white
    # you pick a start node: it become grey
    # if its neighbour is not black

    # make a a_component from while doing dfs

    # be able to check if all peptides are covered



# 4. reduce

    # make a list (minimal_protein_list), and then for each make_protein_list, append the result

    # make_protein_list
    # protein_list = a new empty list
    # While not (all_peptides_covered)
        # find_protein_most_edges
        # append that protein into the protein_List
        # for loop over the list of the protein edges
            # and set those peptide's covered attribute to true


    # find_protein_most_edges
        # most = 0
        # protein_most_edges = the first one in the protein dict
        # for loop over the protein dict
            # check if the protein has its "selected" set to true
                # if it does skip a loop
            # find the length of their list of edges (num_edges)
            # if bigger than most, then set protein_most_edges = this one

    # all_peptides_covered
        # flag = true
        # for loop over the peptide dict
        # if peptide is not flagged as covered, then
        # then flag = false
        # return flag
    # do i need to do this, or it is better that I somehow do this while flagging
    # peptides as covered


    # print out that list (minimal_protein_list)





# problem, currently reorganize_protein_neighbours reorganizes it for all protein neighrbours
