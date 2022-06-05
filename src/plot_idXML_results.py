
import sys
import matplotlib.pyplot as plt
from pyopenms import IdXMLFile


"""this is plotting qvalues and pp distributions of epifany result with
respect to target and decoy, and then modified to plot pp distribution of
protein inference"""

if __name__ == "__main__":
    prot_ids = []
    pep_ids = []
    epifany_file = sys.argv[1]

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    score_type = "None"

    target_epifany_q_values = []
    target_epifany_pp = []
    decoy_epifany_q_values = []
    decoy_epifany_pp = []

    all_proteins = []

    for protein_id in prot_ids:

        score_type = protein_id.getScoreType()

        print(protein_id.getScoreType())

        for hit in protein_id.getHits():

            if hit.getMetaValue("target_decoy") == "target":
                if score_type == "q-value":
                    target_epifany_q_values.append(hit.getScore())
                elif score_type == "Posterior Probability":
                    target_epifany_pp.append(hit.getScore())
                    # this would be for pep
                    # decoy_epifany_pp.append(1 - hit.getScore())

            elif hit.getMetaValue("target_decoy") == "decoy":
                if score_type == "q-value":
                    decoy_epifany_q_values.append(hit.getScore())
                elif score_type == "Posterior Probability":
                    decoy_epifany_pp.append(hit.getScore())
                    # this would be for pep
                    # decoy_epifany_pp.append(1 - hit.getScore())

            all_proteins.append(hit.getAccession())

    print(len(all_proteins))

    # the protein inference gives me negative infinity values
    # change -inf to 1
    # target_epifany_pp_cleaned = [1 if x == float('-inf') else x for x in target_epifany_pp]
    # decoy_epifany_pp_cleaned = [1 if x == float('-inf') else x for x in decoy_epifany_pp]

    # remove -inf
    target_epifany_pp_cleaned = [i for i in target_epifany_pp if i != float('-inf')]
    decoy_epifany_pp_cleaned = [i for i in target_epifany_pp if i != float('-inf')]

    # note since even though the bins are auto, they could be different, if one bin looks
    # thicker than the other one (i.e. the auto for target is not the same as
    # the auto for decoy)
    if score_type == "q-value":
        _ = plt.hist(decoy_epifany_q_values, bins=500, alpha=0.5, label="decoy")
        _ = plt.hist(target_epifany_q_values, bins=500, alpha=0.5, label="target")
        plt.title("epifany qvalue")
        plt.xlabel("Q value")
        plt.ylabel("Number of proteins")
        plt.legend()
        plt.show()
    elif score_type == "Posterior Probability":
        _ = plt.hist(target_epifany_pp_cleaned, bins=100, alpha=0.5, label="target")
        _ = plt.hist(decoy_epifany_pp_cleaned, bins=100, alpha=0.5, label="decoy")
        plt.title("Epifany Posterior Probability")
        plt.xlabel("Posterior Probability")
        plt.ylabel("Number of proteins")
        plt.legend()
        plt.show()

