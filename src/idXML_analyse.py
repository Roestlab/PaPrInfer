
import sys
import matplotlib.pyplot as plt
from pyopenms import IdXMLFile
import collections


"""trying to find number of peptides per protein distribution either
1) how many peptides match to a protein or (from pre inference idXML)
2) how many peptide are explained by inferring this protein (post inference idXML"""

if __name__ == "__main__":
    prot_ids = []
    pep_ids = []
    epifany_file = sys.argv[1]

    IdXMLFile().load(epifany_file, prot_ids,
                     pep_ids)

    score_type = "None"

    protein_to_peptide = {}

    for peptide_id in pep_ids:

        for peptide_hit in peptide_id.getHits():



            peptide_evidence_vector = peptide_hit.getPeptideEvidences()

            for peptide_evidence in peptide_evidence_vector:

                protein_accession = peptide_evidence.getProteinAccession()

                # if key was not in the dict, setdefault return the default value,
                # empty list here. if it was then it returns the value
                protein_to_peptide.setdefault(protein_accession, []).append(peptide_hit.getSeq)


                #             if hit.getMetaValue("target_decoy") == "target":




    # if score_type == "q-value":
    #     _ = plt.hist(decoy_epifany_q_values, bins='auto', alpha=0.5, label="decoy")
    #     _ = plt.hist(target_epifany_q_values, bins='auto', alpha=0.5, label="target")
    #     plt.title("epifany qvalue")
    #     plt.xlabel("Q value")
    #     plt.ylabel("Number of proteins")
    #     plt.legend()
    #     plt.show()
    # elif score_type == "Posterior Probability":
    #     _ = plt.hist(target_epifany_pp_cleaned, bins=100, alpha=0.5, label="target")
    #     _ = plt.hist(decoy_epifany_pp_cleaned, bins=100, alpha=0.5, label="decoy")
    #     plt.title("proteininfernce Posterior Probability")
    #     plt.xlabel("Posterior Probability")
    #     plt.ylabel("Number of proteins")
    #     plt.legend()
    #     plt.show()

