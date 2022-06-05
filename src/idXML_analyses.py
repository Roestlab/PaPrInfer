
import sys
import matplotlib.pyplot as plt
from pyopenms import IdXMLFile


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

    protein_to_peptide_target = {}
    protein_to_peptide_decoy = {}

    for peptide_id in pep_ids:

        for peptide_hit in peptide_id.getHits():

            peptide_evidence_vector = peptide_hit.getPeptideEvidences()

            for peptide_evidence in peptide_evidence_vector:

                protein_accession = peptide_evidence.getProteinAccession()

                if peptide_hit.getMetaValue("target_decoy") == "target":
                    # if key was not in the dict, setdefault return the default value,
                    # which is an empty list in this case. if it was then it returns the value
                    protein_to_peptide_target.setdefault(protein_accession, [])\
                        .append(peptide_hit.getSequence())
                elif peptide_hit.getMetaValue("target_decoy") == "decoy":
                    # if key was not in the dict, setdefault return the default value,
                    # which is an empty list in this case. if it was then it returns the value
                    protein_to_peptide_decoy.setdefault(protein_accession, [])\
                        .append(peptide_hit.getSequence())

    target_ptp_distribution = []
    decoy_ptp_distribution = []
    for key, value in protein_to_peptide_target.items():
        target_ptp_distribution.append(len(value))

    for key, value in protein_to_peptide_decoy.items():
        decoy_ptp_distribution.append(len(value))


    _ = plt.hist(target_ptp_distribution, bins=100, alpha=0.5, label="target")
    _ = plt.hist(decoy_ptp_distribution, bins=100, alpha=0.5, label="decoy")
    plt.title("epifany peptides per protein")
    plt.xlabel("Number of peptide")
    plt.ylabel("Number of proteins")
    plt.legend()
    plt.show()
