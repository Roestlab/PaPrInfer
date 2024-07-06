import sys

import get_proteins_at_threshold
import sys

from pyopenms import IdXMLFile

if __name__ == "__main__":
    prot_ids = []
    pep_ids = []

    IdXMLFile().load(sys.argv[1], prot_ids, pep_ids)

    for protein_id in prot_ids:
        for hit in protein_id.getHits():
            print("Protein hit accession:", hit.getAccession())

    for peptide_id in pep_ids:
        # Peptide identification values
      # PeptideHits
        for hit in peptide_id.getHits():
            print("Peptide hit sequence:", hit.getSequence())

            print(" - Peptide hit score:", hit.getScore())

            print(" - Mapping to proteins:", [ev.getProteinAccession() for ev in hit.getPeptideEvidences()])
