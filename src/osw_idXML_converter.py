# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Thu Feb 11 12:00:41 2021
#
# @author: kren
# """


import sqlite3
import sys
from typing import Dict, List, Tuple

from pyopenms import AASequence, IdXMLFile, PeptideEvidence, PeptideHit, \
    PeptideIdentification, ProteinHit, ProteinIdentification

from ppi import get_all_link_for_peptide, get_all_protein_accession, \
    get_all_protein_id


def begin_connection(db_name: str):
    con = sqlite3.connect(db_name)
    return con


def fill_protein_identification(con, protein_identification,
                                protein_accession_dict: Dict[str, List[str]],
                                context,
                                run_id):
    # Each ProteinIdentification object stores a vector of protein hits
    target_protein_hit_list = []
    decoy_protein_hit_list = []
    protein_id_list = get_all_protein_id(con, context, run_id)
    target_protein_accession_list = []
    decoy_protein_accession_list = []

    for protein_id, decoy in protein_id_list:
        list_of_protein_accession = protein_accession_dict[protein_id]
        if decoy == 0:  # if not decoy, is target
            target_protein_accession_list.extend(list_of_protein_accession)
        elif decoy == 1:  # is decoy
            decoy_protein_accession_list.extend(list_of_protein_accession)

    # strip whitespaces
    stripped_target_protein_accession_list = [s.strip() for s in
                                              target_protein_accession_list]
    stripped_decoy_protein_accession_list = [s.strip() for s in
                                             decoy_protein_accession_list]

    # make them distinct
    stripped_target_protein_accession_list = list(set(stripped_target_protein_accession_list))
    stripped_decoy_protein_accession_list = list(set(stripped_decoy_protein_accession_list))

    # last parameter ask whether or not it is a decoy
    make_protein_hit(target_protein_hit_list,
                     stripped_target_protein_accession_list,
                     False)
    make_protein_hit(decoy_protein_hit_list,
                     stripped_decoy_protein_accession_list,
                     True)

    all_protein_hit_list = []
    all_protein_hit_list.extend(target_protein_hit_list)
    all_protein_hit_list.extend(decoy_protein_hit_list)
    protein_identification.setHits(all_protein_hit_list)


def make_protein_hit(protein_hit_list, protein_accession_list: List[str],
                     is_decoy):
    for protein_accession in protein_accession_list:
        protein_hit = ProteinHit()
        protein_hit.setAccession(protein_accession)
        if is_decoy:
            protein_hit.setMetaValue("target_decoy", b"decoy")
        else:
            protein_hit.setMetaValue("target_decoy", b"target")
        protein_hit_list.append(protein_hit)


def fill_peptide_identification(con,
                                peptide_id_list: List[ProteinIdentification],
                                protein_accession_dict: Dict[str, List[str]],
                                context: str, run_id: int, pep_limit: int):
    c = con.cursor()
    # TODO: here is where threshold PEP happens
    if context == 'global':
        c.execute(
            """SELECT Score_Peptide.PEP, Peptide.unmodified_sequence, 
            Score_Peptide.peptide_id, PEPTIDE.DECOY
            FROM Peptide
            INNER JOIN Score_Peptide ON Score_Peptide.peptide_id = Peptide.id
            WHERE SCORE_PEPTIDE.CONTEXT = 'global'
            AND PEP<=:pep_limit""", {'pep_limit': pep_limit})
    elif context == 'run-specific':
        c.execute(
            """SELECT Score_Peptide.PEP, Peptide.unmodified_sequence, 
            Score_Peptide.peptide_id, PEPTIDE.DECOY
            FROM Peptide
            INNER JOIN Score_Peptide ON Score_Peptide.peptide_id = Peptide.id
            WHERE CONTEXT='run-specific' AND RUN_ID=:run_id
            AND PEP<=:pep_limit""", {'run_id': run_id, 'pep_limit': pep_limit})

    linked_protein_dict = get_all_link_for_peptide(con)

    for row in c.fetchall():
        posterior_error_probability = row[0]
        aa_sequence = row[1]
        pep_id = str(row[2])
        peptide_is_decoy = bool(row[3])

        # if this peptide is not mapped to any protein
        # if pep_id not in linked_protein_dict:
        #     continue

        # Create new peptide identification object and fill basic information
        peptide_identification = PeptideIdentification()

        peptide_identification.setScoreType("pep")
        peptide_identification.setHigherScoreBetter(False)
        # I guess this is just name? (actually I think it has to be
        # the same as the protein identification object)
        peptide_identification.setIdentifier("IdentificationRun1")

        # peptide_identification is the peptide Identification object,
        # pep_id is the peptide identifier

        add_peptide_hit(linked_protein_dict, protein_accession_dict,
                        peptide_identification,
                        posterior_error_probability, aa_sequence, pep_id,
                        peptide_is_decoy)

        peptide_id_list.append(peptide_identification)

    c.close()


def add_peptide_hit(linked_protein_dict: Dict[str, List[Tuple[str, int]]],
                    protein_accession_dict: Dict[str, List[str]],
                    peptide_id: PeptideIdentification, pep: int,
                    aa_sequence: str, pep_id: str, peptide_is_decoy: bool):
    # create a new PeptideHit (best PSM, best score)
    peptide_hit = PeptideHit()
    peptide_hit.setScore(pep)
    peptide_hit.setSequence(AASequence.fromString(aa_sequence))
    if peptide_is_decoy:
        peptide_hit.setMetaValue("target_decoy", b"decoy")
    else:
        peptide_hit.setMetaValue("target_decoy", b"target")

    protein_accession_list = []
    # for all protein sqlite id that this peptide sqlite id maps to
    for pro_id, decoy in linked_protein_dict[pep_id]:
        # we combine the list of accession each protein sqlite id maps to
        protein_accession_list.extend(protein_accession_dict[pro_id])


    # strip whitespace
    stripped_protein_accession_list = [s.strip() for s in
                                       protein_accession_list]

    # then we make the list protein accession distinct
    stripped_protein_accession_list = list(set(stripped_protein_accession_list))

    # with the distinct protein accession, we then construct peptide evidence
    ev_list = []
    for protein_accession in stripped_protein_accession_list:
        ev = PeptideEvidence()
        ev.setProteinAccession(protein_accession)
        # ev.setAABefore(b"UNKNOWN_AA")
        # ev.setAAAfter(b"UNKNOWN_AA")
        ev_list.append(ev)

    # there is only 1 peptide hit per peptide identification
    # set all peptide_evidences to peptide_hit
    peptide_hit.setPeptideEvidences(ev_list)

    # add all peptide_hit to peptide identification
    # there should only be one per peptide
    peptide_id.setHits([peptide_hit])


def store_on_disk(out_file_name, protein_id_list, peptide_id_list):
    # Store the identification data in an idXML file
    # the string for file storing is just path
    IdXMLFile().store(out_file_name, protein_id_list, peptide_id_list)


def main(input_file: str, out_file: str, context: str, run_id: str,
         pep_limit: str) -> None:
    run_id = int(run_id)
    pep_limit = int(pep_limit)

    con = begin_connection(input_file)

    protein_accession_dict = get_all_protein_accession(con)

    # Create new protein identification object corresponding to a single search
    # but targets
    # and fill basic information of protein identification object
    protein_identification = ProteinIdentification()
    protein_identification.setIdentifier("IdentificationRun1")

    # then add protein hits to protein identification object
    fill_protein_identification(con, protein_identification,
                                protein_accession_dict, context, run_id)

    # make protein identification into a list
    protein_id_list = [protein_identification]

    # do the same thing with peptide identification
    peptide_id_list = []
    fill_peptide_identification(con, peptide_id_list, protein_accession_dict,
                                context, run_id, pep_limit)

    # for protein_id in protein_id_list:
    #     for hit in protein_id.getHits():
    #         print("Protein hit accession:", hit.getAccession())
    #         if hit.getAccession() == " tr|H7C1C4|H7C1C4_HUMAN":
    #             print("found it")
    #
    # for peptide_id in peptide_id_list:
    #     for hit in peptide_id.getHits():
    #         print(" - Peptide hit sequence:", hit.getSequence())
    #         print(" - Mapping to proteins:", [ev.getProteinAccession() for ev in
    #                                           hit.getPeptideEvidences()])

    # then finally store on disk
    store_on_disk(out_file, protein_id_list, peptide_id_list)


if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # print("""usage: osw_idXML_converter.py <sql file path> <out file path>""")
    # input_osw_file, out_idXML_file, context, run_id, pep_limit
    # if context is global, run_id can be anything
    # also for osw that is with the msfragger peptide, it only has global
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
