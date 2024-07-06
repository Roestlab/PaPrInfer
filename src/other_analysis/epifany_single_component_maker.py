# this code is largely inspired by
# https://pyopenms.readthedocs.io/en/latest/datastructures_id.html
import sys

from pyopenms import AASequence, IdXMLFile, PeptideEvidence, PeptideHit, \
    PeptideIdentification, ProteinHit, ProteinIdentification


def main(out_file_name):
    # for proteins
    protein_id = ProteinIdentification()
    protein_id.setIdentifier("IdentificationRun1")

    target_protein_accession_list = ['PAL4E', 'PAL4H', 'PAL4D',
                                     'PAL4F', 'PAL4A', 'A0A0H2UH3',
                                     'PAL4C', 'PAL4G', 'PPIA',
                                     'C9J5S7', 'F8WE65', 'A0A7I2V5J5',
                                     'A0A7P0S768', 'E5RIZ5', 'A0A7I2V4V1']
    decoy_protein_accession_list = ['DECOY PAL4E', 'DECOY PAL4H', 'DECOY PAL4D',
                                    'DECOY PAL4F', 'DECOY PAL4A', 'DECOY A0A0H2UH3',
                                    'DECOY PAL4C', 'DECOY PAL4G', 'DECOY PPIA',
                                    'DECOY C9J5S7', 'DECOY F8WE65', 'DECOY A0A7I2V5J5',
                                    'DECOY A0A7P0S768', 'DECOY E5RIZ5', 'DECOY A0A7I2V4V1']

    make_protein_hit_list(protein_id, target_protein_accession_list,
                          decoy_protein_accession_list)

    peptide_id_list = fill_all_peptide_id()

    # make protein identification into a list
    protein_id_list = [protein_id]

    IdXMLFile().store(out_file_name, protein_id_list, peptide_id_list)


def fill_all_peptide_id():
    # for peptides
    peptide_id_list = []
    peptide_id = fill_peptide_id('TEWLDGK',
                                 ['PAL4E', 'PAL4H', 'PAL4D', 'PAL4F', 'PAL4A',
                                  'A0A0H2UH3', 'PAL4C', 'PAL4G', 'PPIA'],
                                 1.0852594434272e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('IIPGFMCQGGDFTR',
                                 ['PAL4E', 'PAL4H', 'PAL4D', 'PAL4F', 'PAL4A',
                                  'A0A0H2UH3', 'PAL4C', 'PPIA', 'C9J5S7',
                                  'F8WE65'],
                                 1.0852594434272e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('HTGSGILSMANAGPNTNGSQFFICTAK',
                                 ['PAL4C', 'PAL4G'],
                                 1.52664048381647e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('HTGPGILSMANAGPNTNGSQFFICTAK',
                                 ['PPIA', 'C9J5S7', 'F8WE65',
                                  'A0A7I2V4V1'],
                                 1.52664048381647e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('ITIADCGQLE',
                                 ['PPIA'],
                                 0.000330313269055314, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('EGMNIVEAMER',
                                 ['PPIA'],
                                 1.52664048381647e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('FEDENFILK',
                                 ['PPIA', 'A0A7I2V4V1', 'C9J5S7', 'F8WE65'],
                                 1.0852594434272e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('VSFELFADKVPK',
                                 ['PPIA', 'A0A7I2V5J5', 'C9J5S7',
                                  'F8WE65'],
                                 1.52664048381647e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('VSFELFADK',
                                 ['A0A7I2V5J5', 'C9J5S7',
                                  'F8WE65', 'PPIA'],
                                 1.0852594434272e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('VNPTVFFDIAVDGEPLGR',
                                 ['PPIA', 'A0A7I2V5J5', 'A0A7P0S768', 'C9J5S7',
                                  'E5RIZ5',
                                  'F8WE65'],
                                 0.000111454024641728, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('CQGGDFTR',
                                 ['A0A7I2V4V1'],
                                 3.12238306472365e-05, False)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('MCQGGDFTR',
                                 ['A0A7I2V4V1'],
                                 1.52664048381647e-05, False)
    peptide_id_list.append(peptide_id)





    peptide_id = fill_peptide_id('GDLWETR',
                                 ['DECOY PAL4E', 'DECOY PAL4H', 'DECOY PAL4D', 'DECOY PAL4F', 'DECOY PAL4A',
                                  'DECOY A0A0H2UH3', 'DECOY PAL4C', 'DECOY PAL4G', 'DECOY PPIA'],
                                 0.600474502860762, True)
    peptide_id_list.append(peptide_id)

    # note there were 2 decoy peptide with this unmodified sequence, but there was
    # only 1 target peptide counterpart, I choose the one whose modified sequence
    # is the pseudo reverse of each other
    peptide_id = fill_peptide_id('TFDGGQCMFGPIIK',
                                 ['DECOY PAL4E', 'DECOY PAL4H', 'DECOY PAL4D', 'DECOY PAL4F', 'DECOY PAL4A',
                                  'DECOY A0A0H2UH3', 'DECOY PAL4C', 'DECOY PPIA', 'DECOY C9J5S7',
                                  'DECOY F8WE65'],
                                 0.68715324112149, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('ATCIFFQSGNTNPGANAMSLIGSGTHR',
                                 ['DECOY PAL4C', 'DECOY PAL4G'],
                                 0.536328078844151, True)
    peptide_id_list.append(peptide_id)
    # same situation as TFDGGQCMFGPIIK
    peptide_id = fill_peptide_id('ATCIFFQSGNTNPGANAMSLIGPGTHR',
                                 ['DECOY PPIA', 'DECOY C9J5S7', 'DECOY F8WE65',
                                  'DECOY A0A7I2V4V1'],
                                 0.407953386259084, True)
    peptide_id_list.append(peptide_id)

    # no decoy peptide
    # peptide_id = fill_peptide_id('ELQGCDAITI',
    #                              ['PPIA'],
    #                              0.000330313269055314, False)
    # peptide_id_list.append(peptide_id)

    # same
    peptide_id = fill_peptide_id('EMAEVINMGEK',
                                 ['DECOY PPIA'],
                                 0.315037722210494, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('LIFNEDEFR',
                                 ['DECOY PPIA', 'DECOY A0A7I2V4V1', 'DECOY C9J5S7', 'DECOY F8WE65'],
                                 0.060715732750426, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('PVKDAFLEFSVR',
                                 ['DECOY PPIA', 'DECOY A0A7I2V5J5', 'DECOY C9J5S7',
                                  'DECOY F8WE65'],
                                 0.623341508906059, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('DAFLEFSVR',
                                 ['DECOY PPIA', 'DECOY C9J5S7',
                                  'DECOY F8WE65'],
                                 0.691996399481143, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('GLPEGDVAIDFFVTPNVK',
                                 ['DECOY PPIA', 'DECOY A0A7I2V5J5', 'DECOY A0A7P0S768', 'DECOY C9J5S7',
                                  'DECOY E5RIZ5', 'DECOY F8WE65'],
                                 0.662020173704249, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('TFDGGQCK',
                                 ['DECOY A0A7I2V4V1'],
                                 0.191864127940947, True)
    peptide_id_list.append(peptide_id)
    peptide_id = fill_peptide_id('TFDGGQCMK',
                                 ['DECOY A0A7I2V4V1'],
                                 0.488750833547983, True)
    peptide_id_list.append(peptide_id)

    return peptide_id_list


def make_protein_hit_list(protein_id, target_protein_accession_list,
                          decoy_protein_accession_list):
    protein_hit_list = list()

    for target_protein in target_protein_accession_list:
        protein_hit = make_protein_hit(target_protein, is_decoy=False)
        protein_hit_list.append(protein_hit)

    for decoy_protein in decoy_protein_accession_list:
        protein_hit = make_protein_hit(decoy_protein, is_decoy=True)
        protein_hit_list.append(protein_hit)

    protein_id.setHits(protein_hit_list)


def make_protein_hit(protein_accession, is_decoy):
    # Each ProteinIdentification object stores a vector of protein hits
    protein_hit = ProteinHit()
    protein_hit.setAccession(protein_accession)
    if is_decoy:
        protein_hit.setMetaValue("target_decoy", b"decoy")
    else:
        protein_hit.setMetaValue("target_decoy", b"target")
    return protein_hit


def fill_peptide_id(peptide_sequence, peptide_evidence_list, pep, is_decoy):
    peptide_id = PeptideIdentification()
    peptide_id.setIdentifier("IdentificationRun1")
    peptide_id.setScoreType("pep")
    peptide_id.setHigherScoreBetter(False)

    peptide_hit = make_peptide_hits(peptide_sequence, pep, is_decoy)
    make_peptide_evidence(peptide_hit, peptide_evidence_list)

    peptide_id.setHits([peptide_hit])

    return peptide_id


def make_peptide_hits(peptide_sequence, pep, is_decoy):
    peptide_hit = PeptideHit()
    peptide_hit.setScore(pep)
    peptide_hit.setSequence(AASequence.fromString(peptide_sequence))
    if is_decoy:
        peptide_hit.setMetaValue("target_decoy", b"decoy")
    else:
        peptide_hit.setMetaValue("target_decoy", b"target")
    return peptide_hit


def make_peptide_evidence(peptide_hit, protein_accession_list):

    ev_list = []
    for protein_accession in protein_accession_list:
        ev = PeptideEvidence()
        ev.setProteinAccession(protein_accession)
        ev_list.append(ev)
    peptide_hit.setPeptideEvidences(ev_list)


if __name__ == "__main__":
    main(sys.argv[1])
