import csv
import sys
import sqlite3


def main(tsv_filename, osw_filename):

    all_peptide = get_all_peptide(tsv_filename)

    modify_osw(osw_filename, all_peptide)


def get_all_peptide(tsv_filename):
    """
    0 is peptide sequence without modification
    10 is razor protein (full protein accession?)
    11 is protein ID (the one matching to osw)
    -1 is the last column, other mapped protein
    :param tsv_filename: the file path of the peptide.tsv output
    from msfragger
    """

    all_peptides = []

    with open(tsv_filename) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            if line[-1] == "":
                all_peptides.append([line[0], line[11], line[10]])
            else:
                all_peptides.append([line[0], line[11], line[10] + ',' + line[-1]])

    return all_peptides

# since the rows in peptide.tsv represent a set of edges
# that is, a peptide map to a group of proteins

# erase protein and mapping table (for only targets)
# from delete mapping target entries
# I can do:
# delete from mapping where protein_id in (select protein_id from mapping
# inner join protein on mapping.protein_id = protein.id where protein.decoy = 0)

# for each line I have the peptide sequence and a list of proteins
# find the matching peptide id in peptide table based on sequence
# assign a random new id to the protein
# add this to the protein and mapping table
# (for protein, just the new id, the list of protein, and decoy as 0)
# (for mapping just matching peptide id, and this new id
# only problem, I wont be able to output pyprophet proteins'a accession
# and now score protein do not match with protein


def modify_osw(osw_filename, all_peptide):
    con = sqlite3.connect(osw_filename)
    c = con.cursor()

    # delete target entries from mapping
    # it check target/decoy based on protein
    c.execute(
        """DELETE FROM PEPTIDE_PROTEIN_MAPPING
        WHERE PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID IN
        (SELECT PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING 
        INNER JOIN PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID
        WHERE PROTEIN.DECOY = 0)"""
    )

    print("mapping target entries deleted")

    # need to delete entries in score protein too
    # it check target/decoy based on protein
    c.execute(
        """DELETE FROM SCORE_PROTEIN
        where SCORE_PROTEIN.PROTEIN_ID IN
        (SELECT PROTEIN_ID FROM SCORE_PROTEIN 
        INNER JOIN PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.ID
        WHERE PROTEIN.DECOY = 0)"""
    )

    print("score target entries deleted")

    # delete target entries from protein
    c.execute(
        """delete from PROTEIN
        where protein.DECOY = 0"""
    )

    print("protein target entries deleted")

    c.execute(
        """SELECT ID
        FROM PROTEIN"""
    )

    decoy_protein_id_list = []
    for row in c.fetchall():
        decoy_protein_id_list.append(row[0])

    # for all rows in tsv
    progress_count = 1
    total_proteins = len(all_peptide)
    for protein in all_peptide:
        pep_seq = protein[0]
        accession = protein[1]
        all_mapped_protein = protein[2]

        # find the matching peptide id based on peptide sequence
        c.execute(
            """SELECT PEPTIDE.ID 
            FROM PEPTIDE
            WHERE PEPTIDE.UNMODIFIED_SEQUENCE=:pep_seq
            AND PEPTIDE.DECOY = 0""", {'pep_seq': pep_seq}
        )
        pep_id_list = c.fetchall()
        # TODO: tsv gives unmodified sequence, but the osw there is entries
        #   where unmodified sequence is the same, but the modified is not
        if len(pep_id_list) == 0:
            print(progress_count, total_proteins, "sequence not in osw")
            progress_count += 1
            continue

        # there should be only 1 row that is returned by the
        # sql query, and there is only 1 column, so we want 0, 0
        pep_id = pep_id_list[0][0]

        # it may used the same pro id as an undeleted decoy protein
        # so if it does, then change the new protein id
        if progress_count in decoy_protein_id_list:
            pro_id = progress_count + total_proteins
        else:
            pro_id = progress_count

        decoy = 0

        # create entry in protein
        c.execute(
            """INSERT INTO PROTEIN(ID, PROTEIN_ACCESSION, DECOY) 
            VALUES(:id, :protein_accession, :decoy)""",
            {'id': pro_id, 'protein_accession': all_mapped_protein, 'decoy': decoy}
        )

        # create entry in mapping
        c.execute(
            """INSERT INTO PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID, PROTEIN_ID) 
            VALUES(:peptide_id, :protein_id)""",
            {'peptide_id': pep_id, 'protein_id': pro_id}
        )

        # then add in score protein
        # do global context, null run id, protein id, 0 score, 0 pvalue, 0 qvalue
        # 0 pep
        context = 'global'
        run_id = 0
        score = 0
        pvalue = 0
        qvalue = 0
        pep = 0

        c.execute(
            """INSERT INTO SCORE_PROTEIN(CONTEXT, RUN_ID, PROTEIN_ID, SCORE,
            PVALUE, QVALUE, PEP) VALUES(:context, :run_id, :protein_id, :score,
            :pvalue, :qvalue, :pep)""",
            {'context': context, 'run_id': run_id, "protein_id": pro_id,
             'score': score, 'pvalue': pvalue, 'qvalue': qvalue, 'pep': pep}
        )

        progress_count += 1
        print(progress_count, total_proteins)

    # now delete peptide that have no matching proteins

    # so all mapping target entries were deleted,
    # so does not have mapping
    # only the ones that were added back has a mapping

    # so peptide with no matching protein should have no mapping
    c.execute(
        """DELETE FROM SCORE_PEPTIDE
        WHERE SCORE_PEPTIDE.PEPTIDE_ID NOT IN
        (SELECT PEPTIDE.ID FROM PEPTIDE
        INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID)"""
    )

    c.execute(
        """DELETE FROM PEPTIDE
        WHERE PEPTIDE.ID NOT IN
        (SELECT PEPTIDE.ID FROM PEPTIDE
        INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID)"""
    )

    con.commit()
    c.close()


if __name__ == "__main__":
    # tsv, then osw
    main(sys.argv[1], sys.argv[2])
