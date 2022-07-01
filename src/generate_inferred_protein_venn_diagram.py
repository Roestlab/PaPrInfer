import sys
from typing import List, Set, Tuple, Union

from matplotlib import pyplot as plt

import get_proteins_at_threshold

from matplotlib_venn import venn2, venn3

"""first of all, there is no pyprophet uniprot 
we need swissprot in all 3 to see what is the difference from pyprophet
we have uniprot of 2 to see how much uniprot differs from idpicker
i think we need a uniprot vs swissprot to see how much better is using uniprot"""


def main(uniprot_epifany_file, uniprot_idpicker_file, swissprot_epifany_file,
         swissprot_idpicker_file, swissprot_pyprophet, threshold) -> None:

    uniprot_epifany_proteins_accession, uniprot_idpicker_proteins_accession, \
    swissprot_epifany_proteins_accession, swissprot_idpicker_proteins_accession, \
    swissprot_pyprophet_proteins_accession = \
        get_all_protein_accession(
            uniprot_epifany_file, uniprot_idpicker_file, swissprot_epifany_file,
            swissprot_idpicker_file, swissprot_pyprophet, threshold
        )

    make_venn_diagram_swissprot(swissprot_epifany_proteins_accession,
                                swissprot_idpicker_proteins_accession,
                                swissprot_pyprophet_proteins_accession)

    make_venn_diagram_uniprot(uniprot_epifany_proteins_accession,
                              uniprot_idpicker_proteins_accession)


def get_all_protein_accession(uniprot_epifany_file, uniprot_idpicker_file,
                              swissprot_epifany_file, swissprot_idpicker_file,
                              swissprot_pyprophet, threshold) \
        -> Tuple[Union[List[int], Set[str]], Set[str], Union[List[int], Set[str]], Set[str], Set[str]]:

    uniprot_epifany_proteins_accession = \
        get_proteins_at_threshold.get_epifany_result(uniprot_epifany_file,
                                                     remove_decoy=True,
                                                     return_qvalue=False,
                                                     threshold=threshold)

    uniprot_idpicker_proteins_accession = \
        get_proteins_at_threshold.get_idpicker_accessions(uniprot_idpicker_file,
                                                          threshold=threshold)

    swissprot_epifany_proteins_accession = \
        get_proteins_at_threshold.get_epifany_result(swissprot_epifany_file,
                                                     remove_decoy=True,
                                                     return_qvalue=False,
                                                     threshold=threshold)

    swissprot_idpicker_proteins_accession = \
        get_proteins_at_threshold.get_idpicker_accessions(swissprot_idpicker_file,
                                                          threshold=threshold)

    swissprot_pyprophet_proteins_accession = \
        get_proteins_at_threshold.get_idpicker_accessions(swissprot_pyprophet,
                                                       threshold=threshold)

    return uniprot_epifany_proteins_accession, uniprot_idpicker_proteins_accession, \
           swissprot_epifany_proteins_accession, swissprot_idpicker_proteins_accession, \
           swissprot_pyprophet_proteins_accession


def make_venn_diagram_swissprot(swissprot_epifany_proteins_accession: Set[str],
                                swissprot_idpicker_proteins_accession: Set[str],
                                swissprot_pyprophet_proteins_accession: Set[str])\
        -> None:
    venn3(
        [swissprot_epifany_proteins_accession,
         swissprot_idpicker_proteins_accession,
         swissprot_pyprophet_proteins_accession],
        ('swissprot epifany',
         'swissprot idpicker',
         'swissprot pyprophet')
    )
    plt.title('venn diagram of all 3 in swissprot')
    plt.savefig('figures and files/swissprot venn.pdf')
    plt.show()


def make_venn_diagram_uniprot(uniprot_epifany_proteins_accession: Set[str],
                              uniprot_idpicker_proteins_accession: Set[str])\
        -> None:

    venn2(
        [uniprot_epifany_proteins_accession,
         uniprot_idpicker_proteins_accession],
        ('uniprot epifany',
         'uniprot idpicker')
    )
    plt.title('venn diagram of all 2 in uniprot')
    plt.savefig('figures and files/uniprot venn.pdf')
    plt.show()


if __name__ == "__main__":
    # uniprot_epifany_file, uniprot_idpicker_file, swissprot_epifany_file,
    #          swissprot_idpicker_file, swissprot_pyprophet, threshold

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
         sys.argv[6])

    # maybe I need one to compare swissprot idpicker
