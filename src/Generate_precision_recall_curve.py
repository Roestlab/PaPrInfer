
import sys
import matplotlib.pyplot as plt
import numpy as np

import compare_idpicker_epifany_pyprophet


def main(epifany_file: str, idpicker_file: str, pyprophet_file: str):

    threshold_list = []
    num_idpicker_protein_list = []
    num_epifany_protein_list = []
    num_pyprophet_protein_list = []

    for threshold in np.arange(0.1, 1.01, 0.1):

        print("threshold", threshold)

        num_epifany_protein, num_idpicker_protein, num_pyprophet_protein = \
            compare_idpicker_epifany_pyprophet.main(epifany_file, idpicker_file, pyprophet_file, str(threshold))

        threshold_list.append(threshold)
        num_idpicker_protein_list.append(num_idpicker_protein)
        num_epifany_protein_list.append(num_epifany_protein)
        num_pyprophet_protein_list.append(num_pyprophet_protein)

    plt.plot(threshold_list, num_idpicker_protein_list, "-o", color='blue', label="Idpicker Protein Groups")
    # plt.plot(threshold_list, num_epifany_protein_list, "-s", color='red', label="Epifany Protein Groups")
    plt.plot(threshold_list, num_pyprophet_protein_list, "-^", color='green', label="Pyprophet Protein (SwissProt)")
    # plt.title("Precision Recall Curve: Swiss-Prot")
    plt.title("Precision Recall Curve: UniProt")
    plt.xlabel("FDR threshold")
    plt.ylabel("Number of proteins")
    # plt.legend(loc='lower right')
    plt.legend(loc='center right')
    plt.show()

    # maybe add vertical line at 0.05

    # add horizontal line at maximum protein?


def zoomed_in_version(epifany_file: str, idpicker_file: str, pyprophet_file: str):
    threshold_list = []
    num_idpicker_protein_list = []
    num_epifany_protein_list = []
    num_pyprophet_protein_list = []

    for threshold in np.arange(0.01, 0.101, 0.01):

        print("threshold", threshold)

        num_epifany_protein, num_idpicker_protein, num_pyprophet_protein = \
            compare_idpicker_epifany_pyprophet.main(epifany_file, idpicker_file,
                                                    pyprophet_file,
                                                    str(threshold))

        threshold_list.append(threshold)
        num_idpicker_protein_list.append(num_idpicker_protein)
        num_epifany_protein_list.append(num_epifany_protein)
        num_pyprophet_protein_list.append(num_pyprophet_protein)

    plt.plot(threshold_list, num_idpicker_protein_list, "-o", color='blue',
             label="Idpicker Protein Groups")
    # plt.plot(threshold_list, num_epifany_protein_list, "-s", color='red',
    #          label="Epifany Protein Groups")
    plt.plot(threshold_list, num_pyprophet_protein_list, "-^", color='green',
             label="Pyprophet Proteins (Swiss-Prot)")
    # plt.title("Precision Recall Curve: Swiss-Prot 0-0.1")
    plt.title("Precision Recall Curve: UniProt 0-0.1")
    plt.xlabel("FDR threshold")
    plt.ylabel("Number of proteins")
    plt.legend(loc='lower right')
    plt.show()


def even_more_zoomed_in_version(epifany_file: str, idpicker_file: str, pyprophet_file: str):
    threshold_list = []
    num_idpicker_protein_list = []
    num_epifany_protein_list = []
    num_pyprophet_protein_list = []

    for threshold in np.arange(0.001, 0.011, 0.001):

        print("threshold", threshold)

        num_epifany_protein, num_idpicker_protein, num_pyprophet_protein = \
            compare_idpicker_epifany_pyprophet.main(epifany_file, idpicker_file,
                                                    pyprophet_file,
                                                    str(threshold))

        threshold_list.append(threshold)
        num_idpicker_protein_list.append(num_idpicker_protein)
        num_epifany_protein_list.append(num_epifany_protein)
        num_pyprophet_protein_list.append(num_pyprophet_protein)

    plt.plot(threshold_list, num_idpicker_protein_list, "-o", color='blue',
             label="Idpicker Protein Groups")
    # plt.plot(threshold_list, num_epifany_protein_list, "-s", color='red',
    #          label="Epifany Protein Groups")
    plt.plot(threshold_list, num_pyprophet_protein_list, "-^", color='green',
             label="Pyprophet Proteins (Swiss-Prot)")
    # plt.title("Precision Recall Curve: Swiss-Prot 0-0.1")
    plt.title("Precision Recall Curve: UniProt 0-0.01")
    plt.xlabel("FDR threshold")
    plt.ylabel("Number of proteins")
    plt.legend(loc='lower right')
    plt.show()


def venn_diagram(epifany_file: str, idpicker_file: str, pyprophet_file: str):

    num_epifany_protein, num_idpicker_protein, num_pyprophet_protein = \
        compare_idpicker_epifany_pyprophet.main(epifany_file, idpicker_file,
                                                pyprophet_file,
                                                str(0.05))

if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # print("""usage: osw_idXML_converter.py <sql file path> <q-value threshold for peptides>""")
    # main(sys.argv[1], sys.argv[2])
    # epifany output file, idpicker output file
    # pyprophet output file is the same as idpicker

    # main(sys.argv[1], sys.argv[2], sys.argv[3])
    # zoomed_in_version(sys.argv[1], sys.argv[2], sys.argv[3])
    # venn_diagram(sys.argv[1], sys.argv[2], sys.argv[3])
    even_more_zoomed_in_version(sys.argv[1], sys.argv[2], sys.argv[3])
