import src.compare_idpicker_epifany
import sys
import matplotlib.pyplot as plt
import numpy as np

from src import compare_idpicker_epifany


def main(epifany_file: str, idpicker_file: str):

    threshold_list = []
    num_idpicker_protein_list = []
    num_epifany_protein_list = []

    for threshold in np.arange(0, 1.05, 0.05):

        print("threshold", threshold)

        num_idpicker_protein, num_epifany_protein = \
            compare_idpicker_epifany.main(epifany_file, idpicker_file, str(threshold))

        threshold_list.append(threshold)
        num_idpicker_protein_list.append(num_idpicker_protein)
        num_epifany_protein_list.append(num_epifany_protein)

    plt.plot(threshold_list, num_idpicker_protein_list, "-o", color='blue', label="Idpicker")
    plt.plot(threshold_list, num_epifany_protein_list, "-s", color='red', label="Epifany")
    plt.title("Precision Recall Curve: Idpicker vs Epifany")
    plt.xlabel("FDR threshold")
    plt.ylabel("Number of proteins")
    plt.legend(loc='upper left')
    plt.show()

    # maybe add vertical line at 0.05

    # add horizontal line at maximum protein?


if __name__ == "__main__":
    # if len(sys.argv) != 3:
    # print("""usage: osw_idXML_converter.py <sql file path> <q-value threshold for peptides>""")
    # main(sys.argv[1], sys.argv[2])
    # epifany output file, idpicker output file
    main(sys.argv[1], sys.argv[2])
