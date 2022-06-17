
import sys
import matplotlib.pyplot as plt
import numpy as np

import compare_idpicker_epifany_pyprophet


def main(epifany_file: str, idpicker_file: str, pyprophet_file: str, include_epifany: str, zoom_degree):

    threshold_list = []
    num_idpicker_protein_list = []
    num_epifany_protein_list = []
    num_pyprophet_protein_list = []

    get_protein_at_threshold(epifany_file, idpicker_file,
                             num_epifany_protein_list,
                             num_idpicker_protein_list,
                             num_pyprophet_protein_list, pyprophet_file,
                             threshold_list, zoom_degree)

    plot_curves(include_epifany, num_epifany_protein_list,
                num_idpicker_protein_list, num_pyprophet_protein_list,
                threshold_list, zoom_degree)


def get_protein_at_threshold(epifany_file, idpicker_file,
                             num_epifany_protein_list,
                             num_idpicker_protein_list,
                             num_pyprophet_protein_list, pyprophet_file,
                             threshold_list, zoom_degree):

    # remember it is start, stop, step
    if zoom_degree == 'normal':
        loop_range = [0.1, 1.01, 0.1]
    elif zoom_degree == 'zoomed':
        loop_range = [0.01, 0.101, 0.01]
    elif zoom_degree == 'zoomed plus':
        loop_range = [0.001, 0.011, 0.001]
    else:
        print("not one of the option, please check usage")
        exit()
        loop_range = [0.1, 1.01, 0.1]

    for threshold in np.arange(loop_range[0], loop_range[1], loop_range[2]):

        print("threshold", threshold)

        num_epifany_protein, num_idpicker_protein, num_pyprophet_protein = \
            compare_idpicker_epifany_pyprophet.main(epifany_file, idpicker_file,
                                                    pyprophet_file,
                                                    str(threshold))

        threshold_list.append(threshold)
        num_idpicker_protein_list.append(num_idpicker_protein)
        num_epifany_protein_list.append(num_epifany_protein)
        num_pyprophet_protein_list.append(num_pyprophet_protein)


def plot_curves(include_epifany, num_epifany_protein_list,
                num_idpicker_protein_list, num_pyprophet_protein_list,
                threshold_list, zoom_degree):
    plt.figure(0)
    plt.plot(threshold_list, num_idpicker_protein_list, "-o", color='blue',
             label="Idpicker Protein Groups")
    plt.plot(threshold_list, num_pyprophet_protein_list, "-^", color='green',
             label="Pyprophet Protein (SwissProt)")
    if include_epifany == 'yes':
        print("epifany is included")
        plt.plot(threshold_list, num_epifany_protein_list, "-s", color='red',
                 label="Epifany Protein Groups")
    # plt.title("Precision Recall Curve: Swiss-Prot")
    plt.title("Precision Recall Curve: UniProt")
    plt.xlabel("FDR threshold")
    plt.ylabel("Number of proteins")
    # plt.legend(loc='lower right')
    plt.legend(loc='center right')

    filename = zoom_degree + '.pdf'
    plt.savefig(filename)
    plt.show()
    # maybe add vertical line at 0.05
    # add horizontal line at maximum protein?

    # convert to numpy array, to list divide
    num_idpicker_protein_list = np.array(num_idpicker_protein_list)
    num_pyprophet_protein_list = np.array(num_pyprophet_protein_list)
    num_epifany_protein_list = np.array(num_epifany_protein_list)

    # percent increase from pyprophet
    percent_increase_idpicker_list = num_idpicker_protein_list / num_pyprophet_protein_list
    percent_increase_epifany_list = num_epifany_protein_list / num_pyprophet_protein_list

    plt.figure(1)
    plt.plot(threshold_list)
    plt.title("Percentage increase of the precise recall curves in uniprot")
    plt.xlabel("FDR threshold")
    plt.ylabel("Percentage")
    # plt.legend(loc='lower right')
    plt.legend(loc='center right')
    filename = 'percent ' + zoom_degree + '.pdf'
    plt.savefig(filename)
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
    # I compare if the fourth argument is 'yes' or not

    # TODO: I should use argparse probably

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], 'normal')
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], 'zoomed')
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], 'zoomed plus')


    # venn_diagram(sys.argv[1], sys.argv[2], sys.argv[3])
