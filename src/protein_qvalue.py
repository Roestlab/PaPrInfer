import sys

from get_proteins_at_threshold import get_epifany_result

from collections import Counter

import matplotlib.pyplot as plt


def main(epifany_file):
    qvalue_list_all = get_epifany_result(epifany_file, str(1), remove_decoy=False, return_qvalue=True)
    qvalue_list_target = get_epifany_result(epifany_file, str(1), remove_decoy=True, return_qvalue=True)

    print(type(qvalue_list_all[0]))

    qvalue_list_all_thresholded = [qvalue for qvalue in qvalue_list_all if qvalue >= 0.2]
    qvalue_list_target_thresholded = [qvalue for qvalue in qvalue_list_target if qvalue >= 0.2]

    plt.hist(qvalue_list_all_thresholded, bins=100)
    plt.hist(qvalue_list_target_thresholded, bins=100)
    plt.savefig('target.pdf')
    plt.show()




if __name__ == "__main__":
    # the epifany output file and a threshold
    main(sys.argv[1])
