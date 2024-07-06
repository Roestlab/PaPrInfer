from os import wait
from time import sleep


def extract_pp_from_file(file_path):
    pp_values = []
    is_decoy = False

    with open(file_path, 'r') as file:
        for line in file:
            if 'DECOY' in line:
                is_decoy = True
            elif 'pp:' in line:
                pp_value = float(line.strip().split()[1])
                pp_values.append((1 - pp_value, is_decoy))
                is_decoy = False  # Reset for the next entry

    return pp_values


def calculate_fdr_count_decoy(pp_values):
    # Sort the pp values in ascending order (lower confidence first)
    sorted_pp = sorted(pp_values)
    num_total = len(sorted_pp)
    fdr_values = []
    num_decoys = 0

    # Calculating FDR
    for i, pp in enumerate(sorted_pp, 1):
        if pp[1] and pp[0] != 0.0:  # If it is a decoy
            num_decoys += 1
        if not pp[1]: # If target
            fdr = num_decoys/ i
            fdr_values.append(fdr)

        # print(fdr, num_decoys, i)
        # sleep(1)

    # Adjust FDR values to ensure they are non-decreasing
    for i in range(len(fdr_values) - 2, -1, -1):
            fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])

    return fdr_values


def calculate_fdr_pep_sum(pp_values):
    # Sort the pp values in ascending order (higher confidence first)
    sorted_pp = sorted(pp_values)
    fdr_values = []
    sum_pep = 0
    num_targets = 0

    # Calculating FDR
    for i, pp in enumerate(sorted_pp, 1):
        if not pp[1]:  # If it's not a decoy
            sum_pep += pp[0]
            fdr = sum_pep / i
            fdr_values.append(fdr)

    # Adjust FDR values to ensure they are non-decreasing
    for i in range(len(fdr_values) - 2, -1, -1):
        fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])

    return fdr_values

# Path to your file
file_path = '/Users/kren/PycharmProjects/PaPrInfer/src/figures and files/protein_group_pep_qvalue.txt'




# Extracting pp values and calculating FDR
pp_values = extract_pp_from_file(file_path)
fdr_values = calculate_fdr_count_decoy(pp_values)



print(fdr_values[100])

# Creating an index list for the number of values (1 to the length of the lists)
index_list = list(range(1, len(fdr_values) + 1))



import matplotlib.pyplot as plt

# Creating the plots
# plt.plot(pp_values, index_list, label='PEP Curve', color='blue')
plt.plot(fdr_values, index_list, label='FDR Curve from PEP', color='red')
# plt.plot(q_values, index_list, label='EPIFANY FDR Curve', color='green')

# Adding labels, title, and legend
plt.ylabel('Number of Values')
plt.xlabel('PEP / FDR')
plt.title('PEP and FDR Curves')
plt.legend()

# Show the plot
plt.show()
