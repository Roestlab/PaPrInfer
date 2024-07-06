import plotly.graph_objects as go

file_path_2 = '/Users/kren/PycharmProjects/PaPrInfer/src/figures and files/protein_group_pep_qvalue.txt'


def read_and_plot_pp_distribution_from_idxml(idxml_file_path):
    pass


def read_and_plot_pp_distribution_from_txt_file(txt_file_path):
    # Initialize lists to store pp values for decoys and targets
    pep_values_decoy = []
    pep_values_target = []

    # Open and read the file
    with open(txt_file_path, 'r') as file:
        # Split by empty line, which separates each protein entry
        entries = file.read().strip().split('\n\n')

    for entry in entries:
        lines = entry.strip().split('\n')
        is_decoy = lines[0].startswith('DECOY')
        for line in lines:
            if line.startswith('pp:'):
                pep_value = 1 - float(line.split(':')[1].strip())
                if is_decoy:
                    pep_values_decoy.append(pep_value)
                else:
                    pep_values_target.append(pep_value)

    # Plotting the distribution of pp values using Plotly
    fig = go.Figure()
    # Add histogram for decoy pp values
    fig.add_trace(go.Histogram(x=pep_values_decoy, name='Decoy PEP',
                               opacity=0.75))
    # Add histogram for target pp values
    fig.add_trace(go.Histogram(x=pep_values_target, name='Target PEP',
                               opacity=0.75))

    # Update the layout of the plot
    fig.update_layout(
        title='Distribution of Posterior Error Probabilities (pep) for Decoys and Targets',
        xaxis_title='Posterior Error Probability (pep)',
        yaxis_title='Count',
        bargap=0.2,  # gap between bars of adjacent location coordinates
        bargroupgap=0.1,  # gap between bars of the same location coordinate
        barmode='overlay'  # bars are drawn on top of one another for comparison
    )
    fig.show()


if __name__ == "__main__":
    read_and_plot_pp_distribution_from_txt_file(file_path_2)
