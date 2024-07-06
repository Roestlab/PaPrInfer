import sqlite3
import sys

import plotly.graph_objects as go


def plot_pep_distribution_protein_group(db_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)

    # SQL query to fetch PEP values and decoy status
    query = """
    SELECT PEP, DECOY
    FROM SCORE_PEPTIDE
    JOIN PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE.ID
    """
    cursor = conn.cursor()
    cursor.execute(query)
    fetched_data = cursor.fetchall()

    # Disconnect from database
    conn.close()

    # Separate PEP values by decoy status
    decoy_pep = [pep for pep, decoy in fetched_data if decoy == 1]
    target_pep = [pep for pep, decoy in fetched_data if decoy == 0]

    # Create a figure with subplots
    fig = go.Figure()

    # Add trace for decoy PEP values
    fig.add_trace(go.Histogram(x=decoy_pep, name='Decoy PEP', opacity=0.75))
    # Add trace for target PEP values
    fig.add_trace(go.Histogram(x=target_pep, name='Target PEP', opacity=0.75))

    # Update the layout of the plot
    fig.update_layout(
        title_text='Distribution of PEP Values for Decoy and Target Peptides',  # title of plot
        xaxis_title_text='PEP',  # x-axis label
        yaxis_title_text='Count',  # y-axis label
        bargap=0.2,  # gap between bars of adjacent location coordinates
        bargroupgap=0.1,  # gap between bars of the same location coordinate
        barmode='overlay'  # bars are drawn on top of one another for comparison
    )

    # Show the figure
    fig.show()

# Example usage:
# db_path = 'path_to_your_database.db'
# plot_pep_distribution(db_path)


if __name__ == "__main__":
    pvalue_list = plot_pep_distribution_protein_group(sys.argv[1])
