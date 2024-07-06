import sys

from pyprophet.stats import pi0est

import sqlite3
import plotly.graph_objects as go


def extract_idpicker_protein_pvalue(db_path):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Execute the query to select the PVALUE column from the SCORE_PROTEIN_GROUP table
    cursor.execute("SELECT PVALUE FROM SCORE_PROTEIN_GROUP")

    # Fetch all the results
    pvalues = cursor.fetchall()

    # Close the connection to the database
    conn.close()

    # Return a list of PVALUEs, extracting them from the tuple format returned by fetchall
    return [pvalue[0] for pvalue in pvalues]


def calculate_q_values(db_path, scaling_factor):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)

    # SQL query to merge tables and order by PEP
    query = """
    SELECT DISTINCT PROTEIN_GROUP.PROTEIN_GROUP_ID, PEP, DECOY, QVALUE
    FROM PROTEIN_GROUP
    INNER JOIN SCORE_PROTEIN_GROUP ON PROTEIN_GROUP.PROTEIN_GROUP_ID = SCORE_PROTEIN_GROUP.PROTEIN_GROUP_ID
    ORDER BY PEP
    """

    cursor = conn.cursor()
    cursor.execute(query)
    fetched_data = cursor.fetchall()

    # Process fetched data
    data = []
    decoy_counts = 0
    total_count = 0
    for protein_group_id, pep, decoy, pre_calculated_qvalue in fetched_data:
        total_count += 1
        if decoy == 1:
            decoy_counts += 1
        if decoy == 0:
            unadjusted_q_value = (decoy_counts / total_count)
            adjusted_q_value = (decoy_counts / total_count) * scaling_factor
            data.append((protein_group_id, adjusted_q_value, pre_calculated_qvalue, unadjusted_q_value))

    # Disconnect from database
    conn.close()

    # Sort data by adjusted q_value
    data_sorted = sorted(data, key=lambda x: x[1])  # sort by q_value

    # Calculate cumulative counts, for calculated
    calculated_q_values = []
    cumulative_counts_calc = []
    current_count = 0
    last_q_value = None
    for _, adjusted_q_value, _, _ in data_sorted:
        if adjusted_q_value != last_q_value:
            last_q_value = adjusted_q_value
            calculated_q_values.append(adjusted_q_value)
            cumulative_counts_calc.append(current_count)
        current_count += 1

    # Sort data by pre-calc q_value
    data_sorted = sorted(data, key=lambda x: x[2])  # sort by q_value

    # Calculate cumulative counts, for pre-calculated
    pre_calculated_q_values = []
    cumulative_counts_pre_calc = []
    current_count = 0
    last_q_value = None
    for _, _, pre_calculated_qvalue, _ in data_sorted:
        if pre_calculated_qvalue != last_q_value:
            last_q_value = pre_calculated_qvalue
            pre_calculated_q_values.append(pre_calculated_qvalue)
            cumulative_counts_pre_calc.append(current_count)
        current_count += 1

    # Sort data by unadjusted q_value
    data_sorted = sorted(data, key=lambda x: x[3])  # sort by q_value

    # Calculate cumulative counts, for unadjusted
    unadjusted_q_values = []
    cumulative_counts_unadjust = []
    current_count = 0
    last_q_value = None
    for _, _, _, unadjusted_q_value in data_sorted:
        if unadjusted_q_value != last_q_value:
            last_q_value = unadjusted_q_value
            unadjusted_q_values.append(unadjusted_q_value)
            cumulative_counts_unadjust.append(current_count)
        current_count += 1

    # Plotting using Plotly
    fig = go.Figure()

    # Add traces for different q-values

    # method 2 (used by epifany, implemented on idpicker data, adjusted to account
    # for underestimation of q-value of this method, stemming from the assumption
    # DIA data's decoy is equal in size to false targets)
    fig.add_trace(go.Scatter(x=calculated_q_values, y=cumulative_counts_calc, mode='lines+markers',
                             name='Q-Value from counting decoy with pi0'))
    # method 1 (default method by idpicker through pyprophet)
    fig.add_trace(go.Scatter(x=pre_calculated_q_values, y=cumulative_counts_pre_calc, mode='lines+markers',
                             name='Q-Value from pvalue with pi0'))
    # method 3 (used by epifany, implemented on idpicker data)
    fig.add_trace(go.Scatter(x=unadjusted_q_values, y=cumulative_counts_unadjust, mode='lines+markers',
                             name='Q-Value from counting decoy with no pi0'))

    # method 1 is the best, followed by 2, then 3

    # Update plot layout
    fig.update_layout(title='Comparison of Calculated vs. Pre-Calculated Q-Values for Idpicker inferred proteins',
                      xaxis_title='Q-Value',
                      yaxis_title='Cumulative Number of Protein Groups',
                      legend_title="Q-Value Type")
    fig.show()


# Example usage:
# db_path = 'path_to_your_database.db'
# fetch_and_plot_both_q_values(db_path)


if __name__ == "__main__":
    pvalue_list = extract_idpicker_protein_pvalue(sys.argv[1])
    pi0 = pi0est(pvalue_list)
    print(pi0)
    calculate_q_values(sys.argv[1], pi0['pi0'])
