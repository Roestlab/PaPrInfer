import sys
import click
import pandas as pd
import sqlite3
import numpy as np

from shutil import copyfile
from pyprophet.report import save_report
from pyprophet.stats import error_statistics, final_err_table, \
    lookup_values_from_error_table, summary_err_table


def statistics_report(data, outfile, context, analyte, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    error_stat, pi0 = error_statistics(data[data.decoy == 0]['score'],
                                       data[data.decoy == 1]['score'],
                                       parametric, pfdr, pi0_lambda, pi0_method,
                                       pi0_smooth_df, pi0_smooth_log_pi0, True,
                                       lfdr_truncate, lfdr_monotone,
                                       lfdr_transformation, lfdr_adj, lfdr_eps)

    stat_table = final_err_table(error_stat)
    summary_table = summary_err_table(error_stat)

    # print summary table
    click.echo("=" * 80)
    click.echo(summary_table)
    click.echo("=" * 80)

    p_values, s_values, peps, q_values \
        = lookup_values_from_error_table(data["score"].values, error_stat)
    data["p_value"] = p_values
    data["s_value"] = s_values
    data["q_value"] = q_values
    data["pep"] = peps

    # export PDF report
    save_report(outfile + "_" + context + "_" + analyte + ".pdf", outfile + ": "
                + context + " " + analyte + "-level error-rate control",
                data[data.decoy==1]["score"], data[data.decoy==0]["score"],
                stat_table["cutoff"], stat_table["svalue"],
                stat_table["qvalue"], data[data.decoy==0]["p_value"], pi0)

    return data


def infer_protein_groups(infile, outfile, context, parametric, pfdr, pi0_lambda,
                   pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate,
                   lfdr_monotone, lfdr_transformation, lfdr_adj, lfdr_eps):

    con = sqlite3.connect(infile)

    data = pd.read_sql_query('''
        SELECT PROTEIN_GROUP.PROTEIN_GROUP_ID AS PROTEIN_GROUP_ID, SCORE, DECOY
        FROM PROTEIN_GROUP
        ''', con)

    data.columns = [col.lower() for col in data.columns]
    con.close()

    print(data[data.decoy == 0]['score'])

    # this is technically speaking, only has 1 run,
    # so that run specific and global and experient wide is the same
    # TODO: change this to be able to accommodate run specfic

    data = statistics_report(data, outfile, context, "protein group", parametric,
                             pfdr, pi0_lambda, pi0_method, pi0_smooth_df,
                             pi0_smooth_log_pi0, lfdr_truncate,
                             lfdr_monotone, lfdr_transformation, lfdr_adj,
                             lfdr_eps)

    # store data in table
    if infile != outfile:
        copyfile(infile, outfile)

    con = sqlite3.connect(outfile)

    c = con.cursor()
    c.execute(
        'SELECT count(name) FROM sqlite_master WHERE type="table" '
        'AND name="SCORE_PROTEIN_GROUP"'
    )
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_PROTEIN_GROUP')
    c.fetchall()

    df = data[['protein_group_id', 'score', 'p_value', 'q_value',
               'pep']]
    df.columns = ['PROTEIN_GROUP_ID', 'SCORE', 'PVALUE',
                  'QVALUE', 'PEP']
    table = "SCORE_PROTEIN_GROUP"
    df.to_sql(table, con, index=False, if_exists='append')

    con.close()


if __name__ == "__main__":
    # sys.argv 1 is just the file name, this function takes the infile name and
    # outfile name, but i set it as the same one anyway
    infer_protein_groups(sys.argv[1], sys.argv[1], 'run-specific', False,
                         False, np.arange(0.05,1.0,0.05), 'bootstrap', 3, False, True,
                         True, 'probit', 1.5, np.power(10.0, -8))
