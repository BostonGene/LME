from lme.utils import update_gene_names
from lme.utils import read_dataset
import pandas as pd
from pathlib import Path


def run_progeny(expression_df, sync_gene_names=True, prog_coeffs=None, **kwargs):
    """
    Runs PROGENy pathway scoring on provided expressions dataframe in python
    :param expression_df: pd.DataFrame; rows - Hugo Gene symbols, columns - samples
    :param prog_coeffs: pd.DataFrame, progeny_genes_coefficients; index - HUGO gene symbols, columns - ['pathway', 'coefficient']
    :returns progeny pathway scores dataframe
    """
    if prog_coeffs is None:
        prog_coeffs = read_dataset(Path(__file__).resolve().parent.joinpath('databases', 'progeny_genes_coefficients.tsv'),
                                   index_col=None)
    if sync_gene_names:
        matching_genes = update_gene_names(genes_old=prog_coeffs['hugo_symbol'],
                                           genes_cur=expression_df.index, **kwargs)

        prog_coeffs = prog_coeffs.assign(hugo_symbol=prog_coeffs.hugo_symbol.map(matching_genes))

    coeffs = pd.pivot_table(prog_coeffs, index=['hugo_symbol'], columns=['pathway'], values='coefficient',
                            aggfunc=sum, fill_value=0)

    return expression_df.reindex(coeffs.index, fill_value=0).T.dot(coeffs).T
