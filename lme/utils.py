import warnings

import numpy as np
import pandas as pd


class GeneSet(object):
    def __init__(self, name, descr, genes):
        self.name = name
        self.descr = descr
        self.genes = set(genes)
        self.genes_ordered = list(genes)

    def __str__(self):
        s = ','.join(self.genes)
        return f'{self.name} ({self.descr}): {s}'

    def __repr__(self):
        return self.__str__()


def read_gene_sets(gmt_file):
    """
    Return dict {geneset_name : GeneSet object}

    :param gmt_file: str, path to .gmt file
    :return: dict
    """
    gene_sets = {}
    with open(gmt_file) as handle:
        for line in handle:
            items = line.strip().split('\t')
            name = items[0].strip()
            description = items[1].strip()
            genes = set([gene.strip() for gene in items[2:]])
            gene_sets[name] = GeneSet(name, description, genes)

    return gene_sets


def ssgsea_score(ranks, genes):
    common_genes = list(set(genes).intersection(set(ranks.index)))
    if not len(common_genes):
        return pd.Series([0] * len(ranks.columns), index=ranks.columns)
    sranks = ranks.loc[common_genes]
    return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (len(ranks.index) - len(common_genes) + 1) / 2


def ssgsea_formula(data, gene_sets, rank_method='max'):
    """
    Return DataFrame with ssgsea scores
    Only overlapping genes will be analyzed

    :param data: pd.DataFrame, DataFrame with samples in columns and variables in rows
    :param gene_sets: dict, keys - processes, values - bioreactor.gsea.GeneSet
    :param rank_method: str, 'min' or 'max'.
    :return: pd.DataFrame, ssgsea scores, index - genesets, columns - patients
    """

    ranks = data.T.rank(method=rank_method, na_option='bottom')

    return pd.DataFrame({gs_name: ssgsea_score(ranks, gene_sets[gs_name].genes)
                         for gs_name in list(gene_sets.keys())})


def median_scale(data, clip=None):
    c_data = (data - data.median()) / data.mad()
    if clip is not None:
        return c_data.clip(-clip, clip)
    return c_data


def read_dataset(file, sep='\t', header=0, index_col=0, comment=None):
    return pd.read_csv(file, sep=sep, header=header, index_col=index_col,
                       na_values=['Na', 'NA', 'NAN'], comment=comment)


def item_series(item, indexed=None):
    """
    Creates a series filled with item with indexes from indexed (if Series-like) or numerical indexes (size=indexed)
    :param item: value for filling
    :param indexed:
    :return:
    """
    if indexed is not None:
        if hasattr(indexed, 'index'):
            return pd.Series([item] * len(indexed), index=indexed.index)
        elif type(indexed) is int and indexed > 0:
            return pd.Series([item] * indexed, index=np.arange(indexed))
    return pd.Series()


def to_common_samples(df_list=()):
    """
    Accepts a list of dataframes. Returns all dataframes with only intersecting indexes
    :param df_list: list of pd.DataFrame
    :return: pd.DataFrame
    """
    cs = set(df_list[0].index)
    for i in range(1, len(df_list)):
        cs = cs.intersection(df_list[i].index)

    if len(cs) < 1:
        warnings.warn('No common samples!')
    return [df_list[i].loc[list(cs)] for i in range(len(df_list))]


def query_genes_by_symbol(genes, verbose=False):
    """

    :param genes:
    :param verbose:
    :return:
    """
    import mygene
    mg = mygene.MyGeneInfo()
    q = mg.querymany(genes, species='human', as_dataframe=True, verbose=verbose, df_index=True,
                     scopes=["symbol"], fields="all")
    try:
        q.dropna(subset=['HGNC'], inplace=True)
        q.dropna(subset=['type_of_gene'], inplace=True)
        q.dropna(subset=['map_location'], inplace=True)
    except Exception:
        warnings.warn('Output lacks map_location or type_of_gene')

    return q


def update_gene_names(genes_old, genes_cur, verbose=False):
    """
    Takes a set of gene names genes_old and matches it with genes_cur.
    For all not found tries to match with known aliases using mygene or GeneStorage.
    All not matched are returned as is.
    Returns a dict with matching rule. No duplicates will be in output
    :param genes_old:
    :param genes_cur:
    :param verbose:
    :return:
    """
    c_genes = set(genes_cur)
    old_genes = set(genes_old)

    missing = set()

    common_genes = c_genes.intersection(old_genes)
    if verbose:
        print('Matched: {}'.format(len(common_genes)))

    converting_genes = old_genes.difference(c_genes)
    rest_genes = c_genes.difference(old_genes)
    match_rule = {cg: cg for cg in common_genes}

    if len(converting_genes):

        if verbose:
            print('Trying to find new names for {} genes in {} known'.format(len(converting_genes),
                                                                             len(rest_genes)))

        qr = query_genes_by_symbol(list(converting_genes), verbose=verbose)
        if hasattr(qr, 'alias'):
            cg_ann = qr.alias.dropna()
        else:
            cg_ann = pd.DataFrame()

        for cg in converting_genes:
            if cg in cg_ann.index:
                if (isinstance(cg_ann.loc[cg], list)) | (isinstance(cg_ann.loc[cg], pd.core.series.Series)):
                    al_set = set(cg_ann[cg])
                else:
                    al_set = set([cg_ann.loc[cg]])

                hits = al_set.intersection(rest_genes)
                if len(hits) == 1:
                    match_rule[cg] = list(hits)[0]
                    rest_genes.remove(match_rule[cg])
                elif len(hits) > 1:
                    warnings.warn('{} hits for gene {}'.format(len(hits), cg))
                    match_rule[cg] = list(hits)[0]
                    rest_genes.remove(match_rule[cg])
                else:
                    missing.add(cg)
                    match_rule[cg] = cg
            else:
                missing.add(cg)
                match_rule[cg] = cg
        if verbose and len(missing):
            print('{} genes were not converted'.format(len(missing)))
    return match_rule
