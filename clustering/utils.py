import pandas as pd


class GeneSet(object):
    def __init__(self, name, descr, genes):
        self.name = name
        self.descr = descr
        self.genes = set(genes)
        self.genes_ordered = list(genes)

    def __str__(self):
        return '{}\t{}\t{}'.format(self.name, self.descr, '\t'.join(self.genes))


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
