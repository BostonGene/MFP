import warnings

import community  # louvain
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm


def gen_graph(similarity_matrix, threshold=0.8):
    """
    Generates a graph from the similarity_matrix (square dataframe). Each sample is a node, similarity - edge weight.
    Edges with weight lower than the threshold are ignored.
    Only nodes with at least 1 edge with weight above the threshold will be present in the final graph
    :param similarity_matrix:
    :param threshold:
    :return:
    """
    G = nx.Graph()
    for col_n in similarity_matrix:
        col = similarity_matrix[col_n].drop(col_n)
        mtr = col[col > threshold]
        for row_n, val in list(mtr.to_dict().items()):
            G.add_edge(col_n, row_n, weight=round(val, 2))
    return G


def louvain_community(correlation_matrix, threshold=1.4, resolution=1, random_state=100, **kwargs):
    """
    Generates a graph from a correlation_matrix with weighted edges (weight<threshold excluded).
    Then detects communities using louvain algorithm (https://github.com/taynaud/python-louvain)
    :param correlation_matrix:
    :param threshold:
    :param resolution:
    :param random_state:
    :param kwargs:
    :return:
    """
    G = gen_graph(correlation_matrix, threshold)
    if len(G.edges()) > 10000000:
        warnings.warn("Too many edges will result in huge computational time")
    return pd.Series(community.best_partition(G, resolution=resolution,
                                              random_state=random_state, **kwargs)) + 1


def dense_clustering(data, threshold=0.4, name='MFP', method='louvain', **kwargs):
    """
    Generates a graph from the table features(cols)*samples(rows).
    Then performs community detection using a selected method leiden|louvain
    :param method:
    :param data:
    :param threshold:
    :param name:
    :return:
    """
    if method == 'louvain':  # others could be implemented like Leiden
        partition = louvain_community(data.T.corr() + 1, threshold + 1, **kwargs)
    else:
        raise Exception('Unknown method')

    return partition.rename(name)


def clustering_profile_metrics(data, threshold_mm=(0.3, 0.65), step=0.01, method='louvain'):
    """
    Iterates threshold in threshold_mm area with step. Calculates cluster separation metrics on each threshold.
    Returns a pd.DataFrame with the metrics
    :param data:
    :param threshold_mm:
    :param step:
    :param method:
    :return:
    """
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    cluster_metrics = {}

    for tr in tqdm(np.round(np.arange(threshold_mm[0], threshold_mm[1], step), 3)):
        clusters_comb = dense_clustering(data, threshold=tr, method=method)
        cluster_metrics[tr] = {
            'ch': calinski_harabasz_score(data.loc[clusters_comb.index], clusters_comb),
            'db': davies_bouldin_score(data.loc[clusters_comb.index], clusters_comb),
            'sc': silhouette_score(data.loc[clusters_comb.index], clusters_comb),
            'N': len(clusters_comb.unique()),
            'perc': clusters_comb,
        }

    return pd.DataFrame(cluster_metrics).T


def clustering_profile_metrics_plot(cluster_metrics, num_clusters_ylim_max=7):
    """
    Plots a dataframe from clustering_profile_metrics
    :param cluster_metrics:
    :param num_clusters_ylim_max:
    :return: axis array
    """
    # necessary for correct x axis sharing
    cluster_metrics.index = [str(x) for x in cluster_metrics.index]

    plots_ratios = [3, 3, 3, 1, 2]
    fig, axs = plt.subplots(len(plots_ratios), 1, figsize=(8, np.sum(plots_ratios)),
                            gridspec_kw={'height_ratios': plots_ratios}, sharex=True)
    for ax in axs:
        ax.tick_params(axis='x', which='minor', length=0)
    af = axs.flat

    ax = cluster_metrics.db.plot(ax=next(af), label='Davies Bouldin', color='#E63D06')
    ax.legend()

    ax = cluster_metrics.ch.plot(ax=next(af), label='Calinski Harabasz', color='#E63D06')
    ax.legend()

    ax = cluster_metrics.sc.plot(ax=next(af), label='Silhouette score', color='#E63D06')
    ax.legend()

    ax = cluster_metrics.N.plot(kind='line', ax=next(af), label='# clusters', color='#000000')
    ax.set_ylim(0, num_clusters_ylim_max)
    ax.legend()

    # display percentage for 10 clusters max
    clusters_perc = pd.DataFrame([x.value_counts() for x in cluster_metrics.perc],
                                 index=cluster_metrics.index).iloc[:, :10]

    clusters_perc.plot(kind='bar', stached=True, ax=next(af), offset=.5)

    ax.set_xticks(ax.get_xticks() - .5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    ax.set_ylabel('Cluster %')
    return ax


def clustering_select_best_tr(data, n_clusters=4, threshold_mm=(0.3, 0.6),
                              step=0.025, method='leiden', num_clusters_ylim_max=7, plot=True):
    """
    Selects the best threshold for n_clusters separation using dense_clustering with selected method
        from threshold_mm with a particular step
    :param data: dataframe with processes (rows - samples, columns - signatures)
    :param n_clusters: desired number of clusters
    :param threshold_mm: range of thresholds
    :param step: step to go through range of thresholds
    :param method: clustering method
    :param num_clusters_ylim_max: set y_lim for plot with number of clusters
    :param plot: whether to plot all matrix
    :return: the threshold to get n_clusters
    """
    cl_scs = clustering_profile_metrics(data, threshold_mm=threshold_mm, step=step, method=method)

    if plot:
        clustering_profile_metrics_plot(cl_scs, num_clusters_ylim_max)
        plt.show()

    cl_scs_filtered = cl_scs[cl_scs.N == n_clusters]

    if not len(cl_scs_filtered):
        raise Exception('No partition with n_clusters = {}'.format(n_clusters))

    cl_scs_filtered.sc += 1 - cl_scs_filtered.sc.min()
    return (cl_scs_filtered.ch / cl_scs_filtered.db / cl_scs_filtered.sc).sort_values().index[-1]
