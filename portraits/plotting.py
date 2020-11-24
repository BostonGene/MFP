import copy

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from portraits.utils import item_series, to_common_samples


def axis_net(x, y, title='', x_len=4, y_len=4, title_y=1, gridspec_kw=None):
    """
    Return an axis iterative for subplots arranged in a net
    :param x: int, number of subplots in a row
    :param y: int, number of subplots in a column
    :param title: str, plot title
    :param x_len: float, width of a subplot in inches
    :param y_len: float, height of a subplot in inches
    :param gridspec_kw: is used to specify axis ner with different rows/cols sizes.
            A dict: height_ratios -> list + width_ratios -> list
    :param title_y: absolute y position for suptitle
    :return: axs.flat, numpy.flatiter object which consists of axes (for further plots)
    """
    if x == y == 1:
        fig, ax = plt.subplots(figsize=(x * x_len, y * y_len))
        af = ax
    else:
        fig, axs = plt.subplots(y, x, figsize=(x * x_len, y * y_len), gridspec_kw=gridspec_kw)
        af = axs.flat

    fig.suptitle(title, y=title_y)
    return af


def lin_colors(factors_vector, cmap='default', sort=True, min_v=0, max_v=1, linspace=True):
    """
    Return dictionary of unique features of "factors_vector" as keys and color hexes as entries
    :param factors_vector: pd.Series
    :param cmap: matplotlib.colors.LinearSegmentedColormap, which colormap to base the returned dictionary on
        default - matplotlib.cmap.hsv with min_v=0, max_v=.8, lighten_color=.9
    :param sort: bool, whether to sort the unique features
    :param min_v: float, for continuous palette - minimum number to choose colors from
    :param max_v: float, for continuous palette - maximum number to choose colors from
    :param linspace: bool, whether to spread the colors from "min_v" to "max_v"
        linspace=False can be used only in discrete cmaps
    :return: dict
    """

    unique_factors = factors_vector.dropna().unique()
    if sort:
        unique_factors = np.sort(unique_factors)

    if cmap == 'default':
        cmap = matplotlib.cm.rainbow
        max_v = .92

    if linspace:
        cmap_colors = cmap(np.linspace(min_v, max_v, len(unique_factors)))
    else:
        cmap_colors = np.array(cmap.colors[:len(unique_factors)])

    return dict(list(zip(unique_factors, [matplotlib.colors.to_hex(x) for x in cmap_colors])))


def axis_net(x, y, title='', x_len=4, y_len=4, title_y=1, gridspec_kw=None):
    """
    Return an axis iterative for subplots arranged in a net
    :param x: int, number of subplots in a row
    :param y: int, number of subplots in a column
    :param title: str, plot title
    :param x_len: float, width of a subplot in inches
    :param y_len: float, height of a subplot in inches
    :param gridspec_kw: is used to specify axis ner with different rows/cols sizes.
            A dict: height_ratios -> list + width_ratios -> list
    :param title_y: absolute y position for suptitle
    :return: axs.flat, numpy.flatiter object which consists of axes (for further plots)
    """
    if x == y == 1:
        fig, ax = plt.subplots(figsize=(x * x_len, y * y_len))
        af = ax
    else:
        fig, axs = plt.subplots(y, x, figsize=(x * x_len, y * y_len), gridspec_kw=gridspec_kw)
        af = axs.flat

    fig.suptitle(title, y=title_y)
    return af


def pca_plot(data, grouping=None, order=(), n_components=2, ax=None, palette=None,
             alpha=1, random_state=42, s=20, figsize=(5, 5), title='',
             legend='in', **kwargs):
    kwargs_scatter = dict()
    kwargs_scatter['linewidth'] = kwargs.pop('linewidth', 0)
    kwargs_scatter['marker'] = kwargs.pop('marker', 'o')
    kwargs_scatter['edgecolor'] = kwargs.pop('edgecolor', 'black')

    if grouping is None:
        grouping = item_series('*', data)

    # Common samples
    c_data, c_grouping = to_common_samples([data, grouping])

    if len(order):
        group_order = copy.copy(order)
    else:
        group_order = np.sort(c_grouping.unique())

    if palette is None:
        cur_palette = lin_colors(c_grouping)
    else:
        cur_palette = copy.copy(palette)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # Get model and transform
    n_components = min(n_components, len(c_data.columns))
    from sklearn.decomposition import PCA
    model = PCA(n_components=n_components, random_state=random_state, **kwargs)

    data_tr = pd.DataFrame(model.fit_transform(c_data), index=c_data.index)

    label_1 = 'PCA 1 component {}% variance explained'.format(int(model.explained_variance_ratio_[0] * 100))
    label_2 = 'PCA 2 component {}% variance explained'.format(int(model.explained_variance_ratio_[1] * 100))

    kwargs_scatter = kwargs_scatter or {}
    for group in group_order:
        samples = list(c_grouping[c_grouping == group].index)
        ax.scatter(data_tr[0][samples], data_tr[1][samples], color=cur_palette[group], s=s, alpha=alpha,
                   label=str(group), **kwargs_scatter)

    if legend == 'out':
        ax.legend(scatterpoints=1, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.1)
    elif legend == 'in':
        ax.legend(scatterpoints=1)

    ax.set_title(title)
    ax.set_xlabel(label_1)
    ax.set_ylabel(label_2)

    return ax


def clustering_heatmap(ds, title='', corr='pearson', method='complete',
                       yl=True, xl=True,
                       cmap=matplotlib.cm.coolwarm, col_colors=None,
                       figsize=None, **kwargs):
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage

    dissimilarity_matrix = 1 - ds.T.corr(method=corr)
    hclust_linkage = linkage(squareform(dissimilarity_matrix), method=method)

    g = sns.clustermap(1 - dissimilarity_matrix, method=method,
                       row_linkage=hclust_linkage, col_linkage=hclust_linkage,
                       cmap=cmap, yticklabels=yl, xticklabels=xl,
                       col_colors=col_colors, figsize=figsize, **kwargs)

    g.fig.suptitle(title)

    return g


def patch_plot(patches, ax=None, order='sort', w=0.25, h=0, legend_right=True,
               show_ticks=False):
    cur_patches = pd.Series(patches)

    if order == 'sort':
        order = list(np.sort(cur_patches.index))

    data = pd.Series([1] * len(order), index=order[::-1])
    if ax is None:
        if h == 0:
            h = 0.3 * len(patches)
        _, ax = plt.subplots(figsize=(w, h))

    data.plot(kind='barh', color=[cur_patches[x] for x in data.index], width=1, ax=ax)
    ax.set_xticks([])
    if legend_right:
        ax.yaxis.tick_right()

    sns.despine(offset={'left': -2}, ax=ax)

    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if not show_ticks:
        ax.tick_params(length=0)

    return ax


def draw_graph(G, ax=None, title='', figsize=(12, 12), v_labels=True, e_labels=True, node_color='r', node_size=30,
               el_fs=5, nl_fs=8):
    """
    Draws a graph.
    :param G:
    :param ax:
    :param title:
    :param figsize:
    :param v_labels:
    :param e_labels:
    :param node_color:
    :param node_size:
    :param el_fs: edge label font size
    :param nl_fs: node label font size
    :return:
    """
    import networkx as nx

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    pos = nx.nx_pydot.graphviz_layout(G, prog="neato")
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_color)
    if v_labels:
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=nl_fs, font_family='sans-serif', font_color='blue')

    nx.draw_networkx_edges(G, pos, ax=ax)
    if e_labels:
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, ax=ax, font_size=el_fs, width=labels, edge_labels=labels)

    ax.set_title(title, fontsize=18)
    return ax
