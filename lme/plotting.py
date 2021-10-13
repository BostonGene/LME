import copy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from lme.utils import item_series, to_common_samples


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


def patch_plot(patches, ax=None, order='sort', w=0.25, h=0, legend_right=True,
               show_ticks=False):
    """
    Plot a palette
    :param patches: dict {value -> color}
    :return:
    """
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


def axis_matras(ys, title='', x_len=8, title_y=1, sharex=True):
    """
    Return an axis iterative for subplots stacked vertically
    :param ys: list, list of lengths by 'y'
    :param title: str, title for plot
    :param x_len: int, length by 'x'
    :param sharex: boolean, images will be shared if True
    :param title_y: absolute y position for suptitle
    :return: axs.flat, numpy.flatiter object which consists of axes (for further plots)
    """
    fig, axs = plt.subplots(len(ys), 1, figsize=(x_len, np.sum(ys)), gridspec_kw={'height_ratios': ys}, sharex=sharex)
    fig.suptitle(title, y=title_y)

    for ax in axs:
        ax.tick_params(axis='x', which='minor', length=0)

    return axs.flat


def line_palette_annotation_plot(val_vector, palette, ax=None, nan_color='#ffffff',
                                 hide_ticks=True, hide_borders=True, **kwargs):
    """
    Draws line annotation plot
    :param val_vector: pd.Series with values
    :param palette: dict, palette for values
    :param ax: ax to plot
    :param nan_color: str, color for np.nan
    :param hide_ticks: bool, whether to plot ticks
    :param hide_borders: bool, whether to plot borders
    :return: ax with plot
    """
    return line_annotation_plot(val_vector.map(palette), ax=ax, nan_color=nan_color,
                                hide_ticks=hide_ticks, hide_borders=hide_borders, **kwargs)


def line_annotation_plot(color_vector, ax=None, nan_color='#ffffff', offset=0, hide_ticks=True, hide_borders=True):
    """

    :param color_vector:
    :param ax:
    :param nan_color:
    :param offset:
    :return:
    """

    if ax is None:
        _, ax = plt.subplots(figsize=(max(len(color_vector) / 15.0, 6), 0.5))

    items_amount = len(color_vector)

    xss = np.arange(items_amount) - offset
    yss = pd.Series([1] * items_amount, index=color_vector.index)

    with sns.axes_style("white"):
        ax.bar(
            xss,
            yss,
            color=color_vector.fillna(nan_color),
            width=1,
            align='edge',
            edgecolor=color_vector.fillna(nan_color),
        )

    ax.set_ylim(0, 1)
    ax.set_xlim(0, items_amount)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.xaxis.label.set_visible(False)
    ax.set_ylabel(color_vector.name, rotation=0, labelpad=10, va='center', ha='right')

    if hide_ticks:
        ax.tick_params(length=0)

    if hide_borders:
        for spine in ['bottom', 'top', 'left', 'right']:
            ax.spines[spine].set_visible(False)

    return ax


def vector_pie_plot(data, ax=None, figsize=(4, 4), title='', palette=None, display_counts=False, order=None):
    """
    Constructs pie plot by provided pd.Series
    :param data: pd.Series
    :param ax: matplotlib axis, axis to plot on
    :param figsize: (float, float), figure size in inches
    :param title: str, plot title
    :param palette: dict, palette for plotting. Keys are unique values from groups, entries are color hexes
    :param display_counts: bool
    :param order: list, order to display groups
    :return: matplotlib axis
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    order = order or list(data.unique())

    c_data = data.value_counts()
    c_data = c_data[[x for x in order if x in c_data.index]]

    if palette is not None:
        c_colors = pd.Series(palette)[c_data.index]
    else:
        c_colors = None

    if display_counts:
        actopcl_rule = lambda p: '{:.0f}'.format(p * sum(c_data.values) / 100)
    else:
        actopcl_rule = '%1.1f%%'

    _, _, text_props = ax.pie(
        c_data, labels=c_data.index, autopct=actopcl_rule, startangle=0, textprops={'fontsize': 14}, colors=c_colors
    )

    for i in text_props:
        i.set_color('#ffffff')
    ax.axis('equal')
    ax.set_title(title)
    ax.set_xlabel(data.name)
    return ax


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
    import matplotlib

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
