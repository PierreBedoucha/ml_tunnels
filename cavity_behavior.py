import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance
import numpy as np
import scipy
import Bio.PDB
import Bio.AlignIO as al
from sklearn import preprocessing
import csv


def facet_scatter(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")

    # # prepending and appending 0.0 values to close the shape for filling
    # x[x.index[-1] + 1] = 0.0
    # y[y.index[-1] + 1] = 0.0
    # b = [0.00]
    # b[1:] = x
    # x = pd.Series(b)
    # c = [0.00]
    # c[1:] = y
    # y = pd.Series(c)

    # plt.scatter(x, y, marker='.')
    ax = plt.gca()
    d0 = scipy.zeros(len(y))
    ax.fill_between(x, y, where=y >= d0, interpolate=True, facecolor='limegreen', alpha=0.4, label="Δ positive")
    ax.fill_between(x, y, where=y <= d0, interpolate=True, facecolor='indianred', alpha=0.4, label="Δ negative")
    # xint = range(min(x), math.ceil(max(x)) + 1)
    # plt.xticks(xint)


def facet_auc(x, z, w1, **kwargs):
    kwargs.pop("color")
    ax = plt.gca()
    z1 = np.array(z)
    w2 = np.array(w1)
    ax.fill_between(x, z, w1, where=z1 >= w2, interpolate=True, facecolor='limegreen', alpha=0.4, label="Δ positive")
    ax.fill_between(x, z, w1, where=z1 <= w2, interpolate=True, facecolor='indianred', alpha=0.4, label="Δ negative")


def facet_line(x, y, z, **kwargs):
    kwargs.pop("color")
    if not z.any():
        sns.lineplot(x=x, y=y, linewidth=2, color='black', label='initial structure')
    elif -15.0 in z.values:
        sns.lineplot(x=x, y=y, linewidth=1.3, label='-15')
    elif 15.0 in z.values:
        sns.lineplot(x=x, y=y, linewidth=1.3, label='+15')
    # xint = range(min(x), math.ceil(max(x)) + 1)
    # plt.xticks(xint)


def subspace_ref(row, df):
    mode = row['mode']
    struct = row['struct']
    offset = row['OFFSET']
    air0 = df.query('mode == @mode and struct == @struct and amplitude == 0 and OFFSET == @offset')[
        "SECTION AREA [$\AA^2$]"]
    if air0.empty:
        airdiff = np.nan
    else:
        airdiff = list(air0)[0]
    return airdiff


if __name__ == '__main__':

    data_path = "./NATs_m7-12_a15_plane-15_25/"

    list_files = [x for x in os.listdir(data_path) if x.endswith(".csv")]
    pd_dict = {}

    for file_csv in list_files:
        # wtr = csv.writer(open(file_csv + '_.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
        #                  escapechar='\\')
        # with open(file_csv) as f1_csv:
        #     for line in f1_csv.readlines():

        with open(os.path.join(data_path, file_csv)) as f:
            lines = f.readlines()

        line_array = lines[0].split(";")
        # line_array[-1] = "OFFSET [$\\AA$]\n"
        line_first_el_array = line_array[0].split(" ")
        line_first_el_array[-1] = "[$\\AA^2$]"
        lines[0] = ' '.join(line_first_el_array) + ';' + ';'.join(line_array[1:])
        with open(os.path.join(data_path, file_csv), "w") as f:
            f.writelines(lines)

        try:
            df = pd.read_csv(os.path.join(data_path + file_csv), sep=';',
                             converters={'SECTION AREA [$\AA^2$]': lambda x: (x.replace(",", ".")),
                                         'TIME STEP': lambda x: (x.replace(",", ".")),
                                         'OFFSET': lambda x: (x.replace(",", "."))})
        except pd.errors.EmptyDataError as err:
            print(err.args[0] + " {0}".format(file_csv))
            continue

        types_dict = {'SECTION AREA [$\AA^2$]': float, 'TIME STEP': float, 'OFFSET': float}
        for col, col_type in types_dict.items():
            df[col] = df[col].astype(col_type)

        grouped = df.groupby(["TIME STEP"])
        train = grouped.apply(lambda x: x.sort_values(["OFFSET"], ascending=True)).reset_index(drop=True)

        # write to csv
        # df_filtered = train.loc[train['TIME STEP'].isin([0.0, 4.0, 6.0, 8.0])]
        df_filtered = train.loc[train['TIME STEP'].isin([1.0, 2.0, 3.0])]
        df_filtered['struct'] = ''.join(file_csv.split('_')[2])
        df_filtered['mode'] = int(file_csv.split('_')[3][1:])
        try:
            df_filtered['amplitude'] = float(file_csv.split('_')[4])
        except ValueError as err:
            d = {1.0: -15.0, 2.0: 0.0, 3.0: 15.0}
            df_filtered['amplitude'] = train['TIME STEP'].map(d)

        # Add ref area list column through the whole df
        df_filtered['NEW'] = df_filtered.apply(lambda r: subspace_ref(r, df_filtered), axis=1)

        # df_filtered.to_csv('filtered_{}.csv'.format('_'.join(file_csv.split('_')[:3])))

        pd_dict['_'.join(file_csv.split('_')[:3])] = df_filtered

        # h = sns.FacetGrid(train, col="OFFSET", palette='seismic', sharey=False, sharex=True, col_wrap=6, height=2,
        #                   aspect=1)
        #
        # h.map(sns.lineplot, "TIME STEP", "SECTION AREA [$\AA^2$]")
        #
        # plt.subplots_adjust(top=0.9)
        # h.fig.suptitle(file_csv)
        # plt.show()

    all = pd.concat(pd_dict.values(), axis=0, ignore_index=True)
    all.to_csv('nats_tunnel_area.csv')


    # res = all.groupby(['struct', 'mode']).apply(
    #     lambda g: g['SECTION AREA [$\AA^2$]'] -
    #               g[g.amplitude == 0.00]["SECTION AREA [$\AA^2$]"].values[0])
    # all['diff'] = res.reset_index(drop=True)

    def subspace(row, df):
        mode = row['mode']
        struct = row['struct']
        offset = row['OFFSET']
        air0 = df.query("mode == @mode and struct == @struct and amplitude == 0 and OFFSET == @offset")[
            "SECTION AREA [$\AA^2$]"]
        if air0.empty:
            airdiff = np.nan
        else:
            airdiff = row['SECTION AREA [$\AA^2$]'] - list(air0)[0]
        return airdiff


    all['diff'] = all.apply(lambda r: subspace(r, all), axis=1)

    col_list = list(set(all['struct'].tolist()))
    col_list.sort()

    tunnel_struct_neg_dict = {
        '3tfy': -4,
        '4kvm': -3,
        '4u9v': -15,
        '5icv': -4,
        '5k18': -6,
        '5wjd': -6
    }

    tunnel_struct_pos_dict = {
        '3tfy': 15,
        '4kvm': 15,
        '4u9v': 9,
        '5icv': 15,
        '5k18': 17,
        '5wjd': 17,
    }

    # for mode in [7, 8, 9, 10, 11, 12]:
    # for mode in [7]:
    # all_mode = all.loc[(all['mode'] == mode)]
    all_mode = all

    # all_mode_0 = all_mode.loc[(all_mode['amplitude'] != 0.00)]
    # all_mode_filter = all_mode.loc[(all_mode['OFFSET'] >= tunnel_struct_dict[all_mode['struct']][0])]
    all_mode_filter_neg = all_mode[all_mode['OFFSET'] >= all_mode['struct'].map(tunnel_struct_neg_dict)].copy()
    all_mode_filter = all_mode_filter_neg[all_mode_filter_neg['OFFSET'] <= all_mode_filter_neg['struct']
        .map(tunnel_struct_pos_dict)].copy()

    # h = sns.FacetGrid(all_mode_filter, col="struct", sharex=False, col_wrap=3, height=2,
    #                   aspect=1, hue="amplitude", col_order=col_list)

    # sns.set(font_scale=1)
    sns.set_style("whitegrid", {'axes.grid': False})
    h = sns.FacetGrid(all_mode_filter, col="struct", row="mode", sharex=False, sharey=False, height=6,
                      aspect=1, hue="amplitude",
                      gridspec_kws={"hspace": 0.6})

    sns.color_palette("BuGn_r")
    # all["OFFSET"] = all["OFFSET"].astype(str)

    # h.map(sns.lineplot, "OFFSET", "SECTION AREA [$\AA^2$]")
    h.map(facet_scatter, "OFFSET", "diff")

    # Series for SectionArea for amplitude 0.0
    # w0 = all_mode[(all_mode.amplitude.eq(0.00))]["SECTION AREA [$\AA^2$]"]

    h.map(facet_line, "OFFSET", "SECTION AREA [$\AA^2$]", "amplitude")
    # h.map(facet_auc, "OFFSET", "SECTION AREA [$\AA^2$]", "NEW")

    # h.map(sns.scatterplot, "OFFSET", "SECTION AREA [$\AA^2$]", data=all_mode)
    # h.map(sns.lineplot, "OFFSET", "diff")
    # h.map(facet_scatter, "OFFSET", "diff")

    # hand, labl = plt.gca().get_legend_handles_labels()
    # plt.legend(labl[:-4], framealpha=0.2)

    # this surpresses the x- and y-labels on each axes of the bottom/leftmost column
    h.set_axis_labels('', '')

    h.set_titles('{col_name}' + ' mode ' + '{row_name}',  weight='bold')

    # # overall ylabel
    # h.fig.text(x=0.02, y=0.5,
    #            verticalalignment='center',  # make sure it's aligned at center vertically
    #            s='SECTION AREA [$\AA^2$]',  # this is the text in the ylabel
    #            size=10,  # customize the fontsize if you will
    #            rotation=90)  # vertical text
    #
    # # overall xlabel
    # h.fig.text(x=0.5, y=0,
    #            horizontalalignment='center',  # make sure it's aligned at center horizontally
    #            s='OFFSET [$\AA$]',  # this is the text in the xlabel
    #            size=10)

    plt.subplots_adjust(top=0.975, bottom=0.025, left=0.025, right=0.975)
    # for ax in h.axes.ravel():
    #     ax.set_xticklabels(ax.get_xticklabels(), fontsize=2)

    # axes = g.axes.flatten()
    # axes[0].set_title("Internal")
    # axes[1].set_title("External")

    # h.fig.suptitle("mode " + str(mode))
    # plt.tight_layout()
    plt.show()
