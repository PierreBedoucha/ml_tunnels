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

    plt.scatter(x, y, marker='.')
    ax = plt.gca()
    d = scipy.zeros(len(y))
    ax.fill_between(x, y, where=y >= d, interpolate=True, facecolor='lightcyan', alpha=0.7)
    ax.fill_between(x, y, where=y <= d, interpolate=True, facecolor='peachpuff', alpha=0.7)


def facet_scatter_pos(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, marker='.')
    ax = plt.gca()
    d = scipy.zeros(len(y))
    ax.fill_betweeny(x, y, where=0 - y > 0 , interpolate=True, facecolor='lightcyan', alpha=0.7)


def facet_scatter_neg(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, marker='.')
    ax = plt.gca()
    d = scipy.zeros(len(y))
    ax.fill_betweeny(x, y, where=0 - y < 0, interpolate=True, facecolor='peachpuff', alpha=0.7)


if __name__ == '__main__':

    data_path = "./NATs_m7-12_a15_plane0-21A/"

    list_files = [x for x in os.listdir(data_path) if x.endswith(".csv")]
    pd_dict = {}

    for file_csv in list_files:
        # wtr = csv.writer(open(file_csv + '_.csv', 'w+'), delimiter=',', lineterminator='\n', quoting=csv.QUOTE_NONE,
        #                  escapechar='\\')
        # with open(file_csv) as f1_csv:
        #     for line in f1_csv.readlines():

        try:
            df = pd.read_csv(os.path.join(data_path + file_csv), sep=';',
                             converters={'SECTION AREA [Ų]': lambda x: (x.replace(",", ".")),
                                         'TIME STEP': lambda x: (x.replace(",", ".")),
                                         'OFFSET': lambda x: (x.replace(",", "."))})
        except pd.errors.EmptyDataError as err:
            continue

        types_dict = {'SECTION AREA [Ų]': float, 'TIME STEP': float, 'OFFSET': float}
        for col, col_type in types_dict.items():
            df[col] = df[col].astype(col_type)

        grouped = df.groupby(["TIME STEP"])
        train = grouped.apply(lambda x: x.sort_values(["OFFSET"], ascending=True)).reset_index(drop=True)

        # write to csv
        # df_filtered = train.loc[train['TIME STEP'].isin([0.0, 4.0, 6.0, 8.0])]
        df_filtered = train.loc[train['TIME STEP'].isin([0.0])]
        df_filtered['struct'] = ''.join(file_csv.split('_')[2])
        df_filtered['mode'] = int(file_csv.split('_')[3][1:])
        df_filtered['amplitude'] = float(file_csv.split('_')[4])
        # df_filtered.to_csv('filtered_{}.csv'.format('_'.join(file_csv.split('_')[:3])))

        pd_dict['_'.join(file_csv.split('_')[:3])] = df_filtered

        # h = sns.FacetGrid(train, col="OFFSET", palette='seismic', sharey=False, sharex=True, col_wrap=6, height=2,
        #                   aspect=1)
        #
        # h.map(sns.lineplot, "TIME STEP", "SECTION AREA [Ų]")
        #
        # plt.subplots_adjust(top=0.9)
        # h.fig.suptitle(file_csv)
        # plt.show()

    all = pd.concat(pd_dict.values(), axis=0, ignore_index=True)
    all.to_csv('nats_tunnel_area.csv')


    # res = all.groupby(['struct', 'mode']).apply(
    #     lambda g: g['SECTION AREA [Ų]'] -
    #               g[g.amplitude == 0.00]["SECTION AREA [Ų]"].values[0])
    # all['diff'] = res.reset_index(drop=True)

    def subspace(row, df):
        mode = row['mode']
        struct = row['struct']
        offset = row['OFFSET']
        air0 = df.query("mode == @mode and struct == @struct and amplitude == 0 and OFFSET == @offset")[
            "SECTION AREA [Ų]"]
        if air0.empty:
            airdiff = np.nan
        else:
            airdiff = row['SECTION AREA [Ų]'] - list(air0)[0]
        return airdiff


    all['diff'] = all.apply(lambda r: subspace(r, all), axis=1)

    col_list = list(set(all['struct'].tolist()))
    col_list.sort()

    for mode in [7, 8, 9, 10, 11, 12]:
        all_mode = all.loc[(all['mode'] == mode)]
        all_mode_0 = all_mode.loc[(all_mode['amplitude'] != 0.00)]

        h = sns.FacetGrid(all_mode_0, col="struct", sharex=True, col_wrap=3, height=2,
                          aspect=1, hue="amplitude", col_order=col_list)
        sns.color_palette("BuGn_r")
        # all["OFFSET"] = all["OFFSET"].astype(str)

        # h.map(sns.lineplot, "OFFSET", "SECTION AREA [Ų]")
        # h.map(sns.scatterplot, "OFFSET", "SECTION AREA [Ų]", data=all_mode)
        h.map(sns.lineplot, "OFFSET", "diff")
        h.map(facet_scatter, "OFFSET", "diff")
        # h.map(facet_scatter_pos, "OFFSET", "diff")
        # h.map(facet_scatter_neg, "OFFSET", "diff")

        plt.legend()

        plt.subplots_adjust(top=0.9)
        h.fig.suptitle(str(mode))
        plt.show()
