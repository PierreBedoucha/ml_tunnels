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
    ax.fill_between(x, y, where=y >= d0, interpolate=True, facecolor='indianred', alpha=0.4)
    ax.fill_between(x, y, where=y <= d0, interpolate=True, facecolor='limegreen', alpha=0.4)


def facet_auc(x, z, w1, **kwargs):
    kwargs.pop("color")
    ax = plt.gca()
    z1 = np.array(z)
    w2 = np.array(w1)
    ax.fill_between(x, z, w1, where=z1 >= w2, interpolate=True, facecolor='indianred', alpha=0.5)
    ax.fill_between(x, z, w1, where=z1 <= w2, interpolate=True, facecolor='limegreen', alpha=0.5)


def facet_line(x, y, z, **kwargs):
    kwargs.pop("color")
    if not z.any():
        sns.lineplot(x=x, y=y, linewidth=3, color='black', label='initial structure')
    elif -15.0 in z.values:
        sns.lineplot(x=x, y=y, linewidth=1.2, set_markerfacecolor="white", label='negative')
    elif 15.0 in z.values:
        sns.lineplot(x=x, y=y, linewidth=1.2, set_markerfacecolor="white", label='positive')


def subspace_ref(row, df):
    mode = row['mode']
    struct = row['struct']
    offset = row['OFFSET']
    air0 = df.query("mode == @mode and struct == @struct and amplitude == 0 and OFFSET == @offset")[
        "SECTION AREA [Ų]"]
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
    # for mode in [7]:
        all_mode = all.loc[(all['mode'] == mode)]
        # all_mode_0 = all_mode.loc[(all_mode['amplitude'] != 0.00)]

        h = sns.FacetGrid(all_mode, col="struct", sharex=True, col_wrap=3, height=2,
                          aspect=1, hue="amplitude", col_order=col_list)
        sns.color_palette("BuGn_r")
        # all["OFFSET"] = all["OFFSET"].astype(str)

        # h.map(sns.lineplot, "OFFSET", "SECTION AREA [Ų]")
        # h.map(facet_scatter, "OFFSET", "diff")

        # Series for SectionArea for amplitude 0.0
        # w0 = all_mode[(all_mode.amplitude.eq(0.00))]["SECTION AREA [Ų]"]

        h.map(facet_line, "OFFSET", "SECTION AREA [Ų]", "amplitude")
        h.map(facet_auc, "OFFSET", "SECTION AREA [Ų]", "NEW")

        # h.map(sns.scatterplot, "OFFSET", "SECTION AREA [Ų]", data=all_mode)
        # h.map(sns.lineplot, "OFFSET", "diff")
        # h.map(facet_scatter, "OFFSET", "diff")

        plt.legend()

        plt.subplots_adjust(top=0.9)
        h.fig.suptitle("mode " + str(mode))
        plt.show()
