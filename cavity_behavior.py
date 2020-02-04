import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance
import numpy as np
import Bio.PDB
import Bio.AlignIO as al
from sklearn import preprocessing
import csv


def facet_scatter(x, y, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, **kwargs)


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
        df_filtered['struct'] = '_'.join(file_csv.split('_')[:3])
        df_filtered.to_csv('filtered_{}.csv'.format('_'.join(file_csv.split('_')[:3])))
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

    col_list = list(set(all['struct'].tolist()))
    col_list.sort()

    h = sns.FacetGrid(all, col="struct", sharey='row', sharex=True, col_wrap=6, height=3,
                      aspect=1, hue="OFFSET", col_order=col_list)
    sns.color_palette("BuGn_r")
    # all["OFFSET"] = all["OFFSET"].astype(str)

    h.map(sns.lineplot, "TIME STEP", "SECTION AREA [Ų]")
    h.map(sns.scatterplot, "TIME STEP", "SECTION AREA [Ų]", data=all)

    plt.legend()

    plt.subplots_adjust(top=0.9)
    h.fig.suptitle("all: Area=f(time step)")
    plt.show()