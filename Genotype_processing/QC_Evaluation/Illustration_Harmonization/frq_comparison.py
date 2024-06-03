#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import seaborn as sns


parser = argparse.ArgumentParser(description='pvp plot')

parser.add_argument('-i', dest='freq', type=str,nargs='+', help='output from compare_ancestry_adjusted_af.py')
parser.add_argument('-n', dest='n', type=int, help='number total of test')
parser.add_argument('-t', dest='title', type=str, help='Title (use quotation for multiple words)')
parser.add_argument('-A', dest='Array' , type=str,nargs='+', help='list arrays (space delimited')
parser.add_argument('-out', dest='out', type=str, help='Output names')

args = parser.parse_args()

df=pd.concat([pd.read_csv(i, sep = "\s+", header=0,index_col='VARIANT') for i in args.freq])
df=df[df["LRT_PVALUE"] != "None"]
arrays=df.filter(like='AF')
arrays.columns=[i.replace("AF_","") for i in arrays.columns]
arrays=arrays[args.Array]
arrays_2=arrays[pd.to_numeric(df["FTEST_PVALUE"])<=(0.05/args.n)]
arrays_1=arrays[pd.to_numeric(df["FTEST_PVALUE"])>(0.05/args.n)]



fig = plt.figure(figsize=(12, 12), dpi = 100)
n = 0
n_row = 0
n_col = 0
rows_ax = dict()
cols_ax = dict()
for label_i in arrays.columns:
    for label_j in arrays.columns:
        if n % len(arrays.columns) == 0:
            n_row += 1
            n_col = 0
        n_col += 1
        n += 1
        print(n, n_row, n_col, label_i, label_j)
        column_i = label_i
        column_j = label_j
        sharey_ax = rows_ax.get(n_row, None)
        sharex_ax = cols_ax.get(n_col, None)
        if sharey_ax is None and sharex_ax is None:
            ax = fig.add_subplot(len(arrays.columns), len(arrays.columns), n)
            rows_ax[n_row] = ax
            cols_ax[n_col] = ax
        elif sharey_ax is not None and sharex_ax is not None:
            ax = fig.add_subplot(len(arrays.columns), len(arrays.columns), n, sharey = sharey_ax, sharex = sharex_ax)
        elif sharey_ax is not None:
            ax = fig.add_subplot(len(arrays.columns), len(arrays.columns), n, sharey = sharey_ax)
            cols_ax[n_col] = ax
        else:
            ax = fig.add_subplot(len(arrays.columns), len(arrays.columns), n, sharex = sharex_ax)
            ax.xaxis.set_visible(False)
            rows_ax[n_row] = ax
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        # ax.scatter(-np.log10(df_results_all_good[column_j]), -np.log10(df_results_all_good[column_i]), s = 1)
        # ax.scatter(-np.log10(df_results_all_bad[column_j]), -np.log10(df_results_all_bad[column_i]), s = 1)
        ax.scatter(arrays_1[column_j], arrays_1[column_i], alpha=0.5,s=20, facecolor='none', edgecolors = 'gray')
        ax.scatter(arrays_2[column_j], arrays_2[column_i], alpha=1,s=60, facecolor='none', edgecolors = 'red')
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(column_j, fontsize=15)
        ax.set_ylabel(column_i, fontsize=15)

for n_row, ax in rows_ax.items():
    ax.yaxis.set_visible(True)
    # ax.set_ylim([0, 5])


for n_col, ax in cols_ax.items():
    ax.xaxis.tick_top()
    ax.xaxis.set_visible(True)
    # ax.set_xlim([0, 5])
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.05,
                    hspace=0.05)

fig.suptitle(args.title, fontsize=16)

plt.savefig(args.out+".png")


