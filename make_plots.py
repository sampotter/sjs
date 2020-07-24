#!/usr/bin/env python


import csv
import glob
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


OUTPUT_PATH = os.path.join('.', 'data')
PLOTS_PATH = os.path.join('.', 'plots')
COLUMNS = ['Nx', 'Emax', 'Erms', 'Emax/umax', 'Emax/urms', 'gEmax', 'gErms', 'CPU']
SLO_FUN = ['1', 'p', 'v', 'm', 'g']
METHODS = ['JMM1', 'JMM2', 'JMM3', 'FMM', 'OLIM']
SLO_FUN_NAMES = {
    '1': 'Constant',
    'p': r'Linear #1',
    'v': r'Linear #2',
    'm': 'Sine',
    'g': 'Sloth'
}


def get_df_key(path):
    basename = os.path.basename(path)
    filename = os.path.splitext(basename)[0]
    method, slo_fun = filename.split('_')
    method = method[:4].upper()
    slo_fun = slo_fun[3]
    return slo_fun, method

def read_df_from_file(path):
    with open(path, 'r') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    return pd.DataFrame(
        np.array([[float(_) for _ in row[:8]] for row in rows]),
        index=np.array([int(row[0])**3 for row in rows]),
        columns=COLUMNS
    )

if __name__ == '__main__':
    data_frames = dict()
    for path in glob.glob(os.path.join(OUTPUT_PATH, '*.csv')):
        key = get_df_key(path)
        df = read_df_from_file(path)
        data_frames[key] = df

    for slo_fun in SLO_FUN:
        print(slo_fun)

        # Make size vs error plots

        fig, axes = plt.subplots(1, 4, figsize=(12, 3), sharex=True, sharey=True)
        axes = axes.flatten()

        for jmm in METHODS:
            df = data_frames[slo_fun, jmm]
            N = df.index
            axes[0].loglog(N, df['Emax'], label=jmm, marker='*')
            axes[1].loglog(N, df['Erms'], label=jmm, marker='*')
            axes[2].loglog(N, df['gEmax'], label=jmm, marker='*')
            axes[3].loglog(N, df['gErms'], label=jmm, marker='*')

        axes[0].set_ylabel(r'$\ell_\infty$ error ($T$)')
        axes[1].set_ylabel(r'RMS error ($T$)')
        axes[2].set_ylabel(r'$\ell_\infty$ error ($\nabla T$)')
        axes[3].set_ylabel(r'RMS error ($\nabla T$)')

        for ax in axes:
            ax.set_xlabel(r'$|\Omega_h|$')
            ax.legend(loc='lower left', fontsize=7)

        fig.suptitle('Slowness: %s' % SLO_FUN_NAMES[slo_fun], fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS_PATH, 'size_vs_errors_%s.pdf' % slo_fun))
        plt.close(fig)

        # Make time vs error plots

        fig, axes = plt.subplots(1, 4, figsize=(12, 3), sharex=True, sharey=True)
        axes = axes.flatten()

        for jmm in METHODS:
            df = data_frames[slo_fun, jmm]
            T = df['CPU']
            axes[0].loglog(T, df['Emax'], label=jmm, marker='*')
            axes[1].loglog(T, df['Erms'], label=jmm, marker='*')
            axes[2].loglog(T, df['gEmax'], label=jmm, marker='*')
            axes[3].loglog(T, df['gErms'], label=jmm, marker='*')

        axes[0].set_ylabel(r'$\ell_\infty$ error ($T$)')
        axes[1].set_ylabel(r'RMS error ($T$)')
        axes[2].set_ylabel(r'$\ell_\infty$ error ($\nabla T$)')
        axes[3].set_ylabel(r'RMS error ($\nabla T$)')

        for ax in axes:
            ax.legend(loc='lower left', fontsize=7)
            ax.set_xlabel(r'CPU Time (s)')

        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS_PATH, 'time_vs_errors_%s.pdf' % slo_fun))
        plt.close(fig)
