#!/usr/bin/env python

import glob
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import os

slo_funs = ['1', 'p', 'v', 'm', 'g']

errors = ['E', 'Ex', 'Ey', 'Exx', 'Exy', 'Eyy']
error_labels = {
    'E': r'$\tau - T$',
    'Ex': r'$\tau_x - T_x$',
    'Ey': r'$\tau_y - T_y$',
    'Exx': r'$\tau_{xx} - T_{xx}$',
    'Exy': r'$\tau_{xy} - T_{xy}$',
    'Eyy': r'$\tau_{yy} - T_{yy}$',
}

results = ['T', 'Tx', 'Ty', 'Txx', 'Txy', 'Tyy']

def get_memory_mapped_results(data_dir):
    mapped = {slo_fun: dict() for slo_fun in slo_funs}
    for slo_fun, error in it.product(slo_funs, errors):
        mapped[slo_fun][error] = dict()
        glob_str = os.path.join(data_dir, slo_fun, '%s_*.bin' % error)
        for E_path in glob.glob(glob_str):
            N = int(os.path.splitext(os.path.basename(E_path))[0].split('_')[1])
            mapped[slo_fun][error][N] = \
                np.memmap(E_path, shape=(N, N), dtype=np.float64)
    return mapped

def get_memory_mapped_errors(data_dir):
    mapped = {slo_fun: dict() for slo_fun in slo_funs}
    for slo_fun, error in it.product(slo_funs, errors):
        mapped[slo_fun][error] = dict()
        glob_str = os.path.join(data_dir, slo_fun, '%s_*.bin' % error)
        for E_path in glob.glob(glob_str):
            N = int(os.path.splitext(os.path.basename(E_path))[0].split('_')[1])
            mapped[slo_fun][error][N] = \
                np.memmap(E_path, shape=(N, N), dtype=np.float64)
    return mapped

def rms(E):
    mask = ~np.isnan(E)
    return np.linalg.norm(E[mask])/np.sqrt(E.size)

def maxerr(E, T):
    mask = ~np.isnan(E)
    return np.linalg.norm(E[mask], ord=np.inf)

if __name__ == '__main__':
    mapped = get_memory_mapped_errors('./data')

    fig = plt.figure(figsize=(14, 3))
    for i, slo_fun in enumerate(slo_funs):
        ax = fig.add_subplot(1, len(slo_funs), i + 1)
        for error in errors:
            N = sorted(list(mapped[slo_fun][error].keys()))
            plt.loglog(
                np.array(N)**3,
                [rms(mapped[slo_fun][error][n]) for n in N],
                label=error_labels[error],
                marker='*'
            )
        ax.set_title('slo_fun = %s' % slo_fun)
        ax.set_ylabel('RMS Error')
        ax.set_xlabel('$|\Omega_h|$')
        ax.legend()
    fig.legend()
    fig.savefig('tmp.pdf')

    # del mapped
