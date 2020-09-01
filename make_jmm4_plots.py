#!/usr/bin/env python

import csv
import glob
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import os

slo_funs = ['1', 'p', 'v', 'm', 'g']

title = {
    '1': 'Constant',
    'p': 'Linear #1',
    'v': 'Linear #2',
    'm': 'Sine',
    'g': 'Sloth'
}

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


def get_R(n, slo_fun):
    if slo_fun == 'v':
        R = n//10
    elif slo_fun == 'g':
        R = n//5
    else:
        R = n//20
    return max(5, R)


def get_i0(n, slo_fun):
    if slo_fun in {'v', 'g'}:
        return 0
    else:
        return n//2


def rms(E, n, i0, R):
    I, J = np.meshgrid(np.arange(n), np.arange(n))
    mask = np.maximum(abs(I - i0), abs(J - i0)) > R
    mask = mask & ~np.isnan(E)
    E_mask = E[mask]
    return np.linalg.norm(E_mask)/np.sqrt(E_mask.size)

def maxerr(E, n, i0, R):
    I, J = np.meshgrid(np.arange(n), np.arange(n))
    mask = np.maximum(abs(I - i0), abs(J - i0)) > R
    mask = mask & ~np.isnan(E)
    E_mask = E[mask]
    return abs(E_mask).max()

if __name__ == '__main__':
    fig, axes = plt.subplots(1, len(slo_funs), figsize=(14, 4))

    vmin, vmax = np.inf, -np.inf

    for i, slo_fun in enumerate(slo_funs):
        print(slo_fun)

        rms_errors = dict()
        max_errors = dict()

        ax = axes[i]
        for error in errors:
            N = []
            data = dict()
            for E_path in glob.glob(os.path.join('.', 'data', slo_fun, '%s_*.bin' % error)):
                n = int(os.path.splitext(os.path.basename(E_path))[0].split('_')[1])
                N.append(n)
                mapped = np.memmap(E_path, shape=(n, n), dtype=np.float64)
                data[n] = mapped.copy()
                del mapped
            N = np.array(sorted(N))
            Erms = np.array([
                rms(data[n], n, get_i0(n, slo_fun), get_R(n, slo_fun))
                for n in N
            ])
            rms_errors[error] = Erms
            Emax = np.array([
                maxerr(data[n], n, get_i0(n, slo_fun), get_R(n, slo_fun))
                for n in N
            ])
            max_errors[error] = Emax
            ax.loglog(
                np.array(N)**3,
                Erms,
                label=error_labels[error],
                marker='*'
            )
            vmin = min(vmin, Erms.min())
            vmax = max(vmax, Erms.max())
        ax.set_title(title[slo_fun])
        if i == 0:
            ax.set_ylabel('RMS Error')
            ax.legend(fontsize=8)
        ax.set_xlabel('$|\Omega_h|$')
        if i > 0:
            ax.set_yticklabels([])

        with open('data/%s/errors.csv' % slo_fun, 'w') as f:
            writer = csv.writer(f)
            for i, n in enumerate(N):
                writer.writerow(
                    [n] + [rms_errors[error][i] for error in errors] +
                    [max_errors[error][i] for error in errors]
                )

    for ax in axes:
        ax.set_yticks(1/np.power(10, np.arange(1, 12)))
        ax.set_ylim(1e-11, 1e-1)

    fig.tight_layout()
    fig.savefig('plots/size_vs_errors_jmm4.pdf')
