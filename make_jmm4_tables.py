#!/usr/bin/env python

import csv
import numpy as np

SLO_FUNS = ['1', 'p', 'v', 'm', 'g']

NAMES = ['Constant', 'Linear \#1', 'Linear \#2', 'Sine', 'Sloth']

if __name__ == '__main__':

    print(r'\begin{tabular}{ccccccc}')

    print(r'  & $\tau - T$ & $\tau_x - T_x$ & $\tau_y - T_y$ & $\tau_{xx} - T_{xx}$ & $\tau_{xy} - T_{xy}$ & $\tau_{yy} - T_{yy}$ \\')

    for name, slo_fun in zip(NAMES, SLO_FUNS):
        print('  ' + name, end='')

        N, E, Ex, Ey, Exx, Exy, Eyy = [], [], [], [], [], [], []
        with open('data/%s/errors.csv' % slo_fun) as f:
            for i, (n, e, ex, ey, exx, exy, eyy) in enumerate(csv.reader(f)):
                if 3 <= i and i <= 7:
                    N.append(float(n))
                    E.append(float(e))
                    Ex.append(float(ex))
                    Ey.append(float(ey))
                    Exx.append(float(exx))
                    Exy.append(float(exy))
                    Eyy.append(float(eyy))
        N = np.array(N)
        E = np.array(E)
        Ex = np.array(Ex)
        Ey = np.array(Ey)
        Exx = np.array(Exx)
        Exy = np.array(Exy)
        Eyy = np.array(Eyy)

        pE = -np.mean((np.log(E[1:]) - np.log(E[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))
        pEx = -np.mean((np.log(Ex[1:]) - np.log(Ex[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))
        pEy = -np.mean((np.log(Ey[1:]) - np.log(Ey[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))
        pExx = -np.mean((np.log(Exx[1:]) - np.log(Exx[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))
        pExy = -np.mean((np.log(Exy[1:]) - np.log(Exy[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))
        pEyy = -np.mean((np.log(Eyy[1:]) - np.log(Eyy[:-1]))/(np.log(N[1:]) - np.log(N[:-1])))

        print(r' & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f & %1.2f \\' % (pE, pEx, pEy, pExx, pExy, pEyy))

    print(r'\end{tabular}')
