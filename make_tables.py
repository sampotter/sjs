import glob
import os


OUTPUT_PATH = os.path.join('.', 'data')
SLO_FUN = ['1', 'p', 'v', 'm', 'g']
JMM = ['JMM1', 'JMM2', 'JMM3']
KEY = ['pmax_T', 'prms_T', 'pmax_dT', 'prms_dT']
JMM_SHORT_NAME = {
    'JMM1': r'\#1',
    'JMM2': r'\#2',
    'JMM3': r'\#3'
}
SLO_FUN_NAMES = {
    '1': 'Constant',
    'p': r'Linear \#1',
    'v': r'Linear \#2',
    'm': 'Sine',
    'g': 'Sloth'
}

def get_ps_from_lines(lines):
    _ = lines[0].split(':')[1].split(',')[0]
    pmax_T = float(_.split('=')[1])
    _ = lines[1].split(':')[1].split(',')[0]
    prms_T = float(_.split('=')[1])
    _ = lines[2].split(':')[1].split(',')[0]
    pmax_dT = float(_.split('=')[1])
    _ = lines[3].split(':')[1].split(',')[0]
    prms_dT = float(_.split('=')[1])
    return pmax_T, prms_T, pmax_dT, prms_dT

def get_slo_fun_and_method(path):
    basename = os.path.basename(path)
    filename = os.path.splitext(basename)[0]
    method, slo_fun = filename.split('_')[:2]
    method = method[:4].upper()
    slo_fun = slo_fun[3]
    return slo_fun, method

if __name__ == '__main__':
    stats = dict()
    for path in glob.glob(os.path.join(OUTPUT_PATH, '*.txt')):
        with open(path, 'r') as f:
            slo_fun, jmm = get_slo_fun_and_method(path)
            lines = f.readlines()
            pmax_T, prms_T, pmax_dT, prms_dT = get_ps_from_lines(lines)
            stats[slo_fun, jmm] = {
                'pmax_T': pmax_T,
                'prms_T': prms_T,
                'pmax_dT': pmax_dT,
                'prms_dT': prms_dT
            }

    print(r'\begin{tabular}{cccccc}')
    print(r'  & JMM & $E_{\mbox{max}} (T)$ & $E_{\mbox{RMS}} (T)$', end='')
    print(r'& $E_{\mbox{max}} (\nabla T)$ & $E_{\mbox{RMS}} (\nabla T)$ \\')
    for slo_fun in SLO_FUN:
        print(r'  \midrule')
        print(r'  \multirow{3}{*}{%s}' % SLO_FUN_NAMES[slo_fun])
        for i, jmm in enumerate(JMM):
            print('   ', end='')
            print(r' & %s' % JMM_SHORT_NAME[jmm], end='')
            for key in KEY:
                print(r' & %1.2f' % stats[slo_fun, jmm][key], end='')
            print(r' \\')
    print(r'\end{tabular}')
