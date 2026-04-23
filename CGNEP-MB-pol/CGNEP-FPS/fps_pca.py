from calorine.calculators import CPUNEP
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from multiprocessing import Pool
import sys

def map_fun(frame):
    return np.mean(calc.get_descriptors(frame), axis=0)

def fps_select(x, n):
    selected = [0]
    dmin = np.sum((x - x[0])**2, axis=1)
    dmin[0] = -1
    for _ in range(n - 1):
        i = np.argmax(dmin)
        selected.append(i)
        d = np.sum((x - x[i])**2, axis=1)
        dmin = np.minimum(dmin, d)
        dmin[selected] = -1
    return selected

if __name__ == '__main__':

    proc_n = 4
    in_xyz = sys.argv[1].split(',')
    in_nep = sys.argv[2]
    n_sel = int(sys.argv[3])

    data_current = []
    frame = [0]
    for xyz in in_xyz:
        atom = read(xyz, index=':', format='extxyz')
        data_current.extend(atom)
        frame.append(len(atom))

    calc = CPUNEP(in_nep)

    with Pool(processes=proc_n) as pool:
        des_current = np.array(pool.map(map_fun, data_current))

    selected_i = fps_select(des_current, n_sel)
    write('selected.xyz', [data_current[i] for i in selected_i], format='extxyz')
    np.savetxt('selected_i', selected_i, fmt='%d')

    reducer = PCA(n_components=2)
    proj_current = reducer.fit_transform(des_current)
    proj_selected = proj_current[selected_i]

    plt.figure(figsize=(8, 7))
    for fi in range(1, len(frame)):
        sta = int(np.sum(frame[:fi]))
        end = int(np.sum(frame[:fi+1]))
        plt.scatter(proj_current[sta:end, 0], proj_current[sta:end, 1],
                    s=40, color=f"C{fi}", label=in_xyz[fi-1])

    plt.scatter(proj_selected[:, 0], proj_selected[:, 1],
                s=10, color='k', label='selected')

    plt.legend()
    plt.axis('off')
    plt.savefig('pca.png', dpi=300, bbox_inches='tight')
