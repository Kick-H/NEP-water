from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from multiprocessing import Pool
import sys

def map_fun(frame):
    return np.mean(calc.get_property('descriptor', frame), axis=0)

if __name__=='__main__':

    nep_dir = ""
    nep_dir = "."
    proc_n=int(sys.argv[1])
    in_xyz = nep_dir+"/"+sys.argv[2]
    in_nep = nep_dir+"/nep.txt"
    min_des = float(sys.argv[3])
    minsel = int(sys.argv[4])
    runame=sys.argv[5]

    data_current=read(in_xyz,index=':',format='extxyz')
    calc = NEP(in_nep)
    # print(calc)

    # select data
    with Pool(processes=proc_n) as pool:
        des_current=np.array(pool.map(map_fun,data_current))
    sampler = FarthestPointSample(min_distance=min_des)
    selected_i = sampler.select(des_current, [], min_select=minsel)
    write(f'selected-{runame}.xyz', [data_current[i] for  i in selected_i],format='extxyz')
    np.savetxt(f"selected-{runame}",selected_i)

    # fit PCA model
    reducer = PCA(n_components=2)
    reducer.fit(des_current)
    # current data
    proj_current = reducer.transform(des_current)
    plt.scatter(proj_current[:,0], proj_current[:,1],label=f"init-{runame}",color="orange")
    # selected data
    proj_selected = reducer.transform(np.array([des_current[i] for i in selected_i]))
    plt.scatter(proj_selected[:,0], proj_selected[:,1],s=1,label="selected",color='blue')

    plt.legend()
    plt.axis('off')
    plt.savefig(f'select-{runame}.png')
