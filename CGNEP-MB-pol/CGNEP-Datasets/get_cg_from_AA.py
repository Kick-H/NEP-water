from ase.io import read, write
from tqdm import tqdm
import numpy as np
from ase import Atoms
from calorine.calculators import CPUNEP
import sys


def compute_cg_ox(pos, force, mass, cell, calc=None):
    nmol = pos.shape[0] // 3
    mol_mass = 15.999 + 1.008 * 2  # 水分子质量
    #cell = np.diag(cell)

    # 预分配数组
    pcg = np.zeros((nmol, 3))
    fcg = np.zeros((nmol, 3))

    if calc != None:
        mol_energy = 0

    for mol in range(nmol):
        atoms = slice(3*mol, 3*mol+3)
        # oxygen_pos = pos[atoms][0]  # 氧原子坐标
        # pcg[mol] = oxygen_pos - cell * np.floor(oxygen_pos / cell)
        pcg[mol] = pos[atoms][0]
        fcg[mol] = force[atoms].sum(axis=0)

        if calc != None:
            sel_list = list([3*mol, 3*mol+1, 3*mol+2])
            mol_atom = atom[sel_list]
            mol_atom.calc = calc
            energy = mol_atom.get_potential_energy()
            mol_energy += energy

    if calc != None:
        return pcg, fcg, np.full(nmol, mol_mass), mol_energy
    else:
        return pcg, fcg, np.full(nmol, mol_mass)

def compute_cg_cm(pos, force, mass, cell):
    natoms = pos.shape[0]
    nmol = natoms // 3
    mol_mass = 15.999 + 1.008 * 2
    fcg = np.zeros((nmol, 3))
    pcg = np.zeros((nmol, 3))
    pms = np.full(nmol, mol_mass)
    cell = np.diag(cell)

    for mol_idx in range(nmol):
        start = 3 * mol_idx
        mol_pos = pos[start:start+3]
        mol_force = force[start:start+3]
        mol_mass_vals = mass[start:start+3]

        ref_pos = mol_pos[0]
        adjusted_pos = [ref_pos.copy()]
        for i in range(1, 3):
            delta = mol_pos[i] - ref_pos
            adjusted_delta = delta - cell * np.round(delta / cell)
            adjusted_pos.append(ref_pos + adjusted_delta)
        adjusted_pos = np.array(adjusted_pos)

        weighted_pos = mol_mass_vals[:, np.newaxis] * adjusted_pos
        pcg[mol_idx] = weighted_pos.sum(axis=0) / mol_mass
        weighted_force = mol_force * mol_mass_vals[:, np.newaxis]
        fcg[mol_idx] = weighted_force.sum(axis=0) / mol_mass

    return pcg, fcg, pms

def get_cg_atom(atom, calc=None):
    # print(dir(atom))
    natoms = atom.get_global_number_of_atoms()
    potential = atom.get_potential_energy()
    try:
        virial = atom.info['virial']
        v_flag = True
    except:
        v_flag = False
    cell = atom.cell
    pbc = atom.get_pbc()
    typ = atom.get_chemical_symbols()

    masses = atom.get_masses()
    pos = atom.get_positions()
    force = atom.get_array('force')
    # pcg, fcg, pms = compute_cg_cm(pos, force, masses, cell)
    if calc != None:
        pcg, fcg, pms, mol_energy = compute_cg_ox(pos, force, masses, cell, calc=calc)
    else:
        pcg, fcg, pms = compute_cg_ox(pos, force, masses, cell)
    tcg = ['W'] * int(natoms/3)
    cg_atoms = Atoms(
        symbols=tcg,
        positions=pcg,
        masses=pms,
    )
    if calc != None:
        cg_atoms.info['energy'] = potential - mol_energy
    else:
        cg_atoms.info['energy'] = potential
    if v_flag:
        cg_atoms.info['virial'] = virial
    cg_atoms.cell = cell
    cg_atoms.pbc = pbc
    cg_atoms.set_array('force', fcg)

    return cg_atoms, v_flag


nep_path = "./nep.txt"

inxyz = sys.argv[1]
outxyz = sys.argv[2]

ref_xyz = tqdm(read(inxyz, ":", format="extxyz"))

cg_atoms = []
cg_novirial = []
for i, atom in enumerate(ref_xyz):
    ref_xyz.set_description(f'Processing {i}')
    # cg_atom, v_flag = get_cg_atom(atom)
    cg_atom, v_flag = get_cg_atom(atom, calc=CPUNEP(nep_path))
    if v_flag:
        cg_atoms.append(cg_atom)
    else:
        cg_novirial.append(cg_atom)

if len(cg_atoms) != 0:
    write(outxyz, cg_atoms, format="extxyz")
if len(cg_novirial) != 0:
    write(f"nov-{outxyz}", cg_novirial, format="extxyz")

