import ase
import numpy as np
from ase import Atoms, Atom
from ase.data import covalent_radii, atomic_numbers
from ase.io import read, write
from ase.neighborlist import build_neighbor_list, natural_cutoffs
from scipy.spatial import distance_matrix


def get_bond_length(atom_1_symbol: str, atom_2_symbol: str, skin: float = 0):
    return covalent_radii[atomic_numbers[atom_1_symbol]] + covalent_radii[atomic_numbers[atom_2_symbol]] + skin


def get_nearest_neighbor(a_molecule: ase.Atoms, atom_index: int, cutoff_mult: float = None, sort_flag: bool = True):
    if cutoff_mult:
        cutoffs = natural_cutoffs(atoms=a_molecule, mult=cutoff_mult)
        neighborList = build_neighbor_list(a_molecule, cutoffs=cutoffs, bothways=True, self_interaction=False)
    else:
        neighborList = build_neighbor_list(a_molecule, bothways=True, self_interaction=False)
    dok_matrix = neighborList.get_connectivity_matrix()
    dok_adj = dok_matrix[atom_index]
    adj_array = dok_adj.tocoo().col
    if sort_flag:
        dist = []
        for ii in adj_array:
            dist.append(np.linalg.norm(a_molecule.positions[atom_index] - a_molecule.positions[ii]))
        map_dict = dict(zip(dist, adj_array))
        sorted_dist = sorted(dist)
        sorted_neighbor_indexes = []
        for a_dist in sorted_dist:
            sorted_neighbor_indexes.append(map_dict[a_dist])
        return sorted_neighbor_indexes
    else:
        return adj_array


def swap_atoms(a_molecule: ase.Atoms, atom_index_1: int, atom_index_2: int):
    a_molecule_copy = a_molecule.copy()
    a_molecule[atom_index_1].position = a_molecule_copy[atom_index_2].position
    a_molecule[atom_index_1].symbol = a_molecule_copy[atom_index_2].symbol
    a_molecule[atom_index_2].position = a_molecule_copy[atom_index_1].position
    a_molecule[atom_index_2].symbol = a_molecule_copy[atom_index_1].symbol
    return a_molecule


def rm_and_sort_atoms(a_molecule: ase.Atoms, atom_index_to_rm: int):
    nearest_neighbor_index = get_nearest_neighbor(a_molecule=a_molecule, atom_index=atom_index_to_rm)[0]
    del a_molecule[atom_index_to_rm]
    if nearest_neighbor_index > atom_index_to_rm:
        nearest_neighbor_index = nearest_neighbor_index - 1
    a_molecule = swap_atoms(a_molecule=a_molecule, atom_index_1=0, atom_index_2=nearest_neighbor_index)
    return a_molecule


def get_addon_pos_by_sample(addon_symbol: str, tgt_atom_index: int, tgt_molecule: Atoms, sample_times: int = 1000,
                            skin: float = 0, cutoff_mult: float = 1.5):
    tgt_atom = tgt_molecule[tgt_atom_index]
    tgt_bond_length = get_bond_length(atom_1_symbol=addon_symbol, atom_2_symbol=tgt_atom.symbol, skin=skin)
    neighbor_indexes = get_nearest_neighbor(a_molecule=tgt_molecule, atom_index=tgt_atom_index, cutoff_mult=cutoff_mult)
    neighbor_positions = tgt_molecule.get_positions()[neighbor_indexes]

    min_rev_distance = 100000000
    for i in range(sample_times):
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        new_atom_position = tgt_atom.position + tgt_bond_length * np.array([x, y, z])
        distances = distance_matrix(new_atom_position.reshape((1, -1)), neighbor_positions).reshape(-1)
        a_distance = (1 / distances).sum()

        if a_distance < min_rev_distance:
            min_rev_d_point = new_atom_position
            min_rev_distance = a_distance

    return min_rev_d_point, neighbor_indexes


def combine_2_mols(molecule_1: Atoms, molecule_2: Atoms, tgt_atom_1_index: int, tgt_atom_2_index: int,
                   sample_times: int = 1000, skin: float = 0, cutoff_mult: float = 1.5,
                   rotation_times: int = 3):
    # translation
    if len(molecule_1) > len(molecule_2):
        main_mol, sub_mol = molecule_1, molecule_2
        tgt_atom_main_index, tgt_atom_addon_index = tgt_atom_1_index, tgt_atom_2_index
    else:
        main_mol, sub_mol = molecule_2, molecule_1
        tgt_atom_main_index, tgt_atom_addon_index = tgt_atom_2_index, tgt_atom_1_index
    addon_atom = sub_mol[tgt_atom_addon_index]
    addon_position, neighbor_indexes = get_addon_pos_by_sample(addon_symbol=addon_atom.symbol,
                                                               tgt_atom_index=tgt_atom_main_index,
                                                               tgt_molecule=main_mol,
                                                               sample_times=sample_times, skin=skin,
                                                               cutoff_mult=cutoff_mult)
    translation_vec = addon_position - addon_atom.position
    sub_mol.translate(translation_vec)
    # rotate along with the nearest neighbors to minimize the fragments repulsion
    main_mol_positions = main_mol.get_positions()
    min_rev_distance = 100000000
    sub_mol_nb_indexes = get_nearest_neighbor(a_molecule=sub_mol, atom_index=tgt_atom_addon_index, cutoff_mult=cutoff_mult)
    for i in range(sample_times):
        sub_mol_copy = sub_mol.copy()
        for a_neighbor_index in sub_mol_nb_indexes[:rotation_times]:
            rotate_angle = np.random.uniform(0, 360)
            rotate_vec = addon_position - sub_mol_copy[a_neighbor_index].position
            sub_mol_copy.rotate(a=rotate_angle, v=rotate_vec, center=addon_position)
            sub_mol_copy_positions = sub_mol_copy.get_positions()
            frag_distances = distance_matrix(sub_mol_copy_positions, main_mol_positions)
            a_frag_distance = (1 / frag_distances).sum()
            if a_frag_distance < min_rev_distance:
                min_rev_distance = a_frag_distance
                final_sub_mol = sub_mol_copy.copy()

    # merge
    for an_atom in final_sub_mol:
        main_mol.append(an_atom)

    return main_mol


def combine_2_mols_with_dummy(molecule_1: Atoms, molecule_2: Atoms, dummy_atom_index_1: int, dummy_atom_index_2: int,
                              sample_times: int = 1000, skin: float = 0, cutoff_mult: float = 1.5,
                              rotation_times: int = 3):
    # remove two dummy atoms and sort the to-be-bonded atom to index 0
    mol_1_without_dummy = rm_and_sort_atoms(a_molecule=molecule_1, atom_index_to_rm=dummy_atom_index_1)
    mol_2_without_dummy = rm_and_sort_atoms(a_molecule=molecule_2, atom_index_to_rm=dummy_atom_index_2)
    combinend_mol = combine_2_mols(molecule_1=mol_1_without_dummy, molecule_2=mol_2_without_dummy,
                                   tgt_atom_1_index=0, tgt_atom_2_index=0, rotation_times=rotation_times,
                                   sample_times=sample_times, skin=skin, cutoff_mult=cutoff_mult)
    return combinend_mol


if __name__ == '__main__':
    from ase.visualize import view
    from ase.build.molecule import molecule

    main_mol = molecule('C7NH5')
    sub_mol = molecule('BDA')

    view(main_mol)
    view(sub_mol)

    final_mol = combine_2_mols_with_dummy(molecule_1=main_mol, molecule_2=sub_mol,
                                          dummy_atom_index_1=11, dummy_atom_index_2=6)


    view(final_mol)
