# CombineMols3D
Combine 2 molecules in 3D by sample.


A detailed description could be found in the article: Exploring exohedral functionalization of fullerene with Automation and Neural Network Potential. 

![overview](.//overview.png)

<center>Overview of the combination process.</center>

## How it works?
In short, it could be summarized into 2 steps:

1. Translation:
   * The new bond length is determined by the summation of the covalent radii of the two atoms. 
   * Set one atom as the central of a sphere, sample positions on this sphere and calculate the repulsion between the new addon to the whole molecule. (The repulsion is measured by the summation of the inverse distance.)
   * Choose the minimum repulsion position and translate the whole addon molecule as a rigid body.
2. Rotation:
   * Set the bonded site as center and the nearest neighbor as rotation axes. Rotate the addon molecule as a rigid body
   * Sample rotation angles to minimize the repulsion between the addon molecule to the main molecule. (Again, the repulsion is measured by the summation of the inverse distance.)
   * Choose the minimum repulsion image and this is the final combination state.


## Install
```
pip install combinemols3d
```

## Usage
1. Prepare your to-be-bonded molecules as `ase.Atoms` object. (Read  `ase` documentation for format manipulation.)
2. Determine which sites to be bonded together. (dummy sites)
3. Combine them by eliminate the dummy sites.

```python
    from ase.visualize import view
    from ase.build.molecule import molecule

    main_mol = molecule('C7NH5')
    sub_mol = molecule('BDA')

    view(main_mol)
    view(sub_mol)

    final_mol = combine_2_mols_with_dummy(molecule_1=main_mol, molecule_2=sub_mol,
                                          dummy_atom_index_1=11, dummy_atom_index_2=6)
    view(final_mol)

```
If you don't want to eliminate the dummy sites, use the `combine_2_mols` instead:
```python
    final_mol = combine_2_mols(molecule_1=main_mol, molecule_2=sub_mol,
                                          tgt_atom_1_index=11, tgt_atom_2_index=6)
```
Other parameters are same for the two functions:

* Conformation related:

  * `skin`: to get a proper bond length for the new bond. Here we adopt the [covalent radii method](https://en.wikipedia.org/wiki/Covalent_radius). A skin parameter is added to finetune the bond length, though it is 0 by default.
    $$
    R(\mathrm{AB})=r(\mathrm{A})+r(\mathrm{B})+skin
    $$

  * `cutoff_mult`: there are several places to detect neighbors for a specific atom. By default, it is `natural cutoff` as follow equation 2. 
    $$
    DetectRange(\mathrm{AB})=r(\mathrm{A})+r(\mathrm{B})+0.3
    $$
    We can enlarge our search scope to multiply the  `natural cutoff` as follow equation 3.
    
    $$
    DetectRange(\mathrm{AB})=cutoff\_mult*(r(\mathrm{A})+r(\mathrm{B})+0.3)
    $$

* Performance related:
  * `sample_times`: how many times to sample addon positions in translation and angles in rotation.
  * `rotation_times`: In the rotation stage, we will rotate the whole addon molecule as a rigid body. The orientations of the addon site to its neighbors are taken as fixed rotation axes. How many neighbors to be considered could be performance critical.


## Note

Contact me: 1660810667@qq.com
