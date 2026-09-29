from scm import plams
import numpy as np
from typing import List, Dict


def get_fragments_data(reader: plams.KFReader) -> Dict:
    """
    Read data about the fragments for the system.

    Args:
        reader: the plams.KFReader object that belongs to the system under study.

    Returns:
        :Dictionary containing information about the fragments:

        - **used_atomic_fragments (bool)** – whether atomic or molecular fragments were used.
        - **fragment_names (List[str])** - a list of the fragment names.
        - **fragment_molecules (Dict[str, plams.Molecule])** - a dictionary containing, for each fragment, 
          the molecule containing the atoms that belong to the fragment.
        - **fmo_fragtype_to_fragname_map (Dict[int, str])** - a dictionary containing, for each FMO fragment type,
          the fragment name that belongs to it.
        - **complex_molecule (plams.Molecule)** - the molecule of the complex.
    """
    def get_atoms() -> List[plams.Atom]:
        '''
        Read in the atoms from a reader in internal order.
        '''
        rkf_coordinates      = np.atleast_1d(reader.read('Geometry', 'xyz')).reshape(-1, 3) * 0.529177
        rkf_atomtypes        = np.atleast_1d(reader.read('Geometry', 'atomtype').split()).astype(str)
        # we need the atomtype indices to figure out which atom type belongs to which atom index
        rkf_atomtype_indices = np.atleast_1d(reader.read('Geometry', 'fragment and atomtype index')).astype(int)
        rkf_atomtype_indices = rkf_atomtype_indices[len(rkf_atomtype_indices)//2:]
        # get the atom types
        atom_types = rkf_atomtypes[rkf_atomtype_indices-1]
        # and build the atoms
        atoms = [plams.Atom(symbol=atom_types[i], coords=rkf_coordinates[i]) for i in range(len(atom_types))]
        return atoms

    ret = {}

    # check if the calculations used atomic fragments
    # if it did it should have the same number of atoms as fragments
    ret['used_atomic_fragments'] = reader.read('Geometry', 'nr of fragments') == reader.read('Geometry', 'nr of atoms')

    # we need to treat systems that use atomic fragments differently from systems that use molecular fragments
    if ret['used_atomic_fragments']:
        rkf_atomtypes        = np.atleast_1d(reader.read('Geometry', 'atomtype').split()).astype(str)
        rkf_atomtype_indices = np.atleast_1d(reader.read('Geometry', 'fragment and atomtype index')).astype(int)
        rkf_atomtype_indices = rkf_atomtype_indices[len(rkf_atomtype_indices)//2:]
        rkf_atomtype_order   = np.atleast_1d(reader.read('Geometry', 'atom order index')).astype(int)
        rkf_atomtype_order   = rkf_atomtype_order[len(rkf_atomtype_order)//2:]
        rkf_napp             = np.atleast_1d(reader.read('Symmetry', 'napp')).astype(int)
        rkf_notyps           = np.atleast_1d(reader.read('Symmetry', 'notyps')).astype(int)
        rkf_fmo_fragments    = np.atleast_1d(reader.read('SFOs', 'fragment')).astype(int)

        notyp_atom_types = rkf_atomtypes[rkf_notyps-1]

        # go through each symmetry type and build a new fragment
        fragment_names = []
        for i, atom_type in enumerate(notyp_atom_types):
            # get the number of fragments that have the same atom type
            is_unique = len([atom_type_ for atom_type_ in notyp_atom_types if atom_type == atom_type_]) == 1
            if is_unique:
                fragment_names.append(str(atom_type))
                continue

            # get the indices of the atoms that belong to this symm. type
            atom_internal_order_indices = np.where(rkf_napp == (i + 1))[0]

            # convert to input order
            atom_input_order_indices = rkf_atomtype_order[atom_internal_order_indices]

            # prepare atom indices to be added to the fragment name
            atom_index_label = ",".join(atom_input_order_indices.astype(str))

            # generate the final fragment name
            fragment_names.append(f'{atom_type}:{atom_index_label}')

        atoms = get_atoms()
        fragment_molecules = {}
        for i, napp in enumerate(rkf_napp):
            fragment_name = fragment_names[napp-1]
            fragment_molecules.setdefault(fragment_name, plams.Molecule())
            fragment_molecules[fragment_name].add_atom(atoms[i])
        
        fmo_unique_fragments = np.unique(rkf_fmo_fragments)
        fmo_fragtype_to_name_map = {fragment_number: fragment_name for fragment_number, fragment_name in zip(fmo_unique_fragments, fragment_names)}

        ret['fragment_names'] = fragment_names
        ret['fragment_molecules'] = fragment_molecules
        ret['fmo_fragtype_to_fragname_map'] = fmo_fragtype_to_name_map
    else:
        rkf_fragtype_indices = np.atleast_1d(reader.read('Geometry', 'fragment and atomtype index')).astype(int)
        rkf_fragtype_indices = rkf_fragtype_indices[:len(rkf_fragtype_indices)//2]
        rkf_fragmenttypes    = np.atleast_1d(reader.read('Geometry', 'fragmenttype').split()).astype(str)
        rkf_fmo_fragments    = np.atleast_1d(reader.read('SFOs', 'fragment')).astype(int)

        atoms = get_atoms()
        fragment_molecules = {}
        for i, fragtype_index in enumerate(rkf_fragtype_indices):
            fragment_name = rkf_fragmenttypes[fragtype_index-1]
            fragment_molecules.setdefault(str(fragment_name), plams.Molecule())
            fragment_molecules[str(fragment_name)].add_atom(atoms[i])

        fmo_unique_fragments = np.unique(rkf_fmo_fragments)
        fmo_fragtype_to_name_map = {fragment_number: str(fragment_name) for fragment_number, fragment_name in zip(fmo_unique_fragments, rkf_fragmenttypes)}

        ret['fragment_names'] = rkf_fragmenttypes
        ret['fragment_molecules'] = fragment_molecules
        ret['fmo_fragtype_to_fragname_map'] = fmo_fragtype_to_name_map

    ret['complex_molecule'] = plams.Molecule()
    [ret['complex_molecule'].add_atom(atom) for atom in get_atoms()]

    return ret
