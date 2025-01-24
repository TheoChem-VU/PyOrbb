import numpy as np
from tcutility import results, timer, ensure_list
from math import sqrt
from scm import plams


def _get_molecules(reader):
    used_regions = reader.read('Geometry', 'nr of fragments') != reader.read('Geometry', 'nr of atoms')
    fragment_indices = np.atleast_1d(reader.read('Geometry', 'fragment and atomtype index'))
    fragment_indices = fragment_indices[len(fragment_indices)//2:]
    fragment_types = np.atleast_1d(reader.read('Geometry', 'fragmenttype').split())
    if not used_regions:
        fragment_uniques = [f'{fragment_types[frag_idx-1]}:{idx+1}' for idx, frag_idx in enumerate(fragment_indices)]
        fragment_uniques = np.array(sorted(set(fragment_uniques), key=lambda fu: int(fu.split(':')[1])))
    else:
        fragment_uniques = np.array(fragment_types)

    coords = np.array(reader.read('Geometry', 'xyz')).reshape(-1, 3) * 0.529177249
    atoms = np.array(reader.read('Geometry', 'atomtype').split())

    order_index = np.array(ensure_list(reader.read('Geometry', 'atom order index'))[:coords.shape[0]]) - 1
    fragment_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[:coords.shape[0]]) - 1
    symbol_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[coords.shape[0]:]) - 1

    coords = coords[order_index]
    atoms = atoms[symbol_index][order_index]
    fragment = fragment_uniques[fragment_index][order_index]

    ret = {'complex': plams.Molecule()}
    [ret['complex'].add_atom(plams.Atom(symbol=atom, coords=coord)) for atom, coord in zip(atoms, coords)]
    for name in fragment_uniques:
        ret[name] = plams.Molecule()

        for atom, frag in zip(ret['complex'], fragment):
            if frag != name:
                continue
            ret[name].add_atom(atom.copy())
    return ret


def _get_calc_info(reader):
    '''
    Function to read useful info about orbitals from kf reader
    '''
    ret = results.Result()

    ret.engine = 'ADF'

    # determine if calculation used relativistic corrections
    # if it did, variable 'escale' will be present in 'SFOs'
    # if it didnt, only variable 'energy' will be present
    ret.relativistic = ('SFOs', 'escale') in reader

    ret.symlabels = reader.read('Symmetry', 'symlab').strip().split()
    ret.symmetry = reader.read('Symmetry', 'grouplabel')
    
    # determine if SFOs are unrestricted or not
    ret.unrestricted_sfos = ('SFOs', 'energy_B') in reader

    # determine if MOs are unrestricted or not
    ret.unrestricted_mos = (ret.symlabels[0], 'eps_B') in reader

    ret.sfo_spins = ['A', 'B'] if ret.unrestricted_sfos else ['AB']
    ret.mo_spins = ['A', 'B'] if ret.unrestricted_mos else ['AB']

    # determine if the calculation used regions or not
    ret.used_regions = reader.read('Geometry', 'nr of fragments') != reader.read('Geometry', 'nr of atoms')
    ret.fragments = reader.read('Geometry', 'fragmenttype').split()
    if not ret.used_regions:
        natom = reader.read('Geometry', 'nr of atoms')
        frags = np.array(ret.fragments)
        atom_order = np.array(reader.read('Geometry', 'atom order index'))
        fragment_index = np.array(reader.read('Geometry', 'fragment and atomtype index')) - 1
        frag_per_atom = frags[fragment_index[natom:]]
        ret.fragments = [f'{frag}:{idx}' for frag, idx in zip(frag_per_atom, atom_order[natom:])]

    return ret


def _square_matrix(S):
    S = np.atleast_1d(S).tolist()
    size = len(S)
    n = int(sqrt(.25 + 2*size) - .5)
    Srows = []
    for i in range(n):
        # start index will be the number of elements before this row
        minidx1 = i * (i+1) // 2
        # stop index will be the number of elements of the next row
        maxidx1 = (i+1) * (i+2) // 2
        Srows.append(S[minidx1:maxidx1])

    # then we go through rows again and add the remaining terms
    Srowsfixed = []
    for i, row in enumerate(Srows):
        Srowsfixed.append(row + [row2[i] for row2 in Srows[i+1:]])
    return Srowsfixed


def _read_data(reader, SCF0_reader=None, output=None):
    ret = results.Result()

    with timer.timer('Orbitals.get_data.read_calc_info'):
        ret.calc_info = _get_calc_info(reader)

    def _read_spin_indep(section, variable, spin, R=reader):
        if spin in ['A', 'AB']:
            spin_suffix = '_A'
        else:
            spin_suffix = '_B'

        if not R:
            return

        if (section, variable + spin_suffix) in R:
            return R.read(section, variable + spin_suffix)

        if (section + spin_suffix, variable) in R:
            return R.read(section + spin_suffix, variable)

        if (section, variable) in R:
            return R.read(section, variable)

    def _compose_vector(data, spins):
        return np.hstack([data[spin] for spin in spins])

    def _read_site_energy(spin, scf0=False):
        R = SCF0_reader if scf0 else reader

        site_energy = _read_spin_indep('SFOs', 'site_energy', sfo_spin, R)
        if site_energy:
            return np.atleast_1d(site_energy)

        site_energy = []
        for symlabel in ret.calc_info.symlabels:
            # This is to correct for symlable being split up into :1, :2 and :3, thus 1E:1 becomes 1E
            if ':' in symlabel and symlabel.split(':')[1].isdigit():
                Fock_symlabel = symlabel.split(':')[0]
            else:
                Fock_symlabel = symlabel

            if ret.calc_info.unrestricted_sfos:
                fock_A =  _read_spin_indep('SFO_Fock_A', Fock_symlabel, sfo_spin, R)
                fock_B =  _read_spin_indep('SFO_Fock_B', Fock_symlabel, sfo_spin, R)
                fmats = [fock_A, fock_B]
            else:
                fock = _read_spin_indep('SFO_Fock', Fock_symlabel, sfo_spin, R)
                fmats = [fock]

            for fmat in fmats:
                if fmat is None:
                    return
                idx = 0
                loop = 2
                while idx <= len(fmat):
                    site_energy.append(fmat[idx])
                    idx += loop
                    loop += 1

        return np.atleast_1d(site_energy)

    def _read_kinetic_energy():
        if not output:
            return

        with open(output) as outp:
            lines = outp.readlines()

        read = False
        store = []
        for line in lines:
            if '----------------------------------' in line:
                continue

            if read and 'Total :' in line:
                break

            if read:
                store.extend(line.strip().split())

            if 'Orbital Kinetic Energies (hartree)' in line:
                read = True
                continue

        Ekin = {}
        curr_irrep = None
        for part in store:
            try:
                float(part)
                is_float = True
            except ValueError:
                is_float = False

            if not is_float and not part in ['(equivalent', 'subspecies)', '----']:
                Ekin[part] = []
                curr_irrep = part
            elif not part in ['(equivalent', 'subspecies)', '----']:
                Ekin[curr_irrep].append(float(part))
            if part == 'subspecies)':
                Ekin[curr_irrep] = Ekin[curr_irrep.split(':')[0] + ':1']

        return Ekin

    def _compose_matrix(data, spins):
        blocks = [np.array(data[symlabel][spin]) for symlabel in ret.calc_info.symlabels for spin in spins]
        shapes = [block.shape for block in blocks]
        total_shape = sum(shape[0] for shape in shapes), sum(shape[1] for shape in shapes)
        out = np.zeros(total_shape)

        current_start_index = 0
        for block in blocks:
            out[current_start_index:current_start_index + block.shape[0], current_start_index:current_start_index + block.shape[1]] = block
            current_start_index += block.shape[0]

        return out

    with timer.timer('Orbitals.get_data.read_molecules'):
        ret.SFOs.fragment_molecules = _get_molecules(reader)

    with timer.timer('Orbitals.get_data.read_sfo_data'):
        ret.SFOs.number = reader.read('SFOs', 'number')
        ret.SFOs.fragtype_to_fragidx = np.atleast_1d(reader.read('Geometry', 'fragmenttype').split())
        ret.SFOs.fragment_types = np.atleast_1d(reader.read('SFOs', 'fragtype').split())
        ret.SFOs.fragment_index = np.atleast_1d(reader.read('SFOs', 'fragment'))


        ret.SFOs.fragorb = np.atleast_1d(reader.read('SFOs', 'fragorb'))
        ret.SFOs.subspecies = np.atleast_1d(reader.read('SFOs', 'subspecies').split())
        ret.SFOs.subspecies_fixed = []
        ret.SFOs.symmetry_index = np.atleast_1d(reader.read('SFOs', 'isfo')) - 1

        subspecies_visited_symm_index = {}
        for subsp, isfo in zip(ret.SFOs.subspecies, ret.SFOs.symmetry_index):
            subspecies_visited_symm_index.setdefault(subsp, [])
            if isfo in subspecies_visited_symm_index[subsp]:
                if ret.calc_info.used_regions:
                    n = int(subsp.split(':')[1])
                    subsp = subsp.split(':')[0] + ':' + str(n + 1)
                    subspecies_visited_symm_index.setdefault(subsp, [])
                else:
                    subsp = str(subspecies_visited_symm_index[subsp].count(isfo)) + subsp.split(':')[0]
                    subspecies_visited_symm_index.setdefault(subsp, [])

            subspecies_visited_symm_index[subsp].append(isfo)
            ret.SFOs.subspecies_fixed.append(subsp)

        ret.SFOs.ifo = np.atleast_1d(reader.read('SFOs', 'ifo')) - 1
        ret.SFOs.spin = [spin for spin in ret.calc_info.sfo_spins for _ in range(ret.SFOs.number)]

        ret.SFOs.fragment_unique = {spin: ret.SFOs.fragment_types for spin in ret.calc_info.sfo_spins}
        if not ret.calc_info.used_regions:
            ret.SFOs.fragment_unique = {spin: [f'{frag_name}:{frag_idx}' for frag_name, frag_idx in zip(ret.SFOs.fragment_types, ret.SFOs.fragment_index)] for spin in ret.calc_info.sfo_spins}

        ret.SFOs.fragment_unique.total = _compose_vector(ret.SFOs.fragment_unique, ret.calc_info.sfo_spins)
        for sfo_spin in ret.calc_info.sfo_spins:
            ret.SFOs.energy[sfo_spin] = np.atleast_1d(_read_spin_indep('SFOs', 'escale', sfo_spin))
            ret.SFOs.occupation[sfo_spin] = np.atleast_1d(_read_spin_indep('SFOs', 'occupation', sfo_spin))
            ret.SFOs.order[sfo_spin] = np.argsort(ret.SFos.energy[sfo_spin])

            ret.SFOs.site_energy[sfo_spin] = _read_site_energy(sfo_spin, False)
            if SCF0_reader:
                ret.SFOs.site_energy_SCF0[sfo_spin] = _read_site_energy(sfo_spin, True)

            for symlabel in ret.calc_info.symlabels:
                energy_by_symlabel = ret.SFOs.energy[sfo_spin][ret.SFOs.subspecies_fixed == symlabel]
                ret.SFOs.order_by_symlabel[symlabel][sfo_spin] = np.argsort(energy_by_symlabel)


        ret.SFOs.energy.total = _compose_vector(ret.SFOs.energy, ret.calc_info.sfo_spins)
        ret.SFOs.occupation.total = _compose_vector(ret.SFOs.occupation, ret.calc_info.sfo_spins)
        ret.SFOs.order.total = np.argsort(ret.SFOs.energy.total)

        ###
        # for frag_idx, energy in zip(ret.SFOs.fragment_index, ret.SFOs.energy.total):
        #     print(frag_idx, energy) 
        ###

        for spin in ret.calc_info.sfo_spins:
            if spin == 'AB':
                ret.SFOs.adf_names[spin] = [f'{index + 1}{symlabel}' for index, symlabel in zip(ret.SFOs.ifo, ret.SFOs.subspecies_fixed)]
            else:
                ret.SFOs.adf_names[spin] = [f'{index + 1}{symlabel}_{spin}' for index, symlabel in zip(ret.SFOs.ifo, ret.SFOs.subspecies_fixed)]

        ret.SFOs.adf_names.total = _compose_vector(ret.SFOs.adf_names, ret.calc_info.sfo_spins)
        ret.SFOs.unique_names = {spin: [f'{frag}({name})' for frag, name in zip(ret.SFOs.fragment_unique[spin], ret.SFOs.adf_names[spin])] for spin in ret.calc_info.sfo_spins}
        ret.SFOs.unique_names.total = _compose_vector(ret.SFOs.unique_names, ret.calc_info.sfo_spins)

    with timer.timer('Orbitals.get_data.read_matrices'):
        ret.MOs.nfrozencores = {symlabel: ncbs for symlabel, ncbs in zip(ret.calc_info.symlabels, ensure_list(reader.read('Symmetry', 'ncbs')))}
        ret.MOs.nfrozencores.total = sum(ret.MOs.nfrozencores.values())
        for symlabel in ret.calc_info.symlabels:
            with timer.timer('Orbitals.get_data.read_matrices.overlap'):
                for sfo_spin in ret.calc_info.sfo_spins:
                    S = _read_spin_indep(symlabel, 'S-CoreSFO', sfo_spin)
                    S = _square_matrix(S)
                    ret.matrices.overlap[symlabel][sfo_spin] = S

            with timer.timer('Orbitals.get_data.read_matrices.fock'):
                for sfo_spin in ret.calc_info.sfo_spins:
                    F = _read_spin_indep('SFO_Fock', symlabel.split(':')[0], sfo_spin)
                    if F:
                        F = _square_matrix(F)
                        ret.matrices.fock[symlabel][sfo_spin] = F

            for mo_spin in ret.calc_info.mo_spins:
                nmo = _read_spin_indep(symlabel, 'nmo', mo_spin)
                ret.MOs.number[symlabel][mo_spin] = nmo
                ret.MOs.energy[symlabel][mo_spin] = np.atleast_1d(_read_spin_indep(symlabel, 'escale', mo_spin))
                occupation = np.atleast_1d(_read_spin_indep(symlabel, 'froc', mo_spin))
                ret.MOs.occupation[symlabel][mo_spin] = occupation

                with timer.timer('Orbitals.get_data.read_matrices.coefficients'):
                    coefficients = np.atleast_2d(_read_spin_indep(symlabel, 'Eig-CoreSFO', mo_spin))
                    coefficients = coefficients.reshape(nmo, -1)
                    ret.matrices.coefficients[symlabel][mo_spin] = coefficients

                with timer.timer('Orbitals.get_data.read_matrices.mulliken_analysis'):
                    if ret.calc_info.sfo_spins == ret.calc_info.mo_spins:
                        S = ret.matrices.overlap[symlabel][mo_spin]
                    else:
                        S = ret.matrices.overlap[symlabel].AB

                    with timer.timer('Orbitals.get_data.read_matrices.mulliken_analysis.contribution'):
                        ret.matrices.mulliken_contribution[symlabel][mo_spin] = coefficients * (coefficients @ S)
                    with timer.timer('Orbitals.get_data.read_matrices.mulliken_analysis.population'):
                        ret.matrices.mulliken_population[symlabel][mo_spin] = np.atleast_2d(occupation).T * ret.matrices.mulliken_contribution[symlabel][mo_spin]

    with timer.timer('Orbitals.get_data.read_mo_data'):
        ret.MOs.energy.total = np.hstack([_compose_vector(ret.MOs.energy[symlabel], ret.calc_info.mo_spins) for symlabel in ret.calc_info.symlabels])
        ret.MOs.occupation.total = np.hstack([_compose_vector(ret.MOs.occupation[symlabel], ret.calc_info.mo_spins) for symlabel in ret.calc_info.symlabels])
        ret.MOs.order.total = np.argsort(ret.MOs.energy.total)
        ret.MOs.number.total = len(ret.MOs.energy.total)
        ret.MOs.spin = [spin for spin in ret.calc_info.mo_spins for _ in range(ret.MOs.number.total)]
        ret.MOs.kinetic_energy = _read_kinetic_energy()

    with timer.timer('Orbitals.get_data.compose_matrices'):
        ret.matrices.overlap.total =                _compose_matrix(ret.matrices.overlap,                ret.calc_info.sfo_spins)
        ret.matrices.coefficients.total =           _compose_matrix(ret.matrices.coefficients,           ret.calc_info.mo_spins)
        ret.matrices.mulliken_contribution.total =  _compose_matrix(ret.matrices.mulliken_contribution,  ret.calc_info.mo_spins)
        ret.matrices.mulliken_population.total =    _compose_matrix(ret.matrices.mulliken_population,    ret.calc_info.mo_spins)

    with timer.timer('Orbitals.get_data.gross_population'):
        ret.SFOs.symlabel = []
        ret.MOs.symlabel = []
        ret.MOs.symmetry_index = []

        for mo_spin in ret.calc_info.mo_spins:
            ret.SFOs.gross_population[mo_spin] = []

        for symlabel in ret.calc_info.symlabels:
            norb = ret.MOs.number[symlabel][ret.calc_info.mo_spins[0]]
            ret.SFOs.symlabel.extend([symlabel] * norb)
            ret.MOs.symlabel.extend([symlabel] * norb)
            ret.MOs.symmetry_index.extend(range(norb))
            for mo_spin in ret.calc_info.mo_spins:
                gp = ret.matrices.mulliken_population[symlabel][mo_spin]
                gp = np.sum(gp, axis=0)[ret.MOs.nfrozencores[symlabel]:]
                ret.SFOs.gross_population[mo_spin].extend(gp.tolist())

    return ret
