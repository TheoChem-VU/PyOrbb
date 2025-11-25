import numpy as np
from math import sqrt
from scm import plams
from pyfmo.nested_dict import NestedDict
from pyfmo.orbitals import fragments

ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else list(x)  # noqa: E731


def _first_principal_numbers(nfrozen_cores):
    aufbau = [
        "S", 
        "S", 
        "P", "P", "P", 
        "S", 
        "P", "P", "P", 
        "D", "D", "D", "D", "D", 
        "S", 
        "P", "P", "P",  
        "D", "D", "D", "D", "D", 
        "F", "F", "F", "F", "F", "F", "F",
        "S",
        "P", "P", "P", 
        "D", "D", "D", "D", "D", 
        "S",
        "F", "F", "F", "F", "F", "F", "F", 
        "P", "P", "P", 
        "D", "D", "D", "D", "D", 
        "S", 
        "P", "P", "P",
    ]

    S = 1
    P = 2
    D = 3
    F = 4
    for typ in aufbau[:nfrozen_cores]:
        if typ == 'S':
            S += 1
        if typ == 'P':
            P += 1/3
        if typ == 'D':
            D += 1/5
        if typ == 'F':
            F += 1/7

    return {'S': round(S), 'P': round(P), 'D': round(D), 'F': round(F)}


def _get_fragoccupations(reader: plams.KFReader) -> dict:
    '''
    Read the fragment occupations from a calculation.
    
    Returns:
        A dictionary containing the alpha and beta spin occupations for each fragment and irrep.
    '''
    inp = reader.read('General', 'engine input')
    if 'fragoccupations' not in inp.lower():
        return {}

    lines = []
    read = False
    for line in inp.splitlines():
        line = line.strip()
        if line.lower() == 'fragoccupations':
            read = True
            continue
        if read and line.lower() == 'end':
            break

        if read:
            lines.append(line)

    lines_lower = [line.lower() for line in lines]
    indices = [0]
    while 'subend' in lines_lower[indices[-1]+1:]:
        index = lines_lower[indices[-1]+1:].index('subend')
        indices.append(index + indices[-1] + 2)

    blocks = []
    for start, end in zip(indices, indices[1:]):
        blocks.append(lines[start:end-1])
    
    data = {}
    for block in blocks:
        frag = block[0]
        data[frag] = {}
        for occ in block[1:]:
            irrep, rest = occ.split(' ', 1)
            a, b = rest.split('//')
            Na = sum(int(float(part)) for part in a.split())
            Nb = sum(int(float(part)) for part in b.split())
            data[frag][irrep] = (Na, Nb)

    return data


def _get_calc_info(reader: plams.KFReader) -> dict:
    '''
    Function to read useful info about orbitals from kf reader
    '''
    ret = NestedDict()

    ret.set('engine', 'ADF')

    # determine if calculation used relativistic corrections
    # if it did, variable 'escale' will be present in 'SFOs'
    # if it didnt, only variable 'energy' will be present
    ret.set('relativistic', ('SFOs', 'escale') in reader)

    ret.set('symlabels', reader.read('Symmetry', 'symlab').strip().split())
    ret.set('symmetry', reader.read('Symmetry', 'grouplabel'))

    # determine if SFOs are unrestricted or not
    ret.set('unrestricted_sfos', ('SFOs', 'energy_B') in reader)
    ret.set('sfo_spins', ['A', 'B'] if ret['unrestricted_sfos'] else ['AB'])

    # determine if MOs are unrestricted or not
    ret.set('unrestricted_mos', (ret['symlabels'][0], 'eps_B') in reader)
    ret.set('mo_spins', ['A', 'B'] if ret['unrestricted_mos'] else ['AB'])

    frag_data = fragments.get_fragments_data(reader)

    # determine if the calculation used regions or not
    ret.set('used_regions', not frag_data['used_atomic_fragments'])
    ret.set('fragments', frag_data['fragment_names'])

    # determine the spin polarization of the complex and fragments
    # if we were given fragoccupations we use those
    spin_pols = _get_fragoccupations(reader)
    for frag in ret['fragments']:
        spin_pols.setdefault(frag, {})

    # otherwise, if unrestricted, we use the occupations
    if ret['unrestricted_sfos']:
        frag_index = np.array(reader.read('SFOs', 'fragment'))
        subspecies = np.array([subsp.split(':')[0] for subsp in reader.read('SFOs', 'subspecies').split()])
        occs_A = np.array(reader.read('SFOs', 'occupation'))
        occs_B = np.array(reader.read('SFOs', 'occupation_B'))
        occ_diff = occs_A - occs_B
        for i, frag in enumerate(ret['fragments'], start=1):
            subspecies_of_frag = subspecies[frag_index == i]
            for subsp in np.unique(subspecies_of_frag):
                spin_pols[frag][subsp] = (sum(occs_A[np.logical_and(frag_index == i, subspecies == subsp)]), sum(occs_B[np.logical_and(frag_index == i, subspecies == subsp)]))

    ret.set('sfo_spinpolarizations', spin_pols)

    # determine if we have access to effective orbital energies
    if ('SFOs', 'site_energy') in reader or 'SFO_Fock_A' in reader or 'SFO_Fock' in reader:
        ret.set('has_site_energy', True)
    else:
        ret.set('has_site_energy', False)

    return ret


def _square_matrix(S: np.ndarray) -> np.ndarray:
    '''
    Convert a flattened lower-echelon type matrix into its square matrix.
    This is useful for reading data from AMS calculations as they often
    store symmetric square matrices in this form, e.g. overlap and Fock matrices.
    '''
    # lists are easier to work with in this case
    S = np.atleast_1d(S).tolist()
    size = len(S)
    n = int(sqrt(.25 + 2*size) - .5)  # number of rows and columns
    Srows = []  # this will hold the first m elements of the row
    for i in range(n):
        # start index will be the number of elements before this row
        min_idx = i * (i+1) // 2
        # stop index will be the number of elements of the next row
        max_idx = (i+1) * (i+2) // 2
        Srows.append(S[min_idx:max_idx])

    # then we go through rows again and add the remaining (n-m) terms
    Srowsfixed = []
    for i, row in enumerate(Srows):
        Srowsfixed.append(row + [row2[i] for row2 in Srows[i+1:]])
    return Srowsfixed


def read_data(reader: plams.KFReader, SCF0_reader: plams.KFReader = None, output: str = None):
    ## HELPER FUNCTIONS
    def _read_spin_indep(section, variable, spin, R=reader):
        # this function is used to read a section%variable from
        # an rkf file that can contain the spin state in its variable name
        if spin in ['A', 'AB']:
            spin_suffix = '_A'
        else:
            spin_suffix = '_B'

        if not R:
            return

        # check for section%variable_{spin}
        if (section, variable + spin_suffix) in R:
            return R.read(section, variable + spin_suffix)

        # check for section_{spin}%variable
        if (section + spin_suffix, variable) in R:
            return R.read(section + spin_suffix, variable)

        # check for default case section%variable
        if (section, variable) in R:
            return R.read(section, variable)

    def _compose_vector(data, spins):
        # concatenate vectors of all spins
        return np.hstack([data[spin] for spin in spins])

    def _read_site_energy(spin, scf0=False):
        # function used to read site energies from rkf file with specific spin
        # also can read from scf0 kfreader if specified
        R = SCF0_reader if scf0 else reader

        site_energy = _read_spin_indep('SFOs', 'site_energy', sfo_spin, R)
        if site_energy:
            return np.atleast_1d(site_energy)

        site_energy = []
        for symlabel in ret['calc_info']['symlabels']:
            # This is to correct for symlable being split up into :1, :2, etc., e.g. 1E:1 becomes 1E
            if ':' in symlabel and symlabel.split(':')[1].isdigit():
                Fock_symlabel = symlabel.split(':')[0]
            else:
                Fock_symlabel = symlabel

            if ret['calc_info']['unrestricted_sfos']:
                fock_A =  _read_spin_indep('SFO_Fock_A', Fock_symlabel, sfo_spin, R)
                fock_B =  _read_spin_indep('SFO_Fock_B', Fock_symlabel, sfo_spin, R)
                fmats = [fock_A, fock_B]
            else:
                fock = _read_spin_indep('SFO_Fock', Fock_symlabel, sfo_spin, R)
                fmats = [fock]

            # the fmat is given as a vector instead of a matrix
            # this algorithm is used to read the correct element on the diagonal
            # fmat is given as a diagonal matrix and not a full matrix
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
        # the kinetic energy is read from the output file
        # and not the rkf file
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

            if not is_float and part not in ['(equivalent', 'subspecies)', '----']:
                Ekin[part] = []
                curr_irrep = part
            elif part not in ['(equivalent', 'subspecies)', '----']:
                Ekin[curr_irrep].append(float(part))
            if part == 'subspecies)':
                Ekin[curr_irrep] = Ekin[curr_irrep.split(':')[0] + ':1']

        return Ekin

    def _compose_matrix(data, spins):
        # form a block-diagonal matrix from its corresponding blocks
        blocks = [np.array(data[symlabel][spin]) for symlabel in ret['calc_info']['symlabels'] for spin in spins]
        shapes = [block.shape for block in blocks]
        total_shape = sum(shape[0] for shape in shapes), sum(shape[1] for shape in shapes)
        out = np.zeros(total_shape)

        current_start_index = 0
        for block in blocks:
            out[current_start_index:current_start_index + block.shape[0], current_start_index:current_start_index + block.shape[1]] = block
            current_start_index += block.shape[0]

        return out

    ## MAIN FUNCTION
    ret = NestedDict()

    ret.set('calc_info', _get_calc_info(reader))

    ret.set('fragment_data', fragments.get_fragments_data(reader))
    molecules = ret['fragment_data']['fragment_molecules']
    molecules['complex'] = ret['fragment_data']['complex_molecule']
    ret.set('molecules', molecules)

    # # if we used atomic basis we always have zero spinpol
    # for frag in ret['fragment_data']['fragment_names']:
    #     ret['calc_info']['sfo_spinpolarizations'][frag] = {}


    ret.set('SFOs', 'number', reader.read('SFOs', 'number'))
    # the name of the fragment
    ret.set('SFOs', 'fragment_index', np.atleast_1d(reader.read('SFOs', 'fragment')))
    # we use the fragment data to map the sfo fragment index to the fragment name
    ret.set('SFOs', 'fragment_types', [ret['fragment_data']['sfo_fragtype_to_fragname_map'][i] for i in ret['SFOs']['fragment_index']])
    # the symmlabel of the SFO
    ret.set('SFOs', 'subspecies', np.atleast_1d(reader.read('SFOs', 'subspecies').split()))
    # index of the SFO in its symmlabel
    ret.set('SFOs', 'symmetry_index', np.atleast_1d(reader.read('SFOs', 'isfo')) - 1)

    # the index of the SFO in its symlabel
    ret.set('SFOs', 'ifo', np.atleast_1d(reader.read('SFOs', 'ifo')) - 1)
    ret.set('SFOs', 'spin', [spin for spin in ret['calc_info']['sfo_spins'] for _ in range(ret['SFOs']['number'])])

    ret.set('SFOs', 'subspecies_fixed', [])

    # some symmetry species can have a subspecies
    # for example, C(3V) symmetry has the E1:1 and E1:2 symmetry species
    # however, ADF only reports for one of the (general E1 label)
    if ret['calc_info']['used_regions']:
        subspecies_visited_symm_index = {}
        for subsp, isfo in zip(ret['SFOs']['subspecies'], ret['SFOs']['symmetry_index']):
            subspecies_visited_symm_index.setdefault(subsp, [])
            if isfo in subspecies_visited_symm_index[subsp]:
                n = int(subsp.split(':')[1])
                subsp = subsp.split(':')[0] + ':' + str(n + 1)
                subspecies_visited_symm_index.setdefault(subsp, [])
            subspecies_visited_symm_index[subsp].append(int(isfo))
            ret['SFOs']['subspecies_fixed'].append(subsp)
    else:
        # if we have atomic fragments we might have SFOs with the same name for the same fragment
        # we should give these unique names as well
        # first get the total number of SFOs that have the same subspecies and fragment
        total_counts = {}
        for subsp, ifo, frag in zip(ret['SFOs']['subspecies'], ret['SFOs']['ifo'], ret['SFOs']['fragment_types']):
            total_counts.setdefault(str(frag), {})
            total_counts[str(frag)].setdefault(str(subsp), {})
            total_counts[str(frag)][str(subsp)].setdefault(int(ifo), 0)
            total_counts[str(frag)][str(subsp)][int(ifo)] += 1

        counts = {}
        # then loop again and set the fixed subspecies
        for subsp, ifo, frag in zip(ret['SFOs']['subspecies'], ret['SFOs']['ifo'], ret['SFOs']['fragment_types']):
            # if there is only one SFO with this subspecies for this fragment
            # we simply set the subspecies as its corrected name
            if total_counts[frag][subsp][ifo] == 1:
                ret['SFOs']['subspecies_fixed'].append(subsp)
                continue
            # otherwise we count which SFO we are at and use that
            # to name the SFO subspecies
            counts.setdefault(str(frag), {})
            counts[str(frag)].setdefault(str(subsp), {})
            counts[str(frag)][str(subsp)].setdefault(int(ifo), 0)
            counts[str(frag)][str(subsp)][int(ifo)] += 1
            ret['SFOs']['subspecies_fixed'].append(f'{subsp}/{counts[frag][subsp][ifo]}')

    # read basic information about the sfos here
    for sfo_spin in ret['calc_info']['sfo_spins']:
        ret.set('SFOs', 'energy', sfo_spin, np.atleast_1d(_read_spin_indep('SFOs', 'escale', sfo_spin)))
        ret.set('SFOs', 'occupation', sfo_spin, np.atleast_1d(_read_spin_indep('SFOs', 'occupation', sfo_spin)))
        # the order in terms of the energy of the SFO
        ret.set('SFOs', 'order', sfo_spin, np.argsort(ret['SFOs']['energy'][sfo_spin]))
        s =  _read_site_energy(sfo_spin, False)
        if s is not None:
            ret.set('SFOs', 'site_energy', sfo_spin, s)
        if SCF0_reader:
            s = _read_site_energy(sfo_spin, True)
            if s is not None:
                ret.set('SFOs', 'site_energy_SCF0', sfo_spin, s)

        for symlabel in ret['calc_info']['symlabels']:
            energy_by_symlabel = ret['SFOs']['energy'][sfo_spin][ret['SFOs']['subspecies_fixed'] == symlabel]
            ret.set('SFOs', 'order_by_symlabel', symlabel, sfo_spin, np.argsort(energy_by_symlabel))

    if 'site_energy' in ret['SFOs']:
        ret.set('SFOs', 'site_energy', 'total', _compose_vector(ret['SFOs']['site_energy'], ret['calc_info']['sfo_spins']))

    if 'site_energy_SCF0' in ret['SFOs']:
        ret.set('SFOs', 'site_energy_SCF0', 'total', _compose_vector(ret['SFOs']['site_energy_SCF0'], ret['calc_info']['sfo_spins']))

    ret.set('SFOs', 'energy', 'total', _compose_vector(ret['SFOs']['energy'], ret['calc_info']['sfo_spins']))
    ret.set('SFOs', 'occupation', 'total', _compose_vector(ret['SFOs']['occupation'], ret['calc_info']['sfo_spins']))
    ret.set('SFOs', 'order', 'total', np.argsort(ret['SFOs']['energy']['total']))

    # construct the names of the SFOs as they would appear in ADFLevels
    for spin in ret['calc_info']['sfo_spins']:
        if spin == 'AB':
            ret.set('SFOs', 'adf_names', spin, [f'{index + 1}{symlabel}' for index, symlabel in zip(ret['SFOs']['ifo'], ret['SFOs']['subspecies_fixed'])])
        else:
            ret.set('SFOs', 'adf_names', spin, [f'{index + 1}{symlabel}_{spin}' for index, symlabel in zip(ret['SFOs']['ifo'], ret['SFOs']['subspecies_fixed'])])

    ret.set('SFOs', 'adf_names', 'total', _compose_vector(ret['SFOs']['adf_names'], ret['calc_info']['sfo_spins']))
    # construct here the unique names, e.g. ``NH3(4E1:1)`` that contains both the SFO orbital name and its fragment name
    # ret.set('SFOs', 'unique_names', {spin: [f'{frag}({name})' for frag, name in zip(ret['SFOs']['fragment_unique'][spin], ret['SFOs']['adf_names'][spin])] for spin in ret['calc_info']['sfo_spins']})
    # ret.set('SFOs', 'unique_names', 'total', _compose_vector(ret['SFOs']['unique_names'], ret['calc_info']['sfo_spins']))

    # read the matrix data such as overlaps, coefficients, etc.
    # we correct for the number of frozen cores later
    ret.set('MOs', 'nfrozencores', {symlabel: ncbs for symlabel, ncbs in zip(ret['calc_info']['symlabels'], ensure_list(reader.read('Symmetry', 'ncbs')))})
    ret.set('MOs', 'nfrozencores', 'total', sum(ret['MOs']['nfrozencores'].values()))

    charge_diff = np.atleast_1d(reader.read('Geometry', 'atomtype total charge')) - np.atleast_1d(reader.read('Geometry', 'atomtype effective charge'))
    charge_diff = {atomtype: charge_diff[i] for i, atomtype in enumerate(reader.read('Geometry', 'atomtype').split())}
    # also get the frozencores per atom
    ret.set('MOs', 'nfrozencores_atom', {atomtype: int(diff/2) for atomtype, diff in charge_diff.items()})

    if not ret['calc_info']['used_regions']:
        atomtypes = np.atleast_1d(reader.read('SFOs', 'fragtype').split())
        first_principal = {atomtype: _first_principal_numbers(ret['MOs']['nfrozencores_atom'][atomtype]) for atomtype in np.unique(atomtypes)}
        for spin in ret['calc_info']['sfo_spins']:
            if spin == 'AB':
                ret.set('SFOs', 'adf_names_fixed_principal', spin, [f'{index + first_principal[atomtype][symlabel.split(":")[0].split("/")[0]]}{symlabel}' for index, atomtype, symlabel in zip(ret['SFOs']['ifo'], atomtypes, ret['SFOs']['subspecies_fixed'])])
            else:
                ret.set('SFOs', 'adf_names_fixed_principal', spin, [f'{index + first_principal[atomtype][symlabel.split(":")[0].split("/")[0]]}{symlabel}_{spin}' for index, atomtype, symlabel in zip(ret['SFOs']['ifo'], atomtypes, ret['SFOs']['subspecies_fixed'])])

    for symlabel in ret['calc_info']['symlabels']:
        for sfo_spin in ret['calc_info']['sfo_spins']:
            S = _read_spin_indep(symlabel, 'S-CoreSFO', sfo_spin)
            S = _square_matrix(S)
            ret.set('matrices', 'overlap', symlabel, sfo_spin, S)

        for sfo_spin in ret['calc_info']['sfo_spins']:
            F = _read_spin_indep('SFO_Fock', symlabel.split(':')[0], sfo_spin)
            if F:
                F = _square_matrix(F)
                ret.set('matrices', 'fock', symlabel, sfo_spin, F)

        for mo_spin in ret['calc_info']['mo_spins']:
            nmo = _read_spin_indep(symlabel, 'nmo', mo_spin)
            ret.set('MOs', 'number', symlabel, mo_spin, nmo)
            ret.set('MOs', 'energy', symlabel, mo_spin, np.atleast_1d(_read_spin_indep(symlabel, 'escale', mo_spin)))
            occupation = np.atleast_1d(_read_spin_indep(symlabel, 'froc', mo_spin))
            ret.set('MOs', 'occupation', symlabel, mo_spin, occupation)

            coefficients = np.atleast_2d(_read_spin_indep(symlabel, 'Eig-CoreSFO', mo_spin))
            coefficients = coefficients.reshape(nmo, -1)
            ret.set('matrices', 'coefficients', symlabel, mo_spin, coefficients)

            # perform mulliken analysis here
            # contribution = C * (C @ S)
            # population   = O * contribution
            # see: https://github.com/TheoChem-VU/PyFMO/issues/28
            if ret['calc_info']['sfo_spins'] == ret['calc_info']['mo_spins']:
                S = ret['matrices']['overlap'][symlabel][mo_spin]
            else:
                S = ret['matrices']['overlap'][symlabel]['AB']

            contr = coefficients * (coefficients @ S)
            contr_normed = (contr.T / np.sum(np.maximum(0, contr), axis=1)).T
            ret.set('matrices', 'mulliken_contribution', symlabel, mo_spin, contr)
            ret.set('matrices', 'mulliken_contribution_normalized', symlabel, mo_spin, contr_normed)
            ret.set('matrices', 'mulliken_population', symlabel, mo_spin, np.atleast_2d(occupation).T * contr)

    ret.set('MOs', 'energy', 'total', np.hstack([_compose_vector(ret['MOs']['energy'][symlabel], ret['calc_info']['mo_spins']) for symlabel in ret['calc_info']['symlabels']]))
    ret.set('MOs', 'occupation', 'total', np.hstack([_compose_vector(ret['MOs']['occupation'][symlabel], ret['calc_info']['mo_spins']) for symlabel in ret['calc_info']['symlabels']]))
    ret.set('MOs', 'order', 'total', np.argsort(ret['MOs']['energy']['total']))
    ret.set('MOs', 'number', 'total', len(ret['MOs']['energy']['total']))
    ret.set('MOs', 'spin', [spin for spin in ret['calc_info']['mo_spins'] for _ in range(ret['MOs']['number']['total'])])
    K = _read_kinetic_energy()
    if K is not None:
        ret.set('MOs', 'kinetic_energy', K)

    ret.set('matrices', 'overlap', 'total',                _compose_matrix(ret['matrices']['overlap'],                ret['calc_info']['sfo_spins']))
    ret.set('matrices', 'coefficients', 'total',           _compose_matrix(ret['matrices']['coefficients'],           ret['calc_info']['mo_spins']))
    ret.set('matrices', 'mulliken_contribution', 'total',  _compose_matrix(ret['matrices']['mulliken_contribution'],  ret['calc_info']['mo_spins']))
    ret.set('matrices', 'mulliken_population', 'total',    _compose_matrix(ret['matrices']['mulliken_population'],    ret['calc_info']['mo_spins']))

    # gross population is the vertical marginal of the Mulliken population matrix
    ret.set('SFOs', 'symlabel', [])
    ret.set('MOs', 'symlabel', [])
    ret.set('MOs', 'symmetry_index', [])

    for mo_spin in ret['calc_info']['mo_spins']:
        ret.set('SFOs', 'gross_population', mo_spin, [])

    for symlabel in ret['calc_info']['symlabels']:
        norb = ret['MOs']['number'][symlabel][ret['calc_info']['mo_spins'][0]]
        ret['SFOs']['symlabel'].extend([symlabel] * norb)
        ret['MOs']['symlabel'].extend([symlabel] * norb)
        ret['MOs']['symmetry_index'].extend(range(norb))
        for mo_spin in ret['calc_info']['mo_spins']:
            gp = ret['matrices']['mulliken_population'][symlabel][mo_spin]
            gp = np.sum(gp, axis=0)[ret['MOs']['nfrozencores'][symlabel]:]
            ret['SFOs']['gross_population'][mo_spin].extend(gp.tolist())

    ret.set('SFOs', 'gross_population', 'total', _compose_vector(ret['SFOs']['gross_population'], ret['calc_info']['mo_spins']))

    return ret
