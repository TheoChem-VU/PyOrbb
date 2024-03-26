from pyfmo.orbitals import info
import numpy as np
from tcutility import ensure_list



def read_SFO_data(reader):  # noqa: N802
    def square_overlaps(S, nmo):
        Srows = []
        for i in range(nmo):
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

    calc_info = info.get_calc_info(reader)

    ret = {}
    # symlabels
    if ('Symmetry', 'symlab') in reader:
        ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    else:
        # if we cannot read the symlabels we default to NOSYM which has only the fully symmetric irrep "A"
        ret['symlabels'] = ['A']

    ret['symlabel_by_sfo'] = [] 

    # number of SFOs per symlabel
    ret['nsfo'] = {}
    total = 0
    for i, symlabel in enumerate(ret['symlabels']):
        norb = int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_A'))))
        total += norb
        ret['nsfo'][symlabel] = norb
        ret['symlabel_by_sfo'].extend(norb * [symlabel])

    ret['nsfo']['total'] = total

    # isfo, index of sfo in symlabel
    if ('SFOs', 'isfo') not in reader:
        ret['isfo'] = range(1, total+1)
    else:
        ret['isfo'] = reader.read('SFOs', 'isfo')

    # fragment, index of fragment sfo belongs to
    ret['fragidx'] = reader.read('SFOs', 'fragment')

    # fragmenttype, name of fragment sfo belongs to
    if ('SFOs', 'fragtype') not in reader:
        fragtypes = [typ.strip() for typ in reader.read('Geometry', 'fragmenttype').split()]
        frag_order = reader.read('Geometry', 'fragment and atomtype index')
        frag_order = frag_order[len(frag_order)//2:]
        ret['fragtypes'] = [fragtypes[frag_order[i-1]-1] for i in ret['fragidx']]
    else:
        ret['fragtypes'] = reader.read('SFOs', 'fragtype').strip().split()

    # get fragment + index unique pairs to construct unique fragment names
    ret['fraguniquenames'] = []
    unqs = set([(name, i) for name, i in zip(ret['fragtypes'], ret['fragidx'])])
    for name, i in zip(ret['fragtypes'], ret['fragidx']):
        # check if fragname is unique
        if len([pair for pair in unqs if pair[0] == name]) == 1:
            ret['fraguniquenames'].append(name)
        else:
            ret['fraguniquenames'].append(f'{name}:{i}')
            
    # fragorb, index of sfo in fragment
    if ('SFOs', 'fragorb') in reader:
        ret['fragorb'] = reader.read('SFOs', 'fragorb')
    else:
        ret['fragorb'] = None

    # subspecies sfo belongs to
    ret['subspecies'] = reader.read('SFOs', 'subspecies').strip().split()

    # ifo, order of sfo
    ret['ifo'] = reader.read('SFOs', 'ifo')

    # if the calculation did not use regions we have to
    # fix the ordering of the atoms (i.e. fragment indices)
    if not calc_info['used_regions']:
        atom_order_index = reader.read('Geometry', 'atom order index')
        atom_order_index = atom_order_index[:len(atom_order_index)//2]
        ret['fragidx'] = [atom_order_index.index(i) + 1 for i in ret['fragidx']]

    # SFO energies
    # and get index if energies are sorted
    energyprefix = 'escale' if calc_info['relativistic'] else 'energy'
    if ('SFOs', energyprefix) in reader:
        if calc_info['unrestricted_sfos']:
            ret['energy'] = {
                'A': reader.read('SFOs', f'{energyprefix}'),
                'B': reader.read('SFOs', f'{energyprefix}_B')
            }
            ret['energyidx'] = {
                'A': np.argsort(ret['energy']['A']),
                'B': np.argsort(ret['energy']['B']),
            }
        else:
            ret['energy'] = {'AB': reader.read('SFOs', f'{energyprefix}')}
            ret['energyidx'] = {'AB': np.argsort(ret['energy']['AB'])}
    else:
        ret['energy'] = {
                'A': None,
                'B': None,
                'AB': None
            }

    # SFO occupations
    if ('SFOs', 'occupation') in reader:
        if calc_info['unrestricted_sfos']:
            ret['occupations'] = {
                'A': np.array(reader.read('SFOs', 'occupation')),
                'B': np.array(reader.read('SFOs', 'occupation_B'))
            }
        else:
            ret['occupations'] = {'AB': np.array(reader.read('SFOs', 'occupation'))}
    else:
        ret['occupations'] = None

    # get indices relative to HOMO for each fragment
    def get_rel_indices(spin):
        relindices = []
        newidx = 0
        homoidx = 0
        for i in range(len(ret['ifo'])-1):
            occ = ret['occupations'][spin][i]
            occ2 = ret['occupations'][spin][i+1]
            if occ > 0 and occ2 == 0:
                homoidx = i
            if (occ == 0 and occ2 > 0):
                relindices.extend(np.arange(newidx - homoidx, i - homoidx + 1))
                newidx = i + 1
            if i == len(ret['ifo'])-2:
                relindices.extend(np.arange(newidx - homoidx, i - homoidx + 2))
        return np.array(relindices)

    if ret['occupations'] is not None:
        if calc_info['unrestricted_sfos']:
            ret['relindices'] = {
                'A': get_rel_indices('A'),
                'B': get_rel_indices('B')
            }
        else:
            ret['relindices'] = {'AB': get_rel_indices('AB')}
    else:
        ret['relindices'] = None

    # overlaps
    ret['overlaps'] = {}
    for symlabel in ret['symlabels']:
        nsfo = ret['nsfo'][symlabel]
        if calc_info['unrestricted_sfos']:
            if (symlabel, 'S-CoreSFO') in reader:
                S_A = ensure_list(reader.read(symlabel, 'S-CoreSFO'))
                S_B = ensure_list(reader.read(symlabel, 'S-CoreSFO_B'))
                ret['overlaps'][symlabel] = {
                    'A': square_overlaps(S_A, nsfo),
                    'B': square_overlaps(S_B, nsfo)
                }
            else:
                ret['overlaps'][symlabel] = {
                    'A': None,
                    'B': None
                }
        else:
            if (symlabel, 'S-CoreSFO') in reader:
                S = ensure_list(reader.read(symlabel, 'S-CoreSFO'))
                ret['overlaps'][symlabel] = {'AB': square_overlaps(S, nsfo)}
            else:
                ret['overlaps'][symlabel] = {'AB': None}

    return ret


def read_MO_data(reader):  # noqa: N802
    calc_info = info.get_calc_info(reader)

    ret = {}
    # number of MOs
    ret['nmo'] = {}
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            print(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_B'))))
            ret['nmo'][symlabel] = {
                'A': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_A')))),
                'B': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_B'))))
            }
        else:
            ret['nmo'][symlabel] = {'AB': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_A'))))}

    # MO energies
    ret['energy'] = {}
    ret['energy']['sorted'] = {'A': [], 'B': [], 'AB': []}
    energyprefix = 'escale' if calc_info['relativistic'] else 'eps'
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['energy'][symlabel] = {
                'A': ensure_list(reader.read(symlabel, f'{energyprefix}_A')),
                'B': ensure_list(reader.read(symlabel, f'{energyprefix}_B'))
            }
            ret['energy']['sorted']['A'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_A')))
            ret['energy']['sorted']['B'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_B')))
        else:
            ret['energy'][symlabel] = {'AB': ensure_list(reader.read(symlabel, f'{energyprefix}_A'))}
            ret['energy']['sorted']['AB'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_A')))
    ret['energy']['sorted']['A'] = sorted(ret['energy']['sorted']['A'])
    ret['energy']['sorted']['B'] = sorted(ret['energy']['sorted']['B'])
    ret['energy']['sorted']['AB'] = sorted(ret['energy']['sorted']['AB'])

    ret['energy_idx'] = {}
    all_idx = {'A': [], 'B': [], 'AB': []}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            energies_A = ensure_list(reader.read(symlabel, f'{energyprefix}_A'))
            idx_A = []
            for energy_A in energies_A:
                idx = ret['energy']['sorted']['A'].index(energy_A)
                while idx in all_idx['A']:
                    idx += 1
                idx_A.append(idx) 
                all_idx['A'].append(idx)

            energies_B = ensure_list(reader.read(symlabel, f'{energyprefix}_B'))
            idx_B = []
            for energy_B in energies_B:
                idx = ret['energy']['sorted']['B'].index(energy_B)
                while idx in all_idx['B']:
                    idx += 1
                idx_B.append(idx)
                all_idx['B'].append(idx)

            ret['energy_idx'][symlabel] = {
                'A': idx_A,
                'B': idx_B
            }
        else:
            energies = ensure_list(reader.read(symlabel, f'{energyprefix}_A'))
            idxs = []
            for energy_B in energies:
                idx = ret['energy']['sorted']['AB'].index(energy_B)
                while idx in all_idx['AB']:
                    idx += 1
                idxs.append(idx)
                all_idx['AB'].append(idx)
            ret['energy_idx'][symlabel] = {'AB': idxs}

    # MO occupations
    ret['occs'] = {}
    ret['noccs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['occs'][symlabel] = {
                'A': ensure_list(reader.read(symlabel, 'froc_A')),
                'B': ensure_list(reader.read(symlabel, 'froc_B'))
            }

            ret['noccs'][symlabel] = {
                'A': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_A')) if occ > 0]),
                'B': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_B')) if occ > 0])
            }
        else:
            ret['occs'][symlabel] = {'AB': ensure_list(reader.read(symlabel, 'froc_A'))}
            ret['noccs'][symlabel] = {'AB': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_A')) if occ > 0]),}
 
    # MO coefficients
    ret['coeffs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            c_A = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_A'))
            c_B = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_B'))
            ret['coeffs'][symlabel] = {
                'A': np.array(c_A).reshape(ret['nmo'][symlabel]['A'], ret['nmo'][symlabel]['A']),
                'B': np.array(c_B).reshape(ret['nmo'][symlabel]['B'], ret['nmo'][symlabel]['B'])
            }
        else:
            c = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_A'))
            ret['coeffs'][symlabel] = {'AB': np.array(c).reshape(ret['nmo'][symlabel]['AB'], ret['nmo'][symlabel]['AB'])}

    # get index of MO in symmetry label
    ret['symmidx'] = {}

    return ret
