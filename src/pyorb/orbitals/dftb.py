from pyorb.orbitals import info
import numpy as np
from scm import plams



def read_SFO_data(reader):  # noqa: N802
    calc_info = info.get_calc_info(reader)

    assert calc_info['engine'] == 'dftb'

    ret = {}

    # number of SFOs and number of fragments
    nsfo = reader.read('FragmentOrbitals', 'nOrbitals')
    ret['nsfo'] = {'A': nsfo, 'total': nsfo}
    ret['nfrags'] = reader.read('FragmentOrbitals', 'nFragments')

    ret['symlabel_by_sfo'] = ['A'] * nsfo
    ret['isfo'] = list(range(1, nsfo+1))

    # read the fragment index the orbital belongs to
    ret['fragidx'] = reader.read('FragmentOrbitals', 'Fragment')

    # read names of the fragments (use element names if we did not do a fragment analysis)
    if calc_info['used_regions']:
        ret['fragtypes'] = [reader.read('FragmentOrbitals', 'BaseNameFragFile').strip().split()[i-1] for i in ret['fragidx']]
        # fragment names are already uniques if regions were used
        ret['fraguniquenames'] = ret['fragtypes']
    else:
        ret['fragtypes'] = [reader.read('FragmentOrbitals', 'FragmentSymbols').strip().split()[i-1] for i in ret['fragidx']]
        # make fragment names unique as {element}:{index}
        ret['fraguniquenames'] = [f'{name}:{i}' for i, name in zip(ret['fragidx'], ret['fragtypes'])]

    # fragorb, index of sfo in the fragment its belongs to
    ret['fragorb'] = reader.read('FragmentOrbitals', 'iFO')
    ret['ifo'] = reader.read('FragmentOrbitals', 'iFO')

    # subspecies sfo belongs to
    ret['subspecies'] = reader.read('FragmentOrbitals', 'SubSpecies').strip().split()

    # SFO energies
    # and get index if energies are sorted
    ret['energy'] = {'AB': reader.read('FragmentOrbitals', 'Energies(1)')}
    ret['energyidx'] = np.argsort(ret['energy'])

    # SFO occupations
    ret['occupations'] = {'AB': np.array(reader.read('FragmentOrbitals', 'Occupations(1)'))}

    # get the indices of the orbitals relative to the HOMO for each fragment
    # loop through each fragment
    ret['relindices'] = {'AB': []}
    for fragidx in range(1, ret['nfrags']+1):
        # first make a list of bools signifying if the fragorb is occupied
        is_occupied = []
        for i, occ in enumerate(ret['occupations']['AB']):
            if ret['fragidx'][i] != fragidx:
                continue
            is_occupied.append(occ > 0)
        # build the indices for each fragorb. The index is determined by the number of occupied orbitals higher than this orb
        # if there are only unoccupied orbitals higher than this orb the index is determined by the number of fragments remaining
        relindices = [-sum(is_occupied[i:]) + 1 if sum(is_occupied[i:]) > 0 else i - sum(is_occupied) + 1 for i in range(len(is_occupied))]
        ret['relindices']['AB'].extend(relindices)

    # overlaps
    ret['overlaps'] = {'A': {'AB': np.array(reader.read('FragmentOrbitals', 'Overlaps(1)')).reshape(nsfo, nsfo)}}

    return ret


def read_MO_data(reader):  # noqa: N802
    calc_info = info.get_calc_info(reader)

    assert calc_info['engine'] == 'dftb'

    ret = {}

    ret['symlabels'] = ['A']
    # read the number of MOs
    nmo = reader.read('Orbitals', 'nOrbitals')
    ret['nmo'] = {'A': {'AB': nmo}}
    ret['energy'] = {'A': {'AB': reader.read('Orbitals', 'Energies(1)')}}
    occs = reader.read('Orbitals', 'Occupations(1)')
    ret['occs'] = {'A': {'AB': occs}}
    ret['noccs'] = {'A': {'AB': sum([occ > 0 for occ in occs])}}
    ret['coeffs'] = {'A': {'AB': np.array(reader.read('Orbitals', 'Coefficients(1)')).reshape(nmo, nmo)}}
    ret['energy_idx'] = {'A': {'AB': list(range(nmo))}}

    return ret

if __name__ == '__main__':
    reader = plams.KFReader('../test/orbitals/rkf/substrate_dftb.rkf')
    reader = plams.KFReader('../test/orbitals/rkf/substrate_cat_complex_dftb.rkf')
    print(read_SFO_data(reader))
    print(read_MO_data(reader))
