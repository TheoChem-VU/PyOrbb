import numpy as np
from TCutility import ensure_list
from scm import plams


def get_calc_info(reader):
    '''
    Function to read useful info about orbitals from kf reader
    '''
    ret = {}
    # check what the engine is first. We can read DFTB and ADF files for now
    if ('General', 'program') in reader:
        engine = reader.read('General', 'program').strip()
    # if program cannot be read from reader it is probably an old version of ADF, so we should default to ADF
    else:
        engine = 'ADF'
    ret['engine'] = engine

    if engine == 'ADF':
        # determine if calculation used relativistic corrections
        # if it did, variable 'escale' will be present in 'SFOs'
        # if it didnt, only variable 'energy' will be present
        ret['relativistic'] = ('SFOs', 'escale') in reader

        # determine if SFOs are unrestricted or not
        ret['unrestricted_sfos'] = ('SFOs', 'energy_B') in reader

        # get the symmetry labels
        if ('Symmetry', 'symlab') in reader:
            ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
        else:
            # if we cannot read the symlabels we default to NOSYM which has only the fully symmetric irrep "A"
            ret['symlabels'] = symmetry.labels['NOSYM']

        # determine if MOs are unrestricted or not
        ret['unrestricted_mos'] = (ret['symlabels'][0], 'eps_B') in reader

        # determine if the calculation used regions or not
        frag_order = reader.read('Geometry', 'fragment and atomtype index')
        frag_order = frag_order[:len(frag_order)//2]
        ret['used_regions'] = max(frag_order) != len(frag_order)

    elif engine == 'dftb':
        ret['relativistic'] = None
        ret['symlabels'] = None
        ret['unrestricted_mos'] = None
        ret['unrestricted_sfos'] = None
        ret['used_regions'] = reader.read('FragmentOrbitals', 'AtomicFragmentOrbitals') == 0

    return ret


def check_rkf(path):
    reader = plams.KFReader(path)
    missing = []
    for required in required_variables(reader):
        if required not in reader:
            missing.append(required)

    for required in missing:
        print(f'Missing {required[0]}: {required[1]}')

    return len(missing) == 0


def required_variables(reader):
    if ('Symmetry', 'symlab') in reader:
        symlabels = reader.read('Symmetry', 'symlab').strip().split()
    elif ('Geometry', 'grouplabel') in reader:
        symlabels = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]
    else:
        symlabels = symmetry.labels['NOSYM']

    ret = []
    ret = ret + [(symlabel, 'Eig-CoreSFO_A') for symlabel in symlabels] + [(symlabel, 'Eig-CoreSFO_B') for symlabel in symlabels]
    ret = ret + [(symlabel, 'S-CoreSFO') for symlabel in symlabels] + [(symlabel, 'S-CoreSFO_B') for symlabel in symlabels]
    ret = ret + [(symlabel, 'escale_A') for symlabel in symlabels] + [(symlabel, 'escale_B') for symlabel in symlabels]
    ret = ret + [(symlabel, 'eps_A') for symlabel in symlabels] + [(symlabel, 'eps_B') for symlabel in symlabels]
    ret = ret + [(symlabel, 'froc_A') for symlabel in symlabels] + [(symlabel, 'froc_B') for symlabel in symlabels]
    ret = ret + [
        ('Symmetry', 'symlab'),
        ('Geometry', 'grouplabel'),
        ('Geometry', 'atom order index'),
        ('SFOs', 'number'),
        ('SFOs', 'isfo'),
        ('SFOs', 'fragment'),
        ('SFOs', 'fragtype'),
        ('SFOs', 'fragorb'),
        ('SFOs', 'subspecies'),
        ('SFOs', 'ifo'),
        ('SFOs', 'escale'),
        ('SFOs', 'escale_B'),
        ('SFOs', 'energy'),
        ('SFOs', 'energy_B'),
        ('SFOs', 'occupation'),
        ('SFOs', 'occupation_B'),
    ]

    return set(ret)


def read_SFO_data(reader):
    program = get_calc_info(reader)['engine']
    if program == 'ADF':
        from TCutility.orbitals import adf
        return adf.read_SFO_data(reader)

    elif program == 'dftb':
        from TCutility.orbitals import dftb
        return dftb.read_SFO_data(reader)


def read_MO_data(reader):
    program = get_calc_info(reader)['engine']
    if program == 'ADF':
        from TCutility.orbitals import adf
        return adf.read_MO_data(reader)

    elif program == 'dftb':
        from TCutility.orbitals import dftb
        return dftb.read_MO_data(reader)
