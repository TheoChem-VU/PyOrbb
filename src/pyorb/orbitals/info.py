from pyorb import orbitals


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



def read_SFO_data(reader):
    program = get_calc_info(reader)['engine']
    if program == 'ADF':
        return orbitals.adf.read_SFO_data(reader)

    elif program == 'dftb':
        return orbitals.dftb.read_SFO_data(reader)


def read_MO_data(reader):
    program = get_calc_info(reader)['engine']
    if program == 'ADF':
        return orbitals.adf.read_MO_data(reader)

    elif program == 'dftb':
        return orbitals.dftb.read_MO_data(reader)
