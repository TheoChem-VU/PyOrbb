from scm import plams
import os
from tcutility import ensure_list, squeeze_list
from pyfmo.orbitals import info

j = os.path.join


class MOs:
    def __init__(self, kfpath=None, reader=None, moleculename=None, **kwargs):
        assert reader or kfpath, 'Please provide a KFReader or path to a kf-file'
        self.reader = reader or plams.KFReader(kfpath)
        self.kfpath = kfpath
        if not self.kfpath:
            self.kfpath = os.path.abspath(self.reader.path)
        else:
            self.kfpath = os.path.abspath(kfpath)

        self.moleculename = moleculename
        self.get_calc_info()
        self.get_mos()

    def __iter__(self):
        return iter(self.mos)

    def _decode_key(self, key):
        '''
        Keys are given in the following format:

            {orbname}[_{spin}]

        Where [_{spin}] is optional
        '''

        # get spin from the key
        index = None
        orbname = None
        spin = None
        if isinstance(key, int):
            index = key
        elif '_' in key:
            orbname, spin = key.split('_')
        else:
            orbname = key

        return {'orbname': orbname, 'spin': spin, 'index': index}

    def get_mo(self, orbname=None, spin=None, index=None):
        ret = []
        for mo in self.mos:
            if orbname is not None and orbname not in [mo.name, mo.relative_name, mo.relname]:
                continue

            if spin is not None and spin != mo.spin:
                continue

            if index is not None and index != mo.index:
                continue

            ret.append(mo)

        return ret

    def __getitem__(self, key):
        if isinstance(key, int):
            return squeeze_list(self.get_mo(index=key))

        if isinstance(key, (tuple, list)):
            ret = []
            for key_ in key:
                ret.extend(ensure_list(self.__getitem__(key_)))    
            return squeeze_list(ret)

        if isinstance(key, str):
            decoded = self._decode_key(key)
            ret = self.get_mo(**decoded)
            return ret[0] if len(ret) == 1 else ret

        if isinstance(key, slice):
            start_decoded = self._decode_key(key.start)
            stop_decoded = self._decode_key(key.stop)
            
            start_sfo = ensure_list(self.get_mo(**start_decoded))
            stop_sfo = ensure_list(self.get_mo(**stop_decoded))

            start_index = min(sfo.index for sfo in start_sfo)
            stop_index = min(sfo.index for sfo in stop_sfo)

            ret = []
            for index in range(start_index, stop_index + 1):
                ret.extend(self.get_mo(index=index))
            return ret

    def get_calc_info(self):
        calc_info = info.get_calc_info(self.reader)
        self.is_unrestricted = calc_info['unrestricted_mos']
        self.uses_molecular_fragments = calc_info['used_regions']
        self.uses_atomic_fragments = not calc_info['used_regions']
        self.is_relativistic = calc_info['relativistic']

        self.spins = ['A', 'B'] if self.is_unrestricted else ['AB']
        self.symlabels = calc_info['symlabels']

    def get_mos(self):
        data = info.read_MO_data(self.reader)

        mo_data = []
        for symlabel in data['symlabels']:
            for spin in self.spins:
                for idx in range(data['nmo'][symlabel][spin]):
                    energy = data['energy'][symlabel][spin][idx]
                    mo_index = data['energy_idx'][symlabel][spin][idx]

                    mo_data.append({
                        'index':                mo_index + 1,
                        'index_in_symlabel':    idx,
                        'name':                 f'{idx + 1}{symlabel}',
                        'moleculename':         self.moleculename,
                        'symmetry':             symlabel,
                        'energy':               energy * 27.21139664,
                        'spin':                 spin,
                        'reader':               self.reader,
                        'kfpath':               self.kfpath,
                        # 'overlaps':             data['overlaps'][symlabel][spin][isfo],
                        'coeffs':               data['coeffs'][symlabel][spin][idx],
                        'occupation':           data['occs'][symlabel][spin][idx],
                        'atomic_fragments':     self.uses_atomic_fragments,
                    })
        self.mos = [MO(**mo_datum) for mo_datum in mo_data]
        
        #sort the MOs based on their energy 
        self.mos = sorted(self.mos, key=lambda mo: mo.energy) 
        
        # determine if there are degenerate MOs (mo.degenerate) and, if so, how many (mo.n_degenerate)
        energies = [round(mo.energy, 6) for mo in self.mos]

        for mo in self.mos:
            degeneracy = energies.count(round(mo.energy, 6))
            mo.n_degenerate = degeneracy
            mo.degenerate = degeneracy > 1


        #Assign correct HOMO-LUMO nomenclature, Regardless of the sequence of the occupied and unoccupied orbitals 
        HOMO_list = [mo for mo in self.mos if mo.occupation > 0] 
        LUMO_list = [mo for mo in self.mos if mo.occupation == 0]
        
        for index_LUMO, mo_LUMO in enumerate(LUMO_list):
             mo_LUMO.relindex = index_LUMO + 1
             mo_LUMO.relname = f'LUMO+{index_LUMO}' if index_LUMO > 0 else 'LUMO'
        
        max_HOMO_index = len(HOMO_list) - 1
        for index_HOMO, mo_HOMO in enumerate(HOMO_list):
            mo_HOMO.relindex = index_HOMO + 1 - len(HOMO_list)
            mo_HOMO.relname = f'HOMO-{max_HOMO_index - index_HOMO}' if max_HOMO_index - index_HOMO > 0 else 'HOMO'


        # Assign name to singly occupied orbitals in the HOMO LUMO nomenclature
        somo_idx = None
        if self.is_unrestricted:
            for idx in range(1, len(self.mos)//2 + 1):
                mos = self.get_mo(index=idx)
                if 0 < sum(mo.occupation for mo in mos) % 2 < 2:
                    somo_idx = idx
                    break
            if somo_idx is not None:
                for idx in range(1, len(self.mos)//2 + 1):
                    relindex = idx - somo_idx
                    mos = self.get_mo(index=idx)
                    if relindex == -1:
                        mos[0].relname = 'HOMO'
                        mos[1].relname = 'HOMO'
                    elif relindex < -1:
                        mos[0].relname = f'HOMO{relindex+1}'
                        mos[1].relname = f'HOMO{relindex+1}'
                    elif relindex == 1:
                        mos[0].relname = 'LUMO'
                        mos[1].relname = 'LUMO'
                    elif relindex > 1:
                        mos[0].relname = f'LUMO+{relindex-1}'
                        mos[1].relname = f'LUMO+{relindex-1}'
                    elif idx == somo_idx:
                        if mos[0].occupation == 1:
                            mos[0].relname = 'SOMO'
                            mos[1].relname = 'SUMO'
                        else:
                            mos[1].relname = 'SOMO'
                            mos[0].relname = 'SUMO'


class MO:
    '''
    lightweight class holding MO data
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return self.full_name

    def __matmul__(self, other):
        return self.get_coeff(other)
    
    def __rmatmul__(self, other):
        return self.__matmul__(other)        

    def get_coeff(self, other):
        if isinstance(other, list):
            return [self.get_coeff(sfo) for sfo in other]

        if self.symmetry != other.symmetry:
            return 0

        if self.spin != other.spin and 'AB' not in [self.spin, other.spin]:
            return 0
            
        return self.coeffs[other.symmetry_type_index]

    def make_name(self, spin=True, relative_name=False, index_name=False, **kwargs):
        name = ''

        if relative_name:
            name += self.relname
        elif index_name:
            name += str(self.index) + self.symmetry
        else:
            name += self.name

        if self.spin != 'AB' and spin:
            name += f'_{self.spin}'

        return name

    @property
    def full_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.name}{spin_part}'

    @property
    def relative_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.relname}{spin_part}'

    @property
    def index_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.index}{self.symmetry}{spin_part}' 

    @property
    def occupied(self):
        return self.occupation > 0

    def cube_file(self, gridsize: str = 'medium'):
        '''
        Generate a cube-file for this MO with a certain grid-size.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
        '''
        from tcutility.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True) as job:
            job.orbital(self)
            job.gridsize(gridsize)

        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_cub_file(job.output_cub_paths[0])

    def draw(self, gridsize: str = 'medium', isovalue: float = 0.03):
        '''
        Generate and draw a cube-file for this MO object.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this MO.

        .. seealso::
            :meth:`MO.cube_file` to generate and return a cube-file for this MO.
        '''
        import tcviewer

        # generate a cube-file or load an existing one
        cub = self.cube_file(gridsize=gridsize)

        # and draw it with a specified isovalue
        with tcviewer.Screen() as scr:
            scr.draw_cub(cub, isovalue, material=tcviewer.materials.orbital_shiny)


if __name__ == '__main__':
    # p = '../test/orbitals/rkf/methyl.rkf'
    p = '../test/orbitals/rkf/substrate_cat_complex_dftb.rkf'
    mos = MOs(kfpath=p)
    for mo in sorted(mos, key=lambda mo: mo.energy):
        print(mo, mo.occupation, mo.relative_name)
