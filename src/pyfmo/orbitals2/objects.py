import pyfmo
from scm import plams
from tcutility import timer, cache
import os
import numpy as np


class OrbitalSelector:
    def __init__(self, orbitals, parent):
        self.orbitals = orbitals
        self.parent = parent

    def __getitem__(self, key):
        if isinstance(key, int):
            return [orb for orb in self.orbitals if orb.index == key]
        if isinstance(key, str):
            return self.get(**self.decode_key(key))

    def decode_key(self, key):
        '''
        Keys are given in the following format:

            {fragname}[:{fragidx}]({orbname}[ {symmetry}])[_{spin}]

        Where [:fragidx] is optional
        '''
        decoded = {
            'index': None,
            'fragment': None,
            'fragidx': None,
            'orbname': None,
            'spin': None,
            'symmetry': None,
        }

        if isinstance(key, int):
            decoded['index'] = key
            return decoded

        # get spin from the key
        decoded['orbname'] = key
        for spin_part in ['_A', '_B', '_AB']:
            if key.endswith(spin_part):
                decoded['spin'] = spin_part[1:]
                decoded['orbname'] = key[:-len(spin_part)]

        # split key into fragment name and orbname 
        if '(' in decoded['orbname']:
            decoded['fragment'], decoded['orbname'] = decoded['orbname'].split('(')
            decoded['orbname'] = decoded['orbname'].strip(')')

        if ' ' in decoded['orbname']:
            decoded['orbname'], decoded['symmetry'] = decoded['orbname'].split()

        # extract fragment index from fragment name if present
        if decoded['fragment'] is not None and ':' in decoded['fragment']:
            decoded['fragment'], decoded['fragidx'] = decoded['fragment'].split(':')
            decoded['fragidx'] = int(decoded['fragidx'])
        return decoded


    def get(self, symmetry=None, spin=None, fragment=None, orbname=None, **kargs):
        orbs = self.orbitals
        # print([orb.subspecies for orb in orbs])
        if symmetry:
            if self.parent.data.calc_info.used_regions:
                orbs = [orb for orb in orbs if orb.symmetry == symmetry]
            else:
                orbs = [orb for orb in orbs if orb.subspecies == symmetry]
        if spin:
            orbs = [orb for orb in orbs if orb.spin == spin]
        if fragment:
            orbs = [orb for orb in orbs if orb.fragment_unique == fragment or orb.fragment == fragment]
        if orbname:
            orbs = [orb for orb in orbs if orb.name == orbname or orb.relative_name == orbname]

        if len(orbs) == 0:
            return None
        if len(orbs) == 1:
            return orbs[0]
        return orbs

    def __len__(self):
        return len(self.orbitals)

    def __iter__(self):
        return iter(sorted(self.orbitals, key=lambda orb: orb.energy))

    @property
    def spins(self):
        return {orb.spin for orb in self.orbitals}

    @property
    def unrestricted(self):
        return all(orb.spin in ['A', 'B'] for orb in self.orbitals)


class SFOs(OrbitalSelector):
    @property
    def fragments(self):
        frags = []
        for sfo in self.orbitals:
            if sfo.fragment_unique not in frags:
                frags.append(sfo.fragment_unique)
        return frags

    def get_fragment_sfos(self, fragment):
        return [sfo for sfo in self.orbitals if sfo.fragment_unique == fragment]


class MOs(OrbitalSelector):
    ...


class Orbital:
    def __init__(self, data, parent):
        self.data = data
        for key, value in data.items():
            setattr(self, key, value)
        self.parent = parent

    def __repr__(self):
        return str(self)

    @property
    @cache.cache
    @timer.timer
    def relative_name(self):
        orbitals = [orb for orb in self.parent.orbitals if orb.spin == self.spin and orb.spin_total_occupation == self.spin_total_occupation]
        if hasattr(self, 'fragment_unique'):
            orbitals = [orb for orb in orbitals if orb.fragment_unique == self.fragment_unique]

        energies = sorted([orb.energy for orb in orbitals])
        order = energies.index(self.energy) + self.degeneracy_index

        if self.doubly_occupied:
            order = len(orbitals) - order - 1
            return f'HOMO-{order}' if order > 0 else 'HOMO'

        if self.singly_occupied:
            if self.occupation == 1:
                order = len(orbitals) - order - 1
                return f'SOMO-{order}' if order > 0 else 'SOMO'
            else:
                return f'SUMO+{order}' if order > 0 else 'SUMO'

        if self.unoccupied:
            return f'LUMO+{order}' if order > 0 else 'LUMO'

    @property
    def doubly_occupied(self):
        return self.spin_total_occupation == 2

    @property
    def singly_occupied(self):
        return self.spin_total_occupation == 1

    @property
    def unoccupied(self):
        return self.spin_total_occupation == 0

    @property
    @cache.cache
    @timer.timer
    def spin_total_occupation(self):
        matching_orbs = [orb for orb in self.parent.orbitals if orb.name == self.name]
        if hasattr(self, 'fragment_unique'):
            matching_orbs = [orb for orb in matching_orbs if orb.fragment_unique == self.fragment_unique]

        return sum(orb.occupation for orb in matching_orbs)

    def cube_file(self, gridsize: str = 'medium', overwrite: bool = False):
        '''
        Generate a cube-file for this SFO with a certain grid-size.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
        '''
        from tcutility.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, overwrite=overwrite) as job:
            # job.orbital(self)
            if isinstance(self, SFO):
                job._sfos.append(self)
            else:
                job._mos.append(self)
            job.settings.ADFFile = self.parent.parent.kfpath
            job.gridsize(gridsize)

        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_cub_file(job.output_cub_paths[0])

    def draw(self, gridsize: str = 'medium', isovalue: float = 0.03, overwrite: bool = False):
        '''
        Generate and draw a cube-file for this SFO object.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
            isovalue: the value with which to generate the isosurface of this SFO.

        .. seealso::
            :meth:`SFO.cube_file` to generate and return a cube-file for this SFO.
        '''
        import tcviewer

        # generate a cube-file or load an existing one
        cub = self.cube_file(gridsize=gridsize, overwrite=overwrite)

        # and draw it with a specified isovalue
        with tcviewer.Screen() as scr:
            scr.draw_cub(cub, isovalue, material=tcviewer.materials.orbital_shiny)

    @property
    def degeneracy_index(self):
        return self.degenerate_orbitals.index(self)

    @property
    def degenerate(self):
        return self.degeneracy > 1

    @property
    def degeneracy(self):
        return len(self.degenerate_orbitals)

    @property
    @cache.cache
    @timer.timer
    def degenerate_orbitals(self):
        return [orb for orb in self.parent.orbitals if orb.energy == self.energy]


class MO(Orbital):
    def __str__(self):
        if self.spin == 'AB':
            return f'{self.name}'
        return f'{self.name}_{self.spin}'


class SFO(Orbital):
    def __str__(self): 
        if self.spin == 'AB':
            return f'{self.fragment_unique}({self.name})'
        return f'{self.fragment_unique}({self.name})_{self.spin}'

    def overlap(self, other):
        assert isinstance(other, SFO)

        if self.spin != other.spin:
            return 0

        if self.symmetry != other.symmetry:
            return 0

        S = self.parent.parent.data.matrices.overlap[self.symmetry][self.spin]
        return S[other.symmetry_index][self.symmetry_index]


    def mulliken_contribution(self, other):
        assert isinstance(other, MO)

        if self.spin != other.spin and other.spin != 'AB':
            return 0

        if self.symmetry != other.symmetry:
            return 0

        c = self.parent.parent.data.matrices.mulliken_contribution[self.symmetry][other.spin]
        return c[other.symmetry_index][self.symmetry_index]


    def coefficient(self, other):
        assert isinstance(other, MO)

        if self.spin != other.spin and other.spin != 'AB':
            return 0

        if self.symmetry != other.symmetry:
            return 0

        c = self.parent.parent.data.matrices.coefficients[self.symmetry][other.spin]
        return c[other.symmetry_index][self.symmetry_index]


    def __matmul__(self, other):
        return self.overlap(other)


    def make_name(self, spin=True, frag_name=False, relative_name=False, index_name=False):
        name = ''

        if frag_name:
            name += self.fragment_unique + '('

        if relative_name:
            name += self.relative_name
        elif index_name:
            name += str(self.index) + self.symmetry
        else:
            name += self.name

        if frag_name:
            name += ')'

        if self.spin != 'AB' and spin:
            name += f'_{self.spin}'

        return name


class Orbitals:
    '''
    Container class that stores information about both MO's and SFO's.
    '''
    def __init__(self, path: str, path_SCF0: str = None, moleculename: str = None):
        r'''
        Two kind of readers are constucted.
        1. path provides the path to a fully converged Fragment analyses calculation with a full SCF. From this, all 
            information regarding the fragment analysis is extracted. This includes the SFO energies of the fully isolated 
            fragments and, if available, the site energies or Fock matrix. From this can return the site energies (diagonal 
            of the Fock matrix).

            The energies taken from this file are the SFO energies of the fully isolated fragments and the site energies 
            (diagonal of the Fock matrix) of the fully relaxed complex.

        2. The path_SCF0 is the pathway to the fragment analysis where SCF is set to zero (SCF=0). This is necessary for 
            reading the site energies (diagonal of the Fock matrix) to obtain the corrected energies of the SFOs. No other 
            information is read from this file.

            The energies extracted from this file are the site energies (diagonal of the Fock matrix) of the two fragments 
            in the field of the second respective fragment. This correction is often considered superior to the SFO energies 
            for the fully isolated fragments.
        '''
        self.reader = plams.KFReader(path)
        self.kfpath = os.path.abspath(path)
        self.SCF0_kfpath = path_SCF0
        self.SCF0_reader = plams.KFReader(path_SCF0) if path_SCF0 else None
        with timer.timer('Orbitals.get_data'):
            self.get_data()
        with timer.timer('Orbitals.gather_sfos'):
            self.gather_sfos()
        with timer.timer('Orbitals.gather_mos'):
            self.gather_mos()

    def get_data(self):
        self.data = pyfmo.orbitals2.adf._read_data(self.reader, SCF0_reader=self.SCF0_reader)

    def gather_sfos(self):
        self.sfos = SFOs([], self)
        sfo_mo_spin_match = self.data.calc_info.unrestricted_mos == self.data.calc_info.unrestricted_sfos
        for sfoi in range(self.data.SFOs.number):
            for spin_idx, sfo_spin in enumerate(self.data.calc_info.sfo_spins):
                if not sfo_mo_spin_match:
                    if self.data.calc_info.unrestricted_mos:
                        gross_pop = self.data.SFOs.gross_population.A[sfoi] + self.data.SFOs.gross_population.B[sfoi]
                        gross_spin = self.data.SFOs.gross_population.A[sfoi] - self.data.SFOs.gross_population.B[sfoi]
                    else:
                        gross_pop = self.data.SFOs.gross_population.AB[sfoi]
                        gross_spin = 0
                else:
                    gross_pop = self.data.SFOs.gross_population[sfo_spin][sfoi]
                    gross_spin = 0

                data = {
                    'index': sfoi + 1,
                    # 'name': f'{self.data.SFOs.ifo[sfoi]}{self.data.SFOs.subspecies[sfoi]}',
                    'name': self.data.SFOs.adf_names[sfo_spin][sfoi],
                    'subspecies': self.data.SFOs.subspecies[sfoi],
                    'symmetry': self.data.SFOs.symlabel[sfoi],
                    'symmetry_index': self.data.SFOs.symmetry_index[sfoi],
                    'index_in_symlabel': self.data.SFOs.symmetry_index[sfoi], # rmove this later
                    'fragment': self.data.calc_info.fragments[self.data.SFOs.fragment_index[sfoi] - 1].split(':')[0],
                    'fragment_unique': self.data.SFOs.fragment_unique.total[sfoi],
                    'spin': sfo_spin,
                    'energy': self.data.SFOs.energy[sfo_spin][sfoi] * 27.2114079527,
                    'occupation': int(self.data.SFOs.occupation[sfo_spin][sfoi]),
                    'occupied': int(self.data.SFOs.occupation[sfo_spin][sfoi]) > 0,
                    'gross_population': gross_pop,
                    'gross_spin': gross_spin,
                }

                data['site_energy'] = np.nan
                if isinstance(self.data.SFOs.site_energy[sfo_spin], np.ndarray):
                    data['site_energy'] = self.data.SFOs.site_energy[sfo_spin][sfoi] * 27.2114079527

                data['site_energy_SCF0'] = np.nan
                if isinstance(self.data.SFOs.site_energy_SCF0[sfo_spin], np.ndarray):
                    data['site_energy_SCF0'] = self.data.SFOs.site_energy_SCF0[sfo_spin][sfoi] * 27.2114079527

                sfo = SFO(data, self.sfos)
                self.sfos.orbitals.append(sfo)

    def gather_mos(self):
        self.mos = MOs([], self)
        for moi in range(self.data.SFOs.number):
            for spin_idx, mo_spin in enumerate(self.data.calc_info.mo_spins):
                symm_idx = self.data.MOs.symmetry_index[moi]
                symlabel = self.data.MOs.symlabel[moi]
                data = {
                    'index': moi + 1,
                    'name': f'{symm_idx+1}{symlabel}',
                    'symmetry': symlabel,
                    'symmetry_index': self.data.MOs.symmetry_index[moi],
                    'index_in_symlabel': self.data.MOs.symmetry_index[moi], # rmove this later
                    'spin': mo_spin,
                    'energy': self.data.MOs.energy[symlabel][mo_spin][symm_idx] * 27.2114079527,
                    'occupation': int(self.data.MOs.occupation[symlabel][mo_spin][symm_idx]),
                    'occupied': int(self.data.MOs.occupation[symlabel][mo_spin][symm_idx]) > 0,
                    # 'gross_population': self.data.SFOs.gross_population[mo_spin][moi],
                }
                sfo = MO(data, self.mos)
                self.mos.orbitals.append(sfo)

    @property
    def fragments(self):
        return self.sfos.fragments

    def write_excel(self, out_file: str = 'pyfmo.xlsx'):
        from pyfmo import write_excel
        
        write_excel.to_excel(self, out_file)


    def write_excel2(self, out_file: str = 'pyfmo2.xlsx'):
        from pyfmo import write_excel2
        
        write_excel2.to_excel(self, out_file)
