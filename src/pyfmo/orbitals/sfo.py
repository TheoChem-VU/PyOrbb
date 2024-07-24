from scm import plams
import numpy as np
import os
# import matplotlib.pyplot as plt
from tcutility import ensure_list, squeeze_list
from pyfmo.orbitals import info, mo
from typing import List

j = os.path.join


class SFOs:
    def __init__(self, kfpath=None, reader=None, path_SCF0=None, **kwargs):
        assert reader or kfpath, 'Please provide a KFReader or path to a kf-file'
        self.reader = reader or plams.KFReader(kfpath)
        self.kfpath = kfpath
        self.path_SCF0 = path_SCF0
        if not self.kfpath:
            self.kfpath = os.path.abspath(self.reader.path)
        else:
            self.kfpath = os.path.abspath(kfpath)
        self.get_calc_info()
        self.get_sfos()

    def __iter__(self):
        return iter(self.sfos)

    def _decode_key(self, key):
        '''
        Keys are given in the following format:

            {fragname}[:{fragidx}]({orbname}[ {symmetry}])[_{spin}]

        Where [:fragidx] is optional
        '''

        if isinstance(key, int):
            return {'index': key}

        # get spin from the key
        spin = None
        orbname = key
        for spin_part in ['_A', '_B', '_AB']:
            if key.endswith(spin_part):
                spin = spin_part[1:]
                orbname = key.removesuffix(spin_part)
                break

        # split key into fragment name and orbname 
        if '(' in orbname:
            fragname, orbname = orbname.split('(')
            orbname = orbname.strip(')')
        else:
            fragname = None

        symmetry = None
        if ' ' in orbname:
            orbname, symmetry = orbname.split()

        # extract fragment index from fragment name if present
        if fragname is not None and ':' in fragname:
            fragname, fragidx = fragname.split(':')
            fragidx = int(fragidx)
        else:
            fragidx = None

        return {'fragname': fragname, 'fragidx': fragidx, 'orbname': orbname, 'spin': spin, 'symmetry': symmetry}

    def get_sfo(self, orbname=None, fragname=None, fragidx=None, spin=None, index=None, symmetry=None):
        ret = []
        for sfo in self.sfos:
            if orbname is not None and orbname not in [sfo.name, sfo.relname, sfo.index_name]:
                if ':' in sfo.name and sfo.name.startswith(orbname):
                    pass
                else:
                    continue

            if fragname is not None and fragname != sfo.fragment:
                continue

            if spin is not None and spin != sfo.spin:
                continue

            if fragidx is not None and fragidx != sfo.fragment_index:
                continue

            if index is not None and index != sfo.index:
                continue

            if symmetry is not None and symmetry != sfo.symmetry:
                continue

            ret.append(sfo)

        return ret

    def __getitem__(self, key):
        # if isinstance(key, int):
        #     return squeeze_list(self.get_sfo(index=key))

        if isinstance(key, (tuple, list)):
            ret = []
            for key_ in key:
                ret.extend(ensure_list(self.__getitem__(key_)))    
            return squeeze_list(ret)

        if isinstance(key, (str, int)):
            if key in self.fragments:
                return squeeze_list(self.get_fragment_sfos(key))

            decoded = self._decode_key(key)
            ret = self.get_sfo(**decoded)
            return squeeze_list(ret)

        if isinstance(key, slice):
            start_sfo = ensure_list(self.__getitem__(key.start) or self.sfos[0])
            stop_sfo = ensure_list(self.__getitem__(key.stop) or self.sfos[-1])

            if key.start:
                start_decoded = self._decode_key(key.start)
                start_decoded.pop('orbname', None)
                start_decoded.pop('index', None)
            elif key.stop:
                start_decoded = self._decode_key(key.stop)
                start_decoded.pop('orbname', None)
                start_decoded.pop('index', None)

            start_index = min(sfo.index for sfo in start_sfo)
            stop_index = max(sfo.index for sfo in stop_sfo)

            ret = []
            for index in range(start_index, stop_index + 1):
                ret.extend(self.get_sfo(index=index, **start_decoded))
            return squeeze_list(ret)

    def get_fragment_sfos(self, fragname):
        ret = [sfo for sfo in self.sfos if sfo.fragment == fragname]
        if len(ret) == 0:
            ret = [sfo for sfo in self.sfos if sfo.fragment.startswith(fragname)]
        return ret

    @property
    def fragments(self):
        return set([sfo.fragment_unique_name for sfo in self.sfos])

    def rename_fragments(self, old, new):
        old = ensure_list(old)
        new = ensure_list(new)
        for sfo in self:
            if sfo.fragment not in old:
                continue
            idx = old.index(sfo.fragment)
            sfo.fragment_unique_name = sfo.fragment_unique_name.replace(sfo.fragment, new[idx])
            sfo.fragment = new[idx]

    def get_calc_info(self):
        calc_info = info.get_calc_info(self.reader)
        self.is_unrestricted = calc_info['unrestricted_sfos']
        self.uses_molecular_fragments = calc_info['used_regions']
        self.uses_atomic_fragments = not calc_info['used_regions']
        self.is_relativistic = calc_info['relativistic']

        self.spins = ['A', 'B'] if self.is_unrestricted else ['AB']
        self.symlabels = calc_info['symlabels']

    def get_sfos(self):
        data = info.read_SFO_data(self.reader, path_SCF0=self.path_SCF0)
        sfo_data = []
        for idx in range(data['nsfo']['total']):
            for spin in self.spins:
                symlabel = data['symlabel_by_sfo'][idx]
                isfo = data['isfo'][idx] - 1
                subspecies = data['subspecies'][idx]
                ifo = data['ifo'][idx]

                sfo_data.append({
                    'index':                idx + 1,
                    'name':                 f'{ifo}{subspecies}',
                    'fragment_index':       data['fragidx'][idx],
                    'fragment':             data['fragtypes'][idx],
                    'fragment_unique_name': data['fraguniquenames'][idx],
                    'fragment_orb_index':   None if data['fragorb'] is None else data['fragorb'][idx],
                    'symmetry_type_index':  isfo,
                    'symmetry':             symlabel, 
                    'energy':               None if data['energy'][spin] is None else data['energy'][spin][idx] * 27.21139664,
                    'spin':                 spin,
                    'reader':               self.reader,
                    'kfpath':               self.kfpath,
                    'overlaps':             None if data['overlaps'][symlabel][spin] is None else data['overlaps'][symlabel][spin][isfo],
                    'occupation':           None if data['occupations'] is None else data['occupations'][spin][idx],
                    'atomic_fragments':     self.uses_atomic_fragments,
                })
                if 'site_energy' in data:
                    sfo_data[-1]['site_energy'] = None if data['site_energy'][spin] is None else data['site_energy'][spin][idx] * 27.21139664
                if 'site_energy_SCF0' in data:
                   sfo_data[-1]['site_energy_SCF0'] = None if data['site_energy_SCF0'][spin] is None else data['site_energy_SCF0'][spin][idx] * 27.21139664
                

        self.sfos = [SFO(**sfo_datum) for sfo_datum in sfo_data]

        #sort the SFOs based on their energy 
        self.sfos = sorted(self.sfos, key=lambda sfo: sfo.energy) 
        
        sfos_frag = {}

        for sfo in self.sfos:
            sfos_frag.setdefault(sfo.fragment_unique_name, [])
            sfos_frag[sfo.fragment_unique_name].append(sfo)
        

        # determine if there are degenerate MOs (mo.degenerate) and, if so, how many (mo.n_degenerate)
        for frag, sfos_list in sfos_frag.items():
            energies = []
            for sfo in sfos_list:
                energies.append(round(sfo.energy, 6))

            for sfo in sfos_list:
                degeneracy = energies.count(round(sfo.energy, 6))
                sfo.n_degenerate = degeneracy
                sfo.degenerate = degeneracy > 1

            HOMO_list = [sfo for sfo in sfos_list if sfo.occupation > 0] 
            LUMO_list = [sfo for sfo in sfos_list if sfo.occupation == 0]
            
            for index_LUMO, sfo_LUMO in enumerate(LUMO_list):
                 sfo_LUMO.relindex = index_LUMO + 1 
                 sfo_LUMO.relname = f'LUMO+{index_LUMO}' if index_LUMO > 0 else 'LUMO'
            
            max_HOMO_index = len(HOMO_list) - 1
            for index_HOMO, sfo_HOMO in enumerate(HOMO_list):
                sfo_HOMO.relindex = index_HOMO + 1 - len(HOMO_list)
                sfo_HOMO.relname = f'HOMO-{max_HOMO_index - index_HOMO}' if max_HOMO_index - index_HOMO > 0 else 'HOMO'


        # for unrestricted molecular fragments we want to assign the SOMO and SUMO and fix the HOMO and LUMO names
        if self.is_unrestricted and self.uses_molecular_fragments:
            # this only applies to unrestricted calculations using molecular fragments
            # SOMO and SUMO are independent for each fragment
            for fragname in self.fragments:
                fragsfos = self.get_fragment_sfos(fragname)
                # first check if fragment is radical or closed shell
                if sum(sfo.occupation for sfo in fragsfos) % 2 == 0:
                    continue
                
                # loop through all sfo's to locate the SOMO, this will be the index where the sum of occupations for a and b spin sfos is 1
                for idx in range(1, len(fragsfos)//2 + 1):
                    sfo_of_idx = [sfo for sfo in fragsfos if sfo.fragment_orb_index == idx]
                    if 0 < sfo_of_idx[0].occupation + sfo_of_idx[1].occupation < 2:
                        somo_idx = idx
                        break

                # then loop through all sfo's again to fix their names
                for idx in range(1, len(fragsfos)//2 + 1):
                    # relindex is 0 for SOMO and SUMO, -1 for HOMO and +1 for LUMO
                    relindex = idx - somo_idx
                    sfo_of_idx = [sfo for sfo in fragsfos if sfo.fragment_orb_index == idx]
                    if relindex == -1:
                        sfo_of_idx[0].relname = 'HOMO'
                        sfo_of_idx[1].relname = 'HOMO'
                    elif relindex < -1:
                        sfo_of_idx[0].relname = f'HOMO{relindex+1}'
                        sfo_of_idx[1].relname = f'HOMO{relindex+1}'
                    elif relindex == 1:
                        sfo_of_idx[0].relname = 'LUMO'
                        sfo_of_idx[1].relname = 'LUMO'
                    elif relindex > 1:
                        sfo_of_idx[0].relname = f'LUMO+{relindex-1}'
                        sfo_of_idx[1].relname = f'LUMO+{relindex-1}'
                    elif idx == somo_idx:
                        if sfo_of_idx[0].occupation == 1:
                            sfo_of_idx[0].relname = 'SOMO'
                            sfo_of_idx[1].relname = 'SUMO'
                        else:
                            sfo_of_idx[1].relname = 'SOMO'
                            sfo_of_idx[0].relname = 'SUMO'


class SFO:
    '''
    lightweight class holding SFO data
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return self.full_name

    def __matmul__(self, other):
        if isinstance(other, mo.MO):
            return other.__matmul__(self)
        return self.get_overlap(other)

    def __rmatmul__(self, other):
        if isinstance(other, mo.MO):
            return other.__rmatmul__(self)
        return self.__matmul__(other)

    def __sub__(self, other):
        if isinstance(other, list):
            return [self.__sub__(sfo) for sfo in other]
        return abs(self.energy - other.energy)

    def __rsub__(self, *args, **kwargs):
        return self.__sub__(*args, **kwargs)

    def get_overlap(self, other: 'SFO') -> float:
        '''
        Return overlap between two SFO objects
        SFOs must have same symmetry and spin, else overlap is 0
        '''
        if isinstance(other, list):
            return [self.get_overlap(sfo) for sfo in other]
        else:
            if self.symmetry != other.symmetry:
                return 0
            if self.spin != other.spin:
                return 0

            return self.overlaps[other.symmetry_type_index]

    def make_name(self, spin=True, frag_name=False, relative_name=False, index_name=False):
        name = ''

        if frag_name:
            name += self.fragment_unique_name + '('

        if relative_name:
            name += self.relname
        elif index_name:
            name += str(self.index) + self.symmetry
        else:
            name += self.name

        if frag_name:
            name += ')'

        if self.spin != 'AB' and spin:
            name += f'_{self.spin}'

        return name

    @property
    def full_name(self):
        return self.make_name(frag_name=True, spin=True)

    @property
    def index_name(self):
        return self.make_name(frag_name=False, spin=True, index_name=True)

    @property
    def relative_name(self):
        return self.make_name(frag_name=False, spin=True, relative_name=True)

    @property
    def occupied(self):
        return self.occupation > 0

    @property
    def singly_occupied(self):
        return self.occupation == 1

    def cube_file(self, gridsize: str = 'medium'):
        '''
        Generate a cube-file for this SFO with a certain grid-size.

        Args:
            gridsize: the size of the grid to generate the cube-file with.
        '''
        from tcutility.job.adf import DensfJob
        from tcintegral import grid

        # start a Densf job to calculate the cube-file. 
        # We want to return the cube-file, so we should wait for it to finish.
        with DensfJob(wait_for_finish=True, use_slurm=False) as job:
            job.orbital(self)
            job.gridsize(gridsize)

        # output_cub_paths returns a list of cube-files generated by the job.
        # we only generate one, so we simply return the first element
        return grid.from_cub_file(job.output_cub_paths[0])

    def draw(self, gridsize: str = 'medium', isovalue: float = 0.03):
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
        cub = self.cube_file(gridsize=gridsize)

        # and draw it with a specified isovalue
        with tcviewer.Screen() as scr:
            scr.draw_cub(cub, isovalue, material=tcviewer.materials.orbital_shiny)



def occ_virt_mask(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO) -> float or np.ndarray:
    '''
    Mask specifying whether a pair of SFOs has one occupied and one virtual orbital
    '''
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append((sfo1.occupied != sfo2.occupied) or sfo1.singly_occupied or sfo2.singly_occupied)
    return np.array(ret).squeeze()


def occ_occ_mask(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO) -> float or np.ndarray:
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(sfo1.occupied and sfo2.occupied)
    return np.array(ret).squeeze()


# @decorators.add_to_func(title=r'$S$', scale=100)
# def overlap(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO) -> float or np.ndarray:
#     ret = []
#     for sfo1 in ensure_list(sfos1):
#         ret.append([])
#         for sfo2 in ensure_list(sfos2):
#             ret[-1].append(abs(sfo1 @ sfo2))
#     return np.array(ret).squeeze()


# @decorators.add_to_func(title=r'$\epsilon$', unit='eV')
# def energy_gap(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO) -> float or np.ndarray:
#     ret = []
#     for sfo1 in ensure_list(sfos1):
#         ret.append([])
#         for sfo2 in ensure_list(sfos2):
#             ret[-1].append(abs(sfo1.energy - sfo2.energy))
#     return np.array(ret).squeeze()


# @decorators.add_to_func(title=r'$\Delta E_{oi}$', scale=1e3, unit=r'10$^3$ eV$^{-1}$')
# def orbint(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO, use_mask: bool = True) -> float or np.ndarray:
#     S = overlap(sfos1, sfos2)
#     dE = energy_gap(sfos1, sfos2)
#     oi = np.array(S**2/dE)
#     if use_mask:
#         mask = occ_virt_mask(sfos1, sfos2)
#         oi[np.logical_not(mask)] = None
#     return oi.squeeze()


# @decorators.add_to_func(title=r'$\Delta E_{Pauli}$', cmap='Reds')
# def pauli(sfos1: List[SFO] or SFO, sfos2: List[SFO] or SFO, use_mask: bool = True) -> float or np.ndarray:
#     S = overlap(sfos1, sfos2)
#     pauli = np.array(S**2)
#     if use_mask:
#         mask = occ_occ_mask(sfos1, sfos2)
#         pauli[np.logical_not(mask)] = None
#     return pauli.squeeze()


# def sort_sfo_pairs(sfos1, sfos2, prop=orbint):
#     '''
#     Sort pairs from sfos1 and sfos2 based on the values prop(sfos1, sfos2)
#     args:
#         sfos1, sfos2: lists of SFO objects
#         prop:         function taking SFO objects or lists of SFO objects
#     return:
#         list of tuples containing (sfo1, sfo2, prop(sfo1, sfo2)) sorted by prop(sfo1, sfo2)
#         here sfo1 and sfo2 are taken from sfos1 and sfos2
#     '''
#     M = prop(sfos1, sfos2)
#     ret = []
#     for i, sfo1 in enumerate(ensure_list(sfos1)):
#         for j, sfo2 in enumerate(ensure_list(sfos2)):
#             if np.isnan(M[i, j]):
#                 continue
#             ret.append((sfo1, sfo2, M[i, j]))

#     ret = sorted(ret, key=lambda pair: pair[-1])
#     return ret


# def plot_sfos_prop(sfos1, sfos2, prop=orbint, propargs=None, propkwargs=None, cmap='Greens', title=None, use_relname=False, use_indexname=False, scale=None):
#     if cmap is None:
#         cmap = 'Greens'
#         if hasattr(prop, 'cmap'):
#             cmap = prop.cmap

#     if scale is None:
#         scale = 1
#         if hasattr(prop, 'scale'):
#             scale = prop.scale

#     prop_name = ''
#     if hasattr(prop, '__name__'):
#         prop_name = prop.__name__

#     if title is None:
#         title = prop_name
#         if hasattr(prop, 'title'):
#             title = prop.title

#     if callable(prop):
#         M = prop(sfos1, sfos2, *propargs, **propkwargs)
#     else:
#         M = prop

#     plotname = sfos1[0].spin + ' ' + sfos1[0].kfpath
#     plt.figure(figsize=(10, 8), num=f'{prop_name} {plotname}')
#     occ_virt_border1 = [i for i in range(1, len(sfos1)) if sfos1[i-1].occupation != sfos1[i].occupation]
#     occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
#     occ_virt_border2 = [i for i in range(1, len(sfos2)) if sfos2[i-1].occupation != sfos2[i].occupation]
#     occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
#     plt.imshow(M, origin='lower', cmap=cmap)
#     # gridlines
#     plt.hlines(y=np.arange(0, len(sfos1))+0.5, xmin=np.full(len(sfos1), -0.5), xmax=np.full(len(sfos1), len(sfos2)-0.5), color="w", linewidth=1.5)
#     plt.vlines(x=np.arange(0, len(sfos2))+0.5, ymin=np.full(len(sfos2), -0.5), ymax=np.full(len(sfos2), len(sfos1)-0.5), color="w", linewidth=1.5)
#     # occ_virt border lines
#     plt.vlines(occ_virt_border2-.5, -.5, len(sfos1)-.5, colors='k', linewidth=2)
#     plt.hlines(occ_virt_border1-.5, -.5, len(sfos2)-.5, colors='k', linewidth=2)
#     # text inside cells
#     for i in range(len(sfos1)):
#         for k in range(len(sfos2)):
#             val = M[i, k]
#             if np.isnan(val):
#                 continue
#             color = 'w' if val > np.nanmax(M) / 2 else 'k'
#             plt.gca().text(k, i, f'{val*scale:.2f}', ha="center", va="center", color=color, fontsize=8)

#     try:
#         psi1 = r'\phi_{' + sfos1[0].fragment_unique_name + r'}'
#         psi2 = r'\phi_{' + sfos2[0].fragment_unique_name + r'}'
#     except:
#         psi1 = 'Frag1'
#         psi2 = 'Frag2'

#     plt.xlabel('$'+psi2+'$', fontsize=16)
#     plt.ylabel('$'+psi1+'$', fontsize=16)
#     yticks = range(len(sfos1))
#     xticks = range(len(sfos2))
#     if use_relname:
#         plt.xticks(xticks, [orb.relative_name for orb in sfos2], rotation=90)
#         plt.yticks(yticks, [orb.relative_name for orb in sfos1], rotation=0)
#     elif use_indexname:
#         plt.xticks(xticks, [orb.index_name for orb in sfos2], rotation=90)
#         plt.yticks(yticks, [orb.index_name for orb in sfos1], rotation=0)
#     else:
#         plt.xticks(xticks, [repr(orb) for orb in sfos2], rotation=90)
#         plt.yticks(yticks, [repr(orb) for orb in sfos1], rotation=0)
#     plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$', fontsize=16)
#     plt.tight_layout()

#     return plot.ShowCaller()


if __name__ == '__main__':
    # p = '../test/orbitals/rkf/BH3NH3.rkf'
    # sfos = SFOs(kfpath=p)

    # sfos_donor = sfos[:'Donor(LUMO+2)']
    # sfos_acceptor = sfos['Acceptor(1A)':'Acceptor(LUMO+2)']
    # sfo_donor_best, sfo_acceptor_best, oi = sort_sfo_pairs(sfos_donor, sfos_acceptor, orbint)[-1]
    # plot_sfos_prop(sfos_donor, sfos_acceptor, orbint, use_relname=True).hold()


    # reader = plams.KFReader('../test/orbitals/rkf/substrate_cat_complex_dftb.rkf')
    # sfos = SFOs(reader=reader)

    # p = '/Users/yumanhordijk/PhD/ychem/calculations2/c1d4ca95a3911eb1f79bf4ef91cc7a88b479d7dc8357860bfdb3e577747ebc3a/transitionstate/EDA/EDA/full/adf.rkf'
    p = '../test/orbitals/rkf/substrate_cat_complex_dftb.rkf'
    sfos = SFOs(kfpath=p)
    for sfo in sfos:
        print(sfo, sfo.relative_name)
    sfos_c = sfos['substrate(HOMO-11)':'substrate(LUMO+3)']
    sfos_h = sfos[:'cat(LUMO+3)']
    # sfos_c_best, sfos_h_best, oi = sort_sfo_pairs(sfos_c, sfos_h, orbint)[-1]
    # plot_sfos_prop(sfos_c, sfos_h, overlap, use_relname=False, use_indexname=True).hold()


    # p = '../test/orbitals/rkf/substrate_cat_complex.rkf'
    # sfos = SFOs(kfpath=p)

    # print(sfos['268A'], sfos['12A'], sfos['268A'] @ sfos['12A'])
