import pyorbb
import tcmu
import os
import numpy as np
from typing import List, Union


def dir2mp4(input_dir: str):
    '''
    Turn files in a folder into mp4 format.
    '''
    import moviepy as mvp

    figure_files = [f for f in os.listdir(input_dir) if f.endswith('.png')]
    frames = [os.path.join(input_dir, f) for f in sorted(figure_files, key=lambda f: int(f.split('.')[0]))]

    # take the file locations of the frames we saved and make a movie
    clips = [mvp.ImageClip(frame).with_duration(1/14) for frame in frames]
    concat_clip = mvp.concatenate_videoclips(clips, method="compose")
    concat_clip.write_videofile(input_dir + '.mp4', fps=14)
    concat_clip.write_gif(input_dir + '.gif', fps=14)


def switch_files(file1, file2):
    os.rename(file1, '.tmp')
    os.rename(file2, file1)
    os.rename('.tmp', file2)


def get_histogram(image_path):
    img = PIL.Image.open(image_path).convert("RGB")
    return img.histogram()


def compare_histograms(hist1, hist2):
    # Simple distance (sum of absolute differences)
    return sum(abs(a - b) for a, b in zip(hist1, hist2))


class MOTracker:
    def __init__(self, orbitals: List[pyorbb.Orbitals] = None, rkf_files: List[os.PathLike] = None):
        # if we are given orbital objects we load them first
        if orbitals is not None:
            self.orbital_objects = orbitals
        # otherwise we load rkf files if they were given
        elif rkf_files is not None:
            self.orbital_objects = [pyorbb.Orbitals(rkf_file) for rkf_file in rkf_files]
        # otherwise we make an empty list and can fill it later
        else:
            self.orbital_objects = []

    def append(self, value: Union[pyorbb.Orbitals, os.PathLike]):
        '''
        Append a new |Orbitals| object to the ``MOTracker``.

        Args:
            value: The value to be appended. Given either as a |Orbitals| object or as the path
                to an ``adf.rkf`` file that will be loaded.
        '''
        if isinstance(value, pyorbb.Orbitals):
            self.orbital_objects.append(value)
        elif isinstance(value, os.PathLike):
            self.orbital_objects.append(pyorbb.Orbitals(value))
        else:
            raise TypeError(f'Appended value must be of pyorbb.Orbitals type or os.PathLike type.')

    def insert(self, pos: int, value: Union[pyorbb.Orbitals, os.PathLike]):
        '''
        Insert a new |Orbitals| object in the ``MOTracker``.

        Args:
            pos: The position where to insert the new |Orbitals| value.
            value: The value to be appended. Given either as a |Orbitals| object or as the path
                to an ``adf.rkf`` file that will be loaded.
        '''
        if isinstance(value, pyorbb.Orbitals):
            self.orbital_objects.insert(pos, value)
        elif isinstance(value, os.PathLike):
            self.orbital_objects.insert(pos, pyorbb.Orbitals(value))
        else:
            raise TypeError(f'Inserted value must be of pyorbb.Orbitals type or os.PathLike type.')

    def remove(self, value: Union[pyorbb.Orbitals, os.PathLike]):
        '''
        Remove a new |Orbitals| object from the ``MOTracker``.

        Args:
            value: The value to be removed. Given either as a |Orbitals| object or as the path
                to an ``adf.rkf`` file that was loaded.
        '''
        if isinstance(value, pyorbb.Orbitals):
            self.orbital_objects.remove(value)
        elif isinstance(value, os.PathLike):
            orbs = [orbs for orbs in self.orbital_objects if os.path.samefile(orbs.kfpath, value)]
            self.orbital_objects.remove(orbs)
        else:
            raise TypeError(f'Removed value must be of pyorbb.Orbitals type or os.PathLike type.')

    def index(self, value: Union[pyorbb.Orbitals, os.PathLike]) -> int:
        '''
        Return the index of an |Orbitals| object in the ``MOTracker``.

        Args:
            value: The value to return the index for. Given either as a |Orbitals| object or as the path
                to an ``adf.rkf`` file that was loaded.
        '''
        if isinstance(value, pyorbb.Orbitals):
            return self.orbital_objects.index(value)
        elif isinstance(value, os.PathLike):
            orbs = [orbs for orbs in self.orbital_objects if os.path.samefile(orbs.kfpath, value)]
            return self.orbital_objects.index(orbs)
        else:
            raise TypeError(f'Indexed value must be of pyorbb.Orbitals type or os.PathLike type.')

    def reverse(self):
        '''
        Reverse the order of the |Orbitals| inside the ``MOTracker``.
        '''
        self.orbital_objects.reverse()

    def energy(self, mo_name: str, start_index: int = 0) -> List[float]:
        '''
        Obtain the energy of an MO along the loaded |Orbitals| objects.
        '''
        mos = self.track_mo(mo_name, start_index=start_index)
        return [mo.energy for mo in mos]

    def geometry(self, *atom_indices: int, relative_to_idx: int = None, **kwargs) -> List[float]:
        '''
        Calculate a geometric parameter of the molecules associated with the loaded |Orbitals| objects.

        Args:
            atom_indices: The indices used to compute the parameter.
            relative_to_idx: If given, calculate the parameter relative to the structure associated with this index.

        .. seealso::
            See `tcmu.geometry.parameter <https://theochem-vu.github.io/TCMU/api/tcmu.html#tcmu.geometry.parameter>`_ 
            for information about the function and arguments used to calculate the geometric parameters.
        
        Examples:
            Obtain the bond distance between atoms 1 and 3:

            .. code-block:: python

                >>> tracker.geometry(1, 3)
                [1.1999995214413095, 1.262068462303791, 1.3241374030662725, ...]
            
            We can obtain the bond-stretching by specifying the ``relative_to_idx`` argument:

            .. code-block:: python

                >>> tracker.geometry(1, 3, relative_to_idx=0)
                [0.0, 0.06206894086248149, 0.12413788162496298, ...]

        .. note::
            Atom counting starts at 1 and the indices used are the internal atom numbers.

        '''
        mols = [orbs.molecule for orbs in self.orbital_objects]
        params = [float(tcmu.geometry.parameter(mol, *[i - 1 for i in atom_indices])) for mol in mols]

        if relative_to_idx is not None:
            return [p - params[relative_to_idx] for p in params]
        return params

    def track_mo(self, mo_name: str, start_index: int = 0) -> List[pyorbb.orbitals.objects.MO]:
        initial_mo = self.orbital_objects[start_index].mos[mo_name]
        return self.__track_mo(initial_mo)

    def __track_mo(self, initial_mo: pyorbb.orbitals.objects.MO) -> List[pyorbb.orbitals.objects.MO]:
        '''
        Track an |MO| called ``mo_name`` across a number of |Orbitals| object starting at the first |Orbitals| index.
        It calculates for each |Orbitals| object the overlap population matrix and compares it to the previous |Orbitals| object.
        This way we can track the |MO| across multiple geometries, even if the order of the |MO|s changes.
        '''
        # first find which Orbital object the initial_mo corresponds to
        initial_orb = None
        for i, orb in enumerate(self.orbital_objects):
            if initial_mo in orb.mos.orbitals:
                initial_orb = orb
                initial_orb_index = i
                break

        if initial_orb is None:
            raise ValueError(f'{initial_mo} is not part of any of the supplied Orbitals objects!')

        def _track(_orbs):
            # this tracks the initial_mo across the specified orbital objects (not including initial_orb)

            # retrieve the data associated with the initial MO
            # we need the overlap matrix S and coefficient vector C
            S_prev = initial_orb.data['matrices']['overlap']['total']
            C_prev = np.array([[fmo.coefficient(initial_mo) for fmo in initial_orb.fmos.orbitals]])

            # we calculate the initial overlap population P
            P_prev = C_prev.T * C_prev * S_prev

            # we then obtain a sequence of MOs that match
            mo_sequence = []
            # starting from the second orbital
            for orb in _orbs:
                # obtain the new overlap and coefficient matrices
                S_curr = orb.data['matrices']['overlap']['total']

                # Cmat is now the coefficient matrix instead of the coefficient
                # vector. This greatly speeds up comparison
                Cmat_curr = orb.data['matrices']['coefficients']['total']

                # we wrangle the matrix into the correct shape so that
                # we can multiply the coefficients in one go
                Cmat_curr = Cmat_curr.reshape(-1, *Cmat_curr.shape)

                # precalculate the products of the coefficient vectors
                Cmat2_curr = (Cmat_curr.T * Cmat_curr)

                # we need to reshape the product matrix so the first index
                # corresponds to the MO indices
                Cmat2_curr = np.moveaxis(Cmat2_curr, 1, 0)

                # calculate the overlap populations
                # the first index corresponds to the MOs
                Pmat_curr = Cmat2_curr * S_curr

                # calculate the differences of the populations
                closeness_curr = np.sum(abs(Pmat_curr - P_prev), axis=(1, 2))

                # get the closest orbital
                closest_idx = np.argmin(closeness_curr)

                # obtain the corresponding MO
                closest = orb.mos.orbitals[closest_idx]

                mo_sequence.append(closest)

                # set the current data to the previous data for next cycle
                P_prev = Pmat_curr[closest_idx]

            return mo_sequence

        # based on which initial orbital was selected we 
        if initial_orb_index == 0:
            return [initial_mo] + _track(self.orbital_objects[1:])
        elif initial_orb_index == (len(self.orbital_objects) - 1):
            return _track(self.orbital_objects[:-1][::-1])[::-1] + [initial_mo]
        else:
            return _track(self.orbital_objects[:initial_orb_index][::-1])[::-1] + [initial_mo] + _track(self.orbital_objects[initial_orb_index+1:][::-1])[::-1]



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    calc_dir = '/Users/yumanhordijk/PhD/Projects/Colleagues/formartin/calcs'
    calc_files = [f for f in os.listdir(calc_dir) if f.startswith('complex')]
    calc_files = list(sorted(calc_files, key=lambda f: int(f.split('.')[1].removeprefix('0'))))
    adf_rkf_files = [os.path.join(calc_dir, f, 'adf.rkf') for f in calc_files]

    T = MOTracker(rkf_files=adf_rkf_files)
    plt.figure(figsize=(3,2))
    # plt.plot(T.geometry(1, 3), T.energy('11SIGMA'), label=r'11$\Sigma$')
    plt.plot(T.geometry(1, 3), T.energy('12SIGMA'), label=r'12$\Sigma$')
    plt.plot(T.geometry(1, 3), T.energy('13SIGMA'), label=r'13$\Sigma$')
    plt.xlabel('Cr---Cr distance / A')
    plt.ylabel('MO energy / eV')
    # plt.legend(frameon=False)
    plt.gca().spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    # plt.show()

    plt.figure(figsize=(3,2))
    # plt.plot(T.geometry(1, 3), [orbs.mos['11SIGMA'].energy for orbs in T.orbital_objects], label=r'11$\Sigma$')
    plt.plot(T.geometry(1, 3), [orbs.mos['12SIGMA'].energy for orbs in T.orbital_objects], label=r'12$\Sigma$')
    plt.plot(T.geometry(1, 3), [orbs.mos['13SIGMA'].energy for orbs in T.orbital_objects], label=r'13$\Sigma$')
    plt.xlabel('Cr---Cr distance / A')
    plt.ylabel('MO energy / eV')
    plt.legend(frameon=False)
    plt.gca().spines[['right', 'top']].set_visible(False)
    plt.tight_layout()
    plt.show()


# # print(adf_rkf_files)
# # adf_rkf_files = list(sorted(adf_rkf_files, key=lambda f: int(f.split('(')[1].split(')')[0])))
# adf_rkf_files = list(sorted(adf_rkf_files, key=lambda f: int(f.split('/')[-2].split('.')[1].removeprefix('0'))))
# # print(adf_rkf_files)
# # exit()

# orbs = []
# for file in adf_rkf_files:
#     orbs.append(pyorbb.Orbitals(file))

# # mols = [orb.molecule for orb in orbs]
# # dists = [mol[1].distance_to(mol[2]) for mol in mols]


# import matplotlib.pyplot as plt
# import tcviewer
# import tcutility
# import shutil
# import PIL

# T = tcutility.geometry.Transform()
# T.rotate(y=0.5*np.pi, )

# with tcviewer.Screen(headless=True) as scr:
#     for mo_name in ['1DELTA:x2-y2', '5PI:y', '11SIGMA', '12SIGMA', '13SIGMA']:
#         mos = track_mo(orbs, orbs[0].mos[mo_name])
#         mo_name = mos[0].name

#         fig_dir_neg = f'figures/neg_{mo_name}'
#         fig_dir_pos = f'figures/pos_{mo_name}'
#         os.makedirs(fig_dir_neg, exist_ok=True)
#         os.makedirs(fig_dir_pos, exist_ok=True)
#         for i, mo in enumerate(mos):
#             with scr.add_molscene() as scene:
#                 scene.transform = T.to_vtkTransform()
#                 scene.draw_molecule(mo.molecule)
#                 scene.draw_dual_isosurface(mo.vtk_file(), isovalue=0.03, opacity=0.35, shininess=1)
#                 scene.screenshot(f'{fig_dir_neg}/{i}.png', enable_transparency=False)

#             with scr.add_molscene() as scene:
#                 scene.transform = T.to_vtkTransform()
#                 scene.draw_molecule(mo.molecule)
#                 scene.draw_dual_isosurface(mo.vtk_file(), isovalue=-0.03, opacity=0.35, shininess=1)
#                 scene.screenshot(f'{fig_dir_pos}/{i}.png', enable_transparency=False)

#             if i > 0:
#                 hist_target = get_histogram(f'{fig_dir_pos}/{i-1}.png')

#                 hist_1 = get_histogram(f'{fig_dir_pos}/{i}.png')
#                 hist_2 = get_histogram(f'{fig_dir_neg}/{i}.png')

#                 # Compare
#                 score1 = compare_histograms(hist_target, hist_1)
#                 score2 = compare_histograms(hist_target, hist_2)

#                 if score1 > score2:
#                     switch_files(f'{fig_dir_pos}/{i}.png', f'{fig_dir_neg}/{i}.png')

#         dir2mp4(fig_dir_pos)
#         dir2mp4(fig_dir_neg)


# from itertools import product

# def plot_overlaps(orbs, sfo_names):
#     plt.figure()

#     # combs = list(product(sfo_names, sfo_names))
#     combs = []
#     for i, name1 in enumerate(sfo_names):
#         for j, name2 in enumerate(sfo_names[i:]):
#             combs.append((name1, name2))
#     print(combs)
#     # combs = [(sfo_names[i],)]
#     # print(combs)
#     overlaps = {}

#     for comb in combs:
#         overlaps[comb] = []
#         for orb in orbs:
#             sfo1, sfo2 = orb.sfos[f'frag1({comb[0]})'], orb.sfos[f'frag2({comb[1]})']
#             print(sfo1, sfo2)
#             overlaps[comb].append(abs(sfo1 @ sfo2))

#     for comb, overlaps in overlaps.items():
#         plt.plot(overlaps, label=f'<{comb[0]} | {comb[1]} >')

#     plt.legend()




# def plot_contr(orbs, mo_name, sfo_names):
#     plt.figure()
#     mos = track_mo(orbs, orbs[0].mos[mo_name])
#     contributions = {sfo_name: [] for sfo_name in sfo_names}
#     total = [0] * len(orbs)
#     for i, (orb, mo) in enumerate(zip(orbs, mos)):
#         for sfo_name in sfo_names:
#             sfo = orb.sfos[sfo_name]
#             contributions[sfo_name].append(sfo.mulliken_contribution(mo)*100*2)

#             total[i] += sfo.mulliken_contribution(mo)*100*2

#     for sfo_name, contr in contributions.items():
#         plt.plot(contr, label=sfo_name)

#     plt.plot(total, label='CumSum', c='k', linestyle='dashed')

#     plt.title(mo_name)
#     plt.legend()
#     plt.ylabel('Contribution / %')


# # mo = orbs[0].mos['11SIGMA']
# # for sfo in orbs[0].sfos:
# #     print(sfo, sfo.mulliken_contribution(mo))
# # exit()

# plot_overlaps(orbs, ['2PI:x', '3PI:x', '4PI:x'])
# plt.savefig('figures/overlap_PI.png', dpi=500)
# plt.close()
# plot_overlaps(orbs, ['6SIGMA', '7SIGMA', '8SIGMA'])
# plt.savefig('figures/overlap_SIGMA.png', dpi=500)
# plt.close()
# plot_overlaps(orbs, ['1DELTA:x2-y2'])
# plt.savefig('figures/overlap_DELTA.png', dpi=500)
# plt.close()
# # plt.show()
# plot_contr(orbs, '5PI:x', ['frag1(2PI:x)', 'frag1(3PI:x)', 'frag1(4PI:x)'])
# plt.savefig('figures/contr_5PI.png', dpi=500)
# plt.close()
# plot_contr(orbs, '11SIGMA', ['frag1(6SIGMA)', 'frag1(7SIGMA)', 'frag1(8SIGMA)'])
# plt.savefig('figures/contr_11SIGMA.png', dpi=500)
# plt.close()
# plot_contr(orbs, '12SIGMA', ['frag1(6SIGMA)', 'frag1(7SIGMA)', 'frag1(8SIGMA)'])
# plt.savefig('figures/contr_12SIGMA.png', dpi=500)
# plt.close()
# plot_contr(orbs, '13SIGMA', ['frag1(6SIGMA)', 'frag1(7SIGMA)', 'frag1(8SIGMA)'])
# plt.savefig('figures/contr_13SIGMA.png', dpi=500)
# plt.close()
# # plt.show()
# # pi5x = track_mo(orbs, orbs.mos[0]['5PI:x'])

# # for orb, mo in zip(orbs, pi5x):
# #     sfos = 


# # sigma11 = track_mo(orbs, orbs.mos[0]['11SIGMA'])
# # sigma12 = track_mo(orbs, orbs.mos[0]['12SIGMA'])
# # sigma13 = track_mo(orbs, orbs.mos[0]['13SIGMA'])
# # delta1 = track_mo(orbs, orbs.mos[0]['1DELTA:x2-y2'])



