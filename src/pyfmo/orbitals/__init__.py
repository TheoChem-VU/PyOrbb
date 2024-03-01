from pyfmo.orbitals import sfo, mo, adf, dftb  # noqa
from scm import plams
import numpy as np
import matplotlib.pyplot as plt
# from yutility import plot, ensure_list
from tcutility import ensure_list, pathfunc
from typing import Union


class Orbitals:
    def __init__(self, path, path_SCF0=None, moleculename=None):
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
        if isinstance(path, (plams.KFReader, plams.KFFile)):
            self.reader = path
        else:
            self.reader = plams.KFReader(path)
        self.mos = mo.MOs(reader=self.reader, moleculename=moleculename)
        self.sfos = sfo.SFOs(reader=self.reader, path_SCF0=path_SCF0)
        self.rename_fragments = self.sfos.rename_fragments

    def mulliken_contribution(self, sfo, mo):
        r'''
        Calculate the Mulliken contribution of a selected SFO to a selected MO.
        The Mulliken contribution originates from the Mulliken population analysis,
        but before scaling the result by the number of electrons occupying the MO.

        The Mulliken contribution of SFO $\mu$ to MO $i$ is given as

            $$\hat{C}_{i\mu} = c_{i\mu}^2 + \sum_{\nu \neq \mu} c_{i\mu} c_{i\nu} S_{\mu\nu}
                             = \sum_\nu c_{i\mu} c_{i\nu} S_{\mu\nu}$$

            Where the index $\nu$ denotes all SFOs, $c_{i\nu}$ is the coefficient 
            of SFO $\nu$ in MO $i$ and $S_{\mu\nu}$ is the overlap between SFOs
            $\mu$ and $\nu$.

        The marginals of the resulting Mulliken contribution matrix $\hat{C}_{i\mu}$ times the 
        occupation is the gross Mulliken population of the orbital.
        '''
        # coefficient of all SFOs contributing to the selected MO
        c_iv = np.array(mo @ self.sfos.sfos)
        # coefficient of the selected SFO
        c_iu = mo @ sfo
        # overlaps between selected SFO and all other SFOs
        S_uv = np.array(sfo @ self.sfos.sfos)
        # calculate the mulliken contribution
        return np.sum(c_iv*c_iu*S_uv)

    def mulliken_population(self, sfo, mo):
        r'''
        Calculate the Mulliken population of a selected SFO to a selected MO.
        The Mulliken population is the product of the occupation of the MO with the 
        Mulliken contribution of the SFO to the MO.

        The Mulliken population of SFO $\mu$ to MO $i$ is then
            
            $$P_{i\mu} = n_i\hat{C}_{i\mu}$$

            Where $n_i$ is the occupation of MO $i$ and $\hat{C}_{i\mu}$ is the 
            Mulliken contribution of SFO $\mu$ to MO $i$.

        The marginals of the Mulliken population matrix $P_{i\mu}$ should equal
        the occupation numbers of the MOs on axis 1 and the 
        '''
        return mo.occupation * self.mulliken_contribution(sfo, mo)

    @property
    def fragments(self):
        return self.sfos.fragments

    @property
    def spins(self):
        return self.sfos.spins

    def get_alike_orbital(self, other_orb: Union[sfo.SFO, mo.MO], kernel='abscosine', plot=False) -> Union[sfo.SFO, mo.MO]:
        '''
        Find an orbital in this collection of orbitals that looks like the given orbital.
        The selection criterium is based on the coefficients that are applied to the SFO or are used for the MO.
        '''
        if isinstance(other_orb, sfo.SFO):
            abscosine = lambda A, B: abs(A) @ abs(B) / np.sqrt((A @ A) * (B @ B))
            cosine = lambda A, B: abs(A @ B) / np.sqrt((A @ A) * (B @ B))

            def kernel(sfo):
                ret = cosine(np.array(sfo.overlaps), np.array(other_orb.overlaps))
                return ret

            most_similar = max(self.sfos, key=kernel)
            if plot:
                plt.plot([kernel(sfo) for sfo in self.sfos])
                plt.show()

            return most_similar

        elif isinstance(other_orb, mo.MO):
            # if kernel == 'dot':
            #     # k_dot = |A @ B|
            #     kernel = lambda mo: abs(mo.coeffs @ other_orb.coeffs)
            # elif kernel == 'rmsd':
            #     # k_rmsd = 1 - sum((A - B)^2)/N
            #     kernel = lambda mo: 1 - sum((mo.coeffs - other_orb.coeffs)**2)/len(other_orb.coeffs)
            # elif kernel == 'cosine':
            #     # k_cosine = |A @ B| / (||A|| ||B||)
            #     kernel = lambda mo: abs(mo.coeffs @ other_orb.coeffs / ((mo.coeffs @ mo.coeffs) * (other_orb.coeffs @ other_orb.coeffs)))
            # elif kernel == 'abscosine':
            #     # k_cosine = |A| @ |B| / (||A|| ||B||)
            #     kernel = lambda mo: abs(mo.coeffs) @ abs(other_orb.coeffs) / ((mo.coeffs @ mo.coeffs) * (other_orb.coeffs @ other_orb.coeffs))

            abscosine = lambda A, B, C, D: abs(A) @ abs(B) / np.sqrt((C @ C) * (D @ D))
            cosine = lambda A, B, C, D: abs(A @ B) / np.sqrt((C @ C) * (D @ D))

            def kernel(mo):
                ret = 0
                for fragment in self.fragments:
                    idx = [i for i, sfo in enumerate(self.sfos) if sfo.fragment == fragment]
                    ret += abscosine(mo.coeffs[idx], other_orb.coeffs[idx], mo.coeffs, other_orb.coeffs)
                return ret

            most_similar = max(self.mos, key=kernel)
            if plot:
                plt.plot([kernel(mo) for mo in self.mos])
                plt.show()

            return most_similar


def sort_orb_pairs(orbs1, orbs2, prop=None):
    '''
    Sort pairs from sfos1 and sfos2 based on the values prop(sfos1, sfos2)
    args:
        sfos1, sfos2: lists of SFO objects
        prop:         function taking SFO objects or lists of SFO objects
    return:
        list of tuples containing (sfo1, sfo2, prop(sfo1, sfo2)) sorted by prop(sfo1, sfo2)
        here sfo1 and sfo2 are taken from sfos1 and sfos2
    '''
    M = prop(orbs1, orbs2)
    ret = []
    for i, orb1 in enumerate(ensure_list(orbs1)):
        for j, orb2 in enumerate(ensure_list(orbs2)):
            if np.isnan(M[i, j]):
                continue
            ret.append((orb1, orb2, M[i, j]))

    ret = sorted(ret, key=lambda pair: pair[-1])
    return ret


# def plot_property(orbs1, orbs2, prop=None, propargs={}, propkwargs={}, cmap='Greens', title=None, unit=None, use_relname=False, use_indexname=False, scale=None, **kwargs):
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

#     if unit is None:
#         if hasattr(prop, 'unit'):
#             unit = prop.unit

#     if unit is not None:
#         unit = '(' + unit + ')'
#     else:
#         unit = ''
 
#     if callable(prop):
#         M = prop(orbs1, orbs2, *propargs, **propkwargs)
#     else:
#         M = prop

#     plotname = orbs1[0].spin + ' ' + orbs1[0].kfpath
#     plt.figure(figsize=kwargs.get('figsize', (10, 8)), label=f'{prop_name} {plotname}')
#     occ_virt_border1 = [i for i in range(1, len(orbs1)) if orbs1[i-1].occupation != orbs1[i].occupation]
#     occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
#     occ_virt_border2 = [i for i in range(1, len(orbs2)) if orbs2[i-1].occupation != orbs2[i].occupation]
#     occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
#     plt.imshow(M, origin='lower', cmap=cmap)
#     # gridlines
#     plt.hlines(y=np.arange(0, len(orbs1))+0.5, xmin=np.full(len(orbs1), -0.5), xmax=np.full(len(orbs1), len(orbs2)-0.5), color="w", linewidth=1.5)
#     plt.vlines(x=np.arange(0, len(orbs2))+0.5, ymin=np.full(len(orbs2), -0.5), ymax=np.full(len(orbs2), len(orbs1)-0.5), color="w", linewidth=1.5)
#     # occ_virt border lines
#     plt.vlines(occ_virt_border2-.5, -.5, len(orbs1)-.5, colors='k', linewidth=2)
#     plt.hlines(occ_virt_border1-.5, -.5, len(orbs2)-.5, colors='k', linewidth=2)
#     # text inside cells
#     for i in range(len(orbs1)):
#         for k in range(len(orbs2)):
#             val = M[i, k]
#             if np.isnan(val):
#                 continue
#             color = 'w' if val > np.nanmax(M) / 2 else 'k'
#             plt.gca().text(k, i, f'{val*scale:.2f}', ha="center", va="center", color=color, fontsize=8)

#     try:
#         psi1 = r'\phi_{' + orbs1[0].fragment_unique_name + r'}'
#         psi2 = r'\phi_{' + orbs2[0].fragment_unique_name + r'}'
#     except:
#         try:
#             psi1 = r'\phi_{' + orbs1[0].moleculename + r'}'
#             psi2 = r'\phi_{' + orbs2[0].moleculename + r'}'
#         except:
#             psi1 = r'\phi_{' + 'Frag1' + r'}'
#             psi2 = r'\phi_{' + 'Frag2' + r'}'

#     plt.xlabel('$'+psi2+'$', fontsize=16)
#     plt.ylabel('$'+psi1+'$', fontsize=16)
#     yticks = range(len(orbs1))
#     xticks = range(len(orbs2))
#     if use_relname:
#         plt.xticks(xticks, [orb.relative_name for orb in orbs2], rotation=90)
#         plt.yticks(yticks, [orb.relative_name for orb in orbs1], rotation=0)
#     elif use_indexname:
#         plt.xticks(xticks, [orb.index_name for orb in orbs2], rotation=90)
#         plt.yticks(yticks, [orb.index_name for orb in orbs1], rotation=0)
#     else:
#         plt.xticks(xticks, [repr(orb) for orb in orbs2], rotation=90)
#         plt.yticks(yticks, [repr(orb) for orb in orbs1], rotation=0)
#     plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$ ' + unit, fontsize=16)
#     plt.tight_layout()

#     return plot.ShowCaller()


if __name__ == '__main__':
    import os 

    dirs = pathfunc.get_subdirectories('../../../test/fixtures/pyfrag')
    files = [os.path.join(d, 'adf.rkf') for d in sorted(dirs, key=lambda d: int(d.split('.')[-1]))]
    orbs = [Orbitals(f) for f in files]

    print(orbs[0].fragments)
    for idx in range(1, 13):
        plt.plot([orb.mos[f'HOMO-{idx}'].energy for orb in orbs], label=f'HOMO-{idx}')
    # plt.plot([orb.mos['HOMO'].energy for orb in orbs], label='HOMO')
    # plt.plot([orb.mos['HOMO-1'].energy for orb in orbs], label='HOMO-1')
    plt.xlabel('IRC Step')
    plt.ylabel(r'$\epsilon$ (eV)')
    plt.title('Without orbital tracking')
    plt.legend()

    def MO_track(start_mo, orbs):
        mos = [start_mo]
        for orb in orbs[1:]:
            best = orb.get_alike_orbital(mos[-1])
            mos.append(best)
        return mos

    plt.figure()
    for idx in range(1, 13):
        # plt.plot([orb.mos[f'HOMO-{idx}'].energy for orb in orbs], label=f'HOMO-{idx}')
        plt.plot([mo.energy for mo in MO_track(orbs[1].mos[f'HOMO-{idx}'], orbs[1:-1])], label=f'HOMO-{idx}')
    # plt.plot([mo.energy for mo in MO_track(orbs[1].mos['HOMO-1'], orbs[1:-1])], label='HOMO-1')
    plt.xlabel('IRC Step')
    plt.ylabel(r'$\epsilon$ (eV)')
    plt.title('With orbital tracking')
    plt.legend()
    plt.show()



    # dirs = pathfunc.get_subdirectories('../../../test/fixtures/pyfrag')
    dirs = pathfunc.match('/Users/yumanhordijk/PhD/MM2024/calculations/IRC/pi_beta_trans_TS1/pi_beta/pi_beta_trans', 'complex.{index}')
    files = [os.path.join(d, 'adf.rkf') for d in sorted(dirs, key=lambda d: int(d.split('.')[-1]))]
    orbs = [Orbitals(f) for f in files]
    for idx in range(1, 13):
        plt.plot([orb.sfos[f'frag1(HOMO-{idx})'].energy for orb in orbs], label=f'HOMO-{idx}')
    # plt.plot([orb.sfos['frag1(HOMO-5)'].energy for orb in orbs], label='HOMO-5')
    # plt.plot([orb.sfos['frag1(HOMO-6)'].energy for orb in orbs], label='HOMO-6')
    # plt.plot([orb.sfos['frag1(HOMO-7)'].energy for orb in orbs], label='HOMO-7')
    plt.xlabel('IRC Step')
    plt.ylabel(r'$\epsilon$ (eV)')
    plt.title('Without orbital tracking')
    plt.legend()

    def SFO_track(start_sfo, orbs):
        sfos = [start_sfo]
        for orb in orbs[1:]:
            best = orb.get_alike_orbital(sfos[-1])
            sfos.append(best)
        return sfos

    plt.figure()
    for idx in range(1, 13):
        plt.plot([mo.energy for mo in SFO_track(orbs[1].sfos[f'frag1(HOMO-{idx})'], orbs[1:-1])], label=f'HOMO-{idx}')
    # plt.plot([mo.energy for mo in SFO_track(orbs[1].sfos['frag1(HOMO-5)'], orbs[1:-1])], label='HOMO-5')
    # plt.plot([mo.energy for mo in SFO_track(orbs[1].sfos['frag1(HOMO-6)'], orbs[1:-1])], label='HOMO-6')
    # plt.plot([mo.energy for mo in SFO_track(orbs[1].sfos['frag1(HOMO-7)'], orbs[1:-1])], label='HOMO-7')
    plt.xlabel('IRC Step')
    plt.ylabel(r'$\epsilon$ (eV)')
    plt.title('With orbital tracking')
    plt.legend()
    plt.show()
