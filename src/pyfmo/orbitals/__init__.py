from pyfmo.orbitals import sfo, mo, adf, dftb  # noqa
from scm import plams
import numpy as np
# import matplotlib.pyplot as plt
# from yutility import plot, ensure_list
from tcutility import ensure_list, cache


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

    @cache.cache
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

    @cache.cache
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

    def write_excel(self, out_file: str = 'pyfmo.xlsx'):
        from pyfmo import write_excel
        
        write_excel.to_excel(self, out_file)


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
    