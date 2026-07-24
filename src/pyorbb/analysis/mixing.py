import pyorbb
import functools
import numpy as np
import itertools as it  # noqa: F401
import matplotlib.pyplot as plt
import os
import networkx as nx  # noqa: F401
import warnings
warnings.filterwarnings('ignore')

ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else list(x)  # noqa: E731

INTERACTION_COLORS = {'OI': '#00FF00', 'PR': '#FF0000', 'Sanitization': '#FF00FF', 'Multiple': '#000000'}


class Mixer2:
    def __init__(self, orbs: pyorbb.Orbitals or str, pr_min_thresh=1e-3, oi_min_thresh=1e-7, oi_max_N=50, pr_max_N=50):
        self.orbs = orbs
        if isinstance(orbs, str):
            self.orbs = pyorbb.Orbitals(str)

        self.fmos = self.orbs.fmos.orbitals
        self.allowed_mos = self.orbs.mos.orbitals
        self.allowed_fmos = self.orbs.fmos.orbitals

        self.data = {}
        self.main_mix = Mixing(self.orbs)
        self.mixes = {'OI': {}, 'PR': {}}

        self.energy_type = 'energy'
        self.set_enable_oi(True)
        self.set_enable_pr(True)
        self.set_allowed_mos(self.orbs.mos.orbitals)
        self.set_allowed_fmos(self.orbs.fmos.orbitals)
        self._prepare()
        self.oi_min_thresh = oi_min_thresh
        self.pr_min_thresh = pr_min_thresh
        self.oi_max_N = oi_max_N
        self.pr_max_N = pr_max_N
        self._get_orbital_interactions()
        self._get_pauli_repulsions()
        self.set_energy_type('energy')

    def find_two_mixing(self, orb1, orb2):
        # self.reset_mixes()
        for mix in self.main_mix.two_mixings:
            if len(mix.fmos) > 2:
                continue
            if len(mix.mos) > 2:
                continue
            if not (orb1 in mix.fmos or orb1 in mix.mos):
                continue
            if not (orb2 in mix.fmos or orb2 in mix.mos):
                continue
            if mix.fraction is None:
                continue
            return mix

    def _prepare(self):
        '''
        Prepare the data used to construct the mixing situations.
        We construct an array ``Ei`` for each FMO energy-type we have in our system.
        The values of the array indicate the strength and type of interaction.
        Positive values indicate destabilizing Pauli repulsions, while 
        negative values indicate stabilizing orbital interactions.
        '''

        # prepare the data we will use
        S = np.array([[fmo1.overlap(fmo2) for fmo1 in self.fmos] for fmo2 in self.fmos])
        o = np.array([fmo.occupation for fmo in self.fmos])
        p = np.array([fmo.gross_population for fmo in self.fmos])

        # get the maximum occupation of an FMO
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_fmos'] else 2

        # make sure all populations are between 0 and max_pop
        p = np.clip(p, 0, max_pop)

        # mangle some data into various matrices
        S2 = S*S  # overlap squared
        dp = (abs(p - o) * abs(p - o).reshape(-1, 1))  # electron gains and losses
        P = p + p.reshape(-1, 1)  # sum of FMO populations
        O = o + o.reshape(-1, 1)  # sum of occupations

        # first max_pop electrons go to the bonding MO
        Pbond = np.clip(P, 0, max_pop)
        # any remaining electrons go to the anti-bonding MO
        Panti = np.clip(P - Pbond, 0, max_pop)

        # the number of electrons involved in pauli repulsion
        Epr = np.clip(O - max_pop, 0, max_pop) * S2
        # remove the upper echelon and diagonal
        # this prevents FMO pair double counting
        # and self-interactions
        Epr = np.tril(Epr, k=-1)

        for energy_type in self.orbs.fmo_energy_types:
            # we calculate the Eoi for each energy type we have available
            e = np.array([getattr(fmo, energy_type) for fmo in self.fmos])
            de = abs(e - e.reshape(-1, 1))  # energy gaps

            # calculate the non-degenerate orbital interaction terms
            Eoi = - (Pbond - Panti) * dp * (S2 / de)
            # for degenerate elements we replace S^2/de with S
            degenerate_mask = np.isclose(de, 0, atol=0.002)
            Eoi[degenerate_mask] = (-Pbond * dp * abs(S))[degenerate_mask]

            # remove upper echelon plus diagonal
            # since the matrix should be symmetric and the diagonal 
            # terms are the self-interactions
            Eoi = np.tril(Eoi, k=-1)
            self.data[energy_type] = (
                Eoi, np.argsort(Eoi, axis=None),
                Epr, np.argsort(-Epr, axis=None),
                )

        self.data['mo_occ'] = np.array([mo.occupied for mo in self.orbs.mos])

    def set_enable_oi(self, val):
        self.enable_oi = val

    def set_enable_pr(self, val):
        self.enable_pr = val

    def set_allowed_mos(self, allowed_mos):
        self.allowed_mos = allowed_mos

    def set_allowed_fmos(self, allowed_fmos):
        self.allowed_fmos = allowed_fmos

    def set_oi_threshold(self, thresh):
        self.oi_threshold = thresh
        self.oi_N = None

    def get_next_oi_threshold(self):
        vals = [-val for val in self.mixes['OI'][self.energy_type].values() if -val <= self.oi_threshold]
        if len(vals) == 0:
            return self.oi_threshold
        return vals[0]

    def get_previous_oi_threshold(self):
        vals = [-val for val in self.mixes['OI'][self.energy_type].values() if -val >= self.oi_threshold]
        if len(vals) == 0:
            return self.oi_threshold
        return vals[-1]
        
    def get_next_pr_threshold(self):
        vals = [val for val in self.mixes['PR'][self.energy_type].values() if val <= self.pr_threshold]
        if len(vals) == 0:
            return self.pr_threshold
        return vals[0]

    def get_previous_pr_threshold(self):
        vals = [val for val in self.mixes['PR'][self.energy_type].values() if val >= self.pr_threshold]
        if len(vals) == 0:
            return self.pr_threshold
        return vals[-1]

    def get_oi_default_threshold(self, fraction=0.7):
        '''
        Get a threshold that makes sure that at least 70% of the OI interactions are included.
        '''
        vals = list(sorted([abs(val) for val in self.mixes['OI'][self.energy_type].values()]))[::-1]
        total = sum(vals)
        fracs = [val/total for val in vals]
        cumsum = np.cumsum(fracs)
        idx = np.where(cumsum >= fraction)[0][0]

        return vals[idx]

    def get_pr_default_threshold(self, fraction=0.1):
        '''
        Get a threshold that makes sure that at least 10% of the PR interactions are included.
        '''
        vals = list(sorted([abs(val) for val in self.mixes['PR'][self.energy_type].values()]))[::-1]
        total = sum(vals)
        fracs = [val/total for val in vals]
        cumsum = np.cumsum(fracs)
        idx = np.where(cumsum >= fraction)[0][0]

        if len(vals) > 2:
            return max(vals[idx], vals[2])
        return vals[idx]

    def get_oi_fraction(self, threshold):
        '''
        Get the fraction of the total OI captured by a specific threshold.
        '''
        vals = list(sorted([abs(val) for val in self.mixes['OI'][self.energy_type].values()]))
        total = sum(vals)
        fracs = [val/total for val in vals if val >= threshold]
        return sum(fracs)

    def get_pr_fraction(self, threshold):
        '''
        Get the fraction of the total OI captured by a specific threshold.
        '''
        vals = list(sorted([abs(val) for val in self.mixes['PR'][self.energy_type].values()]))
        total = sum(vals)
        fracs = [val/total for val in vals if val >= threshold]
        return sum(fracs)

    def set_pr_threshold(self, thresh):
        self.pr_threshold = thresh
        self.pr_N = None

    def set_oi_N(self, N):
        self.oi_threshold = None
        self.oi_N = N

    def set_pr_N(self, N):
        self.pr_threshold = None
        self.pr_N = N

    def set_energy_type(self, typ):
        self.energy_type = typ
        if typ not in self.mixes['OI']:
            self._get_orbital_interactions()
        if typ not in self.mixes['PR']:
            self._get_pauli_repulsions()

    def _get_orbital_interactions(self):
        Eoi, Eoi_order = self.data[self.energy_type][0], self.data[self.energy_type][1]
        self._get_mixes(Eoi, Eoi_order, self.oi_min_thresh, 'OI')

    def _get_pauli_repulsions(self):
        Epr, Epr_order = self.data[self.energy_type][2], self.data[self.energy_type][3]
        self._get_mixes(Epr, Epr_order, self.pr_min_thresh, 'PR')

    def _get_mixes(self, 
                   M, 
                   order, 
                   min_thresh,
                   interaction_type):
        self.mixes[interaction_type][self.energy_type] = {}
        n = 0
        while 1:
            i, j = np.unravel_index(order[n], M.shape)
            v = M[i, j]

            if abs(v) < min_thresh:
                break

            fmo1, fmo2 = self.fmos[i], self.fmos[j]

            # if fmo1 not in self.allowed_fmos or fmo2 not in self.allowed_fmos:
                # continue

            mo1, mo2 = self._get_mos(fmo1, fmo2, interaction_type=interaction_type)
            # if mo1 not in self.allowed_mos or mo2 not in self.allowed_mos:
                # continue

            mix = Mixing(self.orbs, [mo1, mo2], [fmo1, fmo2], fraction=v/np.sum(M), connection_type=interaction_type)
            self.mixes[interaction_type][self.energy_type][mix] = v
            n += 1

            if interaction_type == 'OI':
                if n == self.oi_max_N:
                    break
            else:
                if n == self.pr_max_N:
                    break

    def _get_mos(self, fmo1, fmo2, interaction_type=None):
        fmo1_contr = np.array([fmo1.mulliken_contribution(mo, normalized=True) for mo in self.orbs.mos.orbitals])
        fmo2_contr = np.array([fmo2.mulliken_contribution(mo, normalized=True) for mo in self.orbs.mos.orbitals])

        occ_contrs = abs(fmo1_contr * fmo2_contr) * self.data['mo_occ']
        virt_contrs = abs(fmo2_contr * fmo1_contr) * (1-self.data['mo_occ'])

        # handle orbital interactions
        if interaction_type == 'OI':
            occ_mo = self.orbs.mos.orbitals[argNmax(occ_contrs, 0)]
            virt_mo = self.orbs.mos.orbitals[argNmax(virt_contrs, 0)]
            return occ_mo, virt_mo

        elif interaction_type == 'PR':
            occ_mo1 = self.orbs.mos.orbitals[argNmax(occ_contrs, 0)]
            occ_mo2 = self.orbs.mos.orbitals[argNmax(occ_contrs, 1)]
            return occ_mo1, occ_mo2


    def reset_mixes(self):
        self.main_mix = Mixing(self.orbs, energy_type=self.energy_type)
        if self.enable_oi:
            for mix, strength in self.mixes['OI'][self.energy_type].items():
                if not any(mo in self.allowed_mos for mo in mix.mos):
                    continue
                if any(fmo not in self.allowed_fmos for fmo in mix.fmos):
                    continue
                if self.oi_threshold is not None and abs(strength) >= self.oi_threshold:
                    self.main_mix += mix

        if self.enable_pr:
            for mix, strength in self.mixes['PR'][self.energy_type].items():
                if not any(mo in self.allowed_mos for mo in mix.mos):
                    continue
                if any(fmo not in self.allowed_fmos for fmo in mix.fmos):
                    continue
                if self.pr_threshold is not None and abs(strength) >= self.pr_threshold:
                    self.main_mix += mix

        self.sanitize()


    def sanitize(self):
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_fmos'] else 2

        for spin in ['A', 'B']:
            for symm in self.orbs.mos.symmetry:
                relevant_mos = [mo for mo in self.orbs.mos if mo.spin in (spin, 'AB') and mo.symmetry == symm]
                relevant_mix_mos = [mo for mo in relevant_mos if mo in self.main_mix.mos]
                N_virt_MO = len([mo for mo in relevant_mix_mos if mo.occupation == 0])
                N_occ_MO = len([mo for mo in relevant_mix_mos if mo.occupation > 0])

                # allowed_fmos = [fmo for fmo in self.orbs.fmos]
                relevant_fmos = [fmo for fmo in self.orbs.fmos if fmo.spin in (spin, 'AB') and fmo.symmetry == symm]
                relevant_mix_fmos = [fmo for fmo in relevant_fmos if fmo in self.main_mix.fmos]
                N_elec_FMO = round(sum([fmo.gross_population for fmo in relevant_mix_fmos]))
                N_virt_FMO = len(relevant_mix_fmos) - N_elec_FMO/max_pop
                N_occ_FMO = len(relevant_mix_fmos) - N_virt_FMO

                # checking some requirements
                missing_occ_MOs = N_occ_MO < N_occ_FMO
                missing_occ_FMOs = N_occ_FMO < N_occ_MO
                missing_virt_MOs = N_virt_MO < N_virt_FMO
                missing_virt_FMOs = N_virt_FMO < N_virt_MO

                if not any([missing_occ_MOs, missing_occ_FMOs, missing_virt_MOs, missing_virt_FMOs]):
                    continue

                ## GENERATE CANDIDATE MOs AND FMOs
                candidate_occ_mos = {}
                candidate_virt_mos = {}
                for mo in relevant_mos:
                    if mo in self.main_mix.mos:
                        continue

                    # skip if we don't need the occupied MOs
                    if mo.occupied and not missing_occ_MOs:
                        continue

                    # same for virtual
                    if not mo.occupied and not missing_virt_MOs:
                        continue

                    highest = 0
                    for i in range(len(relevant_mix_fmos)):
                        C1 = relevant_mix_fmos[i].mulliken_contribution(mo, normalized=True)
                        for j in range(i+1, len(relevant_mix_fmos)):
                            C2 = relevant_mix_fmos[j].mulliken_contribution(mo, normalized=True)
                            if abs(C1 * C2) > highest:
                                highest = abs(C1*C2)

                    if mo.occupied:
                        candidate_occ_mos[mo] = highest
                    else:
                        candidate_virt_mos[mo] = highest

                candidate_occ_mos = sorted(candidate_occ_mos.items(), key=lambda r: -r[1])
                candidate_virt_mos = sorted(candidate_virt_mos.items(), key=lambda r: -r[1])


                candidate_occ_fmos = {}
                candidate_virt_fmos = {}
                for fmo in relevant_fmos:
                    if fmo in relevant_mix_fmos:
                        continue

                    # skip if we don't need the occupied FMOs
                    if fmo.occupied and not missing_occ_FMOs:
                        continue

                    # same for virtual
                    if not fmo.occupied and not missing_virt_FMOs:
                        continue

                    highest = 0
                    for i in range(len(relevant_mix_mos)):
                        C1 = fmo.mulliken_contribution(relevant_mix_mos[i], normalized=True)
                        for j in range(i+1, len(relevant_mix_mos)):
                            C2 = fmo.mulliken_contribution(relevant_mix_mos[j], normalized=True)
                            if abs(C1 * C2) > highest:
                                highest = abs(C1*C2)

                    if fmo.occupied:
                        candidate_occ_fmos[fmo] = highest
                    else:
                        candidate_virt_fmos[fmo] = highest

                candidate_occ_fmos = sorted(candidate_occ_fmos.items(), key=lambda r: -r[1])
                candidate_virt_fmos = sorted(candidate_virt_fmos.items(), key=lambda r: -r[1])


                # ADD MOs and FMOs BASED ON UNMET REQUIREMENTS
                if missing_occ_MOs:
                    N_occ_MO_missing = N_occ_FMO - N_occ_MO
                    for i in range(int(N_occ_MO_missing)):
                        self.main_mix.add_mo(candidate_occ_mos[i][0])

                if missing_occ_FMOs:
                    N_occ_FMO_missing = N_occ_MO - N_occ_FMO
                    for i in range(int(N_occ_FMO_missing)):
                        self.main_mix.add_fmo(candidate_occ_fmos[i][0])

                if missing_virt_MOs:
                    N_virt_MO_missing = N_virt_FMO - N_virt_MO
                    for i in range(int(N_virt_MO_missing)):
                        self.main_mix.add_mo(candidate_virt_mos[i][0])

                if missing_virt_FMOs:
                    N_virt_FMO_missing = N_virt_MO - N_virt_FMO
                    for i in range(int(N_virt_FMO_missing)):
                        self.main_mix.add_fmo(candidate_virt_fmos[i][0])

        [delattr(fmo, '_display_occupation') for fmo in self.main_mix.fmos if hasattr(fmo, '_display_occupation')]
        for mix in self.split():
            # check if there are too many electrons in the mos compared to the fmos
            total_mo_occ = sum(mo.occupation for mo in mix.mos)
            total_fmo_occ = sum(fmo.occupation for fmo in mix.fmos)
            diff = total_mo_occ - total_fmo_occ
            # if there are more electrons in the mos we add an electron to the FMOs
            if diff > 0:
                while diff > 0:
                    max_fmo = max([fmo for fmo in mix.fmos if not hasattr(fmo, '_display_occupation')], key=lambda fmo: fmo.gross_population)
                    max_fmo._display_occupation = min(diff, max_pop)
                    diff -= max_fmo._display_occupation
            else:
                while diff < 0:
                    min_fmo = min([fmo for fmo in mix.fmos if not hasattr(fmo, '_display_occupation')], key=lambda fmo: fmo.gross_population)
                    min_fmo._display_occupation = min_fmo.occupation - min(min_fmo.occupation, abs(diff))
                    # print(min_fmo, min_fmo._display_occupation, to_remove, diff)
                    diff += min(min_fmo.occupation, abs(diff))

    def split(self):
        return self.main_mix.split()

    def draw_diagram(self, *args, **kwargs):
        return self.main_mix.draw_diagram(*args, **kwargs)

    @property
    def connections(self):
        return self.main_mix.connections


class Mixer:
    '''
    The main class responsible for generating 2-mixing situations.

    Args:
        orbs: 
    '''
    def __init__(self, orbs: pyorbb.Orbitals, energy_type: str = 'energy'):
        self.orbs = orbs
        self.energy_type = energy_type
        self._prepare()

    def _prepare(self):
        self.fmos = {}
        self.fmos_occ = {}
        self.fmos_vir = {}
        self.fmos_energy = {}
        for i, frag in enumerate(self.orbs.fmos.fragments):
            self.fmos[frag] = [fmo for fmo in self.orbs.fmos if fmo.fragment_unique == frag]
            self.fmos_occ[frag] = np.array([fmo.occupation > 0 for fmo in self.fmos[frag]]).reshape(-1, 1)
            if not self.orbs.fmos.unrestricted:
                self.fmos_vir[frag] = np.array([fmo.occupation < 2 for fmo in self.fmos[frag]]).reshape(-1, 1)
            else:
                self.fmos_vir[frag] = np.array([fmo.occupation < 1 for fmo in self.fmos[frag]]).reshape(-1, 1)

            self.fmos_energy[frag] = np.array([getattr(fmo, self.energy_type) for fmo in self.fmos[frag]]).reshape(-1, 1)

        self.mos = list(self.orbs.mos)
        self.mo_occ = np.array([mo.occupied for mo in self.mos])

        self.S_oi = {}
        self.dE_oi = {}
        self.oi = {}
        self.oi_approx_total = 0

        self.S_pauli = {}
        self.pauli = {}
        self.pauli_approx_total = 0
        for i, frag in enumerate(self.orbs.fmos.fragments):
            for frag2 in self.orbs.fmos.fragments[i+1:]:
                # get data for oi
                occ_virt_mask = np.logical_or(np.logical_and(self.fmos_occ[frag], self.fmos_vir[frag2].T), np.logical_and(self.fmos_vir[frag], self.fmos_occ[frag2].T))

                self.S_oi[(frag, frag2)] = overlap_mat(self.fmos[frag], self.fmos[frag2])
                self.dE_oi[(frag, frag2)] = abs(self.fmos_energy[frag] - self.fmos_energy[frag2].T)
                nogap_mask = self.dE_oi[(frag, frag2)] != 0
                self.dE_oi[(frag, frag2)] += (1 - nogap_mask)

                self.oi[(frag, frag2)] = -self.S_oi[(frag, frag2)]**2 / self.dE_oi[(frag, frag2)] * occ_virt_mask
                self.oi_approx_total += self.oi[(frag, frag2)].sum()

                # get data for pauli
                occ_occ_mask = np.logical_and(self.fmos_occ[frag], self.fmos_occ[frag2].T)
                self.S_pauli[(frag, frag2)] = overlap_mat(self.fmos[frag], self.fmos[frag2])
                self.pauli[(frag, frag2)] = self.S_pauli[(frag, frag2)]**2 * occ_occ_mask
                self.pauli_approx_total += self.pauli[(frag, frag2)].sum()

        self.oi_ref = self.orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503
        self.pauli_ref = self.orbs.reader.read('Energy', 'Pauli Total') * 627.503


    def orbital_interactions(self, fraction_thresh=None, N=None):
        """
        Yield the first ``N`` strongest orbital interactions.
        """
        if fraction_thresh is None and N is None:
            fraction_thresh = .03
        elif N is not None:
            fraction_thresh = None
        ret = []
        for i, frag in enumerate(self.orbs.fmos.fragments):
            for frag2 in self.orbs.fmos.fragments[i+1:]:
                j = 0
                while 1:
                    best = np.unravel_index(argNmax(-self.oi[(frag, frag2)], j), self.oi[(frag, frag2)].shape)
                    best_oi = self.oi[(frag, frag2)][best]
                    fmo1, fmo2  = self.fmos[frag][best[0]], self.fmos[frag2][best[1]]

                    fmo1_contr = np.array([fmo1.mulliken_contribution(mo) for mo in self.mos])
                    fmo2_contr = np.array([fmo2.mulliken_contribution(mo) for mo in self.mos])

                    occ_contrs = abs(fmo1_contr * fmo2_contr) * self.mo_occ
                    virt_contrs = abs(fmo2_contr * fmo1_contr) * (1-self.mo_occ)

                    occ_mo = self.mos[np.argmax(occ_contrs)]
                    virt_mo = self.mos[np.argmax(virt_contrs)]
                    frac = best_oi / self.oi_approx_total
                    if fraction_thresh is None and j == N:
                        break
                    elif fraction_thresh is not None and frac < fraction_thresh:
                        break
                    stab = frac * self.oi_ref

                    typ = {conn: 'OI' for conn in list(it.product([fmo1, fmo2], [occ_mo, virt_mo]))}
                    mix = Mixing(self.orbs, [occ_mo, virt_mo], [fmo1, fmo2], stab, frac, connection_type=typ)
                    ret.append(mix)
                    j += 1

        ret = sorted(ret, key=lambda row: row.strength)
        if N is not None:
            ret = ret[:N]
        return ret


    def pauli_repulsions(self, fraction_thresh=None, N=None):
        """
        Yield the first ``N`` strongest Pauli repulsion interactions.
        """
        if fraction_thresh is None and N is None:
            fraction_thresh = .03
        elif N is not None:
            fraction_thresh = None
        ret = []
        for i, frag in enumerate(self.orbs.fmos.fragments):
            for frag2 in self.orbs.fmos.fragments[i+1:]:
                j = 0
                while 1:
                    best = np.unravel_index(argNmax(self.pauli[(frag, frag2)], j), self.pauli[(frag, frag2)].shape)
                    best_pauli = self.pauli[(frag, frag2)][best]
                    fmo1, fmo2  = self.fmos[frag][best[0]], self.fmos[frag2][best[1]]

                    fmo1_contr = np.array([fmo1.mulliken_contribution(mo) * mo.occupation for mo in self.mos])
                    fmo2_contr = np.array([fmo2.mulliken_contribution(mo) * mo.occupation for mo in self.mos])

                    occ_contrs = abs(fmo1_contr * fmo2_contr)

                    occ_mo1_idx = argNmax(occ_contrs, 0)
                    occ_mo2_idx = argNmax(occ_contrs, 1)

                    occ_mo1 = self.mos[occ_mo1_idx]
                    occ_mo2 = self.mos[occ_mo2_idx]

                    frac = best_pauli / self.pauli_approx_total
                    if fraction_thresh is None and j == N:
                        break
                    elif fraction_thresh is not None and frac < fraction_thresh:
                        break
                    stab = frac * self.pauli_ref

                    typ = {conn: 'PR' for conn in list(it.product([fmo1, fmo2], [occ_mo1, occ_mo2]))}
                    mix = Mixing(self.orbs, [occ_mo1, occ_mo2], [fmo1, fmo2], stab, frac, connection_type=typ)
                    ret.append(mix)
                    j += 1

        ret = sorted(ret, key=lambda row: -row.strength)
        if N is not None:
            ret = ret[:N]
        return ret


class Mixing:
    def __init__(self, orbs, mos=None, fmos=None, strength=None, fraction=None, connections=None, connection_type=None, energy_type='energy'):
        self.orbs = orbs
        self.mos = mos or []
        self.fmos = fmos or []
        self.strength = strength or 0
        self.fraction = fraction or 0
        self.energy_type = energy_type
        self.connections = connections
        self.connection_type = connection_type
        if connections is None:
            self.connections = list(it.product(self.fmos, self.mos))
        if connection_type is None:
            self.connection_type = {conn: 'Multiple' for conn in self.connections}
        if isinstance(connection_type, str):
            self.connection_type = {conn: connection_type for conn in self.connections}

        self.two_mixings = [self]

    def __str__(self):
        s = f'{self.__class__.__name__}('
        s += f'[{", ".join([fmo.make_name(frag_name=True, relative_name=True, spin=True) for fmo in self.fmos])}]'
        s += ' -> '
        s += f'[{", ".join([mo.relative_name for mo in self.mos])}]'
        if self.strength:
            s += f', strength={self.strength:+.1f} kcal/mol'
        if self.fraction:
            s += f', fraction={self.fraction:.1%}' 
        s += f', nelectrons={self.nelectrons()})'
        return s

    @property
    def PR_is_empty(self):
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_fmos'] else 2
        for mix in self.two_mixings:
            if mix is self:
                continue
            if mix.nelectrons() == 2 * max_pop:
                return False
        return True

    @property
    def OI_is_empty(self):
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_fmos'] else 2
        for mix in self.two_mixings:
            if mix is self:
                continue
            if mix.nelectrons() == max_pop:
                return False
        return True

    def add_mo(self, mo, typ='Sanitization', connections=None):
        '''
        Add an MO to this mixing diagram. 
        '''
        self.mos.append(mo)
        if connections is None:
            for fmo in self.fmos:
                if abs(fmo.mulliken_contribution(mo)) > 0.03:
                    self.connections.append((fmo, mo))
                    self.connection_type[(fmo, mo)] = typ
        else:
            self.connections.append(connections)
            for conn in self.connections:
                self.connection_type[conn] = typ


    def add_fmo(self, fmo, typ='Sanitization', connections=None):
        '''
        Add an FMO to this mixing diagram. 
        '''
        self.fmos.append(fmo)
        if connections is None:
            for mo in self.mos:
                if abs(fmo.mulliken_contribution(mo)) > 0.03:
                    self.connections.append((fmo, mo))
                    self.connection_type[(fmo, mo)] = typ
        else:
            self.connections.append(connections)
            for conn in self.connections:
                self.connection_type[conn] = typ


    def draw_diagram(self, ax=None, ylim=None, simple=False, **kwargs):
        if simple:
            pyorbb.plotting.simple_orbital_diagram.draw_interaction(self.fmos, self.mos, self.connections, None, energy_type=self.energy_type, connection_types=self.connection_type, ax=ax, ylim=ylim)
        else:
            pyorbb.plotting.orbital_diagram.draw_interaction(self.fmos, self.mos, self.connections, None, energy_type=self.energy_type, connection_types=self.connection_type, ax=ax, ylim=ylim, **kwargs)

    def draw_fmos(self, overlap=False, screen=None):
        import tcviewer  # noqa: F811

        if screen is None:
            scr = tcviewer.Screen()
            scr.__enter__()
        else:
            scr = screen

        for fmo in self.fmos:
            with scr.add_molscene() as scene:
                cub = fmo.cube_file()
                scene.draw_text(str(fmo))

                if fmo.occupied:
                    colors = ([1, 0, 0], [0, 0, 1])
                else:
                    colors = ([0, 1, 1], [1, .5, 0])

                scene.draw_isosurface(cub, -.03, opacity=.25, color=colors[0])
                scene.draw_isosurface(cub,  .03, opacity=.25, color=colors[1])
                scene.draw_molecule(fmo.molecule)

        if screen is None:
            scr.__exit__()

    def screenshot_fmos(self, outdir='FMO_pictures'):
        import tcviewer  # noqa: F811

        os.makedirs(outdir, exist_ok=True)
        with tcviewer.Screen(headless=True) as scr:
            for fmo in self.fmos:
                with scr.add_molscene() as scene:
                    cub = fmo.cube_file()
                    scene.draw_molecule(fmo.molecule)

                    if fmo.occupied:
                        colors = ([1, 0, 0], [0, 0, 1])
                    else:
                        colors = ([0, 1, 1], [1, .5, 0])
                
                    scene.draw_isosurface(cub, -.03, opacity=.5, color=colors[0])
                    scene.draw_isosurface(cub,  .03, opacity=.5, color=colors[1])
                    scene.screenshot(os.path.join(outdir, f'{str(fmo)}.png'))

    def nelectrons(self):
        return sum(mo.occupation for mo in self.mos)

    def __add__(self, other: 'Mixing'):
        self.fmos.extend([ofmo for ofmo in other.fmos if ofmo not in self.fmos])
        self.mos.extend([omo for omo in other.mos if omo not in self.mos])
        self.connections.extend([oconn for oconn in other.connections if oconn not in self.connections])
        for conn, typ in other.connection_type.items():
            if conn in self.connection_type:
                if typ == self.connection_type[conn]:
                    continue
                else:
                    self.connection_type[conn] = 'Multiple'
            else:
                self.connection_type[conn] = typ

        self.strength = None
        self.fraction = None
        self.two_mixings.extend(other.two_mixings)
        return self

    def fits(self, other: 'Mixing'):
        if len(self.fmos) == 0 and len(self.mos) == 0:
            return True
        if not any(fmo in other.fmos for fmo in self.fmos):
            return False
        if not any(mo in other.mos for mo in self.mos):
            return False
        return True

    def xiaobo_value(self):
        if all(fmo.occupied for fmo in self.fmos):
            maxval = 0
            for fmo1 in self.fmos:
                for fmo2 in self.fmos:
                    if fmo1 == fmo2:
                        continue
                    maxval = max(maxval, abs(fmo1.overlap(fmo2)))
            return maxval

        val = (self.fmos[0] @ self.fmos[1])**2 / abs(self.fmos[0].energy - self.fmos[1].energy)

        for fmo in self.fmos:
            gp = fmo.gross_population
            gp = np.clip(gp, 0, 2)
            excess = abs(fmo.occupation - gp)
            val *= excess

        return val

    def xiaobo_check(self, threshold=0.01):
        if all(fmo.occupied for fmo in self.fmos):
            for fmo1 in self.fmos:
                for fmo2 in self.fmos:
                    if fmo1 == fmo2:
                        continue
                    if abs(fmo1.overlap(fmo2)) < threshold:
                        return False
            return True

        val = (self.fmos[0] @ self.fmos[1])**2 / abs(self.fmos[0].energy - self.fmos[1].energy)

        for fmo in self.fmos:
            gp = fmo.gross_population
            gp = np.clip(gp, 0, 2)
            excess = abs(fmo.occupation - gp)
            val *= excess

        return val > threshold

    @property
    def fragments(self):
        return set(fmo.fragment_unique for fmo in self.fmos)

    @property
    def lowest_contribution(self):
        return min([abs(fmo.mulliken_contribution(mo)) for fmo in self.fmos for mo in self.mos])

    @property
    def irrep(self):
        return self.mos[0].symmetry

    @property
    def spin(self):
        return self.mos[0].spin

    def split(self):
        G = nx.Graph()
        G.add_edges_from(self.connections)
        subGs = [G.subgraph(c) for c in nx.connected_components(G)]
        mixes = []
        for subG in subGs:
            orbs_ = subG.nodes()
            mos = [orb for orb in orbs_ if isinstance(orb, pyorbb.orbitals.objects.MO)]
            fmos = [orb for orb in orbs_ if isinstance(orb, pyorbb.orbitals.objects.FMO)]
            connections = subG.edges()
            connections = [conn[::-1] if isinstance(conn[0], pyorbb.orbitals.objects.MO) else conn for conn in connections]
            mixes.append(Mixing(
                self.orbs,
                mos=mos, 
                fmos=fmos, 
                connections=connections,
                connection_type={conn: self.connection_type[conn] for conn in connections},
                energy_type=self.energy_type))
        return mixes

    def find_closed_interactions(self, orb):
        G = nx.Graph()
        G.add_edges_from(self.connections)
        cycles = [cycle for cycle in nx.algorithms.cycles.simple_cycles(G, length_bound=4) if orb in cycle]
        return cycles
        ret = []
        for two_mixing in self.two_mixings:
            if orb in two_mixing.fmos or orb in two_mixing.mos:
                ret.append([*two_mixing.fmos, *two_mixing.mos])
        return ret


def overlap_mat(fmos1, fmos2):
    ret = []
    for fmo1 in ensure_list(fmos1):
        ret.append([])
        for fmo2 in ensure_list(fmos2):
            ret[-1].append(abs(fmo1 @ fmo2))
    return np.array(ret).squeeze()


def argNmax(arr, N):
    """
    Get the Nth maximum element of an array.
    """
    return np.argsort(-arr, axis=None)[N]


def track_mixing(orbss, mixing):
    out_dir = 'pyfrag_imgs'
    os.makedirs(out_dir, exist_ok=True)
    for i, orbs in enumerate(orbss):
        fmos = [orbs.fmos[str(mix_fmo)] for mix_fmo in mixing.fmos]
        mos = [orbs.mos[str(mix_mo)] for mix_mo in mixing.mos]
        pyorbb.plotting.orbital_diagram.draw_interaction(fmos, mos, it.product(fmos, mos))

        plt.savefig(os.path.join(out_dir, f'{i}.jpg'))
        plt.close()


def _is_bonding(fmo1, fmo2, mo):
    return round(fmo1.coefficient(mo) * fmo2.coefficient(mo) * (fmo1 @ fmo2), 4) >= 0
