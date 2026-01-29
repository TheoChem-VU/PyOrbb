import pyfmo
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
    def __init__(self, orbs: pyfmo.Orbitals or str, pr_min_thresh=1e-3, oi_min_thresh=1e-7):
        self.orbs = orbs
        if isinstance(orbs, str):
            self.orbs = pyfmo.Orbitals(str)

        self.sfos = self.orbs.sfos.orbitals
        self.allowed_mos = self.orbs.mos.orbitals
        self.allowed_sfos = self.orbs.sfos.orbitals

        self.data = {}
        self.main_mix = Mixing(self.orbs)
        self.mixes = {'OI': {}, 'PR': {}}

        self.energy_type = 'energy'
        self.set_enable_oi(True)
        self.set_enable_pr(True)
        self.set_allowed_mos(self.orbs.mos.orbitals)
        self.set_allowed_sfos(self.orbs.sfos.orbitals)
        self._prepare()
        self.oi_min_thresh = oi_min_thresh
        self.pr_min_thresh = pr_min_thresh
        self._get_orbital_interactions()
        self._get_pauli_repulsions()
        self.set_energy_type('energy')

    def _prepare(self):
        '''
        Prepare the data used to construct the mixing situations.
        We construct an array ``Ei`` for each SFO energy-type we have in our system.
        The values of the array indicate the strength and type of interaction.
        Positive values indicate destabilizing Pauli repulsions, while 
        negative values indicate stabilizing orbital interactions.
        '''

        # prepare the data we will use
        S = np.array([[sfo1.overlap(sfo2) for sfo1 in self.sfos] for sfo2 in self.sfos])
        o = np.array([sfo.occupation for sfo in self.sfos])
        p = np.array([sfo.gross_population for sfo in self.sfos])

        # get the maximum occupation of an SFO
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_sfos'] else 2

        # make sure all populations are between 0 and max_pop
        p = np.clip(p, 0, max_pop)

        # mangle some data into various matrices
        S2 = S*S  # overlap squared
        dp = abs(p - o) * abs(p - o).reshape(-1, 1)  # electron gains and losses
        P = p + p.reshape(-1, 1)  # sum of SFO populations
        O = o + o.reshape(-1, 1)  # sum of occupations

        # first max_pop electrons go to the bonding MO
        Noi = np.clip(P, 0, max_pop)
        # any remaining electrons go to the anti-bonding MO
        Npr = np.clip(P - Noi, 0, max_pop)

        # the number of electrons involved in pauli repulsion
        Epr = np.clip(O - max_pop, 0, max_pop) * S2
        # remove the upper echelon and diagonal
        # this prevents SFO pair double counting
        # and self-interactions
        Epr = np.tril(Epr, k=-1)

        for energy_type in self.orbs.sfo_energy_types:
            # we calculate the Eoi for each energy type we have available
            e = np.array([getattr(sfo, energy_type) for sfo in self.sfos])
            de = abs(e - e.reshape(-1, 1))  # energy gap

            # calculate the non-degenerate orbital interaction terms
            Eoi = -(Noi - Npr) * dp * S2 / de
            # for degenerate elements we replace S^2/de with S
            degenerate_mask = np.isclose(de, 0, atol=0.002)
            Eoi[degenerate_mask] = (-Noi * dp * abs(S))[degenerate_mask]

            # remove upper echelon plus diagonal
            # since the matrix should be symmetric and the diagonal 
            # terms are the self-interactions
            Eoi = np.tril(Eoi, k=-1)
            self.data[energy_type] = (
                Eoi, np.argsort(Eoi, axis=None), 
                Epr, np.argsort(-Epr, axis=None)
                )

        self.data['mo_occ'] = np.array([mo.occupied for mo in self.orbs.mos])

    def set_enable_oi(self, val):
        self.enable_oi = val

    def set_enable_pr(self, val):
        self.enable_pr = val

    def set_allowed_mos(self, allowed_mos):
        self.allowed_mos = allowed_mos

    def set_allowed_sfos(self, allowed_sfos):
        self.allowed_sfos = allowed_sfos

    def set_oi_threshold(self, thresh):
        self.oi_threshold = thresh
        self.oi_N = None

    def get_next_oi_threshold(self):
        vals = [-val for val in self.mixes['OI'][self.energy_type].values() if -val < self.oi_threshold]
        if len(vals) == 0:
            return self.oi_threshold
        return vals[0]

    def get_previous_oi_threshold(self):
        vals = [-val for val in self.mixes['OI'][self.energy_type].values() if -val > self.oi_threshold]
        if len(vals) == 0:
            return self.oi_threshold
        return vals[-1]
        
    def get_next_pr_threshold(self):
        vals = [val for val in self.mixes['PR'][self.energy_type].values() if val < self.pr_threshold]
        if len(vals) == 0:
            return self.pr_threshold
        return vals[0]

    def get_previous_pr_threshold(self):
        vals = [val for val in self.mixes['PR'][self.energy_type].values() if val > self.pr_threshold]
        if len(vals) == 0:
            return self.pr_threshold
        return vals[-1]

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
        v = float('inf')

        self.mixes[interaction_type][self.energy_type] = {}
        n = 0
        while 1:
            i, j = np.unravel_index(order[n], M.shape)
            v = M[i, j]

            if abs(v) < min_thresh:
                break

            sfo1, sfo2 = self.sfos[i], self.sfos[j]

            if sfo1 not in self.allowed_sfos or sfo2 not in self.allowed_sfos:
                continue

            mo1, mo2 = self._get_mos(sfo1, sfo2, interaction_type=interaction_type)
            if mo1 not in self.allowed_mos or mo2 not in self.allowed_mos:
                continue

            mix = Mixing(self.orbs, [mo1, mo2], [sfo1, sfo2], connection_type=interaction_type)
            self.mixes[interaction_type][self.energy_type][mix] = v
            n += 1


    def _get_mos(self, sfo1, sfo2, interaction_type=None):
        sfo1_contr = np.array([sfo1.mulliken_contribution(mo, normalized=True) for mo in self.orbs.mos.orbitals])
        sfo2_contr = np.array([sfo2.mulliken_contribution(mo, normalized=True) for mo in self.orbs.mos.orbitals])

        occ_contrs = abs(sfo1_contr * sfo2_contr) * self.data['mo_occ']
        virt_contrs = abs(sfo2_contr * sfo1_contr) * (1-self.data['mo_occ'])

        max_occ = 1 if self.orbs.data['calc_info']['unrestricted_mos'] else 2
        occ_total = sfo1.occupation + sfo2.occupation

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
                if any(sfo not in self.allowed_sfos for sfo in mix.sfos):
                    continue
                if self.oi_threshold is not None and abs(strength) >= self.oi_threshold:
                    self.main_mix += mix

        if self.enable_pr:
            for mix, strength in self.mixes['PR'][self.energy_type].items():
                if not any(mo in self.allowed_mos for mo in mix.mos):
                    continue
                if any(sfo not in self.allowed_sfos for sfo in mix.sfos):
                    continue
                if self.pr_threshold is not None and abs(strength) >= self.pr_threshold:
                    self.main_mix += mix

        self.sanitize()


    def sanitize(self):
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_sfos'] else 2

        for spin in ['A', 'B']:
            for symm in self.orbs.mos.symmetry:
                relevant_mos = [mo for mo in self.orbs.mos if mo.spin in (spin, 'AB') and mo.symmetry == symm]
                relevant_mix_mos = [mo for mo in relevant_mos if mo in self.main_mix.mos]
                N_virt_MO = len([mo for mo in relevant_mix_mos if mo.occupation == 0])
                N_occ_MO = len([mo for mo in relevant_mix_mos if mo.occupation == max_pop])

                # allowed_sfos = [sfo for sfo in self.orbs.sfos]
                relevant_sfos = [sfo for sfo in self.orbs.sfos if sfo.spin in (spin, 'AB') and sfo.symmetry == symm]
                relevant_mix_sfos = [sfo for sfo in relevant_sfos if sfo in self.main_mix.sfos]
                N_elec_SFO = round(sum([sfo.gross_population for sfo in relevant_mix_sfos]))
                N_virt_SFO = len(relevant_mix_sfos) - N_elec_SFO/max_pop
                N_occ_SFO = len(relevant_mix_sfos) - N_virt_SFO

                # checking some requirements
                missing_occ_MOs = N_occ_MO < N_occ_SFO
                missing_occ_SFOs = N_occ_SFO < N_occ_MO
                missing_virt_MOs = N_virt_MO < N_virt_SFO
                missing_virt_SFOs = N_virt_SFO < N_virt_MO

                if not any([missing_occ_MOs, missing_occ_SFOs, missing_virt_MOs, missing_virt_SFOs]):
                    continue

                ## GENERATE CANDIDATE MOs AND SFOs
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
                    for i in range(len(relevant_mix_sfos)):
                        C1 = relevant_mix_sfos[i].mulliken_contribution(mo, normalized=True)
                        for j in range(i+1, len(relevant_mix_sfos)):
                            C2 = relevant_mix_sfos[j].mulliken_contribution(mo, normalized=True)
                            if abs(C1 * C2) > highest:
                                highest = abs(C1*C2)

                    if mo.occupied:
                        candidate_occ_mos[mo] = highest
                    else:
                        candidate_virt_mos[mo] = highest

                candidate_occ_mos = sorted(candidate_occ_mos.items(), key=lambda r: -r[1])
                candidate_virt_mos = sorted(candidate_virt_mos.items(), key=lambda r: -r[1])


                candidate_occ_sfos = {}
                candidate_virt_sfos = {}
                for sfo in relevant_sfos:
                    if sfo in relevant_mix_sfos:
                        continue

                    # skip if we don't need the occupied SFOs
                    if sfo.occupied and not missing_occ_SFOs:
                        continue

                    # same for virtual
                    if not sfo.occupied and not missing_virt_SFOs:
                        continue

                    highest = 0
                    for i in range(len(relevant_mix_mos)):
                        C1 = sfo.mulliken_contribution(relevant_mix_mos[i], normalized=True)
                        for j in range(i+1, len(relevant_mix_mos)):
                            C2 = sfo.mulliken_contribution(relevant_mix_mos[j], normalized=True)
                            if abs(C1 * C2) > highest:
                                highest = abs(C1*C2)

                    if sfo.occupied:
                        candidate_occ_sfos[sfo] = highest
                    else:
                        candidate_virt_sfos[sfo] = highest

                candidate_occ_sfos = sorted(candidate_occ_sfos.items(), key=lambda r: -r[1])
                candidate_virt_sfos = sorted(candidate_virt_sfos.items(), key=lambda r: -r[1])


                # ADD MOs and SFOs BASED ON UNMET REQUIREMENTS
                if missing_occ_MOs:
                    N_occ_MO_missing = N_occ_SFO - N_occ_MO
                    for i in range(int(N_occ_MO_missing)):
                        self.main_mix.add_mo(candidate_occ_mos[i][0])

                if missing_occ_SFOs:
                    N_occ_SFO_missing = N_occ_MO - N_occ_SFO
                    for i in range(int(N_occ_SFO_missing)):
                        self.main_mix.add_sfo(candidate_occ_sfos[i][0])

                if missing_virt_MOs:
                    N_virt_MO_missing = N_virt_SFO - N_virt_MO
                    for i in range(int(N_virt_MO_missing)):
                        self.main_mix.add_mo(candidate_virt_mos[i][0])

                if missing_virt_SFOs:
                    N_virt_SFO_missing = N_virt_MO - N_virt_SFO
                    for i in range(int(N_virt_SFO_missing)):
                        self.main_mix.add_sfo(candidate_virt_sfos[i][0])


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
    def __init__(self, orbs: pyfmo.Orbitals, energy_type: str = 'energy'):
        self.orbs = orbs
        self.energy_type = energy_type
        self._prepare()

    def _prepare(self):
        self.sfos = {}
        self.sfos_occ = {}
        self.sfos_vir = {}
        self.sfos_energy = {}
        for i, frag in enumerate(self.orbs.sfos.fragments):
            self.sfos[frag] = [sfo for sfo in self.orbs.sfos if sfo.fragment_unique == frag]
            self.sfos_occ[frag] = np.array([sfo.occupation > 0 for sfo in self.sfos[frag]]).reshape(-1, 1)
            if not self.orbs.sfos.unrestricted:
                self.sfos_vir[frag] = np.array([sfo.occupation < 2 for sfo in self.sfos[frag]]).reshape(-1, 1)
            else:
                self.sfos_vir[frag] = np.array([sfo.occupation < 1 for sfo in self.sfos[frag]]).reshape(-1, 1)

            self.sfos_energy[frag] = np.array([getattr(sfo, self.energy_type) for sfo in self.sfos[frag]]).reshape(-1, 1)

        self.mos = list(self.orbs.mos)
        self.mo_occ = np.array([mo.occupied for mo in self.mos])

        self.S_oi = {}
        self.dE_oi = {}
        self.oi = {}
        self.oi_approx_total = 0

        self.S_pauli = {}
        self.pauli = {}
        self.pauli_approx_total = 0
        for i, frag in enumerate(self.orbs.sfos.fragments):
            for frag2 in self.orbs.sfos.fragments[i+1:]:
                # get data for oi
                occ_virt_mask = np.logical_or(np.logical_and(self.sfos_occ[frag], self.sfos_vir[frag2].T), np.logical_and(self.sfos_vir[frag], self.sfos_occ[frag2].T))

                self.S_oi[(frag, frag2)] = overlap_mat(self.sfos[frag], self.sfos[frag2])
                self.dE_oi[(frag, frag2)] = abs(self.sfos_energy[frag] - self.sfos_energy[frag2].T)
                nogap_mask = self.dE_oi[(frag, frag2)] != 0
                self.dE_oi[(frag, frag2)] += (1 - nogap_mask)

                self.oi[(frag, frag2)] = -self.S_oi[(frag, frag2)]**2 / self.dE_oi[(frag, frag2)] * occ_virt_mask
                self.oi_approx_total += self.oi[(frag, frag2)].sum()

                # get data for pauli
                occ_occ_mask = np.logical_and(self.sfos_occ[frag], self.sfos_occ[frag2].T)
                self.S_pauli[(frag, frag2)] = overlap_mat(self.sfos[frag], self.sfos[frag2])
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
        for i, frag in enumerate(self.orbs.sfos.fragments):
            for frag2 in self.orbs.sfos.fragments[i+1:]:
                j = 0
                while 1:
                    best = np.unravel_index(argNmax(-self.oi[(frag, frag2)], j), self.oi[(frag, frag2)].shape)
                    best_oi = self.oi[(frag, frag2)][best]
                    sfo1, sfo2  = self.sfos[frag][best[0]], self.sfos[frag2][best[1]]

                    sfo1_contr = np.array([sfo1.mulliken_contribution(mo) for mo in self.mos])
                    sfo2_contr = np.array([sfo2.mulliken_contribution(mo) for mo in self.mos])

                    occ_contrs = abs(sfo1_contr * sfo2_contr) * self.mo_occ
                    virt_contrs = abs(sfo2_contr * sfo1_contr) * (1-self.mo_occ)

                    occ_mo = self.mos[np.argmax(occ_contrs)]
                    virt_mo = self.mos[np.argmax(virt_contrs)]
                    frac = best_oi / self.oi_approx_total
                    if fraction_thresh is None and j == N:
                        break
                    elif fraction_thresh is not None and frac < fraction_thresh:
                        break
                    stab = frac * self.oi_ref

                    typ = {conn: 'OI' for conn in list(it.product([sfo1, sfo2], [occ_mo, virt_mo]))}
                    mix = Mixing(self.orbs, [occ_mo, virt_mo], [sfo1, sfo2], stab, frac, connection_type=typ)
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
        for i, frag in enumerate(self.orbs.sfos.fragments):
            for frag2 in self.orbs.sfos.fragments[i+1:]:
                j = 0
                while 1:
                    best = np.unravel_index(argNmax(self.pauli[(frag, frag2)], j), self.pauli[(frag, frag2)].shape)
                    best_pauli = self.pauli[(frag, frag2)][best]
                    sfo1, sfo2  = self.sfos[frag][best[0]], self.sfos[frag2][best[1]]

                    sfo1_contr = np.array([sfo1.mulliken_contribution(mo) * mo.occupation for mo in self.mos])
                    sfo2_contr = np.array([sfo2.mulliken_contribution(mo) * mo.occupation for mo in self.mos])

                    occ_contrs = abs(sfo1_contr * sfo2_contr)

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

                    typ = {conn: 'PR' for conn in list(it.product([sfo1, sfo2], [occ_mo1, occ_mo2]))}
                    mix = Mixing(self.orbs, [occ_mo1, occ_mo2], [sfo1, sfo2], stab, frac, connection_type=typ)
                    ret.append(mix)
                    j += 1

        ret = sorted(ret, key=lambda row: -row.strength)
        if N is not None:
            ret = ret[:N]
        return ret


class Mixing:
    def __init__(self, orbs, mos=None, sfos=None, strength=None, fraction=None, connections=None, connection_type=None, energy_type='energy'):
        self.orbs = orbs
        self.mos = mos or []
        self.sfos = sfos or []
        self.strength = strength or 0
        self.fraction = fraction or 0
        self.energy_type = energy_type
        self.connections = connections
        self.connection_type = connection_type
        if connections is None:
            self.connections = list(it.product(self.sfos, self.mos))
        if connection_type is None:
            self.connection_type = {conn: 'Multiple' for conn in self.connections}
        if isinstance(connection_type, str):
            self.connection_type = {conn: connection_type for conn in self.connections}

        self.two_mixings = [[self]]

    def __str__(self):
        s = f'{self.__class__.__name__}('
        s += f'[{", ".join([sfo.make_name(frag_name=True, relative_name=True, spin=True) for sfo in self.sfos])}]'
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
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_sfos'] else 2
        for mix in self.two_mixings:
            mix = mix[0]
            if mix is self:
                continue
            if mix.nelectrons() == 2 * max_pop:
                return False
        return True

    @property
    def OI_is_empty(self):
        max_pop = 1 if self.orbs.data['calc_info']['unrestricted_sfos'] else 2
        for mix in self.two_mixings:
            mix = mix[0]
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
            for sfo in self.sfos:
                if abs(sfo.mulliken_contribution(mo)) > 0.03:
                    self.connections.append((sfo, mo))
                    self.connection_type[(sfo, mo)] = typ
        else:
            self.connections.append(connections)
            for conn in self.connections:
                self.connection_type[conn] = typ


    def add_sfo(self, sfo, typ='Sanitization', connections=None):
        '''
        Add an SFO to this mixing diagram. 
        '''
        self.sfos.append(sfo)
        if connections is None:
            for mo in self.mos:
                if abs(sfo.mulliken_contribution(mo)) > 0.03:
                    self.connections.append((sfo, mo))
                    self.connection_type[(sfo, mo)] = typ
        else:
            self.connections.append(connections)
            for conn in self.connections:
                self.connection_type[conn] = typ


    def draw_diagram(self, ax=None, ylim=None, simple=False, **kwargs):
        if simple:
            pyfmo.plotting.simple_orbital_diagram.draw_interaction(self.sfos, self.mos, self.connections, None, energy_type=self.energy_type, connection_types=self.connection_type, ax=ax, ylim=ylim)
        else:
            pyfmo.plotting.orbital_diagram.draw_interaction(self.sfos, self.mos, self.connections, None, energy_type=self.energy_type, connection_types=self.connection_type, ax=ax, ylim=ylim, **kwargs)

    def draw_sfos(self, overlap=False, screen=None):
        import tcviewer  # noqa: F811

        if screen is None:
            scr = tcviewer.Screen()
            scr.__enter__()
        else:
            scr = screen

        for sfo in self.sfos:
            with scr.add_molscene() as scene:
                cub = sfo.cube_file()
                scene.draw_text(str(sfo))

                if sfo.occupied:
                    colors = ([1, 0, 0], [0, 0, 1])
                else:
                    colors = ([0, 1, 1], [1, .5, 0])

                scene.draw_isosurface(cub, -.03, opacity=.25, color=colors[0])
                scene.draw_isosurface(cub,  .03, opacity=.25, color=colors[1])
                scene.draw_molecule(sfo.molecule)

        if screen is None:
            scr.__exit__()

    def screenshot_sfos(self, outdir='SFO_pictures'):
        import tcviewer  # noqa: F811

        os.makedirs(outdir, exist_ok=True)
        with tcviewer.Screen(headless=True) as scr:
            for sfo in self.sfos:
                with scr.add_molscene() as scene:
                    cub = sfo.cube_file()
                    scene.draw_molecule(sfo.molecule)

                    if sfo.occupied:
                        colors = ([1, 0, 0], [0, 0, 1])
                    else:
                        colors = ([0, 1, 1], [1, .5, 0])
                
                    scene.draw_isosurface(cub, -.03, opacity=.5, color=colors[0])
                    scene.draw_isosurface(cub,  .03, opacity=.5, color=colors[1])
                    scene.screenshot(os.path.join(outdir, f'{str(sfo)}.png'))

    def nelectrons(self):
        return sum(mo.occupation for mo in self.mos)

    def __add__(self, other: 'Mixing'):
        self.sfos.extend([osfo for osfo in other.sfos if osfo not in self.sfos])
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
        if len(self.sfos) == 0 and len(self.mos) == 0:
            return True
        if not any(sfo in other.sfos for sfo in self.sfos):
            return False
        if not any(mo in other.mos for mo in self.mos):
            return False
        return True

    def xiaobo_value(self):
        if all(sfo.occupied for sfo in self.sfos):
            maxval = 0
            for sfo1 in self.sfos:
                for sfo2 in self.sfos:
                    if sfo1 == sfo2:
                        continue
                    maxval = max(maxval, abs(sfo1.overlap(sfo2)))
            return maxval

        val = (self.sfos[0] @ self.sfos[1])**2 / abs(self.sfos[0].energy - self.sfos[1].energy)

        for sfo in self.sfos:
            gp = sfo.gross_population
            gp = np.clip(gp, 0, 2)
            excess = abs(sfo.occupation - gp)
            val *= excess

        return val

    def xiaobo_check(self, threshold=0.01):
        if all(sfo.occupied for sfo in self.sfos):
            for sfo1 in self.sfos:
                for sfo2 in self.sfos:
                    if sfo1 == sfo2:
                        continue
                    if abs(sfo1.overlap(sfo2)) < threshold:
                        return False
            return True

        val = (self.sfos[0] @ self.sfos[1])**2 / abs(self.sfos[0].energy - self.sfos[1].energy)

        for sfo in self.sfos:
            gp = sfo.gross_population
            gp = np.clip(gp, 0, 2)
            excess = abs(sfo.occupation - gp)
            val *= excess

        return val > threshold

    @property
    def fragments(self):
        return set(sfo.fragment_unique for sfo in self.sfos)

    @property
    def lowest_contribution(self):
        return min([abs(sfo.mulliken_contribution(mo)) for sfo in self.sfos for mo in self.mos])

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
            mos = [orb for orb in orbs_ if isinstance(orb, pyfmo.orbitals.objects.MO)]
            sfos = [orb for orb in orbs_ if isinstance(orb, pyfmo.orbitals.objects.SFO)]
            connections = subG.edges()
            connections = [conn[::-1] if isinstance(conn[0], pyfmo.orbitals.objects.MO) else conn for conn in connections]
            mixes.append(Mixing(
                self.orbs,
                mos=mos, 
                sfos=sfos, 
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
            if orb in two_mixing.sfos or orb in two_mixing.mos:
                ret.append([*two_mixing.sfos, *two_mixing.mos])
        return ret


def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
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
        sfos = [orbs.sfos[str(mix_sfo)] for mix_sfo in mixing.sfos]
        mos = [orbs.mos[str(mix_mo)] for mix_mo in mixing.mos]
        pyfmo.plotting.orbital_diagram.draw_interaction(sfos, mos, it.product(sfos, mos))

        plt.savefig(os.path.join(out_dir, f'{i}.jpg'))
        plt.close()


def _is_bonding(sfo1, sfo2, mo):
    return round(sfo1.coefficient(mo) * sfo2.coefficient(mo) * (sfo1 @ sfo2), 4) >= 0
