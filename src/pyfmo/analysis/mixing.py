import pyfmo
import tcutility
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import os
import scipy
import networkx as nx


class Mixer:
    def __init__(self, orbs, energy_type='energy'):
    def __init__(self, orbs: pyfmo.Orbitals, energy_type: str = 'energy'):
        self.orbs = orbs
        self.energy_type = energy_type
        self._prepare()

    def _prepare(self):
        C = self.orbs.data.matrices.mulliken_contribution.total
        self.sfos = {}
        self.sfos_occ = {}
        self.sfos_vir = {}
        self.sfos_energy = {}
        max_occ = max(sfo.occupation for sfo in self.orbs.sfos)
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
                # print(self.dE_oi[(frag, frag2)])
                nogap_mask = self.dE_oi[(frag, frag2)] != 0
                # total_mask = occ_virt_mask * nogap_mask
                self.dE_oi[(frag, frag2)] += (1 - nogap_mask)
                # print(self.dE_oi[(frag, frag2)])

                self.oi[(frag, frag2)] = -self.S_oi[(frag, frag2)]**2 / self.dE_oi[(frag, frag2)] * occ_virt_mask
                self.oi_approx_total += self.oi[(frag, frag2)].sum()
                # print(self.dE_oi[(frag, frag2)])

                # get data for pauli
                occ_occ_mask = np.logical_and(self.sfos_occ[frag], self.sfos_occ[frag2].T)
                # plt.imshow(occ_occ_mask)
                # plt.show()
                self.S_pauli[(frag, frag2)] = overlap_mat(self.sfos[frag], self.sfos[frag2])
                self.pauli[(frag, frag2)] = self.S_pauli[(frag, frag2)]**2 * occ_occ_mask
                self.pauli_approx_total += self.pauli[(frag, frag2)].sum()

        self.oi_ref = self.orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503
        self.pauli_ref = self.orbs.reader.read('Energy', 'Pauli Total') * 627.503


    def orbital_interactions(self, fraction_thresh=None, N=None):
        """
        Yield the first ``N`` strongest orbital interactions.
        """
        # print(self.oi_ref, self.oi_approx_total)
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
                    # print(best_oi, self.oi_approx_total)
                    frac = best_oi / self.oi_approx_total
                    # print(frac, fraction_thresh)
                    if fraction_thresh is None and j == N:
                        break
                    elif fraction_thresh is not None and frac < fraction_thresh:
                        break
                    stab = frac * self.oi_ref

                    col = {conn: 'g' for conn in list(it.product([sfo1, sfo2], [occ_mo, virt_mo]))}
                    mix = Mixing(self.orbs, [occ_mo, virt_mo], [sfo1, sfo2], stab, frac, connection_colors=col)
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
                    # plt.imshow(self.pauli[(frag, frag2)])
                    # plt.show()
                    best = np.unravel_index(argNmax(self.pauli[(frag, frag2)], j), self.pauli[(frag, frag2)].shape)
                    best_pauli = self.pauli[(frag, frag2)][best]
                    sfo1, sfo2  = self.sfos[frag][best[0]], self.sfos[frag2][best[1]]

                    sfo1_contr = np.array([sfo1.mulliken_contribution(mo) * mo.occupation for mo in self.mos])
                    sfo2_contr = np.array([sfo2.mulliken_contribution(mo) * mo.occupation for mo in self.mos])

                    occ_contrs = abs(sfo1_contr * sfo2_contr)

                    occ_mo1_idx = argNmax(occ_contrs, 0)
                    occ_mo2_idx = argNmax(occ_contrs, 1)

                    # print(sfo1, sfo2)
                    # print(occ_mo1_idx, occ_mo2_idx)

                    occ_mo1 = self.mos[occ_mo1_idx]
                    occ_mo2 = self.mos[occ_mo2_idx]

                    # print(occ_mo1, occ_mo1.occupation, occ_mo1.occupied)
                    # print(occ_mo2, occ_mo2.occupation, occ_mo2.occupied)
                    frac = best_pauli / self.pauli_approx_total
                    if fraction_thresh is None and j == N:
                        break
                    elif fraction_thresh is not None and frac < fraction_thresh:
                        break
                    stab = frac * self.pauli_ref

                    col = {conn: 'r' for conn in list(it.product([sfo1, sfo2], [occ_mo1, occ_mo2]))}
                    mix = Mixing(self.orbs, [occ_mo1, occ_mo2], [sfo1, sfo2], stab, frac, connection_colors=col)
                    ret.append(mix)
                    j += 1

        ret = sorted(ret, key=lambda row: -row.strength)
        if N is not None:
            ret = ret[:N]
        return ret



class Mixing:
    def __init__(self, orbs, mos=None, sfos=None, strength=None, fraction=None, connections=None, connection_colors=None, energy_type='energy'):
        self.orbs = orbs
        self.mos = mos or []
        self.sfos = sfos or []
        self.strength = strength or 0
        self.fraction = fraction or 0
        self.energy_type = energy_type
        self.connections = connections
        self.connection_colors = connection_colors
        if connections is None:
            self.connections = list(it.product(self.sfos, self.mos))
        if connection_colors is None:
            self.connection_colors = {conn: 'k' for conn in self.connections}
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

    def draw_diagram(self, ax=None, ylim=None):
        if self.strength is not None:
            title = rf'$\Delta E^{{({self.nelectrons()})}}_{{ij}} \approx {self.strength:5.1f}$ kcal/mol (${self.fraction:5.1%}$% of total)'
        else:
            title = None

        pyfmo.plotting.orbital_diagram.draw_interaction(self.sfos, self.mos, self.connections, None, energy_type=self.energy_type, connection_colors=self.connection_colors, ax=ax, ylim=ylim)

    def draw_sfos(self, overlap=False, screen=None):
        import tcviewer

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
        import tcviewer

        mols = list(set(sfo.molecule for sfo in self.sfos))
        centroids = [np.mean(mol, axis=0) for mol in mols]
        # print(centroids)
        mol = mols[0]
        for mol_ in mols[1:]:
            mol = mol + mol_

        coordinates = np.array(mol)

        # for this to work we should first get the centroid of our molecule
        centroid = np.mean(coordinates, axis=0)
        # and get the centered coordiantes
        Xc = coordinates - centroid

        # we then do a singular-value decomposition to obtain
        # the three principle components (Vh) with their eigenvalues (s)
        _, s, Vh = scipy.linalg.svd(Xc)

        # then compute a transformation matrix for generating the correct spheroid
        transform = tcutility.geometry.Transform()
        transform.translate(centroid)
        transform.rotate((np.diag(s/2) @ Vh).T)
        centroids = transform(centroids)
        transform.rotate(tcutility.geometry.vector_align_rotmat(centroids[0] - centroids[1], [0, 1, 0]))
        # mol = transform(mol)

        os.makedirs(outdir, exist_ok=True)
        with tcviewer.Screen(headless=True) as scr:
            for sfo in self.sfos:
                with scr.add_molscene() as scene:
                    scene.transform = transform.to_vtkTransform()
                    cub = sfo.cube_file()
                    scene.draw_molecule(mol)

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
        for conn, col in other.connection_colors.items():
            if conn in self.connection_colors:
                if col == self.connection_colors[conn]:
                    continue
                else:
                    self.connection_colors[conn] = 'k'
            else:
                self.connection_colors[conn] = col

        # self.connection_colors.update(other.connection_colors)
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

        val = 1
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

        val = 1
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
            mos = [orb for orb in orbs_ if isinstance(orb, pyfmo.orbitals2.objects.MO)]
            sfos = [orb for orb in orbs_ if isinstance(orb, pyfmo.orbitals2.objects.SFO)]
            connections = subG.edges()
            connections = [conn[::-1] if isinstance(conn[0], pyfmo.orbitals2.objects.MO) else conn for conn in connections]
            mixes.append(Mixing(
                self.orbs,
                mos=mos, 
                sfos=sfos, 
                connections=connections,
                connection_colors={conn: self.connection_colors[conn] for conn in connections},
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

    def _add_extra_orbital(self, orbitals, occupied, nelec=2, add_MO=True):
        max_contr = -1
        max_contr_orb = None
        max_contr_orb1 = None
        max_contr_orb2 = None

        for orbital in orbitals:
            if orbital.occupied != occupied:
                continue
            if orbital in self.mos or orbital in self.sfos:
                continue
            if orbital.occupation > nelec:
                continue

            if add_MO:
                contr = np.array([sfo.mulliken_contribution(orbital) for sfo in self.sfos])
                orb1 = self.sfos[argNmax(contr, 0)]
                orb2 = self.sfos[argNmax(contr, 1)]
                c1, c2 = contr[argNmax(contr, 0)], contr[argNmax(contr, 1)]
            else:
                contr = np.array([orbital.mulliken_contribution(mo) for mo in self.mos])
                orb1 = self.mos[argNmax(contr, 0)]
                orb2 = self.mos[argNmax(contr, 1)]
                c1, c2 = contr[argNmax(contr, 0)], contr[argNmax(contr, 1)]

            contr1 = np.clip(c1, 0, 1)
            contr2 = np.clip(c2, 0, 1)

            if contr1 * contr2 > max_contr:
                max_contr = contr1 * contr2
                max_contr_orb = orbital
                max_contr_orb1 = orb1
                max_contr_orb2 = orb2

        if add_MO:
            self.mos.append(max_contr_orb)
            self.connections.append([max_contr_orb1, max_contr_orb])
            self.connections.append([max_contr_orb2, max_contr_orb])
            self.connection_colors[max_contr_orb1, max_contr_orb] = 'purple'
            self.connection_colors[max_contr_orb2, max_contr_orb] = 'purple'
        else:
            self.sfos.append(max_contr_orb)
            self.connections.append([max_contr_orb, max_contr_orb1])
            self.connections.append([max_contr_orb, max_contr_orb2])
            self.connection_colors[max_contr_orb, max_contr_orb1] = 'purple'
            self.connection_colors[max_contr_orb, max_contr_orb2] = 'purple'


    def _add_extra_virtual_mo(self, nelec=2):
        max_contr = -1
        max_contr_mo = None
        max_contr_sfo1 = None
        max_contr_sfo2 = None
        for mo in self.orbs.mos:
            if mo.occupied:
                continue
            if mo in self.mos:
                continue
            if mo.occupation > nelec:
                continue

            contr = np.array([sfo.mulliken_contribution(mo) for sfo in self.sfos])
            sfo1 = self.sfos[argNmax(contr, 0)]
            sfo2 = self.sfos[argNmax(contr, 1)]
            contr1 = np.clip(sfo1.mulliken_contribution(mo), 0, 1)
            contr2 = np.clip(sfo2.mulliken_contribution(mo), 0, 1)

            if contr1 * contr2 > max_contr:
                max_contr = contr1 * contr2
                max_contr_mo = mo
                max_contr_sfo1 = sfo1
                max_contr_sfo2 = sfo2

        self.mos.append(max_contr_mo)
        self.connections.append([max_contr_sfo1, max_contr_mo])
        self.connections.append([max_contr_sfo2, max_contr_mo])
        self.connection_colors[max_contr_sfo1, max_contr_mo] = 'purple'
        self.connection_colors[max_contr_sfo2, max_contr_mo] = 'purple'

    def _add_extra_occupied_mo(self, nelec=2):
        max_contr = -1
        max_contr_mo = None
        max_contr_sfo1 = None
        max_contr_sfo2 = None
        for mo in self.orbs.mos:
            if not mo.occupied:
                continue
            if mo in self.mos:
                continue

            print(mo, mo.occupation, nelec, mo.occupation > nelec)
            if mo.occupation > nelec:
                continue
            print('hello')

            contr = np.array([sfo.mulliken_contribution(mo) for sfo in self.sfos])
            sfo1 = self.sfos[argNmax(contr, 0)]
            sfo2 = self.sfos[argNmax(contr, 1)]
            contr1 = np.clip(sfo1.mulliken_contribution(mo), 0, 1)
            contr2 = np.clip(sfo2.mulliken_contribution(mo), 0, 1)
            # print(mo, contr1 * contr2, max_contr, contr1 * contr2 > max_contr)
            if contr1 * contr2 > max_contr:
                max_contr = contr1 * contr2
                max_contr_mo = mo
                max_contr_sfo1 = sfo1
                max_contr_sfo2 = sfo2
        # print(max_contr_mo)
        self.mos.append(max_contr_mo)
        self.connections.append([max_contr_sfo1, max_contr_mo])
        self.connections.append([max_contr_sfo2, max_contr_mo])
        self.connection_colors[max_contr_sfo1, max_contr_mo] = 'purple'
        self.connection_colors[max_contr_sfo2, max_contr_mo] = 'purple'

    def _add_extra_virtual_sfo(self, nelec=2):
        max_contr = -1
        max_contr_sfo = None
        max_contr_mo1 = None
        max_contr_mo2 = None
        virt_mos = [mo for mo in self.mos if not mo.occupied]
        for sfo in self.orbs.sfos:
            if sfo.occupied:
                continue
            if sfo in self.sfos:
                continue
            if sfo.occupation > nelec:
                continue

            contr = np.array([sfo.mulliken_contribution(mo) for mo in virt_mos])
            mo1 = virt_mos[argNmax(contr, 0)]
            mo2 = virt_mos[argNmax(contr, 1)]
            contr1 = np.clip(sfo.mulliken_contribution(mo1), 0, 1)
            contr2 = np.clip(sfo.mulliken_contribution(mo2), 0, 1)

            if contr1 * contr2 > max_contr:
                max_contr = contr1 * contr2
                max_contr_sfo = sfo
                max_contr_mo1 = mo1
                max_contr_mo2 = mo2

        self.sfos.append(max_contr_sfo)
        self.connections.append([max_contr_sfo, max_contr_mo1])
        self.connections.append([max_contr_sfo, max_contr_mo2])
        self.connection_colors[max_contr_sfo, max_contr_mo1] = 'purple'
        self.connection_colors[max_contr_sfo, max_contr_mo2] = 'purple'

    def _add_extra_occupied_sfo(self, nelec=2):
        max_contr = -1
        max_contr_sfo = None
        max_contr_mo1 = None
        max_contr_mo2 = None
        occ_mos = [mo for mo in self.mos if mo.occupied]
        for sfo in self.orbs.sfos:
            if not sfo.occupied:
                continue
            if sfo in self.sfos:
                continue
            if sfo.occupation > nelec:
                continue

            contr = np.array([sfo.mulliken_contribution(mo) for mo in occ_mos])
            mo1 = occ_mos[argNmax(contr, 0)]
            mo2 = occ_mos[argNmax(contr, 1)]
            contr1 = abs(sfo.mulliken_contribution(mo1))
            contr2 = abs(sfo.mulliken_contribution(mo2))

            if contr1 * contr2 > max_contr:
                max_contr = contr1 * contr2
                max_contr_sfo = sfo
                max_contr_mo1 = mo1
                max_contr_mo2 = mo2

        self.sfos.append(max_contr_sfo)
        self.connections.append([max_contr_sfo, max_contr_mo1])
        self.connections.append([max_contr_sfo, max_contr_mo2])
        self.connection_colors[max_contr_sfo, max_contr_mo1] = 'purple'
        self.connection_colors[max_contr_sfo, max_contr_mo2] = 'purple'


    def sanitize(self):
        def excess_elec():
            nsfos_elec = sum([sfo.occupation for sfo in mix.sfos])
            nmos_elec = sum([mo.occupation for mo in mix.mos])
            return int(nsfos_elec - nmos_elec)

        def excess_virt():
            nsfos_virt = len([sfo for sfo in mix.sfos if sfo.occupation == 0])
            nmos_virt = len([mo for mo in mix.mos if mo.occupation == 0])
            return nsfos_virt - nmos_virt

        def excess_half():
            nsfos_half = len([sfo for sfo in mix.sfos if sfo.occupation == 1])
            nmos_half = len([mo for mo in mix.mos if mo.occupation == 1])
            return nsfos_half - nmos_half

        def excess_occ():
            nsfos_occ = len([sfo for sfo in mix.sfos if sfo.occupation == 2])
            nmos_occ = len([mo for mo in mix.mos if mo.occupation == 2])
            return nsfos_occ - nmos_occ

        sub_mixes = self.split()
        for mix in sub_mixes:
            print('ello')
            print(excess_elec(), excess_virt(), excess_half(), excess_occ())

            if excess_virt() + excess_half() + excess_occ() == 0:
                continue
            
            # check the number of SFOs and MOs
            if excess_virt() > 0:
                for i in range(excess_virt()):
                    try:
                        self._add_extra_orbital(self.orbs.mos, add_MO=True, occupied=False)
                    except Exception:
                        pass
            elif excess_virt() < 0:
                for i in range(-excess_virt()):
                    try:
                        self._add_extra_orbital(self.orbs.sfos, add_MO=False, occupied=False)
                    except Exception:
                        pass

            if excess_occ() > 0:
                for i in range(excess_occ()):
                    try:
                        self._add_extra_orbital(self.orbs.mos, add_MO=True, occupied=True)
                    except Exception:
                        pass

            if excess_occ() < 0:
                for i in range(-excess_occ()):
                    try:
                        self._add_extra_orbital(self.orbs.sfos, add_MO=False, occupied=True)
                    except Exception:
                        pass

            print(excess_elec(), excess_virt(), excess_half(), excess_occ())


def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in tcutility.ensure_list(sfos1):
        ret.append([])
        for sfo2 in tcutility.ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


def argNmax(arr, N):
    """
    Get the Nth maximum element of an array.
    """
    return np.argsort(-arr, axis=None)[N]


def track_mixing(orbss, mixing):
    ret = []
    out_dir = 'pyfrag_imgs'
    os.makedirs(out_dir, exist_ok=True)
    for i, orbs in enumerate(orbss):
        sfos = [orbs.sfos[str(mix_sfo)] for mix_sfo in mixing.sfos]
        mos = [orbs.mos[str(mix_mo)] for mix_mo in mixing.mos]
        pyfmo.plotting.orbital_diagram.draw_interaction(sfos, mos, it.product(sfos, mos))

        plt.savefig(os.path.join(out_dir, f'{i}.jpg'))
        plt.close()


def _find_orbs(sfos1, sfos2, mos):
    '''
    Find MO which has contributions from both fragments 
    and the SFOs that have the strongest contribution to it.
    '''
    res = []
    for mo in mos:
        sfos1_contr = sum([abs(sfo.mulliken_contribution(mo)) for sfo in sfos1])
        sfos2_contr = sum([abs(sfo.mulliken_contribution(mo)) for sfo in sfos2])
        res.append(sfos1_contr * sfos2_contr)

    mo = mos[np.argmax(res)]
    sfos1_contr = [abs(sfo.mulliken_contribution(mo)) for sfo in sfos1]
    sfos2_contr = [abs(sfo.mulliken_contribution(mo)) for sfo in sfos2]

    sfo1 = sfos1[np.argmax(sfos1_contr)]
    sfo2 = sfos2[np.argmax(sfos2_contr)]

    strength = abs(sfo1.mulliken_contribution(mo) * sfo2.mulliken_contribution(mo))

    return sfo1, sfo2, mo, strength


def _is_bonding(sfo1, sfo2, mo):
    return round(sfo1.coefficient(mo) * sfo2.coefficient(mo) * (sfo1 @ sfo2), 4) >= 0

@tcutility.cache.cache
def pauli2(orbs, index=0, mo1=None, mo2=None, sfo1=None, sfo2=None, irrep=None, spin=None):
    frag1, frag2 = orbs.fragments[0], orbs.fragments[1]

    sfos1 = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag1) if sfo.occupied]
    sfos2 = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag2) if sfo.occupied]

    if irrep is not None:
        sfos1 = [sfo for sfo in sfos1 if sfo.symmetry == irrep]
        sfos2 = [sfo for sfo in sfos2 if sfo.symmetry == irrep]

    if spin is not None:
        sfos1 = [sfo for sfo in sfos1 if sfo.spin == spin]
        sfos2 = [sfo for sfo in sfos2 if sfo.spin == spin]

    if sfo1 in sfos1:
        sfos1 = [sfo1]
    if sfo1 in sfos2:
        sfos2 = [sfo1]
    if sfo2 in sfos1:
        sfos1 = [sfo2]
    if sfo2 in sfos2:
        sfos2 = [sfo2]

    if mo1 is None and mo2 is None:
        mos = [mo for mo in orbs.mos if mo.occupied]
    else:
        if mo2 is not None:
            mo1, mo2 = mo2, mo1
        mos = [mo1]

    sfo1, sfo2, mo1, strength1 = _find_orbs(sfos1, sfos2, mos)
    if _is_bonding(sfo1, sfo2, mo1):
        mos = [mo for mo in orbs.mos if mo.occupied and mo.energy > mo1.energy]
    else:
        mos = [mo for mo in orbs.mos if mo.occupied and mo.energy < mo1.energy]

    if mo2 is None:
        _, _, mo2, strength2 = _find_orbs([sfo1], [sfo2], mos)
    else:
        _, _, _, strength2 = _find_orbs([sfo1], [sfo2], [mo2])

    return Mixing(
        orbs,
        mos=[mo1, mo2],
        sfos=[sfo1, sfo2],
        connections=[(sfo1, mo1), (sfo2, mo1), (sfo1, mo2), (sfo2, mo2)],
        strength=strength1 * strength2)


def oi2(orbs, index=0, irrep=None):
    frag1, frag2 = orbs.fragments[0], orbs.fragments[1]

    frag1_sfos_occ = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag1) if sfo.occupied]
    frag2_sfos_occ = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag2) if sfo.occupied]

    frag1_sfos_vir = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag1) if not sfo.occupied]
    frag2_sfos_vir = [sfo for sfo in orbs.sfos.get_fragment_sfos(frag2) if not sfo.occupied]

    if irrep is not None:
        frag1_sfos_occ = [sfo for sfo in frag1_sfos_occ if sfo.symmetry == irrep]
        frag2_sfos_occ = [sfo for sfo in frag2_sfos_occ if sfo.symmetry == irrep]
        frag1_sfos_vir = [sfo for sfo in frag1_sfos_vir if sfo.symmetry == irrep]
        frag2_sfos_vir = [sfo for sfo in frag2_sfos_vir if sfo.symmetry == irrep]

    # first select an occupied orbital that has the most even contribution from both frags
    res = []
    mos = []
    direction = []
    for i, mo in enumerate(orbs.mos):
        if not mo.occupied:
            continue

        frag1_contr_occ = sum([sfo.mulliken_contribution(mo) for sfo in frag1_sfos_occ])
        frag2_contr_occ = sum([sfo.mulliken_contribution(mo) for sfo in frag2_sfos_occ])
        frag1_contr_vir = sum([sfo.mulliken_contribution(mo) for sfo in frag1_sfos_vir])
        frag2_contr_vir = sum([sfo.mulliken_contribution(mo) for sfo in frag2_sfos_vir])

        strength1 = frag1_contr_occ * frag2_contr_vir
        strength2 = frag1_contr_vir * frag2_contr_occ
        if strength1 > strength2:
            res.append(strength1)
            direction.append(0)
        else:
            res.append(strength2)
            direction.append(1)
        mos.append(mo)

    best_mo = mos[np.argmax(res)]
    if direction == 0:
        frag1_contr = [abs(sfo.mulliken_contribution(best_mo)) for sfo in frag1_sfos_occ]
        frag2_contr = [abs(sfo.mulliken_contribution(best_mo)) for sfo in frag2_sfos_vir]
        best_sfo1 = frag1_sfos_occ[np.argmax(frag1_contr)]
        best_sfo2 = frag2_sfos_vir[np.argmax(frag2_contr)]
    else:
        frag1_contr = [abs(sfo.mulliken_contribution(best_mo)) for sfo in frag1_sfos_vir]
        frag2_contr = [abs(sfo.mulliken_contribution(best_mo)) for sfo in frag2_sfos_occ]
        best_sfo1 = frag1_sfos_vir[np.argmax(frag1_contr)]
        best_sfo2 = frag2_sfos_occ[np.argmax(frag2_contr)]

    other_mos = [mo for mo in orbs.mos if not mo.occupied]
    other_mixes = []
    for i, mo in enumerate(other_mos):
        if mo.occupied:
            continue

        frag1_contr = best_sfo1.mulliken_contribution(mo)
        frag2_contr = best_sfo2.mulliken_contribution(mo)

        other_mixes.append(frag1_contr * frag2_contr)

    best_other_mo = other_mos[np.argmax(other_mixes)]

    return Mixing(
        orbs, 
        mos=[best_mo, best_other_mo], 
        sfos=[best_sfo1, best_sfo2], 
        connections=[(best_sfo1, best_mo), (best_sfo2, best_mo), (best_sfo1, best_other_mo), (best_sfo2, best_other_mo)])



# if __name__ == '__main__':
#     orbs = pyfmo.orbitals2.objects.Orbitals('../../../calculations/PyOrb_testing_2022/HydrogenBond/GuanineCytosine.results/adf.rkf')
    
#     C = orbs.data.matrices.mulliken_contribution.total
#     vals, vecs = np.linalg.eig(C.T)
#     print(vecs)
#     plt.imshow(vecs.real)
#     plt.show()
#     sfo1 = orbs.sfos['Cytosine(24AA)']
#     sfo2 = orbs.sfos['Guanine(33AA)']

#     losses = []
#     for mo in orbs.mos:
#         loss = sfo1.mulliken_contribution(mo)
#         losses.append(loss)

#     gains = []
#     for i, mo in enumerate(orbs.mos):
#         gain = sfo2.mulliken_contribution(mo) * losses[i]
#         if gain > 1e-5:
#             print(mo, round(gain, 5))
#         gains.append(gain)

#     losses = np.array(losses)
#     gains = np.array(gains)
#     print(losses @ gains)

#     print(sum(losses))
#     print(sum(gains))
#     # plt.plot(losses)
#     # plt.plot(gains)
#     plt.plot(gains)
#     plt.show()



# # exit()
if __name__ == '__main__':
    import networkx as nx
    import tcviewer
    import matplotlib.pyplot as plt
    p = '../../../calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.results/'
    # p = '../../../calculations/PyOrb_testing_2022/TransitionState/DielsAlder.results/'
    # p = '../../../calculations/PyOrb_testing_2022/CoordinationBondFeCO4CO/FeCO4CO.results/'
    # p = '../../../calculations/PyOrb_testing_2022/CoordinationBondFeCO4CH4/FeCO4CH4.results/'
    # p = '../../../calculations/PyOrb_testing_2022/HydrogenBond/GuanineCytosine.results/'
    # p = '../../../calculations/PyOrb_testing_2022/HeterolyticBond/complex/'
    # p = '../../../calculations/PyOrb_testing_2022/HomolyticBond/complex/'
    # p = '../../../calculations/PyOrb_testing_2022/ChemicalBond/frag.results/'
    orbs = pyfmo.orbitals2.objects.Orbitals(p + 'adf.rkf')
    # res = tcutility.results.read(p)

    frag1_sfos = orbs.sfos.get_fragment_sfos('NH3')
    frag2_sfos = orbs.sfos.get_fragment_sfos('BH3')

    print(frag1_sfos)
    occ1 = np.array([sfo.occupied for sfo in frag1_sfos]).reshape(-1, 1)
    occ2 = np.array([sfo.occupied for sfo in frag2_sfos]).reshape(-1, 1)
    E1 = np.array([sfo.energy for sfo in frag1_sfos]).reshape(-1, 1)
    E2 = np.array([sfo.energy for sfo in frag2_sfos]).reshape(-1, 1).T
    occ_virt_mask = np.logical_xor(occ1.T, occ2)
    plt.imshow(occ_virt_mask)
    plt.show()
    S = np.array([[abs(sfo1 @ sfo2) for sfo1 in frag1_sfos] for sfo2 in frag2_sfos])
    dE = np.array([[abs(sfo1.energy - sfo2.energy) for sfo1 in frag1_sfos] for sfo2 in frag2_sfos])
    
    plt.imshow(S**2/dE * occ_virt_mask)
    plt.show()

    stab = []
    K = 1.75
    for sfo1 in frag1_sfos:
        stab.append([])
        for sfo2 in frag2_sfos:
            S = abs(sfo1 @ sfo2)
            e1, e2 = list(sorted([sfo1.energy, sfo2.energy]))
            print(e1, e2)
            H = K * S * (e1 + e2) / 2
            stab[-1].append((H - e1*S)**2 / (e1 - e2))
            # stab[-1].append(-S**2 * ((K/2 - 1) * e1 + K/2*e2)**2/(e1-e2))

    stab = np.array(stab)
    # stab = S**2 * ((K/2 - 1) * E1 + K/2*E2)**2/(E1-E2)
    plt.imshow(stab * occ_virt_mask)
    plt.show()


    # frag1_orbs = pyfmo.orbitals2.objects.Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.NH3.results/adf.rkf')
    # frag2_orbs = pyfmo.orbitals2.objects.Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/DonorAcceptor/NH3BH3.BH3.results/adf.rkf')
    # mo = orbs.mos['8A1']
    # print(mo)
    # contr = np.array([sfo.mulliken_contribution(mo) for sfo in orbs.sfos])
    # energ = np.array([sfo.energy for sfo in orbs.sfos])

    # print(mo.energy)
    # print(sum(contr * energ))
    # print(mo.energy - sum(contr * energ))

    # mo = orbs.mos['5A1']
    # print(mo)
    # contr = np.array([sfo.mulliken_contribution(mo) for sfo in orbs.sfos])
    # energ = np.array([sfo.energy for sfo in orbs.sfos])

    # print(mo.energy)
    # print(sum(contr * energ))
    # print(mo.energy - sum(contr * energ))
#     # mix = oi2(orbs)
#     # # mix = pauli2(orbs, sfo1=orbs.sfos['Guanine(24AA)'], sfo2=orbs.sfos['Cytosine(23AA)'])
#     # mix += pauli2(orbs)
#     # mix.draw_diagram()
#     # plt.title('New way')
#     # # plt.show()
#     # # exit()


#     # # main_mix = None
#     # main_G = nx.Graph()
#     # with tcviewer.Screen() as scr:
#     #     mixer = Mixer(orbs)
#     #     mixes = mixer.orbital_interactions(N=2)

#     #     for mix in mixes:
#     #         if main_mix is None:
#     #             main_mix = mix
#     #             continue
#     #         main_mix += mix
#     #         continue
#     #         print(mix)
#     #         # plt.figure()
#     #         # mix.draw_diagram()
#     #         # plt.show()
#     #         mix.draw_sfos(screen=scr)

#     #         # for sfo in mix.sfos:
#     #         #     main_G.add_node(sfo)
#     #         # for mo in mix.mos:
#     #         #     main_G.add_node(mo)

#     #         for sfo in mix.sfos:
#     #             main_G.add_node(sfo, pos=(orbs.fragments.index(sfo.fragment_unique) * 2 - 1, sfo.energy), type='sfo')
#     #             for mo in mix.mos:
#     #                 main_G.add_node(mo, pos=(0, mo.energy), type='mo')
#     #                 main_G.add_edge(sfo, mo)

#     #     mixes = mixer.pauli_repulsions(N=2)
#     #     for mix in mixes:
#     #         print(mix)
#     #         main_mix += mix
#     #         continue
#     #         # plt.figure()
#     #         # mix.draw_diagram()
#     #         # plt.show()
#     #         mix.draw_sfos(screen=scr)

#     #         # for sfo in mix.sfos:
#     #         #     main_G.add_node(sfo)
#     #         # for mo in mix.mos:
#     #         #     main_G.add_node(mo)

#     #         for sfo in mix.sfos:
#     #             main_G.add_node(sfo, pos=(orbs.fragments.index(sfo.fragment_unique) * 2 - 1, sfo.energy), type='sfo')
#     #             for mo in mix.mos:
#     #                 main_G.add_node(mo, pos=(0, mo.energy), type='mo')
#     #                 main_G.add_edge(sfo, mo)
#     # main_mix.draw_diagram()
#     # plt.show()
#     # components = [main_G.subgraph(H) for H in nx.connected_components(main_G)]
#     # for component in components:
#     #     nodes = component.nodes(data='type')
#     #     sfos = [node[0] for node in nodes if node[1] == 'sfo']
#     #     mos = [node[0] for node in nodes if node[1] == 'mo']
#     #     print(sfos, mos)
#     #     pyfmo.plotting.orbital_diagram.draw_interaction(sfos, mos, it.product(sfos, mos), orbs, energy_type='energy')
#     #     plt.show()
#     from matplotlib.widgets import Slider, CheckButtons
#     from matplotlib.gridspec import GridSpec


#     mixer = Mixer(orbs, energy_type='energy')

#     def draw_diagram(arg):
#         plt.cla()
#         plt.title('Old way')
#         mix = Mixing(orbs)
#         if oi_b.get_status()[0]:
#             for mix_ in mixer.orbital_interactions(N=10):
#                 # print(mix_, mix_.nelectrons())
#                 if mix_.xiaobo_check(oi_s.val):
#                     mix += mix_

#         if pauli_b.get_status()[0]:
#             for mix_ in mixer.pauli_repulsions(N=1):
#                 if mix_.xiaobo_check(pauli_s.val):
#                     mix += mix_
#         mix.draw_diagram()
#         plt.gcf().canvas.draw_idle()

#     # plt.subplots()
#     # mix = mixer.orbital_interactions(N=1)[0]
#     # mix += mixer.pauli_repulsions(N=1)[0]
#     # mixes.extend()
#     # plt.figure()

#     gs = GridSpec(nrows=4, ncols=2, height_ratios=[1, .05, .05, .05], width_ratios=[.1, .7])

#     oi_bax = plt.gcf().add_subplot(gs[2, 0])
#     pauli_bax = plt.gcf().add_subplot(gs[3, 0])

#     oi_bax.axis('off')
#     pauli_bax.axis('off')

#     oi_sax = plt.gcf().add_subplot(gs[2, 1])
#     pauli_sax = plt.gcf().add_subplot(gs[3, 1])

#     oi_b = CheckButtons(oi_bax, labels=['Show'], actives=[True])
#     pauli_b = CheckButtons(pauli_bax, labels=['Show'], actives=[True])
    
#     oi_s = Slider(oi_sax, 'OI', 0.001, mixer.orbital_interactions(N=1)[0].lowest_contribution, valinit=.025, facecolor='g')
#     pauli_s = Slider(pauli_sax, 'Pauli', 0.001, mixer.pauli_repulsions(N=1)[0].lowest_contribution, valinit=0.2, facecolor='r')
    
#     oi_b.on_clicked(draw_diagram)
#     pauli_b.on_clicked(draw_diagram)

#     oi_s.on_changed(draw_diagram)
#     pauli_s.on_changed(draw_diagram)

#     plt.gcf().add_subplot(gs[0, :])
#     draw_diagram(mixer)

#     # plt.tight_layout()
#     plt.show()

#     exit()
#     # combined_mixes = [Mixing(orbs)]
#     # for mix in mixes:
#     #     print(mix)
#         # for cmix in combined_mixes:
#             # print(mix)
#             # if cmix.fits(mix):
#     #         if True:
#     #             cmix += mix
#     #             break
#     #         else:
#     #             combined_mixes.append(mix)

#     # print(combined_mixes)
#     # for cmix in combined_mixes:
#     #     # cmix.screenshot_sfos(outdir='/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/ChemicalBond/frag.results/orbs')
#     #     # cmix.draw_sfos()
#     #     # plt.figure()
#     #     cmix.draw_diagram()
#     # plt.show()
