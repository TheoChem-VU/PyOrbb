from PySide6 import QtWidgets, QtCore
import pyfmo


class OrbitalSelectionState:
    def __init__(self, orbs):
        self.orbs = orbs
        self.reset()

    def reset(self):
        # to indicate MOs we use the mos object as the key
        self.orbitals = {self.orbs.mos: {}}
        self.spins = {self.orbs.mos: {}}
        self.irreps = {self.orbs.mos: {}}

        for mo in self.orbs.mos:
            self.orbitals[self.orbs.mos][mo] = None
            self.spins[self.orbs.mos].setdefault(mo.spin, True)
            self.irreps[self.orbs.mos].setdefault(mo.symmetry, True)

        for frag in self.orbs.fragments:
            self.orbitals[frag] = {}
            self.spins[frag] = {}
            self.irreps[frag] = {}
            for sfo in self.orbs.sfos.filter(fragment=frag):
                self.orbitals[frag][sfo] = None
                self.spins[frag].setdefault(sfo.spin, True)
                self.irreps[frag].setdefault(sfo.subspecies, True)

    def mo_states(self):
        ret = {}
        for mo in self.orbs.mos:
            if self.orbitals[self.orbs.mos][mo] is False:
                ret[mo] = False
                continue

            if not self.spins[self.orbs.mos][mo.spin]:
                ret[mo] = False
                continue

            if not self.irreps[self.orbs.mos][mo.symmetry]:
                ret[mo] = False
                continue

            ret[mo] = True

        return ret

    def sfo_states(self, frag):
        ret = {}
        for sfo in self.orbs.sfos.filter(fragment=frag):
            if self.orbitals[frag][sfo] is False:
                ret[sfo] = False
                continue

            if not self.spins[frag][sfo.spin]:
                ret[sfo] = False
                continue

            if not self.irreps[frag][sfo.subspecies]:
                ret[sfo] = False
                continue

            ret[sfo] = True

        return ret

    def is_enabled(self, orbital: 'MO or SFO'):
        if isinstance(orbital, pyfmo.orbitals.objects.MO):
            return self.mo_states()[orbital]
        else:
            return self.sfo_states(orbital.fragment_unique)[orbital]

    def set_state(self, orbital: 'MO or SFO', state):
        if isinstance(orbital, pyfmo.orbitals.objects.MO):
            self.orbitals[self.orbs.mos][orbital] = state
        else:
            self.orbitals[orbital.fragment_unique][orbital] = state


class OrbitalSelectionTable(QtWidgets.QTableWidget):
    def __init__(self, parent, orbitals: "List[MO] or List[SFO]"):
        self.parent = parent
        self.orbital_checkboxes = {}
        self.orbitals = orbitals
        is_MO = isinstance(orbitals[0], pyfmo.orbitals.objects.MO)

        column_widths = {
            'Orbital': 150,
            'Spin': 40,
            'Occ.': 40,
            'Gross Pop.': 70,
            'Symm.': 50,
            'Subsp.': 50,
            'Energy (reg.)': 80,
            'Energy (eff.)': 80,
            'Rel. Name': 130,
            'Rel. Name (Symm.)': 130,
            'Rel. Name (Subsp.)': 130,
        }

        headers = list(column_widths)

        if is_MO:
            headers.remove('Subsp.')
            headers.remove('Gross Pop.')
            headers.remove('Rel. Name (Subsp.)')

        if not self.parent.parent.orbs.data['calc_info']['has_site_energy'] or is_MO:
            headers.remove('Energy (eff.)')

        super().__init__(len(orbitals), len(headers), parent=parent)

        self.verticalHeader().setDefaultSectionSize(35)
        self.verticalHeader().setVisible(False)
        self.setHorizontalHeaderLabels(headers)
        for i, header in enumerate(headers):
            self.setColumnWidth(i, column_widths[header])

        for i, orb in enumerate(orbitals):
            orbital_frame = QtWidgets.QFrame()
            layout_ = QtWidgets.QHBoxLayout()
            orbital_frame.setLayout(layout_)
            
            self.orbital_checkboxes[orb] = QtWidgets.QCheckBox()
            self.orbital_checkboxes[orb].setChecked(self.parent.parent.state.is_enabled(orb))
            self.orbital_checkboxes[orb].checkStateChanged.connect(self.multi_select)

            layout_.addWidget(self.orbital_checkboxes[orb])
            layout_.addWidget(QtWidgets.QLabel(pyfmo.generate_label(orb, mode='html')))

            self.setCellWidget(i, headers.index('Orbital'), orbital_frame)

            for j, header in enumerate(headers):
                if header == 'Spin':
                    spin = {'A': '<i>α</i>', 'B': '<i>β</i>', 'AB': '<i>αβ</i>'}[orb.spin]
                    label = QtWidgets.QLabel(spin)
                elif header == 'Occ.':
                    label = QtWidgets.QLabel(str(round(orb.occupation, 3)))
                elif header == 'Gross Pop.':
                    label = QtWidgets.QLabel(f'{orb.gross_population:.3f}')
                elif header == 'Symm.':
                    label = QtWidgets.QLabel(pyfmo.translate_irrep_label(orb.symmetry, mode='html'))
                elif header == 'Subsp.':
                    label = QtWidgets.QLabel(pyfmo.translate_irrep_label(orb.subspecies, mode='html'))
                elif header == 'Rel. Name':
                    label = QtWidgets.QLabel(orb.relative_name)
                elif header == 'Rel. Name (Symm.)':
                    label = QtWidgets.QLabel(orb.symmetry_relative_name)
                elif header == 'Rel. Name (Subsp.)':
                    label = QtWidgets.QLabel(orb.subspecies_relative_name)
                elif header == 'Energy (reg.)':
                    label = QtWidgets.QLabel(f'{orb.energy: .2f}')
                elif header == 'Energy (eff.)':
                    label = QtWidgets.QLabel(f'{orb.site_energy: .2f}')
                else:
                    continue

                label.setAlignment(QtCore.Qt.AlignCenter)
                label.setTextFormat(QtCore.Qt.RichText)
                self.setCellWidget(i, j, label)

    def multi_select(self, check_state):
        selected_rows = []
        for selected_range in self.selectedRanges():
            selected_rows.extend(range(selected_range.topRow(), selected_range.bottomRow() + 1))
        
        for row in selected_rows:
            btn = list(self.orbital_checkboxes.values())[row]
            btn.setCheckState(check_state)

    def select_all_button_handler(self, check_state):
        if check_state == QtCore.Qt.CheckState.Unchecked:
            for orb in self.orbitals:
                self.parent.parent.state.set_state(orb, False)

        if check_state == QtCore.Qt.CheckState.Checked:
            for orb in self.orbitals:
                self.parent.parent.state.set_state(orb, None)

        self.reset_checkboxes()

    def reset_checkboxes(self):
        for orb, checkbox in self.orbital_checkboxes.items():
            checkbox.setChecked(self.parent.parent.state.is_enabled(orb))



class OrbitalSelectionTab(QtWidgets.QFrame):
    def __init__(self, parent, orbitals):
        super().__init__(parent=parent)
        self.parent = parent

        layout = QtWidgets.QGridLayout()
        self.setLayout(layout)

        # each tabs has a table
        table = OrbitalSelectionTable(self, orbitals)
        # a select all checkbox
        select_all_button = QtWidgets.QCheckBox('Select All')
        select_all_button.setChecked(True)
        select_all_button.checkStateChanged.connect(table.select_all_button_handler)

        # and a spins and irreps selection button


        layout.addWidget(select_all_button, 0, 0, 1, 2)
        layout.addWidget(table, 1, 0, 1, 2)



class OrbitalSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent, orbs):
        super().__init__(parent=parent)
        self.parent = parent
        self.orbs = orbs

        # state keeps track of allowed and disallowed orbitals
        # also spins and irreps
        self.state = OrbitalSelectionState(orbs)

        layout = QtWidgets.QGridLayout(self)
        self.setLayout(layout)

        # title of the dialogue
        layout.addWidget(QtWidgets.QLabel('Select allowed orbitals:\n'), 0, 0, 1, 3)

        # all buttons and tables go into the tabs
        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs, 1, 0, 1, 3)
        # add a tab for the MOs
        self.tabs.addTab(OrbitalSelectionTab(self, orbs.mos.orbitals), self.parent.parent._ICONS['mo'], 'Complex')
        # and for each fragment
        for fragment in orbs.fragments:
            self.tabs.addTab(OrbitalSelectionTab(self, orbs.sfos.filter(fragment=fragment)), self.parent.parent._ICONS['sfo'], fragment)

        # some standard buttons
        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        reset_btn = QtWidgets.QPushButton('Reset')
        reset_btn.clicked.connect(self.reset)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)
        layout.addWidget(reset_btn, 2, 2, 1, 1)


    def reset(self, tab=None):
        ...


    def apply(self):
        for column, col_btns in self._btns.items():
            for orb, btn in col_btns.items():
                self.parent._orb_selection[column][orb] = btn.isChecked()

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

    def accept(self):
        self.apply()
        self.parent._update_plot()
        self.hide()
