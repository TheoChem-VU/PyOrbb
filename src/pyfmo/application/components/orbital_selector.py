from PySide6 import QtWidgets, QtCore
from pyfmo.application.components import rich_widgets, latex_renderer
import pyfmo
import dictfunc
import functools


class OrbitalSelectionState:
    def __init__(self, orbs):
        self.orbs = orbs
        self.reset()

    def rename(self, old, new):
        self.orbitals[new] = self.orbitals.pop(old)
        self.spins[new] = self.spins.pop(old)
        self.irreps[new] = self.irreps.pop(old)
        self.select_all_buttons[new] = self.select_all_buttons.pop(old)

    def reset(self):
        # to indicate MOs we use the mos object as the key
        self.orbitals = {self.orbs.mos: {}}
        self.spins = {self.orbs.mos: {}}
        self.irreps = {self.orbs.mos: {}}
        self.select_all_buttons = {self.orbs.mos: True}

        for mo in self.orbs.mos:
            self.orbitals[self.orbs.mos][mo] = True
            self.spins[self.orbs.mos].setdefault(mo.spin, True)
            self.irreps[self.orbs.mos].setdefault(mo.symmetry, True)

        for frag in self.orbs.fragments:
            self.orbitals[frag] = {}
            self.spins[frag] = {}
            self.irreps[frag] = {}
            for sfo in self.orbs.sfos.filter(fragment=frag):
                self.orbitals[frag][sfo] = True
                self.spins[frag].setdefault(sfo.spin, True)
                self.irreps[frag].setdefault(sfo.subspecies, True)
            self.select_all_buttons[frag] = True

    def disable_all(self, system=None):
        # to indicate MOs we use the mos object as the key
        if system is self.orbs.mos:
            for mo in self.orbs.mos:
                self.orbitals[self.orbs.mos][mo] = False

        for frag in self.orbs.fragments:
            if system != frag:
                continue

            for sfo in self.orbs.sfos.filter(fragment=frag):
                self.orbitals[frag][sfo] = False

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

    def allowed_mos(self):
        ret = []
        for mo in self.orbs.mos:
            if self.orbitals[self.orbs.mos][mo] is False:
                continue
            ret.append(mo)
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

    def allowed_sfos(self, frag=None):
        ret = []

        for frag_, orbitals in self.orbitals.items():
            if frag is not None and frag_ != frag:
                continue

            if frag_ == self.orbs.mos:
                continue

            for orb, state in orbitals.items():
                if state is False:
                    continue

                ret.append(orb)
        return ret

    def is_enabled(self, orbital):
        if isinstance(orbital, pyfmo.orbitals.objects.MO):
            return self.mo_states()[orbital]
        else:
            return self.sfo_states(orbital.fragment)[orbital]

    def set_state(self, orbital, state):
        if isinstance(orbital, pyfmo.orbitals.objects.MO):
            self.orbitals[self.orbs.mos][orbital] = state
        else:
            self.orbitals[orbital.fragment][orbital] = state

    def copy(self):
        new = OrbitalSelectionState(self.orbs)
        new.orbitals = dictfunc.list_to_dict(dictfunc.dict_to_list(self.orbitals))
        new.spins = dictfunc.list_to_dict(dictfunc.dict_to_list(self.spins))
        new.irreps = dictfunc.list_to_dict(dictfunc.dict_to_list(self.irreps))
        new.select_all_buttons = dictfunc.list_to_dict(dictfunc.dict_to_list(self.select_all_buttons))
        return new

class TableFloatItem(QtWidgets.QTableWidgetItem):
    def __init__(self, value, precision=3):
        super().__init__()
        self.setData(QtCore.Qt.DisplayRole, round(value, precision))
        # self.setText(f'{value:> .{precision}f}')
        # self.setTextAlignment(QtCore.Qt.AlignHCenter)


class OrbitalSelectionTable(QtWidgets.QTableWidget):
    def __init__(self, parent, orbitals):
        self.parent = parent
        self.orbital_checkboxes = {}
        self.orbitals = orbitals
        is_MO = isinstance(orbitals[0], pyfmo.orbitals.objects.MO)

        column_widths = {
            'Orbital': 150 + 15,
            'Spin': 40 + 15,
            'Occ.': 40 + 15,
            'Gross Pop.': 70 + 15,
            'Symm.': 50 + 15,
            'Subsp.': 50 + 15,
            'Energy (reg.)': 80 + 15,
            'Energy (eff.)': 80 + 15,
            'Rel. Name': 130 + 15,
            'Rel. Name (Symm.)': 130 + 15,
            'Rel. Name (Subsp.)': 130 + 15,
        }

        column_tooltips = {
            'Orbital': None,
            'Spin': 'Spin state of this orbital',
            'Occ.': 'Occupation',
            'Gross Pop.': 'Mulliken gross population',
            'Symm.': 'Irreducible representation in the symmetry of the complex',
            'Subsp.': 'Irreducible representation in the symmetry of the fragment',
            'Energy (reg.)': 'Regular orbital energy',
            'Energy (eff.)': 'Effective orbital energy',
            'Rel. Name': 'Relative name over all orbitals',
            'Rel. Name (Symm.)': 'Relative name over all orbitals with the same complex symmetry',
            'Rel. Name (Subsp.)': 'Relative name over all orbitals with the same fragment symmetry',
        }

        headers = list(column_widths)

        if is_MO:
            headers.remove('Subsp.')
            headers.remove('Gross Pop.')
            headers.remove('Rel. Name (Subsp.)')

        if not self.parent.parent.orbs.data['calc_info']['has_site_energy'] or is_MO:
            headers.remove('Energy (eff.)')

        super().__init__(len(orbitals), len(headers), parent=parent)

        self.setItemDelegate(rich_widgets.HTMLDelegate())
        self.verticalHeader().setDefaultSectionSize(35)
        self.verticalHeader().setVisible(False)
        self.setHorizontalHeaderLabels(headers)
        for i, header in enumerate(headers):
            self.setColumnWidth(i, column_widths[header])

        for i, orb in enumerate(orbitals):
            orbital_frame = QtWidgets.QFrame()
            layout_ = QtWidgets.QHBoxLayout()
            orbital_frame.setLayout(layout_)
            orbital_frame.setToolTip(f'ADF-name: {orb}')
            
            self.orbital_checkboxes[orb] = QtWidgets.QCheckBox()
            self.orbital_checkboxes[orb].setToolTip("Check to include this orbital in the analysis")
            self.orbital_checkboxes[orb].setChecked(self.parent.parent.state.is_enabled(orb))
            self.orbital_checkboxes[orb].checkStateChanged.connect(self.multi_select)

            draw_button = QtWidgets.QPushButton()
            draw_button.setIcon(self.parent.parent.parent.parent._ICONS["draw"])
            draw_button.setIconSize(QtCore.QSize(8, 8))
            draw_button.setFixedSize(QtCore.QSize(16, 16))
            draw_button.setToolTip("Click to draw this orbital")
            draw_button.clicked.connect(functools.partial(self.parent.parent.parent.plot.draw_orbital, orb=orb, draw_type='single'))
            layout_.addWidget(draw_button, stretch=0)
            layout_.addWidget(self.orbital_checkboxes[orb])
            layout_.addWidget(QtWidgets.QLabel(pyfmo.generate_label(orb, mode='html')))
            layout_.addStretch()
            self.setCellWidget(i, headers.index('Orbital'), orbital_frame)

            for j, header in enumerate(headers):
                if header == 'Spin':
                    spin = {'A': '<i>α</i>', 'B': '<i>β</i>', 'AB': '<i>αβ</i>'}[orb.spin]
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(spin)
                elif header == 'Occ.':
                    item = TableFloatItem(orb.occupation, 3)
                elif header == 'Gross Pop.':
                    item = TableFloatItem(orb.gross_population, 3)
                elif header == 'Symm.':
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(pyfmo.translate_irrep_label(orb.symmetry, mode='html'))
                elif header == 'Subsp.':
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(pyfmo.translate_irrep_label(orb.subspecies, mode='html'))
                elif header == 'Rel. Name':
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(orb.relative_name)
                elif header == 'Rel. Name (Symm.)':
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(orb.symmetry_relative_name)
                elif header == 'Rel. Name (Subsp.)':
                    item = QtWidgets.QTableWidgetItem()
                    item.setText(orb.subspecies_relative_name)
                elif header == 'Energy (reg.)':
                    item = TableFloatItem(float(orb.energy), 2)
                elif header == 'Energy (eff.)':
                    item = TableFloatItem(float(orb.site_energy), 2)
                else:
                    continue

                item.setToolTip(column_tooltips[header])
                item.setTextAlignment(QtCore.Qt.AlignCenter)
                self.setItem(i, j, item)

        self.setSortingEnabled(True)
        self.sortItems(headers.index("Energy (reg.)"), QtCore.Qt.AscendingOrder)
        self.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers);

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
            self.parent.parent.state.select_all_buttons[self.parent.state_key] = False

        if check_state == QtCore.Qt.CheckState.Checked:
            for orb in self.orbitals:
                self.parent.parent.state.set_state(orb, None)

            self.parent.parent.state.select_all_buttons[self.parent.state_key] = True

        self.reset_checkboxes()

    def reset_checkboxes(self):
        for orb, checkbox in self.orbital_checkboxes.items():
            checkbox.setChecked(self.parent.parent.state.is_enabled(orb))



class OrbitalSelectionTab(QtWidgets.QFrame):
    def __init__(self, parent, orbitals, state_key):
        super().__init__(parent=parent)
        self.parent = parent

        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)

        self.state_key = state_key

        # each tabs has a table
        self.table = OrbitalSelectionTable(self, orbitals)
        # a select all checkbox
        self.select_all_button = QtWidgets.QCheckBox('Select All')
        self.select_all_button.setChecked(True)
        self.select_all_button.checkStateChanged.connect(self.table.select_all_button_handler)

        # and a spins and irreps selection button
        self.spin_selection_dialog = SpinSelectionDialog(self)
        spin_select_button = QtWidgets.QPushButton('Spins')
        spin_select_button.clicked.connect(self.spin_selection_dialog.open)
        self.irrep_selection_dialog = IrrepSelectionDialog(self)
        irrep_select_button = QtWidgets.QPushButton('Irreps')
        irrep_select_button.clicked.connect(self.irrep_selection_dialog.open)

        selection_frame = QtWidgets.QFrame()
        selection_frame_layout = QtWidgets.QHBoxLayout()
        selection_frame.setLayout(selection_frame_layout)
        selection_frame_layout.addWidget(self.select_all_button)
        selection_frame_layout.addWidget(spin_select_button)
        selection_frame_layout.addWidget(irrep_select_button)
        selection_frame_layout.addStretch()

        layout.addWidget(selection_frame, stretch=0)
        layout.addWidget(self.table)

    def update_state(self):
        for orb, checkbox in self.table.orbital_checkboxes.items():
            self.parent.state.orbitals[self.state_key][orb] = checkbox.isChecked()

    def reset(self):
        for orb, state in self.parent.state.orbitals[self.state_key].items():
            self.table.orbital_checkboxes[orb].setChecked(state)

        self.select_all_button.setChecked(self.parent.state.select_all_buttons[self.state_key])


class SpinSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent=parent)
        self.parent = parent
        self.state_key = parent.state_key

        layout = QtWidgets.QGridLayout()
        self.setLayout(layout)

        spin_state = self.parent.parent.state.spins[self.state_key]
        self.checkboxes = {}
        for i, spin in enumerate(spin_state):
            label = {'A': 'α', 'B': 'β', 'AB': 'αβ'}[spin]
            self.checkboxes[spin] = QtWidgets.QCheckBox(label)
            layout.addWidget(self.checkboxes[spin], i, 0, 1, 2)

        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        layout.addWidget(save_btn, i + 1, 0, 1, 1)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(cancel_btn, i + 1, 1, 1, 1)

        self.setup()

    def setup(self):
        spin_state = self.parent.parent.state.spins[self.state_key]
        for spin, checkbox in self.checkboxes.items():
            checkbox.setChecked(spin_state[spin])

    def open(self):
        self.setup()
        super().open()

    def accept(self):
        # propagate state to parent table's checkboxes
        for spin, checkbox in self.checkboxes.items():
            self.parent.parent.state.spins[self.state_key][spin] = checkbox.isChecked()

            for orb, table_checkbox in self.parent.table.orbital_checkboxes.items():
                if orb.spin != spin:
                    continue
                table_checkbox.setChecked(checkbox.isChecked())

        self.hide()


class TabButton(QtWidgets.QFrame):
    def __init__(self, parent, text):
        super().__init__(parent)
        self.parent = parent
        self.text = text

        self.layout = QtWidgets.QHBoxLayout(self)
        self.setLayout(self.layout)

        pixmap = latex_renderer.convert_to_QPixMap(text, darkmode=self.parent.parent.parent.isDarkMode)
        self._lbl = QtWidgets.QLabel(self)
        self._lbl.setPixmap(pixmap)
        self.layout.addWidget(self._lbl)
        self.layout.setContentsMargins(0, 0, 0, 0)


class IrrepSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super().__init__(parent=parent)
        self.parent = parent
        self.state_key = parent.state_key

        layout = QtWidgets.QGridLayout()
        self.setLayout(layout)

        irrep_state = self.parent.parent.state.irreps[self.state_key]
        self.checkboxes = {}
        for i, irrep in enumerate(irrep_state):
            label = pyfmo.translate_irrep_label(irrep, mode='html')
            self.checkboxes[irrep] = rich_widgets.HTMLCheckBox(label)
            layout.addWidget(self.checkboxes[irrep], i, 0, 1, 2)

        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        layout.addWidget(save_btn, i + 1, 0, 1, 1)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(cancel_btn, i + 1, 1, 1, 1)

        self.setup()

    def setup(self):
        irrep_state = self.parent.parent.state.irreps[self.state_key]
        for irrep, checkbox in self.checkboxes.items():
            checkbox.setChecked(irrep_state[irrep])

    def open(self):
        self.setup()
        super().open()

    def accept(self):
        # propagate state to parent table's checkboxes
        for irrep, checkbox in self.checkboxes.items():
            self.parent.parent.state.irreps[self.state_key][irrep] = checkbox.isChecked()

            for orb, table_checkbox in self.parent.table.orbital_checkboxes.items():
                if isinstance(orb, pyfmo.orbitals.objects.MO):
                    if orb.symmetry != irrep:
                        continue
                else:
                    if orb.subspecies != irrep:
                        continue
                table_checkbox.setChecked(checkbox.isChecked())

        self.hide()


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
        layout.addWidget(QtWidgets.QLabel('<b>Select allowed orbitals:</b>\n'), 0, 0, 1, 3)

        # all buttons and tables go into the tabs
        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs, 1, 0, 1, 3)
        self.tab_indices = {}
        # add a tab for the MOs
        self.tab_stor = {'Complex': OrbitalSelectionTab(self, orbs.mos.orbitals, orbs.mos)}

        self.addTab(self.tab_stor['Complex'], self.parent.parent._ICONS['mo'], "Complex")
        self.tab_indices["Complex"] = 0

        # and for each fragment
        for i, fragment in enumerate(orbs.fragments):
            self.tab_indices[fragment] = i + 1
            self.tab_stor[fragment] = OrbitalSelectionTab(self, orbs.sfos.filter(fragment=fragment), fragment)
            self.addTab(self.tab_stor[fragment], self.parent.parent._ICONS['sfo'], fragment)

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

    def addTab(self, widget, icon, text):
        index = self.tabs.addTab(widget, icon, "")
        btn = TabButton(self, text)
        self.tabs.tabBar().setTabButton(index, QtWidgets.QTabBar.RightSide, btn)

    def rename(self, old_name, new_name):
        index = self.tab_indices.pop(old_name)
        new_btn = TabButton(self, new_name)
        self.tabs.tabBar().setTabButton(index, QtWidgets.QTabBar.RightSide, new_btn)
        self.tab_indices[new_name] = index
        self.tab_stor[new_name] = self.tab_stor.pop(old_name)
        if isinstance(self.tab_stor[new_name].state_key, str):
            self.tab_stor[new_name].state_key = new_name
            self.state.rename(old_name, new_name)

    def reset(self, tab=None):
        self.state.reset()
        for system, tab in self.tab_stor.items():
            tab.reset()

    def apply(self):
        for system, tab in self.tab_stor.items():
            tab.update_state()

    def exec(self, *args):
        self.__old_state = self.state.copy()

        for system, tab in self.tab_stor.items():
            tab.reset()

        super().exec()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

    def accept(self):
        self.apply()
        self.parent._update_plot()
        self.hide()

    def reject(self):
        self.state = self.__old_state
        self.hide()
