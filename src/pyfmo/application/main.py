from PySide6 import QtWidgets, QtCore, QtGui
import pyfmo
from .components import (
    orbital_selector,
    rich_widgets,
    latex_renderer,
    action_widget,
    settings,
    editable_tabs,
    column_dragger,
    spoilers,
    carousel,
    shadow,
    theme_switcher,
    )
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backend_tools import Cursors
from matplotlib.figure import Figure
import numpy as np
import os
from math import floor, ceil
from functools import partial
import pyperclip
import platform
import traceback
from scm import plams
import webbrowser
import tcviewer


slider_resolution = 500

def mol2xyz(mol):
    s = ''
    s += str(len(mol.atoms)) + '\n\n'
    for atom in mol:
        s += f'{atom.symbol:2} {atom.coords[0]} {atom.coords[1]} {atom.coords[2]}\n'
    return s


def _determine_formal_charges(orbs):
    # build up the effective charges of the atoms
    # this takes into account the atom number and number of frozen core electrons
    atomtypes = orbs.reader.read('Geometry', 'atomtype').split()
    eff_charges = orbs.reader.read('Geometry', 'atomtype effective charge')

    if isinstance(eff_charges, float):
        eff_charges = [eff_charges]

    if isinstance(atomtypes, float):
        atomtypes = [atomtypes]

    atomtype_charges = {typ: charge for typ, charge in zip(atomtypes, eff_charges)}

    # calculate the charges for the fragments and the complex
    charges = {}
    for frag in orbs.fragments:
        sfos = orbs.sfos.filter(fragment=frag)
        # we need the atoms in the molecule
        mol = sfos[0].molecule
        expected_Nelectrons = sum(atomtype_charges[atom.symbol] for atom in mol)
        actual_Nelectrons = round(sum(sfo.occupation for sfo in sfos))
        charges[frag] = expected_Nelectrons - actual_Nelectrons
        
    charges['Complex'] = sum(charges.values())
    return charges

def _determine_vdd_charges(orbs):
    # build up the effective charges of the atoms
    # this takes into account the atom number and number of frozen core electrons
    atomtypes = np.atleast_1d(orbs.reader.read('Geometry', 'atomtype').split())
    vdd_charges = np.atleast_1d(orbs.reader.read('Properties', 'AtomCharge_SCF Voronoi')) - np.atleast_1d(orbs.reader.read('Properties', 'AtomCharge_initial Voronoi'))

    def get_atom_indices(mol):
        complex_mol = orbs.data['molecules']['complex']
        atomtype_order = np.atleast_1d(orbs.reader.read('Geometry', 'atom order index'))
        atomtype_order = atomtype_order[len(atomtype_order)//2:]
        indices = []
        for atom in mol:
            for i, atom_ref in enumerate(complex_mol):
                if atom.coords == atom_ref.coords:
                    indices.append(atomtype_order[i] - 1)
        return indices

    # calculate the charges for the fragments and the complex
    charges = {}
    for frag in orbs.fragments:
        mol = orbs.data['molecules'][frag]
        indices = get_atom_indices(mol)
        charges[frag] = sum(vdd_charges[i] for i in indices)
        
    charges['Complex'] = sum(charges.values())
    return charges



class ScrollLabel(QtWidgets.QScrollArea):
    # constructor
    def __init__(self, text='', *args, **kwargs):
        QtWidgets.QScrollArea.__init__(self, *args, **kwargs)

        # making widget resizable
        self.setWidgetResizable(True)
        self.horizontalScrollBar().setEnabled(False)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        # making qwidget object
        content = QtWidgets.QWidget(self)
        self.setWidget(content)

        # vertical box layout
        lay = QtWidgets.QVBoxLayout(content)

        # creating label
        self.label = QtWidgets.QLabel(content)

        # setting alignment to the text
        self.label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.label.setStyleSheet('padding: 10px; font-size: 10px;')
        self.label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

        # adding label to the layout
        lay.addWidget(self.label)

        self.setText(text)

    # the setText method
    def setText(self, text):
        # setting text to the label
        self.label.setText(text)


class ETypeDialog(QtWidgets.QDialog):
    def __init__(self, parent, possibilities, selected):
        super().__init__(parent=parent)
        self.parent = parent
        self.possibilities = possibilities
        self.selected = selected
        layout = QtWidgets.QGridLayout(self)

        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('<b>Select energy type for SFOs:</b>\n'), 0, 0, 1, 2)
        rbtn_layout = QtWidgets.QGridLayout()
        rbtn_group = QtWidgets.QButtonGroup()
        layout.addLayout(rbtn_layout, 1, 0, 1, 2)

        display_names = {
            'energy': 'Regular',
            'approx_site_energy': 'Effective (approximate)',
            'site_energy': 'Effective',
            'site_energy_SCF0': 'Effective (initial density)',
        }

        self.rbuttons = {}
        for i, pos in enumerate(self.possibilities):
            self.rbuttons[pos] = QtWidgets.QRadioButton(display_names[pos])
            if pos == self.selected:
                self.rbuttons[pos].setChecked(True)
            rbtn_layout.addWidget(self.rbuttons[pos], i, 0, 1, 1)
        
        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        shadow.apply(save_btn)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        shadow.apply(cancel_btn)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

    def accept(self):
        for key, val in self.rbuttons.items():
            if val.isChecked():
                self.parent._energytype_selection = key
        self.parent._update_plot()
        self.hide()


class FragRenameDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        layout = QtWidgets.QGridLayout(self)
        self.setLayout(layout)
        title = QtWidgets.QLabel('<b>Rename orbital column name</b>')
        layout.addWidget(title, 0, 0, 1, 2)
        self._frag_rename_textedit = QtWidgets.QLineEdit(self)
        layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        shadow.apply(save_btn)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        shadow.apply(cancel_btn)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, default=''):
        self._text = default
        self._frag_rename_textedit.setText(default)
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        self.accepted.connect(self.get_text)
        loop.exec()

        return self._text

    def get_text(self):
        self._text = self._frag_rename_textedit.text()


class YAxisDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        layout = QtWidgets.QGridLayout(self)
        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('<b>Choose Y-axis limits</b>'), 0, 0, 1, 2)

        label = QtWidgets.QLabel("Min:")
        layout.addWidget(label, 1, 0, 1, 1)
        self._ylim_low_textedit = QtWidgets.QLineEdit(self)
        layout.addWidget(self._ylim_low_textedit, 1, 1, 1, 1)

        label = QtWidgets.QLabel("Max:")
        layout.addWidget(label, 2, 0, 1, 1)
        self._ylim_high_textedit = QtWidgets.QLineEdit(self)
        layout.addWidget(self._ylim_high_textedit, 2, 1, 1, 1)

        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        shadow.apply(save_btn)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        shadow.apply(cancel_btn)
        layout.addWidget(save_btn, 4, 0, 1, 1)
        layout.addWidget(cancel_btn, 4, 1, 1, 1)

        _reset_btn = QtWidgets.QPushButton("Reset Y-axis")
        _reset_btn.clicked.connect(self.reset)
        _reset_btn.clicked.connect(self.reject)
        shadow.apply(_reset_btn)
        layout.addWidget(_reset_btn, 3, 0, 1, 1)

    def reset(self):
        self._tup = None

    def open(self, default):
        self._tup = default
        self._ylim_low_textedit.setText(str(round(default[0], 3)))
        self._ylim_high_textedit.setText(str(round(default[1], 3)))
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        self.accepted.connect(self.get_tup)
        loop.exec()

        return self._tup

    def get_tup(self):
        self._tup = float(self._ylim_low_textedit.text()), float(self._ylim_high_textedit.text())


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=9, height=6.5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.__dragger = column_dragger.Dragger(self.fig, move_callback=self.set_xtick_order, release_callback=self.set_xtick_order)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(top=1, right=1, bottom=0.1, left=0.12)
        self.parent = parent
        super().__init__(self.fig)

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click)
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self._selected_orbitals = []
        self._frag_rename_dialog = FragRenameDialog(self)
        self._yaxis_dialog = YAxisDialog(self)
        self._already_unfaded = True
        self.previous_mouse_pos = None
        self._add_warning = False

    def add_warning(self):
        self._add_warning = True

    def set_xtick_order(self, order):
        self.parent._xtick_order = {tick.get_text(): float(pos) for tick, pos in order.items()}
        self.parent._update_plot()

    def draw_orbital(self, orb=None, draw_type='single'):
        if hasattr(orb, '__len__'):
            kfpath = orb[0].parent.parent.kfpath
        else:
            kfpath = orb.parent.parent.kfpath

        if ' ' in kfpath:
            QtWidgets.QMessageBox.critical(self, 'Error', 'The adf.rkf path contains a space. We will not be able to run Densf properly.\nPlease move the file to a different location.')
            return

        import tcviewer

        amsloc = self.parent.parent.settings_dialog.get("Densf", "General", "amsbin")
        if platform.system() == "Windows":
            preambles = [f'set AMSHOME={os.path.split(self.parent.parent.settings_dialog.get("Densf", "General", "amsbin"))[0]}', f'set AMSBIN={os.path.split(self.parent.parent.settings_dialog.get("Densf", "General", "amsbin"))[0]}/bin']
        else:
            amsbashrc = os.path.join(amsloc, 'Contents', 'Resources', 'amshome', 'amsbashrc.sh')
            preambles = [f'source {amsbashrc}']

        if self.parent.tcviewer_screen is None or self.parent.tcviewer_screen.isclosed:
            self.parent.tcviewer_screen = tcviewer.screen._ScreenWindow()
            self.parent.tcviewer_screen.setWindowIcon(self.parent.parent._ICONS['pyorbb'])
            self.parent.tcviewer_screen.setWindowTitle('PyOrbb Viewer')
            self.parent.tcviewer_screen.show()

            app = QtWidgets.QApplication.instance()
            app.windows.append(self.parent.tcviewer_screen)

        gridsize = self.parent.parent.settings_dialog.get("Densf", "General", "gridsize")
        with self.parent.tcviewer_screen.add_molscene() as scene:
            if draw_type == 'single':
                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb.occupied else ([1, .5, 0], [0, 1, 1])
                scene.draw_molecule(orb.molecule)
                try:
                    data = orb.vtk_file(gridsize=gridsize, preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                scene.draw_dual_isosurface(data, colorm=c1, colorp=c2)
                scene.draw_text(str(orb))

            if draw_type == 'sum':
                mol = orb[0].molecule + orb[1].molecule
                overlap_sign = np.sign(orb[0] @ orb[1])
                if overlap_sign == 0:
                    overlap_sign = 1
                scene.draw_molecule(mol)

                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb[0].occupied else ([1, .5, 0], [0, 1, 1])
                try:
                    cub = orb[0].vtk_file(gridsize=gridsize, preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                scene.draw_dual_isosurface(cub, colorm=c1, colorp=c2)

                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb[1].occupied else ([1, .5, 0], [0, 1, 1])
                try:
                    cub = orb[1].vtk_file(gridsize=gridsize, preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                cub.values *= overlap_sign
                scene.draw_dual_isosurface(cub, colorm=c1, colorp=c2)
                
                scene.draw_text(str(orb[0]) + ' & ' + str(orb[1]))

            if draw_type == 'overlap':
                c1, c2 = [0, 1, 0], [1, 0, 1]
                mol = orb[0].molecule + orb[1].molecule
                scene.draw_molecule(mol)
                try:
                    cub1 = orb[0].vtk_file(gridsize=gridsize, preambles=preambles, grid_around_mol=mol)
                    cub2 = orb[1].vtk_file(gridsize=gridsize, preambles=preambles, grid_around_mol=mol)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                cub1.values *= cub2.values
                scene.draw_dual_isosurface(cub1, colorm=c1, colorp=c2, isovalue=0.03**2)
                
                scene.draw_text(str(orb[0]) + ' * ' + str(orb[1]))

    def draw_molecule(self, mol=None):
        print('before loaded tcviewer')
        import tcviewer
        print('loaded tcviewer')
        # get or make a new viewer
        if self.parent.tcviewer_screen is None or self.parent.tcviewer_screen.isclosed:
            print('making a screen')
            self.parent.tcviewer_screen = tcviewer.screen._ScreenWindow()
            self.parent.tcviewer_screen.setWindowIcon(self.parent.parent._ICONS['pyorbb'])
            self.parent.tcviewer_screen.setWindowTitle('PyOrbb Viewer')
            self.parent.tcviewer_screen.show()

            app = QtWidgets.QApplication.instance()
            app.windows.append(self.parent.tcviewer_screen)

        print('made a screen')
        with self.parent.tcviewer_screen.add_molscene() as scene:
            print('drawing mol')
            scene.draw_molecule(mol)

    def _set_orbital_info_box(self):
        self.parent.orbital_info_box.empty()
        for i, orb in enumerate(self._selected_orbitals):
            s = ''
            if isinstance(orb, pyfmo.orbitals.objects.SFO):
                icon = self.parent.parent._ICONS['sfo']
                submixes = self.parent.main_mix.split()
                submix = [submix for submix in submixes if orb in submix.sfos][0]
                s += f'Name         {pyfmo.generate_label(orb, mode="html", use_formatting=False)} ({orb.relative_name})'
                s += f'\nSymm.        {pyfmo.translate_irrep_label(orb.symmetry, mode="html", use_formatting=False)} ({orb.symmetry_relative_name})'
                s += f'\nSubsp.       {pyfmo.translate_irrep_label(orb.subspecies, mode="html", use_formatting=False)} ({orb.subspecies_relative_name})'
                s += f'\nFragment     {orb.fragment}'
                s += f'\nEnergy      {getattr(orb, self.parent._energytype_selection): .2f} eV'
                s += f'\nOccupation  {orb.occupation: .2f}'
                s += f'\nPop.        {orb.gross_population: .3f}'
                s += f'\nSpin-pop.   {orb.gross_spin: .3f}'
                s += f'\nSpin         {orb.spin}'
                s += f'\nIrrep        {orb.symmetry}'

                s += '\n\nSFO                      S   dE (eV)'
                s += '\n─────────────────── ────── ─────────'
                for sfo2 in sorted(submix.sfos, key=lambda sfo_: -abs(orb @ sfo_)):
                    if sfo2 == orb:
                        continue
                    if sfo2.fragment == orb.fragment:
                        continue
                    s += f'\n{str(sfo2):19.19} {orb @ sfo2: 5.3f} {abs(getattr(orb, self.parent._energytype_selection) - getattr(sfo2, self.parent._energytype_selection)): 8.2f}'
                
                s += '\n\nMO                     Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for mo in sorted(submix.mos, key=lambda mo: -abs(orb.mulliken_contribution(mo))):
                    s += f'\n{str(mo):19.19} {orb.mulliken_contribution(mo): 8.2%} {orb.coefficient(mo): 7.4f}'
                title = f"{orb.fragment}({pyfmo.generate_label(orb, mode='latex')})"

            if isinstance(orb, pyfmo.orbitals.objects.MO):
                icon = self.parent.parent._ICONS['mo']
                submixes = self.parent.main_mix.split()
                submix = [submix for submix in submixes if orb in submix.mos][0]
                s += f'Name         {pyfmo.generate_label(orb, mode="html", use_formatting=False)} ({orb.relative_name})'
                s += f'\nSymm.        {pyfmo.translate_irrep_label(orb.symmetry, mode="html", use_formatting=False)} ({orb.symmetry_relative_name})'
                s += f'\nEnergy      {orb.energy: .2f} eV'
                s += f'\nOccupation  {orb.occupation: .2f}'
                s += f'\nSpin         {orb.spin}'
                s += f'\nIrrep        {orb.symmetry}'

                s += '\n\nSFO                    Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for sfo in sorted(submix.sfos, key=lambda sfo: -abs(sfo.mulliken_contribution(orb))):
                    s += f'\n{str(sfo):19.19} {sfo.mulliken_contribution(orb): 8.2%} {sfo.coefficient(orb): 7.4f}'
                title = pyfmo.generate_label(orb, mode='latex')

            if isinstance(orb, tuple) and isinstance(orb[0], pyfmo.orbitals.objects.SFO) and isinstance(orb[1], pyfmo.orbitals.objects.MO):
                sfo, mo = orb
                icon = self.parent.parent._ICONS['contribution']
                connected_sfos = [conn[0] for conn in self.parent.main_mix.connections if conn[1] == mo and conn[0].fragment != sfo.fragment]
                s += 'SFO'
                s += f'\n  Name     {pyfmo.generate_label(sfo, mode="html", use_formatting=False)} ({sfo.relative_name})'
                s += f'\n  Symm.    {pyfmo.translate_irrep_label(sfo.symmetry, mode="html", use_formatting=False)} {sfo.symmetry_relative_name}\n'
                s += '\nMO'
                s += f'\n  Name     {pyfmo.generate_label(mo, mode="html", use_formatting=False)} ({mo.relative_name})'
                s += f'\n  Symm.    {pyfmo.translate_irrep_label(mo.symmetry, mode="html", use_formatting=False)} {mo.symmetry_relative_name}\n'
                s += f'\nContr.    {sfo.mulliken_contribution(mo): .2%}'
                s += f'\nCoeff.    {sfo.coefficient(mo): .6f}'
                s += f'\nSpin       {sfo.spin}'
                s += f'\nIrrep      {sfo.symmetry}'
                s += '\n\nSecond SFO          Bonding?'
                s += '\n─────────────────── ────────'
                for sfo2 in connected_sfos:
                    is_bonding = ((sfo @ sfo2) * sfo.coefficient(mo) * sfo2.coefficient(mo)) >= 0
                    s += f'\n{str(sfo2):19.19} {"   Yes  " if is_bonding else "    No    "}'
                orb = f'{sfo} ⇒ {mo}'
                title = f"{sfo.fragment}({pyfmo.generate_label(sfo, mode='latex')}) ⇒ {pyfmo.generate_label(mo, mode='latex')}"

            if isinstance(orb, tuple) and isinstance(orb[0], pyfmo.orbitals.objects.SFO) and isinstance(orb[1], pyfmo.orbitals.objects.SFO):
                sfo1, sfo2 = orb
                self.parent._set_mix_settings()
                mix = self.parent.main_mix.find_two_mixing(sfo1, sfo2)
                int_type = list(mix.connection_type.values())[0]
                mo1, mo2 = mix.mos
                icon = self.parent.parent._ICONS['mix']

                if mo2.energy > mo1.energy:
                    mo1, mo2 = mo2, mo1

                name_sfo1 = f'{pyfmo.generate_label(sfo1, mode="html", use_formatting=False)} ({sfo1.relative_name})'
                name_sfo2 = f'{pyfmo.generate_label(sfo2, mode="html", use_formatting=False)} ({sfo2.relative_name})'
                name_mo1 = f'{pyfmo.generate_label(mo1, mode="html", use_formatting=False)} ({mo1.relative_name})'
                name_mo2 = f'{pyfmo.generate_label(mo2, mode="html", use_formatting=False)} ({mo2.relative_name})'
                if int_type == 'PR':
                    s += '\nPauli Repulsive Interaction\n───────────────────────────\n'
                else:
                    s += '\nOrbital Interaction\n───────────────────\n'
                s += '\n' + name_mo1.center(23)
                s += f'\n        ╱       ╲'
                s += f'\n       ╱         ╲'
                s += f'\n      ╱           ╲'
                s += '\n' + name_sfo1 + ' ' * (18 - len(name_sfo1)) + name_sfo2
                s += f'\n      ╲           ╱'
                s += f'\n       ╲         ╱'
                s += f'\n        ╲       ╱'
                s += '\n' + name_mo2.center(23)
                s += f'\n\n𝛙i          {pyfmo.generate_label(sfo1, mode="html", use_formatting=False)}'
                s += f'\n𝛙j          {pyfmo.generate_label(sfo2, mode="html", use_formatting=False)}'
                s += f'\nΨk          {pyfmo.generate_label(mo1, mode="html", use_formatting=False)}'
                s += f'\nΨl          {pyfmo.generate_label(mo2, mode="html", use_formatting=False)}'

                S = sfo1 @ sfo2
                s += f'\n\nSij         {S: 5.3f}'
                max_pop = 1 if self.parent.orbs.data['calc_info']['unrestricted_sfos'] else 2
                if int_type == 'OI':
                    de = abs(getattr(sfo1, self.parent._energytype_selection) - getattr(sfo2, self.parent._energytype_selection))
                    s += f'\n|εi-εj|      {de:.2f} eV'
                    s += f'\nS^2/|εi-εj|  {S**2 / de:.4f} eV⁻¹'
                    dpi = -(sfo1.occupation - sfo1.gross_population)
                    s += f'\nΔpi         {dpi: 5.3f} e⁻'
                    dpj = -(sfo2.occupation - sfo2.gross_population)
                    s += f'\nΔpj         {dpj: 5.3f} e⁻'
                    ptot = sfo1.gross_population + sfo2.gross_population
                    poi = min(max_pop, ptot)
                    ppr = max(0, ptot - poi)

                    s += f'\npoi          {poi:5.3f} e⁻'
                    s += f'\nppr          {ppr:5.3f} e⁻'

                    s += f'\n\nRoi         {-abs(dpi*dpj)*(poi - ppr) * S**2/de:5.2e} ({mix.fraction:.2%})'
                else:
                    oi = sfo1.occupation
                    oj = sfo2.occupation
                    s += f'\nOpr          {max(oi+oj - max_pop, 0):5.3f} e⁻'
                    s += f'\n\nRpr          {max(oi+oj - max_pop, 0) * S**2:5.2e} ({mix.fraction:.2%})'

                s += f'\n\nCik         {sfo1.mulliken_contribution(mo1): 5.3%}'
                s += f'\nCil         {sfo1.mulliken_contribution(mo2): 5.3%}'
                s += f'\nCjk         {sfo2.mulliken_contribution(mo1): 5.3%}'
                s += f'\nCjl         {sfo2.mulliken_contribution(mo2): 5.3%}'
                s += f'\nCik⋅Cil⋅Cjk⋅Cjl {sfo1.mulliken_contribution(mo1)*sfo1.mulliken_contribution(mo2)*sfo2.mulliken_contribution(mo1)*sfo2.mulliken_contribution(mo2): 5.3%}'

                title = f"{sfo1.fragment}({pyfmo.generate_label(sfo1, mode='latex')}) ± {sfo2.fragment}({pyfmo.generate_label(sfo2, mode='latex')})"

            label = QtWidgets.QLabel(s)
            label.setStyleSheet('font: 10px "IBM Plex Mono";')
            self.parent.orbital_info_box.addSpoiler(title, label, icon)

    def on_plot_click(self, event):
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)

        if event.dblclick:
            if any(artist.contains(event)[0] for artist in self.axes.get_yticklabels()):
                self.parent.ylim = self._yaxis_dialog.open(self.axes.get_ylim())
                self.parent._update_plot()

            # self.parent.new_tick_labels = []
            for artist in self.axes.get_xticklabels():
                if not artist.contains(event)[0]:
                    # self.parent.new_tick_labels.append(artist.get_text())
                    continue
                new_txt = self._frag_rename_dialog.open(artist.get_text())

                self.parent.system_info_box.renameSpoiler(artist.get_text(), new_txt)
                self.parent._orb_selection_dialog.rename(artist.get_text(), new_txt)
                self.parent.orbs.rename_fragment(artist.get_text(), new_txt)
                self.parent._xtick_order[new_txt] = self.parent._xtick_order.pop(artist.get_text())

                if artist.is_MO:
                    self.parent.parent.settings_dialog.set('Plot', 'Levels', 'mo_column_name', new_txt)

            self.parent._update_plot()
            self.fig.canvas.draw_idle()
            self.__dragger.mouse_held = False

        # Iterating over each data member plotted
        lines = self.axes.get_children()
        lines = sorted(lines, key=lambda line: -line.zorder)

        shift_is_held = QtGui.QGuiApplication.instance().keyboardModifiers() == QtCore.Qt.KeyboardModifier.ShiftModifier
        if not shift_is_held:
            self._selected_orbitals = []

        for curve in lines:
            gid = curve.get_gid()
            if gid is None:
                continue 

            # Searching which data member corresponds to current mouse position
            if not curve.contains(event)[0]:
                continue

            s = ''
            if gid.startswith('MO_'):
                mo = self.parent.orbs.mos.orbitals[int(gid[3:])]
                if isinstance(mo, list):
                    mo = mo[0]

                if mo not in self._selected_orbitals:
                    self._selected_orbitals.append(mo)
                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False

                self.fig.canvas.draw_idle()

                break

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[4:])]
                if isinstance(sfo, list):
                    sfo = sfo[0]

                if sfo not in self._selected_orbitals:
                    self._selected_orbitals.append(sfo)

                to_add = []
                for orb in self._selected_orbitals:
                    if not isinstance(orb, pyfmo.orbitals.objects.SFO):
                        continue
                    if orb is sfo:
                        continue

                    mix = self.parent.main_mix.find_two_mixing(sfo, orb)
                    if mix is None:
                        continue

                    to_add.extend([(sfo, orb)])
                    to_add.extend(mix.mos)
                self._selected_orbitals.extend(to_add)

                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False
                self.fig.canvas.draw_idle()

                break

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[4:].split('->')[0].strip())]
                mo = self.parent.orbs.mos.orbitals[int(gid[4:].split('->')[1].strip())]
                
                if (sfo, mo) not in self._selected_orbitals:
                    self._selected_orbitals.append((sfo, mo))

                if sfo not in self._selected_orbitals:
                    self._selected_orbitals.append(sfo)

                if mo not in self._selected_orbitals:
                    self._selected_orbitals.append(mo)

                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False
                self.fig.canvas.draw_idle()

                break

        else:
            if not self._already_unfaded:
                self._already_unfaded = True
                self._unfade()
                self.fig.canvas.set_cursor(Cursors.POINTER)
                self.fig.canvas.draw_idle()

        if event.button == 3:
            self.parent.ylim = None
            self.previous_mouse_pos = None
            self.parent._update_plot()

        self._set_orbital_info_box()

        # set up the menu for the pushbutton
        menu = QtWidgets.QMenu(self)
        for orb in self._selected_orbitals:
            if isinstance(orb, tuple):
                continue

            if isinstance(orb, pyfmo.orbitals.objects.SFO):
                icon = self.parent.parent._ICONS['sfo']
                orb_lab = f'{orb.fragment}({pyfmo.generate_label(orb, mode="latex")})'
            else:
                icon = self.parent.parent._ICONS['mo']
                orb_lab = pyfmo.generate_label(orb, mode='latex')

            func = partial(self.draw_orbital, orb=orb, draw_type='single')
            widg = action_widget.DrawAction(self, orb_lab, icon, func)
            menu.addAction(widg)

        # add the overlap actions
        for i, orb in enumerate(self._selected_orbitals):
            if isinstance(orb, pyfmo.orbitals.objects.MO):
                continue

            if isinstance(orb, tuple):
                continue

            for orb2 in self._selected_orbitals[i+1:]:
                if isinstance(orb2, pyfmo.orbitals.objects.MO):
                    continue

                if isinstance(orb2, tuple):
                    continue

                if orb.fragment == orb2.fragment:
                    continue

                icon = self.parent.parent._ICONS['overlap']

                func = partial(self.draw_orbital, orb=[orb, orb2], draw_type='overlap')
                orb_lab1 = pyfmo.generate_label(orb, mode='latex')
                orb_lab2 = pyfmo.generate_label(orb2, mode='latex')
                widg = action_widget.DrawAction(self, f'{orb.fragment}({orb_lab1}) * {orb2.fragment}({orb_lab2})', icon, func)
                menu.addAction(widg)

                # action = QtGui.QAction(f'{orb} * {orb2}', self)
                # action.setIconVisibleInMenu(True)
                # action.setIcon(icon)
                # action.triggered.connect()
                # menu.addAction(action)

        # add the sum actions
        for i, orb in enumerate(self._selected_orbitals):
            if isinstance(orb, pyfmo.orbitals.objects.MO):
                continue

            if isinstance(orb, tuple):
                continue

            for orb2 in self._selected_orbitals[i+1:]:
                if isinstance(orb2, pyfmo.orbitals.objects.MO):
                    continue

                if isinstance(orb2, tuple):
                    continue

                if orb.fragment == orb2.fragment:
                    continue

                icon = self.parent.parent._ICONS['sum']
                # action = QtGui.QAction(f'{orb}, {orb2}', self)
                # action.setIconVisibleInMenu(True)
                # action.setIcon(icon)
                # action.triggered.connect(partial(self.draw_orbital, orb=[orb, orb2], draw_type='sum'))
                # menu.addAction(action)

                func = partial(self.draw_orbital, orb=[orb, orb2], draw_type='sum')
                orb_lab1 = pyfmo.generate_label(orb, mode='latex')
                orb_lab2 = pyfmo.generate_label(orb2, mode='latex')
                widg = action_widget.DrawAction(self, f'{orb.fragment}({orb_lab1}), {orb2.fragment}({orb_lab2})', icon, func)
                menu.addAction(widg)


        self.parent.orbital_draw_button.setMenu(menu)
        self.parent.orbital_draw_button.setEnabled(len(menu.actions()) > 0)
        self.parent.orbital_filter_button.setEnabled(len(self._selected_orbitals) > 0 or self.parent.orbital_filter_button.isChecked())


    def on_plot_hover(self, event):
        # Iterating over each data member plotted
        lines = self.axes.get_children()
        lines = sorted(lines, key=lambda line: -line.zorder)
        for curve in lines:
            gid = curve.get_gid()
            if gid is None:
                continue 

            # Searching which data member corresponds to current mouse position
            if not curve.contains(event)[0]:
                continue

            if gid.startswith('MO_') or gid.startswith('SFO_') or gid.startswith('MIX_'):
                self.fig.canvas.set_cursor(Cursors.HAND)
                break
        else:
            self.fig.canvas.set_cursor(Cursors.POINTER)

        # we have to use QT's mousebutton detection systems instead of pyplot as
        # they can be inconsistent
        if QtWidgets.QApplication.instance().mouseButtons() == QtCore.Qt.LeftButton:
            if self.parent.ylim is None:
                self.parent.ylim = self.axes.get_ylim()

            mouse_pos = event.ydata
            if self.previous_mouse_pos is not None and mouse_pos is not None:
                dy = mouse_pos - self.previous_mouse_pos
                self.parent.ylim = self.parent.ylim[0] - dy, self.parent.ylim[1] - dy
                # self.parent._update_plot()
                self.axes.set_ylim(self.parent.ylim)
                self.axes.redraw_in_frame()
                self.fig.canvas.draw_idle()

                self.previous_mouse_pos = mouse_pos - dy
            else:
                self.previous_mouse_pos = mouse_pos

        elif event.button is None:
            self.previous_mouse_pos = None


    def _unfade(self):
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)

        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue            

            if not hasattr(artist, 'orig_color'):
                try:
                    artist.orig_color = artist.get_color()
                except AttributeError:
                    artist.orig_color = artist.get_fc()
            if not hasattr(artist, 'orig_alpha'):
                artist.orig_alpha = artist.get_alpha() or 1

            artist.set_color('white')
            artist.set_alpha(1)
            self.axes.draw_artist(artist)
            self.fig.canvas.draw_idle()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha)
            self.axes.draw_artist(artist)

            self.fig.canvas.draw_idle()

    def _fade_unrelated_ints(self, orbs):
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)
        faded_artists = []
        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue

            if not hasattr(artist, 'orig_color'):
                try:
                    artist.orig_color = artist.get_color()
                except AttributeError:
                    artist.orig_color = artist.get_fc()

            if not hasattr(artist, 'orig_alpha'):
                artist.orig_alpha = artist.get_alpha() or 1

            if gid.startswith('MO_'):
                mo = self.parent.orbs.mos.orbitals[int(gid[3:])]
                if mo in orbs:
                    continue

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[4:])]
                if sfo in orbs:
                    continue

            if gid.startswith('ARROWMO_'):
                mo = self.parent.orbs.mos.orbitals[int(gid[8:])]
                if mo in orbs:
                    continue

            if gid.startswith('ARROWSFO_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[9:])]
                if sfo in orbs:
                    continue

            if gid.startswith('TEXTMO_'):
                mo = self.parent.orbs.mos.orbitals[int(gid[7:])]
                if mo in orbs:
                    continue

            if gid.startswith('TEXTSFO_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[8:])]
                if sfo in orbs:
                    continue

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos.orbitals[int(gid[4:].split('->')[0].strip())]
                mo = self.parent.orbs.mos.orbitals[int(gid[4:].split('->')[1].strip())]
                if mo in orbs and sfo in orbs:
                    continue

            faded_artists.append(artist)

        for artist in faded_artists:
            artist.set_color('white')
            artist.set_alpha(1)
            self.axes.draw_artist(artist)
            self.fig.canvas.draw_idle()

            artist.set_color(artist.orig_color)
            artist.set_alpha(0.05)
            self.axes.draw_artist(artist)
            self.fig.canvas.draw_idle()

        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue
            if artist in faded_artists:
                continue

            artist.set_color('white')
            artist.set_alpha(1)
            self.axes.draw_artist(artist)
            self.fig.canvas.draw_idle()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha)
            self.axes.draw_artist(artist)

        self.fig.canvas.draw_idle()

    def on_scroll(self, event):
        if self.parent.ylim is None:
            self.parent.ylim = self.axes.get_ylim()
        dy = self.parent.ylim[1] - self.parent.ylim[0]
        mousey = event.ydata
        f = (mousey - self.parent.ylim[0]) / dy
        change = event.step * dy * 0.001
        new_dy = dy + change
        self.parent.ylim = (mousey - new_dy * f, mousey + new_dy * (1 - f))
        self.axes.set_ylim(self.parent.ylim)
        self.axes.redraw_in_frame()
        self.fig.canvas.draw_idle()



class SaveFileDialog(QtWidgets.QFileDialog):
    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.setFileMode(QtWidgets.QFileDialog.AnyFile)
        self.setFilter(QtCore.QDir.Filter.Files)

    def open(self, title, filters, slot=QtCore.SLOT("get_file_from_dialog()")):
        self.setWindowTitle(title)
        self.setNameFilters(filters)
        super().open(self.parent, slot=slot)
        return self.selectedFiles()[0]


class CopyLabel(QtWidgets.QFrame):
    def __init__(self, text, copy_text=None, icon=None):
        super().__init__()

        self.copy_text = copy_text
        if copy_text is None:
            self.copy_text = text

        layout = QtWidgets.QHBoxLayout()
        self.setLayout(layout)
        label = QtWidgets.QLabel(text)
        label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        label.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        layout.addWidget(label)
        icon = QtGui.QIcon(os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons', 'copy.png'))
        copy_button = QtWidgets.QPushButton(icon, '')
        copy_button.setToolTip('Copy')
        copy_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        copy_button.clicked.connect(self.copy)

        layout.addWidget(copy_button)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def copy(self):
        pyperclip.copy(self.copy_text)


class MoleculeLabel(QtWidgets.QFrame):
    def __init__(self, parent, text, molecule, copy_text=None):
        super().__init__()
        self.parent = parent
        self.molecule = molecule
        self.copy_text = copy_text
        if copy_text is None:
            self.copy_text = text

        layout = QtWidgets.QHBoxLayout()
        self.setLayout(layout)
        label = QtWidgets.QLabel(text)
        label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        label.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        layout.addWidget(label)
        copy_icon = QtGui.QIcon(os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons', 'copy.png'))
        copy_button = QtWidgets.QPushButton(copy_icon, '')
        copy_button.setToolTip('Copy')
        copy_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        copy_button.clicked.connect(self.copy)

        layout.addWidget(copy_button)

        draw_icon = QtGui.QIcon(os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons', 'draw.png'))
        draw_button = QtWidgets.QPushButton(draw_icon, '')
        draw_button.setToolTip('Draw')
        draw_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        draw_button.clicked.connect(self.draw)

        layout.addWidget(draw_button)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def copy(self):
        pyperclip.copy(self.copy_text)

    def draw(self):
        print(self.parent.plot.draw_molecule(self.molecule))


class ValueSlider(QtWidgets.QSlider):
    '''
    Taken from https://stackoverflow.com/questions/56694631/how-to-reverse-the-ticklabels-for-a-qslider-element-and-how-to-make-it-clickable
    by user "S. Nick".
    '''
    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            val = self.pixelPosToRangeValue(event.pos())
            self.setValue(val)
        super().mousePressEvent(event) 

    def pixelPosToRangeValue(self, pos):
        opt = QtWidgets.QStyleOptionSlider()
        self.initStyleOption(opt)
        gr = self.style().subControlRect(QtWidgets.QStyle.CC_Slider, opt, QtWidgets.QStyle.SC_SliderGroove, self)
        sr = self.style().subControlRect(QtWidgets.QStyle.CC_Slider, opt, QtWidgets.QStyle.SC_SliderHandle, self)

        if self.orientation() == QtCore.Qt.Horizontal:
            sliderLength = sr.width()
            sliderMin = gr.x()
            sliderMax = gr.right() - sliderLength + 1
        else:
            sliderLength = sr.height()
            sliderMin = gr.y()
            sliderMax = gr.bottom() - sliderLength + 1;
        pr = pos - sr.center() + sr.topLeft()
        p = pr.x() if self.orientation() == QtCore.Qt.Horizontal else pr.y()
        return QtWidgets.QStyle.sliderValueFromPosition(self.minimum(), self.maximum(), p - sliderMin,
                                               sliderMax - sliderMin, opt.upsideDown)



class AnalysisWindow(QtWidgets.QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.setAcceptDrops(True)
        self.open_rkf_filedialog = QtWidgets.QFileDialog(self)
        self.open_rkf_filedialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        self.open_rkf_filedialog.setWindowTitle('Select adf.rkf file')
        self.open_rkf_filedialog.setFilter(QtCore.QDir.Filter.Files)
        self.open_rkf_filedialog.setNameFilters({"RKF file (*.rkf)", "Any file (*)"})

        self._loaded_analysis = False

        self.errordialog = QtWidgets.QErrorMessage(self)
        self.central_layout = QtWidgets.QVBoxLayout(self)

        self.new_tick_labels = None
        self.ylim = None
        self._xtick_order = None

        self.tcviewer_screen = None

    # The following three methods set up dragging and dropping for the app
    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dragMoveEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        """
        Drop files directly onto the widget
        File locations are stored in fname
        :param e:
        :return:
        """
        if e.mimeData().hasUrls:
            e.setDropAction(QtCore.Qt.CopyAction)
            e.accept()
            for url in e.mimeData().urls():
                fname = str(url.toLocalFile())
                if not self._loaded_analysis:
                    self.load_analysis(fname)
                    self.parent.tabs.setTabText(self.parent.tabs.currentIndex(), os.path.split(fname)[1])
                    self.parent.tabs.setTabToolTip(self.parent.tabs.currentIndex(), fname)
                    return
                else:
                    window = self.parent._add_analysis_tab(tabname=os.path.split(fname)[1])
                    window.load_analysis(fname)

        else:
            e.ignore()

    def open_filedialog(self):
        self.open_rkf_filedialog.open(self, QtCore.SLOT("get_file_from_dialog()"))

    def open_sheets_filedialog(self):
        self.new_sheets_filedialog.open(self, QtCore.SLOT("get_sheets_file_from_dialog()"))

    def open_figure_filedialog(self):
        self.new_figure_filedialog.open(self, QtCore.SLOT("get_figure_file_from_dialog()"))

    def open_errorialog(self, message):
        self.errordialog.open(self, QtCore.SLOT("get_file_from_dialog()"))

    @QtCore.Slot()
    def get_file_from_dialog(self):
        file = self.open_rkf_filedialog.selectedFiles()[0]
        self.load_analysis(file)

    @QtCore.Slot()
    def get_figure_file_from_dialog(self):
        file = self.new_figure_filedialog.selectedFiles()[0]

    def _update_plot(self):
        settings = self.parent.settings_dialog.get_flat_state()
        # set the colors
        self._set_mix_settings()
        self._draw_diagram(
            ylim=self.ylim,
            xtick_order=self._xtick_order,
            **settings
            )

    def _set_mix_settings(self):
        allowed_mos = self.allowed_mos_override if self.allowed_mos_override is not None else self._orb_selection_dialog.state.allowed_mos()
        allowed_sfos = self.allowed_sfos_override if self.allowed_sfos_override is not None else self._orb_selection_dialog.state.allowed_sfos()
        self.main_mix.set_oi_threshold(10**(self.slider_OI.value()/slider_resolution))
        self.main_mix.set_pr_threshold(self.slider_PR.value()/slider_resolution/1000)
        self.main_mix.set_enable_oi(self.cbox_OI.isChecked())
        self.main_mix.set_enable_pr(self.cbox_PR.isChecked())
        self.main_mix.set_allowed_mos(allowed_mos)
        self.main_mix.set_allowed_sfos(allowed_sfos)
        self.main_mix.set_energy_type(self._energytype_selection)
        self.main_mix.reset_mixes()

    def _draw_diagram(self, *args, ylim=None, **kwargs):
        ax = self.plot.axes
        fig = self.plot.fig

        ax.clear()
        ax.yaxis.set_major_formatter('{x: 3.0f}')
        self.main_mix.draw_diagram(ax=ax, ylim=ylim, highlighted_orbitals=self.filtered_orbitals, use_darkmode=False, **kwargs)
        if self.new_tick_labels is not None:
            self.plot.axes.set_xticklabels(self.new_tick_labels)

        self.OI_is_empty_label.setVisible(self.main_mix.main_mix.OI_is_empty)
        self.PR_is_empty_label.setVisible(self.main_mix.main_mix.PR_is_empty)

        if self.plot._add_warning:
            ax.text(0.9, 0.9, '⚠︎', transform=fig.transFigure, fontsize=30, c='r')


        props = dict(edgecolor='white', facecolor='white', alpha=1)  # bbox features
        fig.canvas.draw_idle()

    def _set_orbital_filter(self):
        selected_systems = set([self.orbs.mos if isinstance(orb, pyfmo.orbitals.objects.MO) else orb.fragment for orb in self.plot._selected_orbitals if not isinstance(orb, tuple)])

        # self.orbital_filter_button.checked() = not self.orbital_filter_button._is_checked

        if self.orbital_filter_button.isChecked():
            if self.orbs.mos in selected_systems:
                self.allowed_mos_override = [orb for orb in self.plot._selected_orbitals if isinstance(orb, pyfmo.orbitals.objects.MO)]
            else:
                self.allowed_mos_override = self.orbs.mos.orbitals

            self.allowed_sfos_override = []
            for fragment in self.orbs.fragments:
                if fragment in selected_systems:
                    self.allowed_sfos_override.extend([orb for orb in self.plot._selected_orbitals if isinstance(orb, pyfmo.orbitals.objects.SFO) and orb.fragment == fragment])
                else:
                    self.allowed_sfos_override.extend(self.orbs.sfos.filter(fragment=fragment))
            self.filtered_orbitals = self.plot._selected_orbitals
        else:
            self.allowed_mos_override = None
            self.allowed_sfos_override = None
            self.filtered_orbitals = None

        self._update_plot()

        # self.orbital_filter_button.setEnabled(len(self.plot._selected_orbitals) > 0 or self.orbital_filter_button._is_checked)

    def load_analysis(self, file):
        try:
            self.orbs = pyfmo.Orbitals(file)
        except Exception as e:
            import traceback

            reason = ""
            try:
                reader = plams.KFReader(file)
                ident = reader.read('General', 'program')
                if ident != 'ADF':
                    reason = 'Not an adf.rkf file.'
            except (IndexError, plams.core.errors.FileError):
                reason = "Not an RKF file."

            self.errordialog.showMessage(f'Could not load orbital data from file:\n{reason}\n\n{file}')
            traceback.print_exc(e)
            return

        self.main_mix = pyfmo.analysis.mixing.Mixer2(self.orbs, pr_min_thresh=0.001**2, oi_min_thresh=0.00000001)

        self._energytype_selection = 'energy'

        self.setAcceptDrops(True)
        if hasattr(self, "_new_page_frame"):
            self._new_page_frame.hide()
            self.central_layout.removeWidget(self._new_page_frame)

        self._energytype_selection_dialog = ETypeDialog(self, self.orbs.sfo_energy_types, self._energytype_selection)

        self._orb_selection = {'Complex': {}}
        for frag in self.orbs.fragments:
            self._orb_selection[frag] = {}
            for sfo in self.orbs.sfos.filter(fragment=frag):
                self._orb_selection[frag][sfo] = True

        for mo in self.orbs.mos:
            self._orb_selection['Complex'][mo] = True

        self._analysis_page_frame = QtWidgets.QFrame(self)
        self.central_layout.addWidget(self._analysis_page_frame)
        layout = QtWidgets.QGridLayout(self._analysis_page_frame)

        plot_container_layout = QtWidgets.QVBoxLayout()
        plot_container = QtWidgets.QSplitter()
        plot_container.setStyleSheet("QSplitter::handle { border: 2px black; }")
        plot_container.setHandleWidth(10)
        plot_container.setOpaqueResize(True)
        plot_container.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self.plot = MplCanvas(self)
        plot_frame = QtWidgets.QFrame()
        shadow.apply(plot_container)
        plot_frame_layout = QtWidgets.QVBoxLayout()
        plot_frame.setLayout(plot_frame_layout)
        plot_frame_layout.addWidget(self.plot)
        self.plot.setMinimumSize(300, 300)

        self.plot.setFocus()
        # if self.parent.isDarkMode:
            # plot_frame.setStyleSheet('QFrame{padding: 0px; margin: 0px; border: 1px solid darkgray; border-radius: 5px; background-color: white;} QDialog{padding: 0px; margin: 0px; border: 1px solid darkgray; border-radius: 5px; background-color: none;} QLabel{padding: 0px; margin: 0px; border: none; border-radius: 5px; background-color: none;}')
        plot_frame.setObjectName('plot_frame')
        plot_container.addWidget(plot_frame)
        layout.addWidget(plot_container, 0, 0, 1, 2)

        self.info_tabs = QtWidgets.QTabWidget()
        shadow.apply(self.info_tabs)

        # self.info_tabs.setFixedSize(300, 500)
        orbital_info_frame = QtWidgets.QFrame()

        orbital_info_layout = QtWidgets.QGridLayout()
        orbital_info_frame.setLayout(orbital_info_layout)
        orbital_info_frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self.orbital_info_box = spoilers.Spoilers(self)

        orbital_info_layout.addWidget(QtWidgets.QLabel('<i>Use</i> <b>Shift + Click</b> <i>to select multiple orbitals!</i>'), 0, 0, 1, 2)
        orbital_info_layout.addWidget(self.orbital_info_box, 1, 0, 1, 2)
        
        self.info_tabs.addTab(orbital_info_frame, 'Orbitals')

        self.system_info_box = spoilers.Spoilers(self)
        self.system_info_box.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.system_info_box.addSpoiler('General', self._get_general_system_info(), self.parent._ICONS['info'])
        self.system_info_box.addSpoiler('Complex', self._get_complex_system_info(), self.parent._ICONS['mo'])
        for frag in self.orbs.fragments:
            self.system_info_box.addSpoiler(frag, self._get_fragment_system_info(frag), self.parent._ICONS['sfo'])

        
        self.info_tabs.addTab(self.system_info_box, 'System')
        self.info_tabs.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        plot_container.addWidget(self.info_tabs)

        slider_layout = QtWidgets.QGridLayout()
        slider_box = QtWidgets.QFrame()
        shadow.apply(slider_box)
        slider_box.setObjectName('sliderbox')
        slider_box.setLayout(slider_layout)
        layout.addWidget(slider_box, 1, 0)

        has_OI = len(self.main_mix.mixes['OI']['energy']) > 0

        self.OI_is_empty_label = QtWidgets.QLabel('⚠️')
        self.OI_is_empty_label.setToolTip('Could not find any Orbital Interactions for these settings.')
        slider_layout.addWidget(self.OI_is_empty_label, 0, 0)
        sp_retain = self.OI_is_empty_label.sizePolicy()
        sp_retain.setRetainSizeWhenHidden(True)
        self.OI_is_empty_label.setSizePolicy(sp_retain)

        self.cbox_OI = QtWidgets.QCheckBox('Show OI')
        self.cbox_OI.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.cbox_OI.setChecked(True)
        slider_layout.addWidget(self.cbox_OI, 0, 1)
        self.cbox_OI.checkStateChanged.connect(self._update_plot)

        label_OI = QtWidgets.QLabel('τ<sub>oi</sub> =')
        label_OI.setToolTip('The threshold value for Orbital Interactions')
        label_OI.setStyleSheet('QLabel{ font: 12pt}')
        slider_layout.addWidget(label_OI, 0, 2)

        inc_oi_btn = QtWidgets.QPushButton('<')
        inc_oi_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        shadow.apply(inc_oi_btn)
        inc_oi_btn.setToolTip('Add next Orbital Interaction')
        inc_oi_btn.clicked.connect(self._set_next_oi_slider)
        slider_layout.addWidget(inc_oi_btn, 0, 4)

        self.slider_OI = ValueSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        self.slider_OI.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        self.slider_OI.setObjectName('slider_OI')
        if not has_OI:
            slider_OI_max = 1
        else:
            slider_OI_max = abs(min(min(v.values()) for v in self.main_mix.mixes['OI'].values()))

        self.slider_OI.setMinimum(np.log10(0.00000001) * slider_resolution)
        self.slider_OI.setMaximum(floor(np.log10(slider_OI_max) * slider_resolution))
        self.slider_OI.setSliderPosition(np.log10(slider_OI_max/1.5) * slider_resolution)
        self.main_mix.set_oi_threshold(slider_OI_max/1.5)
        slider_layout.addWidget(self.slider_OI, 0, 5)

        dec_oi_btn = QtWidgets.QPushButton('>')
        dec_oi_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        dec_oi_btn.setToolTip('Remove weakest Orbital Interaction')
        dec_oi_btn.clicked.connect(self._set_previous_oi_slider)
        slider_layout.addWidget(dec_oi_btn, 0, 6)
        shadow.apply(dec_oi_btn)

        label_value_OI = QtWidgets.QLabel(f'{10**(self.slider_OI.value()/slider_resolution):.2E}')
        label_value_OI.setToolTip('The threshold value for Orbital Interactions')
        label_value_OI.setStyleSheet('font: 10px "IBM Plex Mono"')
        slider_layout.addWidget(label_value_OI, 0, 3)
        self.slider_OI.valueChanged.connect(lambda value: (self._update_plot(), label_value_OI.setText(f'{10**(value/slider_resolution):.2E}')))

        if not has_OI:
            self.cbox_OI.setEnabled(False)
            label_OI.setEnabled(False)
            label_value_OI.setEnabled(False)
            inc_oi_btn.setEnabled(False)
            self.slider_OI.setEnabled(False)
            dec_oi_btn.setEnabled(False)

        has_PR = len(self.main_mix.mixes['PR']['energy']) > 0

        self.PR_is_empty_label = QtWidgets.QLabel('⚠️')
        if has_PR:
            self.PR_is_empty_label.setToolTip('Could not find any Pauli Repulsions for these settings.')
        else:
            self.PR_is_empty_label.setToolTip('This system has no Pauli Repulsive interactions.')

        slider_layout.addWidget(self.PR_is_empty_label, 1, 0)

        self.cbox_PR = QtWidgets.QCheckBox('Show PR')
        self.cbox_PR.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.cbox_PR.setChecked(True)
        slider_layout.addWidget(self.cbox_PR, 1, 1)
        self.cbox_PR.checkStateChanged.connect(self._update_plot)

        label_PR = QtWidgets.QLabel('τ<sub>pr</sub> =')
        label_PR.setToolTip('The threshold value for Pauli Repulsions')
        label_PR.setStyleSheet('QLabel{ font: 12pt}')
        slider_layout.addWidget(label_PR, 1, 2)

        self.slider_PR = ValueSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        self.slider_PR.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        # slider_PR_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in self.pauli_mixes.values())
        if not has_PR:
            slider_PR_max = 0.001**2
        else:
            # slider_PR_max = max(self.main_mix.mixes['PR'].values())
            slider_PR_max = abs(max(max(v.values()) for v in self.main_mix.mixes['PR'].values()))

        inc_pr_btn = QtWidgets.QPushButton('<')
        inc_pr_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        inc_pr_btn.setToolTip('Show next Pauli Repulsion')
        inc_pr_btn.clicked.connect(self._set_next_pr_slider)
        shadow.apply(inc_pr_btn)
        slider_layout.addWidget(inc_pr_btn, 1, 4)
        self.slider_PR.setMinimum(0.001**2 * 1000 * slider_resolution)
        self.slider_PR.setMaximum(slider_PR_max * 1000 * slider_resolution)
        self.slider_PR.setSliderPosition(slider_PR_max/1.5 * 1000 * slider_resolution)
        self.slider_PR.setObjectName('slider_PR')
        self.main_mix.set_pr_threshold(slider_PR_max/1.1)
        slider_layout.addWidget(self.slider_PR, 1, 5)

        dec_pr_btn = QtWidgets.QPushButton('>')
        dec_pr_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        dec_pr_btn.setToolTip('Remove weakest Pauli Repulsion')
        dec_pr_btn.clicked.connect(self._set_previous_pr_slider)
        shadow.apply(dec_pr_btn)
        slider_layout.addWidget(dec_pr_btn, 1, 6)
        label_value_PR = QtWidgets.QLabel(f'{self.slider_PR.value()/slider_resolution:.3f}')
        label_value_PR.setStyleSheet('font: 10px "IBM Plex Mono"')
        label_value_PR.setToolTip('The threshold value for Pauli Repulsions')
        slider_layout.addWidget(label_value_PR, 1, 3)
        self.slider_PR.valueChanged.connect(lambda value: (self._update_plot(), label_value_PR.setText(f'{value/slider_resolution/1000:.4f}')))
        
        if not has_PR:
            self.cbox_PR.setEnabled(False)
            label_PR.setEnabled(False)
            label_value_PR.setEnabled(False)
            inc_pr_btn.setEnabled(False)
            self.slider_PR.setEnabled(False)
            dec_pr_btn.setEnabled(False)

        slider_layout.setColumnStretch(0, 0)
        slider_layout.setColumnStretch(1, 0)
        slider_layout.setColumnStretch(2, 0)
        slider_layout.setColumnStretch(3, 0)
        slider_layout.setColumnStretch(4, 0)
        slider_layout.setColumnStretch(5, 1)
        slider_layout.setColumnStretch(6, 0)

        selector_box = QtWidgets.QFrame()
        # selector_box
        selector_layout = QtWidgets.QGridLayout()
        selector_box.setLayout(selector_layout)
        layout.addWidget(selector_box, 1, 1)

        etype_btn = QtWidgets.QPushButton('Energy Type')
        etype_btn.clicked.connect(self._energytype_selection_dialog.exec)
        selector_layout.addWidget(etype_btn, 1, 0)
        shadow.apply(etype_btn)

        # self._orb_selection_dialog = OrbitalSelectionDialog(self, self._orb_selection)
        self._orb_selection_dialog = orbital_selector.OrbitalSelectionDialog(self, self.orbs)

        orb_btn = QtWidgets.QPushButton('Orbitals')
        orb_btn.clicked.connect(self._orb_selection_dialog.exec)
        selector_layout.addWidget(orb_btn, 1, 1)
        shadow.apply(orb_btn)

        self.orbital_draw_button = QtWidgets.QPushButton()
        self.orbital_draw_button.setEnabled(False)
        # self.orbital_draw_button.setStyleSheet('QPushButton::menu-indicator { image: none; }')

        menu = QtWidgets.QMenu(self)
        self.orbital_draw_button.setMenu(menu)
        self.orbital_draw_button.setText('Draw Orbitals')
        shadow.apply(self.orbital_draw_button)
        self.orbital_filter_button = QtWidgets.QPushButton('Filter')
        self.orbital_filter_button.setCheckable(True)
        self.orbital_filter_button.clicked.connect(self._set_orbital_filter)
        self.orbital_filter_button.setEnabled(False)
        shadow.apply(self.orbital_filter_button)
        self.allowed_mos_override = None
        self.allowed_sfos_override = None
        self.filtered_orbitals = None
        selector_layout.addWidget(self.orbital_draw_button, 0, 0)
        selector_layout.addWidget(self.orbital_filter_button, 0, 1)

        orb_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        etype_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.orbital_draw_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.orbital_filter_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        # we have to set the size of the buttons manually
        max_width = max(
            self.orbital_filter_button.sizeHint().width(),
            self.orbital_draw_button.sizeHint().width(),
            etype_btn.sizeHint().width(),
            orb_btn.sizeHint().width())

        selector_layout.setColumnMinimumWidth(0, max_width)
        selector_layout.setColumnMinimumWidth(1, max_width)

        misc_box = QtWidgets.QFrame()
        misc_box_layout = QtWidgets.QGridLayout()
        misc_box.setLayout(misc_box_layout)
        layout.addWidget(misc_box, 2, 0)

        make_sheet_btn = QtWidgets.QPushButton('Generate Sheets')
        shadow.apply(make_sheet_btn)
        save_fig_btn = QtWidgets.QPushButton('Save Figure')
        shadow.apply(save_fig_btn)

        make_sheet_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        save_fig_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        make_sheet_btn.clicked.connect(self.get_sheets_save_file)
        save_fig_btn.clicked.connect(self.get_figure_save_file)
        misc_box_layout.addWidget(make_sheet_btn, 0, 0, 1, 1)
        misc_box_layout.addWidget(save_fig_btn, 0, 1, 1, 1)

        self.notice_tab = spoilers.Spoilers(self)

        notices_frame = QtWidgets.QFrame()
        notices_layout = QtWidgets.QVBoxLayout()
        notices_frame.setLayout(notices_layout)
        notices_frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        more_info_label = QtWidgets.QLabel('<a href=\"http://google.com/\">Click here for more information about notices</a>')
        more_info_label.setTextFormat(QtCore.Qt.RichText)
        more_info_label.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        more_info_label.setOpenExternalLinks(True)

        notices_layout.addWidget(more_info_label)
        notices_layout.addWidget(self.notice_tab)
        
        self.notice_tab_idx = self.info_tabs.addTab(notices_frame, 'Notices')

        for title, text in self.orbs.notices['warning']:
            self.add_warning_notice(title, text)
        for title, text in self.orbs.notices['error']:
            self.add_error_notice(title, text)
        for title, text in self.orbs.notices['info']:
            self.add_info_notice(title, text)

        self._loaded_analysis = True
        self.parent.settings_dialog.settingsChanged.connect(self._update_plot)

        msg_title = None
        msg_text = 'The provided calculation triggered '

        if len(self.orbs.notices['warning']) > 0 and len(self.orbs.notices['error']) > 0:
            msg_text = f'The provided calculation triggered {len(self.orbs.notices["warning"])} warning(s) and {len(self.orbs.notices["error"])} error(s).\n\nPlease read the Notices carefully!'
            msg_title = 'Error'
            self.info_tabs.setCurrentIndex(2)
            QtWidgets.QMessageBox.critical(self, msg_title, msg_text)
        elif len(self.orbs.notices['warning']) > 0:
            msg_text = f'The provided calculation triggered {len(self.orbs.notices["warning"])} warning(s).\n\nPlease read the Notices carefully!'
            msg_title = 'Warning'
            self.info_tabs.setCurrentIndex(2)
            QtWidgets.QMessageBox.warning(self, msg_title, msg_text)
        elif len(self.orbs.notices['error']) > 0:
            msg_text = f'The provided calculation triggered {len(self.orbs.notices["error"])} error(s).\n\nPlease read the Notices carefully!'
            msg_title = 'Error'
            self.info_tabs.setCurrentIndex(2)
            QtWidgets.QMessageBox.critical(self, msg_title, msg_text)

        self._update_plot()


    def add_info_notice(self, title, text):
        label = QtWidgets.QLabel(text)
        label.setStyleSheet('font: 10px "IBM Plex Mono";')
        self.notice_tab.addSpoiler(title, label, icon=self.parent._ICONS['info'])
        # self._reset_notice_bar_tabbutton()

    def add_warning_notice(self, title, text):
        label = QtWidgets.QLabel(text)
        label.setStyleSheet('font: 10px "IBM Plex Mono";')
        self.notice_tab.addSpoiler(title, label, icon=self.parent._ICONS['warning'])
        # self._reset_notice_bar_tabbutton()

    def add_error_notice(self, title, text):
        label = QtWidgets.QLabel(text)
        label.setStyleSheet('font: 10px "IBM Plex Mono";')
        self.notice_tab.addSpoiler(title, label, icon=self.parent._ICONS['error'])
        # self._reset_notice_bar_tabbutton()

    def _reset_notice_bar_tabbutton(self):
        n_info = len(self.orbs.notices['info'])
        n_warning = len(self.orbs.notices['warning'])
        n_error = len(self.orbs.notices['error'])
        pixmap = latex_renderer.multicolor_QPixMap((r'$\text{ⓘ}$', f'{n_info} ', r'$\text{⚠}$', f'{n_warning} ', r'$\text{×}$', f'{n_error}'), ('#2ec4b6', 'k', '#ff9f1c', 'k', '#e71d36', 'k'))
        label = QtWidgets.QLabel()
        label.setPixmap(pixmap)
        self.info_tabs.tabBar().setTabButton(self.notice_tab_idx, QtWidgets.QTabBar.LeftSide, label)

    def _get_general_system_info(self):
        frame = QtWidgets.QFrame()
        frame.setStyleSheet('QLabel{padding: 2px; font-size: 10pt} QPushButton{icon-size: 10px;}')
        layout = QtWidgets.QGridLayout()
        frame.setLayout(layout)

        row = 0
        layout.addWidget(CopyLabel('<b>adf.rkf path</b>', self.orbs.reader.path), row, 0, 1, 2)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Symmetry</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(self.orbs.data['calc_info']['symmetry'].strip())
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        # number of fragments
        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Nº Fragments</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(str(len(self.orbs.fragments)))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        # EDA terms
        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>E</i><sub>int</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'Bond Energy') * 627.503:.2f} kcal mol<sup>–1</sup>")
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>V</i><sub>elstat</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'elstat') * 627.503:.2f} kcal mol<sup>–1</sup>")
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>E</i><sub>Pauli</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'Pauli Total') * 627.503:.2f} kcal mol<sup>–1</sup>")
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>E</i><sub>oi</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503:.2f} kcal mol<sup>–1</sup>")
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        irreps = []
        for lab in self.orbs.reader.read('Symmetry', 'symlab').split():
            l = lab.split(':')[0]
            if l not in irreps:
                irreps.append(l)

        if len(irreps) > 1:
            for lab in irreps:
                row += 1
                layout.addWidget(QtWidgets.QLabel(f'<b>Δ<i>E</i><sub>oi</sub></b>({pyfmo.translate_irrep_label(lab, mode="html")})'), row, 0, 1, 1)
                label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', f'Orb.Int. {lab}') * 627.503:.2f} kcal mol<sup>–1</sup>")
                label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>E</i><sub>disp</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'Dispersion Energy') * 627.503:.2f} kcal mol<sup>–1</sup>")
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        return frame

    def _get_complex_system_info(self):
        frame = QtWidgets.QFrame()
        frame.setStyleSheet('QLabel{padding: 2px; font: 10pt} QPushButton{icon-size: 10px;}')
        layout = QtWidgets.QGridLayout()
        frame.setLayout(layout)

        charges = _determine_formal_charges(self.orbs)
        unrestricted_mos = self.orbs.data['calc_info']['unrestricted_mos']

        all_spin_pols = self.orbs.data['calc_info']['sfo_spinpolarizations']
        total_spin_pols = 0
        for frag, frag_spin_pols in all_spin_pols.items():
            for irrep, spin_pols in frag_spin_pols.items():
                total_spin_pols += spin_pols[0] - spin_pols[1]

        row = 0
        layout.addWidget(MoleculeLabel(self, '<b>Geometry (xyz)</b>', self.orbs.molecule, mol2xyz(self.orbs.molecule)), row, 0, 1, 2)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Charge</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{round(charges["Complex"]):+d}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Restricted</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(str(not unrestricted_mos))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Spin-Polarization</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{round(total_spin_pols):+d}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Nº Electrons</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(str(sum(mo.occupation for mo in self.orbs.mos)))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        return frame

    def _get_fragment_system_info(self, frag):
        frame = QtWidgets.QFrame()
        frame.setStyleSheet('QLabel{padding: 2px; font: 10pt;} QPushButton{icon-size: 10px;}')
        layout = QtWidgets.QGridLayout()
        frame.setLayout(layout)

        sfos = self.orbs.sfos.filter(fragment=frag)
        charges = _determine_formal_charges(self.orbs)
        vdd_charges = _determine_vdd_charges(self.orbs)
        unrestricted_sfos = self.orbs.data['calc_info']['unrestricted_sfos']

        row = 0
        layout.addWidget(MoleculeLabel(self, '<b>Geometry (xyz)</b>', sfos[0].molecule, mol2xyz(sfos[0].molecule)), row, 0, 1, 2)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Restricted</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(str(not unrestricted_sfos))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        # print(self.orbs.data['calc_info']['sfo_spinpolarizations'])
        frag_spin_pols = self.orbs.data['calc_info']['sfo_spinpolarizations'][frag]
        # print(frag_spin_pols)
        total_spin_pols = 0
        for irrep, spin_pols in frag_spin_pols.items():
            total_spin_pols += spin_pols[0] - spin_pols[1]
            # print(irrep, spin_pols)

            row += 1
            layout.addWidget(QtWidgets.QLabel(f'<b>Spin-Polarization ({pyfmo.translate_irrep_label(irrep, mode="html")})</b>'), row, 0, 1, 1)
            label = QtWidgets.QLabel(f'{round(spin_pols[0] - spin_pols[1]):+d}')
            label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Spin-Polarization (Total)</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{round(total_spin_pols):+d}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Nº Electrons</b>'), row, 0, 1, 1)
        occ = sum(sfo.occupation for sfo in sfos)
        label = QtWidgets.QLabel(str(round(occ)))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Formal Charge</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{round(charges[frag]):+d}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Mulliken Charge</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{occ - sum(sfo.gross_population for sfo in sfos):+.3f}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>VDD Charge</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{vdd_charges[frag]:+.3f}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        return frame

    def _set_next_pr_slider(self):
        next_val = self.main_mix.get_next_pr_threshold()
        self.slider_PR.setSliderPosition(floor(next_val * 1000 * slider_resolution))

    def _set_previous_pr_slider(self):
        next_val = self.main_mix.get_previous_pr_threshold()
        self.slider_PR.setSliderPosition(ceil(next_val * 1000 * slider_resolution))

    def _set_next_oi_slider(self):
        next_val = self.main_mix.get_next_oi_threshold()
        self.slider_OI.setSliderPosition(floor(np.log10(next_val) * slider_resolution))

    def _set_previous_oi_slider(self):
        next_val = self.main_mix.get_previous_oi_threshold()
        self.slider_OI.setSliderPosition(ceil(np.log10(next_val) * slider_resolution))

    def get_sheets_save_file(self):
        print('hello')
        d = os.path.join(os.path.split(self.orbs.kfpath)[0], 'pyorbb.xlsx')
        filename, v = QtWidgets.QFileDialog.getSaveFileName(None, "Save File")
        
        print('after', repr(filename), repr(v))
        if not filename.strip():
            return

        self.orbs.write_excel(filename)
        os.system(f'open {filename}')

    def get_figure_save_file(self):
        d = os.path.join(os.path.split(self.orbs.kfpath)[0], 'pyorbb.png')
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", dir=d, filter="PNG file (*.png);;Any file (*)")
        
        if not filename.strip():
            return

        self.plot.fig.savefig(filename, dpi=600)

    def update_icons(self):
        app = QtWidgets.QApplication.instance()
        darkmode = app.isDarkMode
        if darkmode:
            self.dropfile_label.setPixmap(app._ICONS['dropfile_dark'].pixmap(50, 50))
            self.open_file_button.setIcon(app._ICONS['folder_dark'])
        else:
            self.dropfile_label.setPixmap(app._ICONS['dropfile'].pixmap(50, 50))
            self.open_file_button.setIcon(app._ICONS['folder'])

    def setup_new(self):
        self._new_page_frame = QtWidgets.QFrame(self)
        new_page_layout = QtWidgets.QVBoxLayout(self._new_page_frame)
        self.central_layout.addWidget(self._new_page_frame, QtCore.Qt.AlignCenter)

        outlined_area = QtWidgets.QFrame()
        layout = QtWidgets.QVBoxLayout()
        drop_area = QtWidgets.QFrame()
        drop_area.setObjectName('drop_area')
        shadow.apply(drop_area)
        layout.addWidget(drop_area)

        new_page_layout.addWidget(drop_area)
        # outlined_area
        outlined_area_layout = QtWidgets.QVBoxLayout(outlined_area)
        drop_area_layout = QtWidgets.QGridLayout(drop_area)

        # image
        pixmap = QtWidgets.QApplication.instance()._ICONS['dropfile'].pixmap(50, 50)
        self.dropfile_label = QtWidgets.QLabel()
        self.dropfile_label.setPixmap(pixmap)
        outlined_area_layout.addWidget(self.dropfile_label, alignment=QtCore.Qt.AlignCenter)

        # Label
        label = QtWidgets.QLabel('<font size="8">Drop a File here</font><br><font size="6"><i>or</i></font>')
        label.setAlignment(QtCore.Qt.AlignCenter)
        label.setTextFormat(QtCore.Qt.RichText)
        outlined_area_layout.addWidget(label)

        # Button
        self.open_file_button = QtWidgets.QPushButton(QtWidgets.QApplication.instance()._ICONS['folder'], ' Select a File')
        shadow.apply(self.open_file_button, radius=30)
        self.open_file_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.open_file_button.setFlat(True)
        self.open_file_button.clicked.connect(self.open_filedialog)

        outlined_area.setObjectName('Outline')        
        outlined_area_layout.addWidget(self.open_file_button, alignment=QtCore.Qt.AlignCenter)
        drop_area_layout.addWidget(outlined_area, 0, 0, 1, 2)

        open_article_btn = QtWidgets.QPushButton('Open the PyOrbb Article')
        open_article_btn.setStyleSheet("margin-left: 20px; margin-right: 20px")
        shadow.apply(open_article_btn, radius=30)
        open_article_btn.clicked.connect(lambda: webbrowser.open('https://github.com/TheoChem-VU/PyFMO'))
        open_article_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        drop_area_layout.addWidget(open_article_btn, 1, 0)

        open_docs_btn = QtWidgets.QPushButton('See the PyOrbb Documentation')
        open_docs_btn.setStyleSheet("margin-left: 20px; margin-right: 20px")
        shadow.apply(open_docs_btn, radius=30)
        open_docs_btn.clicked.connect(lambda: webbrowser.open('https://theochem-vu.github.io/PyFMO/'))
        open_docs_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        drop_area_layout.addWidget(open_docs_btn, 1, 1)

        self.recent_publish_carousel = carousel.PublicationCarousel(self)
        self.recent_publish_carousel.setObjectName('carousel')
        new_page_layout.addWidget(self.recent_publish_carousel)

        # add a button for making an issue on github
        issue_btn = QtWidgets.QPushButton('?')
        issue_btn.clicked.connect(lambda: webbrowser.open('https://github.com/TheoChem-VU/PyFMO/issues/new'))
        issue_btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        issue_btn.setToolTip('Open an issue on GitHub')
        issue_btn.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        issue_btn.setFixedSize(30, 30)
        issue_btn.setStyleSheet('QPushButton{border-radius: 15px; font-weight: bold}')
        shadow.apply(issue_btn)
        new_page_layout.addWidget(issue_btn)


        self.update_icons()



class PyOrbbWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.theme_switcher = theme_switcher.ThemeSwitcher(self)

        ICON_FOLDER = os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons')
        self._ICONS = {file.removesuffix('.png'): QtGui.QIcon(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}
        self._PIXMAPS = {file.removesuffix('.png'): QtGui.QPixmap(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}
        
        self.resize(1030 + 22 + 12, 698 + 52)
        self.layout = QtWidgets.QGridLayout()
        grid_widget = QtWidgets.QWidget()

        self.windows = []

        grid_widget.setLayout(self.layout)
        self.setCentralWidget(grid_widget)

        self.setWindowTitle("PyOrbb Analysis Tool")

        self.tabs = editable_tabs.WindowTabs(self)

        add_tab_button = QtWidgets.QPushButton('+')
        add_tab_button.resize(50, 50)
        add_tab_button.clicked.connect(self._add_analysis_tab)
        add_tab_button.setStyleSheet("""
            QPushButton {
                font-size: 20px;
                border-radius: 16px;
                padding: 6px;
                background-color: transparent;
                border: none;
            }
            """)
        add_tab_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        tab_shape = QtWidgets.QTabBar.Shape.RoundedNorth
        self.tabs.tabBar().setShape(tab_shape)
        self.tabs.tabCloseRequested.connect(self.close_window_on_last_tab)
        self.tabs.setCornerWidget(add_tab_button, QtCore.Qt.TopLeftCorner)

        self.layout.addWidget(self.tabs, 0, 0, 1, 1)

        # File menu
        menuBar = self.menuBar()
        fileMenu = menuBar.addMenu("File")
        action = fileMenu.addAction("New PyOrbb window")
        action.triggered.connect(QtWidgets.QApplication.instance().add_window)

        action = fileMenu.addAction("New PyOrbb Viewer window")
        action.triggered.connect(QtWidgets.QApplication.instance().open_empty_viewer)

        plot_menu = menuBar.addMenu("Plot")

        quit = QtGui.QAction("&Quit", self)
        quit.setShortcut("Ctrl+Q")
        fileMenu.addAction(quit)

        # Edit menu
        editMenu = menuBar.addMenu("Edit")
        editMenu.addAction("Copy")
        editMenu.addAction("Paste")

        # Help menu
        preferenceMenu = menuBar.addMenu('Preferences')

        settings_action = preferenceMenu.addAction("Open Settings")
        settings_action.triggered.connect(self._open_settings)
        self.setWindowIcon(self._ICONS["pyorbb"])

        self.settings_dialog = settings.SettingsDialog(self)

        self.settings_dialog.settingsChanged.connect(self.set_theme)

        # self.tabs.palette().setColor(QtGui.QPalette.ColorRole.Window, 'red')
        # self.tabs.palette().setColor(QtGui.QPalette.ColorRole.Button, 'red')
        # self.tabs.palette().setColor(QtGui.QPalette.ColorRole.Base, 'red')
        self._add_analysis_tab()
        self.tabs.setCurrentIndex(0)

        self.set_theme()

    def set_theme(self):
        settings_theme = self.settings_dialog.get("PyOrbb", "Color Scheme", "theme_mode")
        if settings_theme == 0:
            self.theme_switcher.set_light_theme()
        elif settings_theme == 1:
            self.theme_switcher.set_dark_theme()
        else:
            self.theme_switcher.set_auto_theme()
        self.changeEvent(QtCore.QEvent(QtCore.QEvent.Type.ThemeChange))
        self.set_style()

    def set_style(self):
        shadow.update_style()
        for window in self.windows:
            window.recent_publish_carousel.update_background_color()
            if window._loaded_analysis:
                window.orbital_info_box.themechange()
                window.system_info_box.themechange()
                window.notice_tab.themechange()
                window._orb_selection_dialog.themechange()
                # window._set_orbital_filter_button_icon()
            else:
                window.update_icons()

        if self.isDarkMode:
            with open(os.path.split(__file__)[0] + '/style_dark.qss') as style:
                self.setStyleSheet(style.read())
        else:
            with open(os.path.split(__file__)[0] + '/style.qss') as style:
                self.setStyleSheet(style.read())

    def changeEvent(self, event):
        if event.type() == QtCore.QEvent.Type.ThemeChange:
            self.set_style()

        super().changeEvent(event)

    def _open_settings(self):
        self.settings_dialog.exec()

    def close_window_on_last_tab(self, index):
        if self.tabs.count() == 0:
            QtWidgets.QApplication.instance().remove_window(self)

    @property
    def isDarkMode(self):
        settings_theme = self.settings_dialog.get("PyOrbb", "Color Scheme", "theme_mode")
        if settings_theme == 0:
            return False
        elif settings_theme == 1:
            return True
        else:
            return QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark

    def _add_analysis_tab(self, object=None, tabname='new', tabtooltip=None):
        window = AnalysisWindow(self)
        window.setup_new()
        idx = self.tabs.addTab(window, tabname)
        self.tabs.setTabToolTip(idx, tabtooltip)
        self.tabs.setCurrentIndex(idx)
        self.windows.append(window)
        return window


class PyOrbbApp(QtWidgets.QApplication):
    def __post_init__(self):
        fontpath = os.path.split(__file__)[0] + '/../fonts/Inter/Inter-VariableFont_opsz,wght.ttf'
        QtGui.QFontDatabase.addApplicationFont(fontpath)
        fontpath = os.path.split(__file__)[0] + '/../fonts/ibm_plex_mono/IBMPlexMono-Regular.ttf'
        QtGui.QFontDatabase.addApplicationFont(fontpath)
        self.setStyle('Fusion')
        self.windows = []

        ICON_FOLDER = os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons')
        self._ICONS = {file.removesuffix('.png'): QtGui.QIcon(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}
        self._PIXMAPS = {file.removesuffix('.png'): QtGui.QPixmap(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}

    def __enter__(self):
        self.__post_init__()
        return self

    def __exit__(self, *args):
        self.add_window()
        self.exec()
        self.shutdown()

    def add_window(self):
        win = PyOrbbWindow()
        self.windows.append(win)
        win.show()
        return win

    def remove_window(self, window):
        self.windows.remove(window)
        window.close()

    @property
    def isDarkMode(self):
        if len(self.windows) == 0:
            return QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark
        else:
            return self.windows[0].isDarkMode

    def open_empty_viewer(self):
        import tcviewer
        tcviewer_screen = tcviewer.screen._ScreenWindow()
        tcviewer_screen.setWindowIcon(self._ICONS['pyorbb'])
        tcviewer_screen.setWindowTitle('PyOrbb Viewer')
        tcviewer_screen.show()
        self.windows.append(tcviewer_screen)
        return tcviewer_screen

