from PySide6 import QtWidgets, QtCore, QtGui
import pyfmo
from .components import orbital_selector, rich_widgets
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backend_tools import Cursors
from matplotlib.figure import Figure
import numpy as np
import os
import shutil
import platformdirs
import json
from math import floor, ceil
import tcutility
from functools import partial
import pyperclip
import platform
import traceback

slider_resolution = 500


def mol2xyz(mol):
    s = ''
    s += str(len(mol.atoms)) + '\n\n'
    for atom in mol:
        s += f'{atom.symbol:2} {atom.coords[0]} {atom.coords[1]} {atom.coords[2]}\n'
    return s

def load_setting(key, default=None):
    d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
    if not os.path.exists(os.path.join(d, 'settings.json')):
        return

    with open(os.path.join(d, 'settings.json')) as jf:
        data = json.loads(jf.read())

    return data.get(key, None)

def save_setting(key, value):
    d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
    if not os.path.exists(os.path.join(d, 'settings.json')):
        data = {}
    else:
        with open(os.path.join(d, 'settings.json')) as jf:
            data = json.loads(jf.read())

    data[key] = value
    with open(os.path.join(d, 'settings.json'), 'w+') as jf:
        jf.write(json.dumps(data))

def default_setting(key, value):
    if load_setting(key) is not None:
        return
    save_setting(key, value)

default_setting('amsbin', '$AMSBIN')


def _determine_charges(orbs):
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
        actual_Nelectrons = int(sum(sfo.occupation for sfo in sfos))
        charges[frag] = expected_Nelectrons - actual_Nelectrons
        
    charges['Complex'] = sum(charges.values())
    return charges


class Spoilers(QtWidgets.QScrollArea):
    def __init__(self, parent=None):
        # making widget resizable
        super().__init__(parent=parent)
        self.parent = parent
        self.setWidgetResizable(True)
        self.horizontalScrollBar().setEnabled(False)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        # making qwidget object
        content = QtWidgets.QWidget(self)
        self.setWidget(content)

        # vertical box layout
        self.layout = QtWidgets.QVBoxLayout(content)
        self.layout.addStretch(1)
        if self.parent.parent.isDarkMode:
            self.setStyleSheet('background-color: rgb(50, 50, 50);')
        else:
            self.setStyleSheet('background-color: white;')

    def addSpoiler(self, title, widget, icon=None):
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(widget)
        spoiler = Spoiler(self, title, icon=icon)
        spoiler.setContentLayout(layout)
        self.layout.insertWidget(self.layout.count() - 1, spoiler)

    def empty(self):
        for i in range(self.layout.count()):
            widget = self.layout.takeAt(0)

            if widget is None:
                break

            widget = widget.widget()
            if widget is None:
                break

            widget.setParent(None)
            widget.destroy()

        self.layout.addStretch(1)


class Spoiler(QtWidgets.QWidget):
    def __init__(self, parent=None, title='', animationDuration=100, icon=None):
        """
        References:
            # Adapted from c++ version
            http://stackoverflow.com/questions/32476006/how-to-make-an-expandable-collapsable-section-widget-in-qt
        """
        super().__init__(parent=parent)

        self.animationDuration = animationDuration
        self.toggleAnimation = QtCore.QParallelAnimationGroup()
        self.contentArea = QtWidgets.QScrollArea(self)
        # self.headerLine = QtWidgets.QFrame(self)
        # toggleLayout = QtWidgets.QHBoxLayout()
        # toggleFrame = 
        self.toggleButton = rich_widgets.HTMLToolButton(self)
        # self.toggleButton = QtWidgets.QToolButton(self)
        # self.toggleButton.setStyle(rich_widgets.HTMLStyle())
        self.mainLayout = QtWidgets.QVBoxLayout()
        titleLayout = QtWidgets.QHBoxLayout()
        titleFrame = QtWidgets.QFrame()
        titleFrame.setLayout(titleLayout)

        toggleButton = self.toggleButton
        toggleButton.setStyleSheet("QToolButton { border: none; font-weight: bold; font-size: 20px; text-align: left top}")

        toggleButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        toggleButton.setArrowType(QtCore.Qt.RightArrow)
        toggleButton.setText(f'{title}')
        toggleButton.setCheckable(True)
        toggleButton.setChecked(False)

        self.contentArea.setStyleSheet("QScrollArea { background-color: white; border: none; }")
        self.contentArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        # start out collapsed
        self.contentArea.setMaximumHeight(0)
        self.contentArea.setMinimumHeight(0)
        # let the entire widget grow and shrink with its content
        toggleAnimation = self.toggleAnimation
        toggleAnimation.addAnimation(QtCore.QPropertyAnimation(self, b"minimumHeight"))
        toggleAnimation.addAnimation(QtCore.QPropertyAnimation(self, b"maximumHeight"))
        toggleAnimation.addAnimation(QtCore.QPropertyAnimation(self.contentArea, b"maximumHeight"))
        # don't waste space
        mainLayout = self.mainLayout
        mainLayout.setSpacing(0)
        mainLayout.setContentsMargins(0, 0, 0, 0)
        titleLayout.setSpacing(0)
        titleLayout.setContentsMargins(0, 0, 0, 0)
        if icon is not None:
            icon_lab = QtWidgets.QLabel()
            icon_lab.setPixmap(icon.pixmap(20, 20))
            titleLayout.addWidget(icon_lab)
        titleLayout.addWidget(self.toggleButton)
        titleLayout.addStretch()

        mainLayout.addWidget(titleFrame)
        mainLayout.addWidget(self.contentArea)
        # mainLayout.addStretch()
        self.setLayout(self.mainLayout)

        def start_animation(checked):
            arrow_type = QtCore.Qt.DownArrow if checked else QtCore.Qt.RightArrow
            direction = QtCore.QAbstractAnimation.Forward if checked else QtCore.QAbstractAnimation.Backward
            toggleButton.setArrowType(arrow_type)
            self.toggleAnimation.setDirection(direction)
            self.toggleAnimation.start()

        self.toggleButton.clicked.connect(start_animation)

    def setContentLayout(self, contentLayout):
        # Not sure if this is equivalent to self.contentArea.destroy()
        self.contentArea.destroy()
        self.contentArea.setLayout(contentLayout)
        collapsedHeight = self.sizeHint().height() - self.contentArea.maximumHeight()
        contentHeight = contentLayout.sizeHint().height()
        for i in range(self.toggleAnimation.animationCount()-1):
            spoilerAnimation = self.toggleAnimation.animationAt(i)
            spoilerAnimation.setDuration(self.animationDuration)
            spoilerAnimation.setStartValue(collapsedHeight)
            spoilerAnimation.setEndValue(collapsedHeight + contentHeight)
        contentAnimation = self.toggleAnimation.animationAt(self.toggleAnimation.animationCount() - 1)
        contentAnimation.setDuration(self.animationDuration)
        contentAnimation.setStartValue(0)
        contentAnimation.setEndValue(contentHeight)


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
        self.label.setStyleSheet('padding: 10px; font: 10px "IBM Plex Mono"')
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
            'site_energy': 'Effective',
            'site_energy_SCF0': 'Effective (initial density)',
        }

        self.rbuttons = {}
        for i, pos in enumerate(self.possibilities):
            self.rbuttons[pos] = QtWidgets.QRadioButton(display_names[pos])
            if pos == self.selected:
                self.rbuttons[pos].setChecked(True)
            rbtn_layout.addWidget(self.rbuttons[pos], i, 0, 1, 1)

        if any(charge != 0 for charge in _determine_charges(self.parent.orbs).values()):
            rbtn_layout.addWidget(QtWidgets.QLabel(f'\n<i><b>Note:</b>\nEffective energies are recommended for charged fragments!</i>'), i+1, 0, 1, 0)

        # layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
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
        if self.parent.parent.parent.isDarkMode:
            self.setStyleSheet("""
                QDialog {
                    background-color: none;
                    border: 1px solid darkgray;
                }
                QPushButton {
                    font-size: 12px;
                    border: 1px solid darkgray;
                    border-radius: 5px;
                    padding: 2px;
                    margin: 0px; 
                    background-color: #3A3A3A;
                }
                QPushButton:hover {
                    background-color: ##656565;
                    }
                QLabel: {
                    border: 0px;
                }
                """)
        else:
            self.setStyleSheet("""
                QDialog {
                    background-color: none;
                    border: 1px solid lightgray;
                }
                QPushButton {
                    font-size: 12px;
                    border: 1px solid lightgray;
                    border-radius: 5px;
                    padding: 2px;
                    margin: 0px; 
                    background-color: white;
                }
                QPushButton:hover {
                    background-color: #f0f0f0;
                    }
                QLabel: {
                    border: none;
                }
                """)
        self._frag_rename_textedit = QtWidgets.QLineEdit(self)
        layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('Save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
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
        layout.addWidget(QtWidgets.QLabel('Choose Y-axis limits'), 0, 0, 1, 2)

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
        cancel_btn = QtWidgets.QPushButton('Cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(save_btn, 4, 0, 1, 1)
        layout.addWidget(cancel_btn, 4, 1, 1, 1)

        _reset_btn = QtWidgets.QPushButton("Reset Y-axis")
        _reset_btn.clicked.connect(self.reset)
        _reset_btn.clicked.connect(self.reject)
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
        self.axes = self.fig.add_subplot(111)
        self.parent = parent
        super().__init__(self.fig)

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click)
        self._selected_orbitals = []
        self._frag_rename_dialog = FragRenameDialog(self)
        self._yaxis_dialog = YAxisDialog(self)
        self._already_unfaded = True


    def draw_orbital(self, orb=None, draw_type='single'):
        import tcviewer

        if platform.system() == "Windows":
            preambles = [f'set AMSHOME={os.path.split(self.parent.parent._amsbin_loc)[0]}', f'set AMSBIN={os.path.split(self.parent.parent._amsbin_loc)[0]}/bin']
        else:
            preambles = [f'source {os.path.join(os.path.split(self.parent.parent._amsbin_loc)[0], "amsbashrc.sh")}']

        if self.parent.tcviewer_screen is None or self.parent.tcviewer_screen.isclosed:
            self.parent.tcviewer_screen = tcviewer.screen._ScreenWindow()
            self.parent.tcviewer_screen.setWindowIcon(self.parent.parent._ICONS['pyorbb'])
            self.parent.tcviewer_screen.__enter__()
            self.parent.tcviewer_screen.setWindowTitle('PyOrbb Viewer')
            self.parent.tcviewer_screen.show()

        with self.parent.tcviewer_screen.add_molscene() as scene:
            if draw_type == 'single':
                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb.occupied else ([1, .5, 0], [0, 1, 1])
                scene.draw_molecule(orb.molecule)
                try:
                    data = orb.vtk_file(gridsize='medium', preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                scene.draw_dual_isosurface(data, colorm=c1, colorp=c2)
                scene.draw_text(str(orb))

            if draw_type == 'sum':
                mol = orb[0].molecule + orb[1].molecule
                scene.draw_molecule(mol)

                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb[0].occupied else ([1, .5, 0], [0, 1, 1])
                try:
                    cub = orb[0].vtk_file(preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                scene.draw_dual_isosurface(cub, colorm=c1, colorp=c2)

                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb[1].occupied else ([1, .5, 0], [0, 1, 1])
                try:
                    cub = orb[1].vtk_file(preambles=preambles)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                scene.draw_dual_isosurface(cub, colorm=c1, colorp=c2)
                
                scene.draw_text(str(orb[0]) + ' & ' + str(orb[1]))

            if draw_type == 'overlap':
                c1, c2 = [0, 1, 0], [1, 0, 1]
                mol = orb[0].molecule + orb[1].molecule
                scene.draw_molecule(mol)
                try:
                    cub1 = orb[0].vtk_file(preambles=preambles, grid_around_mol=mol)
                    cub2 = orb[1].vtk_file(preambles=preambles, grid_around_mol=mol)
                except Exception as e:
                    print("".join(traceback.format_exception(type(e), e, e.__traceback__)))
                    QtWidgets.QMessageBox.critical(self, 'Error', 'There was an issue with running densf.\nUse preferences > Set AMS Path to set the AMS installation path.')
                cub1.values *= cub2.values
                scene.draw_dual_isosurface(cub1, colorm=c1, colorp=c2, isovalue=0.03**2)
                
                scene.draw_text(str(orb[0]) + ' * ' + str(orb[1]))



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
                s += f'\nFragment     {orb.fragment_unique}'
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
                    if sfo2.fragment_unique == orb.fragment_unique:
                        continue
                    s += f'\n{str(sfo2):19.19} {orb @ sfo2: 5.3f} {abs(getattr(orb, self.parent._energytype_selection) - getattr(sfo2, self.parent._energytype_selection)): 8.2f}'
                
                s += '\n\nMO                     Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for mo in sorted(submix.mos, key=lambda mo: -abs(orb.mulliken_contribution(mo))):
                    s += f'\n{str(mo):19.19} {orb.mulliken_contribution(mo): 8.2%} {orb.coefficient(mo): 7.4f}'
                title = f"{orb.fragment_unique}({pyfmo.generate_label(orb, mode='html')})"

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
                title = pyfmo.generate_label(orb, mode='html')

            if isinstance(orb, tuple):
                sfo, mo = orb
                icon = self.parent.parent._ICONS['mix']
                connected_sfos = [conn[0] for conn in self.parent.main_mix.connections if conn[1] == mo and conn[0].fragment_unique != sfo.fragment_unique]
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
                title = f"{sfo.fragment_unique}({pyfmo.generate_label(sfo, mode='html')}) ⇒ {pyfmo.generate_label(mo, mode='html')}"

            label = QtWidgets.QLabel(s)
            label.setStyleSheet('padding: 3px; font: 10px "IBM Plex Mono"')
            self.parent.orbital_info_box.addSpoiler(title, label, icon)
        # self.parent.orbital_info_box.layout.addStretch(1)


    def on_plot_click(self, event):
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)

        if event.dblclick:
            if any(artist.contains(event)[0] for artist in self.axes.get_yticklabels()):
                self.parent.ylim = self._yaxis_dialog.open(self.axes.get_ylim())
                self.parent._update_plot()

            self.parent.new_tick_labels = []
            for artist in self.axes.get_xticklabels():
                if not artist.contains(event)[0]:
                    self.parent.new_tick_labels.append(artist.get_text())
                    continue
                new_txt = self._frag_rename_dialog.open(artist.get_text())
                self.parent.new_tick_labels.append(new_txt)

            self.axes.set_xticklabels(self.parent.new_tick_labels)
            self.fig.canvas.draw_idle()

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
                mo = self.parent.orbs.mos[gid[3:]]
                if mo not in self._selected_orbitals:
                    self._selected_orbitals.append(mo)
                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False

                self.selected_orbital = mo
                self.fig.canvas.draw_idle()

                break

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos[gid[4:]]
                if sfo not in self._selected_orbitals:
                    self._selected_orbitals.append(sfo)
                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False
                self.selected_orbital = sfo
                self.fig.canvas.draw_idle()

                break

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = self.parent.orbs.mos[gid[4:].split('->')[1].strip()]
                
                if (sfo, mo) not in self._selected_orbitals:
                    self._selected_orbitals.append((sfo, mo))

                if sfo not in self._selected_orbitals:
                    self._selected_orbitals.append(sfo)

                if mo not in self._selected_orbitals:
                    self._selected_orbitals.append(mo)

                self._fade_unrelated_ints(self._selected_orbitals)
                self._already_unfaded = False
                self.selected_orbital = None
                self.fig.canvas.draw_idle()

                break

        else:
            if not self._already_unfaded:
                self._already_unfaded = True
                self._unfade()
                self.selected_orbital = None
                self.fig.canvas.set_cursor(Cursors.POINTER)
                self.fig.canvas.draw_idle()

        self._set_orbital_info_box()

        # set up the menu for the pushbutton
        menu = QtWidgets.QMenu(self)
        for orb in self._selected_orbitals:
            if isinstance(orb, pyfmo.orbitals.objects.SFO):
                icon = self.parent.parent._ICONS['sfo']
            else:
                icon = self.parent.parent._ICONS['mo']

            action = QtGui.QAction(str(orb), self)
            action.setIconVisibleInMenu(True)
            action.setIcon(icon)
            action.triggered.connect(partial(self.draw_orbital, orb=orb, draw_type='single'))
            menu.addAction(action)

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

                if orb.fragment_unique == orb2.fragment_unique:
                    continue

                icon = self.parent.parent._ICONS['overlap']
                action = QtGui.QAction(f'{orb} * {orb2}', self)
                action.setIconVisibleInMenu(True)
                action.setIcon(icon)
                action.triggered.connect(partial(self.draw_orbital, orb=[orb, orb2], draw_type='overlap'))
                menu.addAction(action)

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

                if orb.fragment_unique == orb2.fragment_unique:
                    continue

                icon = self.parent.parent._ICONS['sum']
                action = QtGui.QAction(f'{orb}, {orb2}', self)
                action.setIconVisibleInMenu(True)
                action.setIcon(icon)
                action.triggered.connect(partial(self.draw_orbital, orb=[orb, orb2], draw_type='sum'))
                menu.addAction(action)


        self.parent.orbital_draw_button.setMenu(menu)
        self.parent.orbital_draw_button.setEnabled(len(menu.actions()) > 0)
        self.parent.orbital_filter_button.setEnabled(len(self._selected_orbitals) > 0 or self.parent.orbital_filter_button._is_checked)


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
                mo = self.parent.orbs.mos[gid[3:]]
                if mo in orbs:
                    continue

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos[gid[4:]]
                if sfo in orbs:
                    continue

            if gid.startswith('ARROWMO_'):
                mo = self.parent.orbs.mos[gid[8:]]
                if mo in orbs:
                    continue

            if gid.startswith('ARROWSFO_'):
                sfo = self.parent.orbs.sfos[gid[9:]]
                if sfo in orbs:
                    continue

            if gid.startswith('TEXTMO_'):
                mo = self.parent.orbs.mos[gid[7:]]
                if mo in orbs:
                    continue

            if gid.startswith('TEXTSFO_'):
                sfo = self.parent.orbs.sfos[gid[8:]]
                if sfo in orbs:
                    continue

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = self.parent.orbs.mos[gid[4:].split('->')[1].strip()]
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


class SaveFileDialog(QtWidgets.QFileDialog):
    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.setFileMode(QtWidgets.QFileDialog.AnyFile)
        self.setFilter(QtCore.QDir.Filter.Files)

    def open(self, title, filters, slot=QtCore.SLOT("get_file_from_dialog()")):
        self.setWindowTitle(title)
        self.setNameFilters(filters)
        super().open(self.parent, slot=QtCore.SLOT("get_file_from_dialog()"))
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
        copy_button.setStyleSheet("""
            QPushButton {
                font-size: 12px;
                border: 1px solid lightgray;
                border-radius: 5px;
                padding: 3px;
                margin: 0px; 
                background-color: white;
            }
            QPushButton:hover {
                background-color: #f0f0f0;
                }
            """)

        layout.addWidget(copy_button)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

    def copy(self):
        pyperclip.copy(self.copy_text)


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

        self.errordialog = QtWidgets.QErrorMessage(self)
        self.central_layout = QtWidgets.QVBoxLayout(self)

        self.new_tick_labels = None
        self.ylim = None

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
                self.load_analysis(fname)
                return
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
        self._draw_diagram(
            oi_thresh=10**(self.slider_OI.value()/slider_resolution),
            pauli_thresh=self.slider_PR.value()/slider_resolution/1000,
            energy_type=self._energytype_selection,
            ylim=self.ylim,
            )

    def _draw_diagram(self, *args, oi_thresh=None, pauli_thresh=None, ylim=None, energy_type=None):
        ax = self.plot.axes
        fig = self.plot.fig
        ax.clear()
        ax.yaxis.set_major_formatter('{x: 3.0f}')
        allowed_mos = self.allowed_mos_override if self.allowed_mos_override is not None else self._orb_selection_dialog.state.allowed_mos()
        allowed_sfos = self.allowed_sfos_override if self.allowed_sfos_override is not None else self._orb_selection_dialog.state.allowed_sfos()
        self.main_mix.set_oi_threshold(oi_thresh)
        self.main_mix.set_pr_threshold(pauli_thresh)
        self.main_mix.set_enable_oi(self.cbox_OI.isChecked())
        self.main_mix.set_enable_pr(self.cbox_PR.isChecked())
        self.main_mix.set_allowed_mos(allowed_mos)
        self.main_mix.set_allowed_sfos(allowed_sfos)
        self.main_mix.set_energy_type(energy_type)
        self.main_mix.reset_mixes()

        self.main_mix.draw_diagram(ax=ax, ylim=ylim, highlighted_orbitals=self.filtered_orbitals, use_darkmode=self.parent.isDarkMode)
        if self.new_tick_labels is not None:
            self.plot.axes.set_xticklabels(self.new_tick_labels)

        self.OI_is_empty_label.setVisible(self.main_mix.main_mix.OI_is_empty)
        self.PR_is_empty_label.setVisible(self.main_mix.main_mix.PR_is_empty)

        props = dict(edgecolor='white', facecolor='white', alpha=1)  # bbox features
        fig.canvas.draw_idle()

    def _set_orbital_filter(self):
        selected_systems = set([self.orbs.mos if isinstance(orb, pyfmo.orbitals.objects.MO) else orb.fragment_unique for orb in self.plot._selected_orbitals if not isinstance(orb, tuple)])

        self.orbital_filter_button._is_checked = not self.orbital_filter_button._is_checked

        if self.orbital_filter_button._is_checked:
            if self.orbs.mos in selected_systems:
                self.allowed_mos_override = [orb for orb in self.plot._selected_orbitals if isinstance(orb, pyfmo.orbitals.objects.MO)]
            else:
                self.allowed_mos_override = self.orbs.mos.orbitals

            self.allowed_sfos_override = []
            for fragment in self.orbs.fragments:
                if fragment in selected_systems:
                    self.allowed_sfos_override.extend([orb for orb in self.plot._selected_orbitals if isinstance(orb, pyfmo.orbitals.objects.SFO) and orb.fragment_unique == fragment])
                else:
                    self.allowed_sfos_override.extend(self.orbs.sfos.filter(fragment=fragment))
            self.filtered_orbitals = self.plot._selected_orbitals
        else:
            self.allowed_mos_override = None
            self.allowed_sfos_override = None
            self.filtered_orbitals = None

        self._update_plot()

        if self.parent.isDarkMode:
            if self.orbital_filter_button._is_checked:
                self.orbital_filter_button.setIcon(self.parent._ICONS['checked_dark'])
            else:
                self.orbital_filter_button.setIcon(self.parent._ICONS['unchecked_dark'])
        else:
            if self.orbital_filter_button._is_checked:
                self.orbital_filter_button.setIcon(self.parent._ICONS['checked'])
            else:
                self.orbital_filter_button.setIcon(self.parent._ICONS['unchecked'])

        self.orbital_filter_button.setEnabled(len(self.plot._selected_orbitals) > 0 or self.orbital_filter_button._is_checked)


    def load_analysis(self, file):
        try:
            self.orbs = pyfmo.Orbitals(file)
        except Exception as e:
            import traceback
            self.errordialog.showMessage(f'Could not load orbital data from file:\n\n{file}')
            traceback.print_exc(e)
            return

        self.main_mix = pyfmo.analysis.mixing.Mixer2(self.orbs, pr_min_thresh=0.001**2, oi_min_thresh=0.00000001)

        self._energytype_selection = 'energy'

        self.setAcceptDrops(False)
        if hasattr(self, "_new_page_frame"):
            self._new_page_frame.hide()
            self.central_layout.removeWidget(self._new_page_frame)

        is_charged = any(charge != 0 for charge in _determine_charges(self.orbs).values())
        if is_charged and 'site_energy' not in self.orbs.sfo_energy_types:
            QtWidgets.QMessageBox.warning(self, "Warning", "WARNING\nYou have charged fragments but the effective energies are not available!\n\n Rerun your calculation with SFOSiteEnergies or FMatSFO enabled.");
        
        if is_charged and 'site_energy' in self.orbs.sfo_energy_types:
            self._energytype_selection = 'site_energy'

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
        plot_container = QtWidgets.QFrame()
        plot_container.setLayout(plot_container_layout)
        plot_container.setFixedSize(700, 500)

        self.plot = MplCanvas(self)
        self.plot.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.plot.setFocus()
        plot_container.setStyleSheet('padding: 0px; margin: 0px; border: 1px solid lightgray; border-radius: 5px; background-color: none;')
        plot_container_layout.addWidget(self.plot, 0)
        layout.addWidget(plot_container, 0, 0, 1, 1, QtCore.Qt.AlignCenter)


        self.info_tabs = QtWidgets.QTabWidget()
        if self.parent.isDarkMode:
            self.info_tabs.setStyleSheet("""
                QTabWidget { 
                    border-radius: 5px; 
                    border: 1px solid darkgray
                    } 
                QTabWidget::pane { 
                    border: 1px solid darkgray; 
                    background-color: rgb(70, 70, 70);
                    border-radius: 5px; 
                    border-top-left-radius: 0px;
                    } 
                QTabWidget::tab-bar {
                    background-color: darkgray; 
                    border: 0px;
                    }
                """)
        else:
            self.info_tabs.setStyleSheet("""
                QTabWidget { 
                    border-radius: 5px; 
                    border: 1px solid lightgray
                    } 
                QTabWidget::pane { 
                    border: 1px solid lightgray; 
                    background-color: white;
                    border-radius: 5px; 
                    border-top-left-radius: 0px;
                    } 
                QTabWidget::tab-bar {
                    background-color: lightgray; 
                    border: 0px;
                    }
                """)
        if self.parent.isDarkMode:
            self.info_tabs.tabBar().setStyleSheet("""
                border-radius: 5px; 
                border: 1px solid darkgray; 
                background-color: #3A3A3A'
                """)
        else:
            self.info_tabs.tabBar().setStyleSheet("""
                border-radius: 5px; 
                border: 1px solid lightgray; 
                background-color: white'
                """)
        

        # self.info_tabs.setFixedSize(300, 500)
        orbital_info_frame = QtWidgets.QFrame()
        if self.parent.isDarkMode:
            orbital_info_frame.setStyleSheet("""
                QPushButton {
                    font-size: 12px;
                    border: 1px solid lightgray;
                    border-radius: 5px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: #3A3A3A;
                }
                QPushButton:hover {
                    background-color: ##656565;
                    }
                """)
        else:
            orbital_info_frame.setStyleSheet("""
                QPushButton {
                    font-size: 12px;
                    border: 1px solid lightgray;
                    border-radius: 5px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: white;
                }
                QPushButton:hover {
                    background-color: #f0f0f0;
                    }
                """)
        orbital_info_layout = QtWidgets.QGridLayout()
        orbital_info_frame.setLayout(orbital_info_layout)
        orbital_info_frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self.orbital_info_box = Spoilers(self)
        orbital_info_layout.addWidget(self.orbital_info_box, 1, 0, 1, 2)

        orbital_info_layout.addWidget(QtWidgets.QLabel('<i>Use</i> <b>Shift + Click</b> <i>to select multiple orbitals!</i>'), 0, 0, 1, 2)
        orbital_info_layout.addWidget(self.orbital_info_box, 1, 0, 1, 2)
        
        self.info_tabs.addTab(orbital_info_frame, 'Orbitals')

        system_info_box = Spoilers(self)
        system_info_box.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        system_info_box.addSpoiler('General', self._get_general_system_info(), self.parent._ICONS['info'])
        system_info_box.addSpoiler('Complex', self._get_complex_system_info(), self.parent._ICONS['mo'])
        for frag in self.orbs.fragments:
            system_info_box.addSpoiler(frag, self._get_fragment_system_info(frag), self.parent._ICONS['sfo'])

        # system_info_box.addSpoiler(self._get_system_info_txt())
        
        self.info_tabs.addTab(system_info_box, 'System')
        self.info_tabs.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        layout.addWidget(self.info_tabs, 0, 1)

        slider_layout = QtWidgets.QGridLayout()
        slider_box = QtWidgets.QFrame()
        slider_box.setObjectName('sliderbox')
        slider_box.setLayout(slider_layout)
        layout.addWidget(slider_box, 1, 0)

        self.OI_is_empty_label = QtWidgets.QLabel('⚠️')
        self.OI_is_empty_label.setToolTip('Could not find any Orbital Interactions for these settings.')
        slider_layout.addWidget(self.OI_is_empty_label, 0, 0)
        sp_retain = self.OI_is_empty_label.sizePolicy()
        sp_retain.setRetainSizeWhenHidden(True)
        self.OI_is_empty_label.setSizePolicy(sp_retain)

        self.cbox_OI = QtWidgets.QCheckBox('Show OI')
        self.cbox_OI.setChecked(True)
        slider_layout.addWidget(self.cbox_OI, 0, 1)
        self.cbox_OI.checkStateChanged.connect(self._update_plot)

        label_OI = QtWidgets.QLabel('τ<sub>oi</sub> =')
        label_OI.setStyleSheet('QLabel{ font: 12pt}')
        slider_layout.addWidget(label_OI, 0, 2)

        inc_oi_btn = QtWidgets.QPushButton('<')
        inc_oi_btn.clicked.connect(self._set_next_oi_slider)
        slider_layout.addWidget(inc_oi_btn, 0, 4)
        if self.parent.isDarkMode:
            slider_box.setStyleSheet("""
                QWidget#sliderbox{
                    padding: 0px; 
                    margin: 0px; 
                    border: 1px solid darkgray; 
                    border-radius: 5px; 
                    background-color: #3A3A3A;
                    }
                QPushButton {
                    font-size: 12px;
                    border: 1px solid darkgray;
                    border-radius: 13px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: #3A3A3A;
                }
                QPushButton:hover {
                    background-color: #656565;
                    }
                """)
        else:
            slider_box.setStyleSheet("""
                QWidget#sliderbox{
                    padding: 0px; 
                    margin: 0px; 
                    border: 1px solid lightgray; 
                    border-radius: 5px; 
                    background-color: white;
                    }
                QPushButton {
                    font-size: 12px;
                    border: 1px solid lightgray;
                    border-radius: 13px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: white;
                }
                QPushButton:hover {
                    background-color: #f0f0f0;
                    }
                """)
        self.slider_OI = QtWidgets.QSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        slider_OI_max = abs(min(min(v.values()) for v in self.main_mix.mixes['OI'].values()))
        self.slider_OI.setMinimum(np.log10(0.00000001) * slider_resolution)
        self.slider_OI.setMaximum(floor(np.log10(slider_OI_max) * slider_resolution))
        self.slider_OI.setSliderPosition(np.log10(slider_OI_max/1.5) * slider_resolution)
        self.main_mix.set_oi_threshold(slider_OI_max/1.5)
        slider_layout.addWidget(self.slider_OI, 0, 5)

        dec_oi_btn = QtWidgets.QPushButton('>')
        dec_oi_btn.clicked.connect(self._set_previous_oi_slider)
        slider_layout.addWidget(dec_oi_btn, 0, 6)

        label_value_OI = QtWidgets.QLabel(f'{10**(self.slider_OI.value()/slider_resolution):.2E}')
        label_value_OI.setStyleSheet('font: 10px "IBM Plex Mono"')
        slider_layout.addWidget(label_value_OI, 0, 3)
        self.slider_OI.valueChanged.connect(lambda value: (self._update_plot(), label_value_OI.setText(f'{10**(value/slider_resolution):.2E}')))

        self.PR_is_empty_label = QtWidgets.QLabel('⚠️')
        self.PR_is_empty_label.setToolTip('Could not find any Pauli Repulsions for these settings.')
        slider_layout.addWidget(self.PR_is_empty_label, 1, 0)

        self.cbox_PR = QtWidgets.QCheckBox('Show PR')
        self.cbox_PR.setChecked(True)
        slider_layout.addWidget(self.cbox_PR, 1, 1)
        self.cbox_PR.checkStateChanged.connect(self._update_plot)

        label_PR = QtWidgets.QLabel('τ<sub>pr</sub> =')
        label_PR.setStyleSheet('QLabel{ font: 12pt}')
        slider_layout.addWidget(label_PR, 1, 2)

        self.slider_PR = QtWidgets.QSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        # slider_PR_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in self.pauli_mixes.values())
        if len(self.main_mix.mixes['PR']) == 0:
            slider_PR_max = 0.001**2
        else:
            # slider_PR_max = max(self.main_mix.mixes['PR'].values())
            slider_PR_max = abs(max(max(v.values()) for v in self.main_mix.mixes['PR'].values()))

        inc_pr_btn = QtWidgets.QPushButton('<')
        inc_pr_btn.clicked.connect(self._set_next_pr_slider)
        slider_layout.addWidget(inc_pr_btn, 1, 4)
        self.slider_PR.setMinimum(0.001**2 * 1000 * slider_resolution)
        self.slider_PR.setMaximum(slider_PR_max * 1000 * slider_resolution)
        self.slider_PR.setSliderPosition(slider_PR_max/1.5 * 1000 * slider_resolution)
        self.main_mix.set_pr_threshold(slider_PR_max/1.1)
        slider_layout.addWidget(self.slider_PR, 1, 5)

        dec_pr_btn = QtWidgets.QPushButton('>')
        dec_pr_btn.clicked.connect(self._set_previous_pr_slider)
        slider_layout.addWidget(dec_pr_btn, 1, 6)
        label_value_PR = QtWidgets.QLabel(f'{self.slider_PR.value()/slider_resolution:.3f}')
        label_value_PR.setStyleSheet('font: 10px "IBM Plex Mono"')
        slider_layout.addWidget(label_value_PR, 1, 3)
        self.slider_PR.valueChanged.connect(lambda value: (self._update_plot(), label_value_PR.setText(f'{value/slider_resolution/1000:.4f}')))
        
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
        etype_btn.clicked.connect(self._energytype_selection_dialog.open)
        selector_layout.addWidget(etype_btn, 1, 0)

        # self._orb_selection_dialog = OrbitalSelectionDialog(self, self._orb_selection)
        self._orb_selection_dialog = orbital_selector.OrbitalSelectionDialog(self, self.orbs)

        orb_btn = QtWidgets.QPushButton('Orbitals')
        orb_btn.clicked.connect(self._orb_selection_dialog.open)
        selector_layout.addWidget(orb_btn, 1, 1)

        self.orbital_draw_button = QtWidgets.QPushButton()
        self.orbital_draw_button.setStyleSheet('QPushButton::menu-indicator { image: none; }')

        menu = QtWidgets.QMenu(self)
        self.orbital_draw_button.setMenu(menu)
        self.orbital_draw_button.setText('Draw')

        self.orbital_filter_button = QtWidgets.QPushButton('Filter')
        if self.parent.isDarkMode:
            self.orbital_filter_button.setIcon(self.parent._ICONS['unchecked_dark'])
        else:
            self.orbital_filter_button.setIcon(self.parent._ICONS['unchecked'])
        self.orbital_filter_button.clicked.connect(self._set_orbital_filter)
        self.orbital_filter_button.setEnabled(False)
        self.orbital_filter_button._is_checked = False
        self.allowed_mos_override = None
        self.allowed_sfos_override = None
        self.filtered_orbitals = None
        selector_layout.addWidget(self.orbital_draw_button, 0, 0)
        selector_layout.addWidget(self.orbital_filter_button, 0, 1)

        misc_box = QtWidgets.QFrame()
        misc_box_layout = QtWidgets.QGridLayout()
        misc_box.setLayout(misc_box_layout)
        layout.addWidget(misc_box, 2, 0)

        make_sheet_btn = QtWidgets.QPushButton('Generate Sheets')
        save_fig_btn = QtWidgets.QPushButton('Save Figure')
        make_sheet_btn.clicked.connect(self.get_sheets_save_file)
        save_fig_btn.clicked.connect(self.get_figure_save_file)
        misc_box_layout.addWidget(make_sheet_btn, 0, 0, 1, 1)
        misc_box_layout.addWidget(save_fig_btn, 0, 1, 1, 1)

        layout.setColumnStretch(0, 0)
        layout.setColumnStretch(1, 1)

        self._update_plot()


    def _get_general_system_info(self):
        frame = QtWidgets.QFrame()
        frame.setStyleSheet('QLabel{padding: 2px; font: 10pt} QPushButton{icon-size: 10px;}')
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

        # # level of theory
        # row += 1
        # layout.addWidget(QtWidgets.QLabel('<b>Level</b>'), row, 0, 1, 1)
        # layout.addWidget(QtWidgets.QLabel('BLYP-D3(BJ)/TZ2P(None)'), row, 1, 1, 1)

        # EDA terms
        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>E</i><sub>int</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'Bond Energy') * 627.503:.2f} kcal mol<sup>–1</sup>")
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

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Δ<i>V</i><sub>elstat</sub></b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f"{self.orbs.reader.read('Energy', 'elstat') * 627.503:.2f} kcal mol<sup>–1</sup>")
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

        charges = _determine_charges(self.orbs)
        unrestricted_mos = self.orbs.data['calc_info']['unrestricted_mos']

        all_spin_pols = self.orbs.data['calc_info']['sfo_spinpolarizations']
        total_spin_pols = 0
        for frag, frag_spin_pols in all_spin_pols.items():
            for irrep, spin_pols in frag_spin_pols.items():
                total_spin_pols += spin_pols[0] - spin_pols[1]

        row = 0
        layout.addWidget(CopyLabel('<b>Geometry (xyz)</b>', mol2xyz(self.orbs.molecule)), row, 0, 1, 2)

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
        charges = _determine_charges(self.orbs)
        unrestricted_sfos = self.orbs.data['calc_info']['unrestricted_sfos']

        row = 0
        layout.addWidget(CopyLabel('<b>Geometry (xyz)</b>', mol2xyz(sfos[0].molecule)), row, 0, 1, 2)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Charge</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(f'{round(charges[frag]):+d}')
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        row += 1
        layout.addWidget(QtWidgets.QLabel('<b>Restricted</b>'), row, 0, 1, 1)
        label = QtWidgets.QLabel(str(not unrestricted_sfos))
        label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        layout.addWidget(label, row, 1, 1, 1)

        frag_spin_pols = self.orbs.data['calc_info']['sfo_spinpolarizations'][frag]
        total_spin_pols = 0
        for irrep, spin_pols in frag_spin_pols.items():
            total_spin_pols += spin_pols[0] - spin_pols[1]

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
        label = QtWidgets.QLabel(str(round(sum(sfo.occupation for sfo in sfos))))
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
        d = os.path.join(os.path.split(self.orbs.kfpath)[0], 'pyorbb.xlsx')
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save File", dir=d, filter="XLSX file (*.xlsx);;Any file (*)")
        
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

    def setup_new(self):
        self._new_page_frame = QtWidgets.QFrame(self)
        self.central_layout.addWidget(self._new_page_frame, QtCore.Qt.AlignCenter)
        layout = QtWidgets.QVBoxLayout(self._new_page_frame)

        # Label
        label = QtWidgets.QLabel('<font size="50">Drop a File here</font>')
        label.setAlignment(QtCore.Qt.AlignCenter)
        label.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(label)

        # second label
        label = QtWidgets.QLabel('<font size="50" style="font-size:300%;color:grey;"><i>or</i></font>')
        label.setAlignment(QtCore.Qt.AlignCenter)
        label.setTextFormat(QtCore.Qt.RichText)
        layout.addWidget(label)

        # Button
        button = QtWidgets.QPushButton('Select a File')
        button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        button.setFlat(True)
        button.clicked.connect(self.open_filedialog)
        label.setTextFormat(QtCore.Qt.RichText)

        self._new_page_frame.setObjectName('ParentWidget')
        if self.parent.isDarkMode:
            self._new_page_frame.setStyleSheet("""
                QWidget#ParentWidget {
                    border: 3px dashed gray;
                    margin: 4px;
                    padding: 80px;
                    border-radius: 10px;
                    background-color: transparent;
                }
                QPushButton:hover {
                    background-color: #656565;
                }
                QPushButton {
                    font-size: 20px;
                    border: 1px solid gray;
                    border-radius: 10px;
                    padding: 8px;
                }
            """)
        else:
            self._new_page_frame.setStyleSheet("""
                QWidget#ParentWidget {
                    border: 3px dashed gray;
                    margin: 4px;
                    padding: 80px;
                    border-radius: 10px;
                    background-color: transparent;
                }
                QPushButton:hover {
                    background-color: #f0f0f0;
                }
                QPushButton {
                    font-size: 20px;
                    border: 1px solid gray;
                    border-radius: 10px;
                    padding: 8px;
                }
            """)
        layout.addWidget(button, alignment=QtCore.Qt.AlignCenter)


class PyOrbbApp(QtWidgets.QApplication):
    def __post_init__(self):
        fontpath = os.path.split(__file__)[0] + '/../cli_scripts/ibm_plex_mono/IBMPlexMono-Regular.ttf'
        QtGui.QFontDatabase.addApplicationFont(fontpath)

        self.window = QtWidgets.QMainWindow()
        self.window.resize(1030 + 22 + 12, 698 + 52)
        self.window.layout = QtWidgets.QGridLayout()
        grid_widget = QtWidgets.QWidget()
        if self.isDarkMode:
            grid_widget.setStyleSheet('''
                QPushButton:hover {
                    background-color: #656565;
                    }
                QPushButton {
                    border: 1px solid darkgray;
                    border-radius: 5px;
                    padding: 4px;
                    margin: 0px; 
                    background-color: #3A3A3A;
                    }
                ''')
        else:
            grid_widget.setStyleSheet('''
                QPushButton:hover {
                    background-color: #f0f0f0;
                    }
                QPushButton {
                    border: 1px solid lightgray;
                    border-radius: 5px;
                    padding: 4px;
                    margin: 0px; 
                    background-color: white;
                    }
                ''')

        grid_widget.setLayout(self.window.layout)
        self.window.setCentralWidget(grid_widget)

        # try to get the densf path from the environment
        # this will be None if it could not be found
        self._amsbin_loc = load_setting('amsbin')
        if 'AMSBIN' not in os.environ:
            os.environ['AMSBIN'] = self._amsbin_loc

        self.window.setWindowTitle("PyOrbb Analysis Tool")

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        add_tab_button = QtWidgets.QPushButton('+')
        # add_tab_button.setFlat(True)
        add_tab_button.resize(50, 50)
        add_tab_button.clicked.connect(self._add_analysis_tab)
        add_tab_button.setStyleSheet("""
            QPushButton {
                font-size: 20px;
                border-radius: 16px;
                padding: 6px;
                background-color: none;
                border: none;
            }
            """)
        add_tab_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        tab_shape = QtWidgets.QTabBar.Shape.RoundedNorth
        self.tabs.tabBar().setShape(tab_shape)
        self.tabs.setCornerWidget(add_tab_button, QtCore.Qt.TopLeftCorner)

        self.window.layout.addWidget(self.tabs, 0, 0, 1, 1)
        self.tabs.tabCloseRequested.connect(self.tabs.removeTab)
        self.tabs.tabBarDoubleClicked.connect(self._edit_tab_title)

        # File menu
        menuBar = self.window.menuBar()
        fileMenu = menuBar.addMenu("File")
        fileMenu.addAction("New")

        save = QtGui.QAction("Save",self)
        save.setShortcut("Ctrl+S")
        fileMenu.addAction(save)

        quit = QtGui.QAction("&Quit", self)
        quit.setShortcut("Ctrl+Q")
        fileMenu.addAction(quit)

        # Edit menu
        editMenu = menuBar.addMenu("Edit")
        editMenu.addAction("Copy")
        editMenu.addAction("Paste")

        # Help menu
        preferenceMenu = menuBar.addMenu("Preferences")
        open_ams_path = QtGui.QAction("Set AMS Path",self)
        open_ams_path.triggered.connect(self.AMS_loc_dialogue)
        preferenceMenu.addAction(open_ams_path)
        
        self.setStyle('Fusion')
        self._add_analysis_tab()

        ICON_FOLDER = os.path.join(os.path.split(__file__)[0], '..', 'application', 'icons')
        self._ICONS = {file.removesuffix('.png'): QtGui.QIcon(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}
        self._PIXMAPS = {file.removesuffix('.png'): QtGui.QPixmap(os.path.join(ICON_FOLDER, file)) for file in os.listdir(ICON_FOLDER)}
        self.window.setWindowIcon(self._ICONS["pyorbb"])


    @property
    def isDarkMode(self):
        return QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark


    def AMS_loc_dialogue(self):
        '''
        opens a filedialog allowing the user to select the path to AMSHOME
        '''
        if platform.system() == 'Darwin':
            d = "/Applications" if os.path.exists("/Applications") else os.getcwd()
            path = QtWidgets.QFileDialog.getOpenFileName(self.window, caption='Select AMS application', dir=d, filter="*.app")[0]
            self._amsbin_loc = os.path.join(path, 'Contents', 'Resources', 'amshome', 'bin')
        else:
            path = QtWidgets.QFileDialog.getExistingDirectory(caption='Select AMS install location', dir=os.getcwd())
            self._amsbin_loc = os.path.join(path, "bin")

        # if it is an app we get the amsbin
        save_setting('amsbin', self._amsbin_loc)
        os.environ['AMSBIN'] = self._amsbin_loc

    def _add_analysis_tab(self, object=None, tabname='new'):
        window = AnalysisWindow(self)
        window.setup_new()
        idx = self.tabs.addTab(window, tabname)
        self.tabs.setCurrentIndex(idx)

    def _edit_tab_title(self, index):
        def text_change_handler(arg):
            self.tabs.tabBar().setTabText(index, lineedit.text())
            rect = self.tabs.tabBar().tabRect(index)
            rect.adjust(34, 0.5, 1, 0.5)
            lineedit.setGeometry(rect)

        lineedit = QtWidgets.QLineEdit(parent=self.tabs)
        lineedit.textChanged.connect(text_change_handler)
        lineedit.editingFinished.connect(lambda: lineedit.hide())
        lineedit.setText(self.tabs.tabText(index))
        lineedit.setStyleSheet("border: 0px; background-color: rgba(0,0,0,0);")
        lineedit.show()
        text_change_handler("")
        lineedit.selectAll()

    def __enter__(self):
        self.__post_init__()
        return self

    def __exit__(self, *args):
        self.window.show()
        self.exec()
        self.shutdown()

