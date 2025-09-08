from PySide6 import *
import pyfmo
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backend_tools import Cursors
from matplotlib.figure import Figure
import numpy as np
import os

slider_resolution = 500


def _detect_charged_fragments(orbs):
    for frag in orbs.fragments:
        sfos = orbs.sfos.filter(fragment=frag)
        mol = sfos[0].molecule
        expected_Nelectrons = sum(atom.atnum for atom in mol)
        actual_Nelectrons = int(sum(sfo.occupation for sfo in sfos))
        if expected_Nelectrons != actual_Nelectrons:
            return True
    return False


class SpinSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent, state):
        super().__init__()
        self.parent = parent
        self.state = state
        layout = QtWidgets.QGridLayout(self)

        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('Select allowed spin-states:\n'), 0, 0, 1, 2)
        cbox_layout = QtWidgets.QGridLayout()
        layout.addLayout(cbox_layout, 1, 0, 1, 2)

        self.cboxes = {}
        for i, (key, val) in enumerate(state.items()):
            self.cboxes[key] = QtWidgets.QCheckBox()
            self.cboxes[key].setChecked(val)
            cbox_layout.addWidget(QtWidgets.QLabel(key), i, 0, 1, 1)
            cbox_layout.addWidget(self.cboxes[key], i, 1, 1, 1)

        # layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

        for key, val in self.cboxes.items():
            self.state[key] = val.isChecked()

        self.parent._update_plot()


class ETypeDialog(QtWidgets.QDialog):
    def __init__(self, parent, possibilities, selected):
        super().__init__()
        self.parent = parent
        self.possibilities = possibilities
        self.selected = selected
        layout = QtWidgets.QGridLayout(self)

        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('Select energy type for SFOs:\n'), 0, 0, 1, 2)
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

        if _detect_charged_fragments(self.parent.orbs):
            rbtn_layout.addWidget(QtWidgets.QLabel(f'\nNote:\nEffective energies are recommended for charge systems!'), i+1, 0, 1, 0)

        # layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()
        for key, val in self.rbuttons.items():
            if val.isChecked():
                self.parent._energytype_selection = key
        self.parent._update_plot()


class SymmSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent, state):
        super().__init__()
        self.parent = parent
        self.state = state
        layout = QtWidgets.QGridLayout(self)

        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('Select allowed irreps:\n'), 0, 0, 1, 2)

        tabs = QtWidgets.QTabWidget()
        layout.addWidget(tabs, 1, 0, 1, 2)

        self.cboxes = {}
        for column, col_state in state.items():
            self.cboxes[column] = {}
            col_frame = QtWidgets.QFrame()
            tabs.addTab(col_frame, column)

            cbox_layout = QtWidgets.QGridLayout()
            col_frame.setLayout(cbox_layout)
            for i, (key, val) in enumerate(col_state.items()):
                self.cboxes[column][key] = QtWidgets.QCheckBox()
                self.cboxes[column][key].setChecked(val)
                cbox_layout.addWidget(QtWidgets.QLabel(key), i, 0, 1, 1)
                cbox_layout.addWidget(self.cboxes[column][key], i, 1, 1, 1)

        # layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

        for column, col_cboxes in self.cboxes.items():
            for key, val in col_cboxes.items():
                self.state[column][key] = val.isChecked()

        self.parent._update_plot()
        # state = {key: val.isChecked() for key, val in self.cboxes.items()}
        # return state


class OrbitalSelectionDialog(QtWidgets.QDialog):
    def __init__(self, parent, state):
        super().__init__()
        self.parent = parent
        self.state = state
        layout = QtWidgets.QGridLayout(self)

        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('Select allowed orbitals:\n'), 0, 0, 1, 2)

        tabs = QtWidgets.QTabWidget()
        layout.addWidget(tabs, 1, 0, 1, 2)

        self.tables = {}
        for column, col_state in state.items():
            # col_frame = QtWidgets.QFrame()
            self.tables[column] = QtWidgets.QTableWidget(len(col_state), 3)
            self.tables[column].setHorizontalHeaderLabels(['Orbital', 'Spin', 'Relative Name'])
            self.tables[column].verticalHeader().setVisible(False)
            tabs.addTab(self.tables[column], column)

            # cbox_layout = QtWidgets.QGridLayout()
            # col_frame.setLayout(cbox_layout)
            for i, (key, val) in enumerate(col_state.items()):
                frame = QtWidgets.QFrame()
                layout_ = QtWidgets.QHBoxLayout()
                frame.setLayout(layout_)
                layout_.addWidget(QtWidgets.QCheckBox())
                layout_.addWidget(QtWidgets.QLabel(key.name))

                self.tables[column].setCellWidget(i, 0, frame)
                self.tables[column].setCellWidget(i, 1, QtWidgets.QLabel(key.spin))
                self.tables[column].setCellWidget(i, 2, QtWidgets.QLabel(key.relative_name))
                ...
                # self.tables[column][key] = QtWidgets.QCheckBox()
                # self.tables[column][key].setChecked(val)
                # cbox_layout.addWidget(QtWidgets.QLabel(key.name), i, 0, 1, 1)
                # cbox_layout.addWidget(self.tables[column][key], i, 1, 1, 1)

        # layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('cancel')
        cancel_btn.clicked.connect(self.reject)
        layout.addWidget(save_btn, 2, 0, 1, 1)
        layout.addWidget(cancel_btn, 2, 1, 1, 1)

    def open(self, *args):
        super().open()

        # wait until the dialog is done
        loop = QtCore.QEventLoop()
        self.finished.connect(loop.quit)
        loop.exec()

        for column, col_tables in self.tables.items():
            for key, val in col_tables.items():
                self.state[column][key] = val.isChecked()

        self.parent._update_plot()
        # state = {key: val.isChecked() for key, val in self.cboxes.items()}
        # return state




class FragRenameDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__()
        layout = QtWidgets.QGridLayout(self)
        self.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel('Rename orbital column name'), 0, 0, 1, 2)
        self._frag_rename_textedit = QtWidgets.QLineEdit(self)
        layout.addWidget(self._frag_rename_textedit, 1, 0, 1, 2)
        save_btn = QtWidgets.QPushButton('save')
        save_btn.clicked.connect(self.accept)
        cancel_btn = QtWidgets.QPushButton('cancel')
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



class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=9, height=6.5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.parent = parent
        super().__init__(self.fig)

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        self._frag_rename_dialog = FragRenameDialog(self)
        self._already_unfaded = True

    def on_click(self, event):

        # global screen
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)

        if event.dblclick:
            self.parent.new_tick_labels = []
            for artist in self.axes.get_xticklabels():
                if not artist.contains(event)[0]:
                    self.parent.new_tick_labels.append(artist.get_text())
                    continue
                new_txt = self._frag_rename_dialog.open(artist.get_text())
                self.parent.new_tick_labels.append(new_txt)
            self.axes.set_xticklabels(self.parent.new_tick_labels)
            self.fig.canvas.draw_idle()

        for artist in artists:
            gid = artist.get_gid()
            if gid is None:
                continue

            # Searching which data member corresponds to current mouse position
            if not artist.contains(event)[0]:
                continue

            if gid.startswith('MO_'):
                orb = self.parent.orbs.mos[gid[3:]]

            elif gid.startswith('SFO_'):
                orb = self.parent.orbs.sfos[gid[4:]]
            else:
                continue


            import tcviewer

            # orb.draw()
            if self.parent.tcviewer_screen is None or self.parent.tcviewer_screen.isclosed:
                self.parent.tcviewer_screen = tcviewer.screen._ScreenWindow()
                self.parent.tcviewer_screen.__enter__()
                self.parent.tcviewer_screen.show()

            with self.parent.tcviewer_screen.add_molscene() as scene:
            # scr.draw_cub(cub, isovalue, material=tcviewer.materials.orbital_shiny)            
                c1, c2 = ([1, 0, 0], [0, 0, 1]) if orb.occupied else ([1, .5, 0], [0, 1, 1])
                scene.draw_molecule(orb.molecule)
                scene.draw_dual_isosurface(orb.cube_file(), colorm=c1, colorp=c2)
                # scene.draw_isosurface(orb.cube_file(), -0.03, c1, opacity=.3)
                # scene.draw_isosurface(orb.cube_file(),  0.03, c2, opacity=.3)
                scene.draw_text(str(orb))


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

            if gid.startswith('MO_') or gid.startswith('SFO_'):
                self.fig.canvas.set_cursor(Cursors.HAND)

            s = '\n'
            if gid.startswith('MO_'):
                mo = self.parent.orbs.mos[gid[3:]]
                submixes = self.parent.main_mix.split()
                submix = [submix for submix in submixes if mo in submix.mos][0]
                s += 'MO'.ljust(35)
                s += f'\n   {mo}'
                s += f'\n   {mo.relative_name}'
                s += f'\n   {mo.symmetry} {mo.irrep_relative_name}\n'
                s += f'\nEnergy     {mo.energy:.2f} eV'
                s += f'\nOccupation {mo.occupation}'
                s += f'\nSpin       {mo.spin}'
                s += f'\nIrrep      {mo.symmetry}'

                s += '\n\nSFO                    Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for sfo in sorted(submix.sfos, key=lambda sfo: -abs(sfo.mulliken_contribution(mo))):
                    s += f'\n{str(sfo):19.19} {sfo.mulliken_contribution(mo): 8.2%} {sfo.coefficient(mo): 7.4f}'

                s += '\n' * (40 - len(s.splitlines()))
                self._fade_unrelated_ints(mo)
                self._already_unfaded = False

                self.parent.info_box.setText(s)
                self.fig.canvas.draw_idle()

                break

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos[gid[4:]]
                submixes = self.parent.main_mix.split()
                submix = [submix for submix in submixes if sfo in submix.sfos][0]
                s += 'SFO'.ljust(35)
                s += f'\n   {sfo}'
                s += f'\n   {sfo.relative_name}'
                s += f'\n   {sfo.symmetry} {sfo.irrep_relative_name}\n'
                s += f'\nFragment   {sfo.fragment_unique}'
                s += f'\nEnergy    {getattr(sfo, "energy"): .2f} eV'
                # s += f'\nEnergy    {getattr(sfo, self.parent.orbs.sfo_energy_types[etype_b.index_selected]): .2f} eV'
                s += f'\nOccupation {sfo.occupation:.2f}'
                s += f'\nPop.      {sfo.gross_population: .3f}'
                s += f'\nSpin-pop. {sfo.gross_spin: .3f}'
                s += f'\nSpin       {sfo.spin}'
                s += f'\nIrrep      {sfo.symmetry}'

                s += '\n\nSFO                      S   dE (eV)'
                s += '\n─────────────────── ────── ─────────'
                for sfo2 in sorted(submix.sfos, key=lambda sfo_: -abs(sfo @ sfo_)):
                    if sfo2 == sfo:
                        continue
                    if sfo2.fragment_unique == sfo.fragment_unique:
                        continue
                    s += f'\n{str(sfo2):19.19} {sfo @ sfo2: 5.3f} {abs(sfo.energy - sfo2.energy): 8.2f}'
                
                s += '\n\nMO                     Contr   Coeff'
                s += '\n─────────────────── ──────── ───────'
                for mo in sorted(submix.mos, key=lambda mo: -abs(sfo.mulliken_contribution(mo))):
                    s += f'\n{str(mo):19.19} {sfo.mulliken_contribution(mo): 8.2%} {sfo.coefficient(mo): 7.4f}'

                s += '\n' * (40 - len(s.splitlines()))
                self._fade_unrelated_ints(sfo)
                self._already_unfaded = False
                self.parent.info_box.setText(s)
                self.fig.canvas.draw_idle()

                break

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = self.parent.orbs.mos[gid[4:].split('->')[1].strip()]
                connected_sfos = [conn[0] for conn in self.parent.main_mix.connections if conn[1] == mo and conn[0].fragment_unique != sfo.fragment_unique]
                s += 'SFO'.ljust(35)
                s += f'\n   {sfo}'
                s += f'\n   {sfo.relative_name}'
                s += f'\n   {sfo.symmetry} {sfo.irrep_relative_name}\n'
                s += '\nMO'
                s += f'\n   {mo}'
                s += f'\n   {mo.relative_name}'
                s += f'\n   {mo.symmetry} {mo.irrep_relative_name}\n'
                s += f'\nContr.    {sfo.mulliken_contribution(mo): .2%}'
                s += f'\nCoeff.    {sfo.coefficient(mo): .6f}'
                s += f'\nSpin       {sfo.spin}'
                s += f'\nIrrep      {sfo.symmetry}'
                s += '\n\nSecond SFO          Bonding?'
                s += '\n─────────────────── ────────'
                for sfo2 in connected_sfos:
                    is_bonding = ((sfo @ sfo2) * sfo.coefficient(mo) * sfo2.coefficient(mo)) >= 0
                    s += f'\n{str(sfo2):19.19} {"   Yes  " if is_bonding else "    No    "}'

                s += '\n' * (40 - len(s.splitlines()))
                self._fade_unrelated_ints(mo)
                self._already_unfaded = False
                self.parent.info_box.setText(s)
                self.fig.canvas.draw_idle()

                break

        else:
            if not self._already_unfaded:
                self._already_unfaded = True
                self._unfade()

                # self.parent.info_box.setText((' '*35 + '\n')*40)
                self.fig.canvas.set_cursor(Cursors.POINTER)
                self.fig.canvas.draw_idle()

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

    def _fade_unrelated_ints(self, orb):
        artists = self.axes.get_children()
        artists = sorted(artists, key=lambda artist: artist.zorder)
        submixes = self.parent.main_mix.split()
        submix = [submix for submix in submixes if orb in submix.sfos or orb in submix.mos][0]
        connections = submix.find_closed_interactions(orb)
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
                if any(mo in conn for conn in connections):
                    continue

            if gid.startswith('SFO_'):
                sfo = self.parent.orbs.sfos[gid[4:]]
                if any(sfo in conn for conn in connections):
                    continue

            if gid.startswith('ARROWMO_'):
                mo = self.parent.orbs.mos[gid[8:]]
                if any(mo in conn for conn in connections):
                    continue

            if gid.startswith('ARROWSFO_'):
                sfo = self.parent.orbs.sfos[gid[9:]]
                if any(sfo in conn for conn in connections):
                    continue

            if gid.startswith('TEXTMO_'):
                mo = self.parent.orbs.mos[gid[7:]]
                if any(mo in conn for conn in connections):
                    continue

            if gid.startswith('TEXTSFO_'):
                sfo = self.parent.orbs.sfos[gid[8:]]
                if any(sfo in conn for conn in connections):
                    continue

            if gid.startswith('MIX_'):
                sfo = self.parent.orbs.sfos[gid[4:].split('->')[0].strip()]
                mo = self.parent.orbs.mos[gid[4:].split('->')[1].strip()]
                if any(mo in conn for conn in connections) and any(sfo in conn for conn in connections):
                    continue
            faded_artists.append(artist)

        for artist in faded_artists:
            artist.set_color('white')
            artist.set_alpha(1)
            self.axes.draw_artist(artist)
            self.fig.canvas.draw_idle()

            artist.set_color(artist.orig_color)
            artist.set_alpha(artist.orig_alpha * 0)
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
        super().open(parent, slot=QtCore.SLOT("get_file_from_dialog()"))
        return self.selectedFiles()[0]


class AnalysisWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setAcceptDrops(True)
        self.open_rkf_filedialog = QtWidgets.QFileDialog(self)
        self.open_rkf_filedialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        self.open_rkf_filedialog.setWindowTitle('Select adf.rkf file')
        self.open_rkf_filedialog.setFilter(QtCore.QDir.Filter.Files)
        self.open_rkf_filedialog.setNameFilters({"RKF file (*.rkf)", "Any file (*)"})

        self.errordialog = QtWidgets.QErrorMessage(self)
        self.warningdialog = QtWidgets.QMessageBox(self)
        self.warningdialog.setIcon(QtWidgets.QMessageBox.Warning);
        self.central_layout = QtWidgets.QVBoxLayout(self)

        self.new_tick_labels = None

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
        # self.load_analysis(file)
        print(file)

    def _update_plot(self):
        self._draw_diagram(
            allowed_spins=self._spin_selection,
            allowed_irreps=self._symmetry_selection,
            oi_thresh=10**(self.slider_OI.value()/slider_resolution),
            pauli_thresh=self.slider_PR.value()/slider_resolution,
            energy_type=self._energytype_selection,
            )


    def _draw_diagram(self, *args, allowed_spins=None, allowed_irreps=None, oi_thresh=None, pauli_thresh=None, ylim=None, energy_type=None):
        ax = self.plot.axes
        fig = self.plot.fig
        ax.clear()
        ax.yaxis.set_major_formatter('{x: 3.0f}')

        # global self.main_mix
        self.main_mix = pyfmo.analysis.mixing.Mixing(self.orbs, energy_type=energy_type)
        if self.cbox_OI.isChecked():
            for mix_ in self.oi_mixes[energy_type]:
                if not all(allowed_irreps['mos'][mo.symmetry] for mo in mix_.mos):
                    continue
                if not all(allowed_irreps[sfo.fragment_unique][sfo.subspecies] for sfo in mix_.sfos):
                    continue
                if not all(allowed_spins[mo.spin] for mo in mix_.mos):
                    continue

                if mix_.xiaobo_check(oi_thresh):
                    self.main_mix += mix_

        if self.cbox_PR.isChecked():
            for mix_ in self.pauli_mixes[energy_type]:
                if not all(allowed_irreps['mos'][mo.symmetry] for mo in mix_.mos):
                    continue
                if not all(allowed_irreps[sfo.fragment_unique][sfo.subspecies] for sfo in mix_.sfos):
                    continue
                if not all(allowed_spins[mo.spin] for mo in mix_.mos):
                    continue
                if mix_.xiaobo_check(pauli_thresh):
                    self.main_mix += mix_

        self.main_mix.sanitize()
        self.main_mix.draw_diagram(ax=ax, ylim=ylim)
        if self.new_tick_labels is not None:
            self.plot.axes.set_xticklabels(self.new_tick_labels)

        props = dict(edgecolor='white', facecolor='white', alpha=1)  # bbox features
        fig.canvas.draw_idle()


    def load_analysis(self, file):
        try:
            self.orbs = pyfmo.Orbitals(file)
        except:
            self.errordialog.showMessage(f'Could not load orbital data from file:\n\n{file}')
            return

        self._spin_selection = {}
        self._symmetry_selection = {}
        # self._symmetry_selection['mos'] = {}
        self._energytype_selection = 'energy'

        for spin in self.orbs.sfos.spins:
            self._spin_selection[spin] = True

        self._spin_selection_dialog = SpinSelectionDialog(self, self._spin_selection)

        for frag in self.orbs.fragments:
            self._symmetry_selection[frag] = {}
            for symm in sorted(set([sfo.subspecies for sfo in self.orbs.sfos.filter(fragment=frag)])):
                self._symmetry_selection[frag][symm] = True

        self._symmetry_selection['mos'] = {}
        for symm in sorted(set([mo.symmetry for mo in self.orbs.mos])):
                self._symmetry_selection['mos'][symm] = True

        self._symmetry_selection_dialog = SymmSelectionDialog(self, self._symmetry_selection)

        self.setAcceptDrops(False)
        if hasattr(self, "_new_page_frame"):
            self._new_page_frame.hide()
            self.central_layout.removeWidget(self._new_page_frame)

        is_charged = _detect_charged_fragments(self.orbs)
        if is_charged and 'site_energy' not in self.orbs.sfo_energy_types:
            self.warningdialog.showMessage("WARNING\nYou have charged fragments but the effective energies are not available!\n\n Rerun your calculation with SFOSiteEnergies or FMatSFO enabled.");
        
        if is_charged and 'site_energy' in self.orbs.sfo_energy_types:
            self._energytype_selection = 'site_energy'


        self._energytype_selection_dialog = ETypeDialog(self, self.orbs.sfo_energy_types, self._energytype_selection)

        self._orb_selection = {}
        for frag in self.orbs.fragments:
            self._orb_selection[frag] = {}
            for sfo in self.orbs.sfos.filter(fragment=frag):
                self._orb_selection[frag][sfo] = True

        self._orb_selection['mos'] = {}
        for mo in self.orbs.mos:
            self._orb_selection['mos'][mo] = True

        self._orb_selection_dialog = OrbitalSelectionDialog(self, self._orb_selection)

        # load a mixer object for each energy type we have available
        mixers = {etype: pyfmo.analysis.mixing.Mixer(self.orbs, energy_type=etype) for etype in self.orbs.sfo_energy_types}
        # and for each mixer generate 20 OI and PR interactions
        self.oi_mixes = {etype: mixer.orbital_interactions(N=40) for etype, mixer in mixers.items()}
        self.pauli_mixes = {etype: mixer.pauli_repulsions(N=40) for etype, mixer in mixers.items()}

        self._analysis_page_frame = QtWidgets.QFrame(self)
        self.central_layout.addWidget(self._analysis_page_frame)
        layout = QtWidgets.QGridLayout(self._analysis_page_frame)
        plot_container_layout = QtWidgets.QVBoxLayout()
        plot_container = QtWidgets.QFrame()
        plot_container.setLayout(plot_container_layout)
        plot_container.setFixedSize(700, 500)

        self.plot = MplCanvas(self)
        plot_container.setStyleSheet('padding: 0px; margin: 0px; border: 2px solid lightgray; border-radius: 5px; background-color: white;')
        plot_container_layout.addWidget(self.plot, 0)
        layout.addWidget(plot_container, 0, 0, 1, 1, QtCore.Qt.AlignCenter)

        self.info_box = QtWidgets.QLabel('')
        self.info_box.setFixedSize(288, 500)
        self.info_box.setStyleSheet('padding: 10px; font: 12pt "IBM Plex Mono"; border-radius: 5px; background-color: white; border: 2px solid lightgray;')
        layout.addWidget(self.info_box, 0, 1, QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)

        slider_layout = QtWidgets.QGridLayout()
        slider_box = QtWidgets.QFrame()
        slider_box.setObjectName('sliderbox')
        slider_box.setStyleSheet('QWidget#sliderbox{padding: 0px; margin: 0px; border: 2px solid lightgray; border-radius: 5px; background-color: white;}')
        slider_box.setLayout(slider_layout)
        layout.addWidget(slider_box, 1, 0)

        self.cbox_OI = QtWidgets.QCheckBox('Show')
        self.cbox_OI.setChecked(True)
        slider_layout.addWidget(self.cbox_OI, 0, 0, QtCore.Qt.AlignCenter)
        self.cbox_OI.checkStateChanged.connect(self._update_plot)

        label_OI = QtWidgets.QLabel('OI')
        slider_layout.addWidget(label_OI, 0, 1, QtCore.Qt.AlignCenter)

        self.slider_OI = QtWidgets.QSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        slider_OI_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in self.oi_mixes.values())
        self.slider_OI.setMinimum(np.log10(0.00000001) * slider_resolution)
        self.slider_OI.setMaximum(np.log10(slider_OI_max) * slider_resolution)
        self.slider_OI.setSliderPosition(np.log10(slider_OI_max/1.5) * slider_resolution)
        slider_layout.addWidget(self.slider_OI, 0, 2, QtCore.Qt.AlignCenter)

        label_value_OI = QtWidgets.QLabel(f'{10**(self.slider_OI.value()/slider_resolution):.2e}')
        slider_layout.addWidget(label_value_OI, 0, 3, QtCore.Qt.AlignCenter)
        self.slider_OI.valueChanged.connect(lambda value: (self._update_plot(), label_value_OI.setText(f'{10**(value/slider_resolution):.2e}')))

        self.cbox_PR = QtWidgets.QCheckBox('Show')
        self.cbox_PR.setChecked(True)
        slider_layout.addWidget(self.cbox_PR, 1, 0, QtCore.Qt.AlignCenter)
        self.cbox_PR.checkStateChanged.connect(self._update_plot)

        label_PR = QtWidgets.QLabel('PR')
        slider_layout.addWidget(label_PR, 1, 1, QtCore.Qt.AlignCenter)

        self.slider_PR = QtWidgets.QSlider(QtCore.Qt.Horizontal, self._analysis_page_frame)
        slider_PR_max = max(max(mix.xiaobo_value() for mix in mixes) for mixes in self.pauli_mixes.values())
        self.slider_PR.setMinimum(0.001 * slider_resolution)
        self.slider_PR.setMaximum(slider_PR_max * slider_resolution)
        self.slider_PR.setSliderPosition(slider_PR_max/1.5 * slider_resolution)
        slider_layout.addWidget(self.slider_PR, 1, 2, QtCore.Qt.AlignCenter)

        label_value_PR = QtWidgets.QLabel(f'{self.slider_PR.value()/slider_resolution:.3f}')
        slider_layout.addWidget(label_value_PR, 1, 3, QtCore.Qt.AlignCenter)
        self.slider_PR.valueChanged.connect(lambda value: (self._update_plot(), label_value_PR.setText(f'{value/slider_resolution:.2f}')))


        selector_box = QtWidgets.QFrame()
        selector_box.setStyleSheet('''
            QPushButton:hover {
                background-color: #f0f0f0;
                }
            QPushButton {
                border: 1px solid lightgray;
                border-radius: 5px;
                padding: 8px;
                margin: 0px; 
                background-color: white;
                }
            ''')
        selector_layout = QtWidgets.QGridLayout()
        selector_box.setLayout(selector_layout)
        layout.addWidget(selector_box, 1, 1)

        spin_btn = QtWidgets.QPushButton('Spin')
        spin_btn.clicked.connect(self._spin_selection_dialog.open)
        selector_layout.addWidget(spin_btn, 0, 0)

        symm_btn = QtWidgets.QPushButton('Irreps')
        symm_btn.clicked.connect(self._symmetry_selection_dialog.open)
        selector_layout.addWidget(symm_btn, 0, 1)

        etype_btn = QtWidgets.QPushButton('Energy Type')
        etype_btn.clicked.connect(self._energytype_selection_dialog.open)
        selector_layout.addWidget(etype_btn, 1, 0)

        orb_btn = QtWidgets.QPushButton('Orbitals')
        orb_btn.clicked.connect(self._orb_selection_dialog.open)
        selector_layout.addWidget(orb_btn, 1, 1)

        misc_box = QtWidgets.QFrame()

        misc_box.setStyleSheet('''
            QPushButton:hover {
                background-color: #f0f0f0;
                }
            QPushButton {
                border: 1px solid lightgray;
                border-radius: 5px;
                padding: 8px;
                margin: 0px; 
                background-color: white;
                }
            ''')
        misc_box_layout = QtWidgets.QGridLayout()
        misc_box.setLayout(misc_box_layout)
        layout.addWidget(misc_box, 2, 0)

        make_sheet_btn = QtWidgets.QPushButton('Generate Sheets')
        save_fig_btn = QtWidgets.QPushButton('Save Figure')
        make_sheet_btn.clicked.connect(self.get_sheets_save_file)
        save_fig_btn.clicked.connect(self.get_figure_save_file)
        misc_box_layout.addWidget(make_sheet_btn, 0, 0, 1, 1)
        misc_box_layout.addWidget(save_fig_btn, 0, 1, 1, 1)
        self._update_plot()

        print(self.size())

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
        _id = QtGui.QFontDatabase.addApplicationFont(fontpath)

        self.window = QtWidgets.QMainWindow()
        self.window.resize(1030 + 22, 698 + 52)
        self.window.layout = QtWidgets.QGridLayout()
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(self.window.layout)
        self.window.setCentralWidget(grid_widget)

        self.window.setWindowTitle("PyOrbb Analysis Tool")

        self.tabs = QtWidgets.QTabWidget()
        # self.tabs.setStyleSheet("""
        #     QTabBar::tab:selected {
        #         background-color: white;
        #         border: 1px lightgray solid;
        #         padding: 9px;
        #     }""")
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
            }
            """)
        add_tab_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))

        tab_shape = QtWidgets.QTabBar.Shape.RoundedNorth
        self.tabs.tabBar().setShape(tab_shape)
        self.tabs.setCornerWidget(add_tab_button, QtCore.Qt.TopLeftCorner)

        self.window.layout.addWidget(self.tabs, 0, 0, 1, 1)
        self.tabs.tabCloseRequested.connect(self.tabs.removeTab)
        self.tabs.tabBarDoubleClicked.connect(self._edit_tab_title)

        self.setStyle('Fusion')
        self._add_analysis_tab()

    def _add_analysis_tab(self, object=None, tabname='new'):
        window = AnalysisWindow()
        window.setup_new()
        idx = self.tabs.addTab(window, tabname)
        self.tabs.setCurrentIndex(idx)

    def _edit_tab_title(self, index):
        def text_change_handler(arg):
            self.tabs.tabBar().setTabText(index, lineedit.text())
            print(lineedit.text(), arg)
            rect = self.tabs.tabBar().tabRect(index)
            rect.adjust(34, 0.5, 1, 0.5)
            lineedit.setGeometry(rect)

        lineedit = QtWidgets.QLineEdit(parent=self.tabs)
        lineedit.textChanged.connect(text_change_handler)
        lineedit.editingFinished.connect(lambda: lineedit.hide())
        lineedit.setText(self.tabs.tabText(index))
        # rect = self.tabs.tabBar().tabRect(index)
        # rect.adjust(34, 0.5, 1, 0.5)
        # lineedit.setGeometry(rect)
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


if __name__ == '__main__':
    with PyOrbbApp():
        ...
