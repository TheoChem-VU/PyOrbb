from PySide6 import QtWidgets, QtCore, QtGui
import os
import platformdirs
import json
import platform
import pyfmo



class SettingSelectionWidget(QtWidgets.QFrame):
    default = None
    def setDefault(self, val):
        self.default = val


class LineEditFileDialogWidget(QtWidgets.QLineEdit):
    def __init__(self, parent=None, filetype="filename", **filedialog_settings):
        super(LineEditFileDialogWidget, self).__init__(parent)
        self.setReadOnly(True)

        icon = QtWidgets.QApplication.style().standardIcon(QtWidgets.QStyle.SP_DirIcon)
        action = self.addAction(icon, QtWidgets.QLineEdit.TrailingPosition)
        self.filedialog_settings = filedialog_settings
        self.filetype = filetype
        action.triggered.connect(self.select_file)

    def select_file(self):
        if self.filetype == 'filename':
            path = QtWidgets.QFileDialog(self).getOpenFileName(**self.filedialog_settings)[0]
        elif self.filetype == 'existingdirectory':
            path = QtWidgets.QFileDialog(self).getExistingDirectory(**self.filedialog_settings)
        if path == "":
            return
        self.setText(path)


class Path(SettingSelectionWidget):
    def __init__(self, parent, filetype=None):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        if platform.system() == 'Darwin':
            d = "/Applications" if os.path.exists("/Applications") else os.getcwd()
            self._filelineedit = LineEditFileDialogWidget(self, filetype="filename", caption='Select AMS application', dir=d, filter="*.app")
        else:
            self._filelineedit = LineEditFileDialogWidget(self, filetype="existingdirectory", caption='Select AMS install location', dir=os.getcwd())
        self.layout.addWidget(self._filelineedit)

    def setValue(self, val):
        self._filelineedit.setText(val)

    def value(self):
        return self._filelineedit.text()

    def reset(self):
        self.setValue(self.default)


class LineEditColorDialogWidget(QtWidgets.QLineEdit):
    def __init__(self, parent=None):
        super(LineEditColorDialogWidget, self).__init__(parent)

        icon = QtWidgets.QApplication.instance()._ICONS['rgb']
        action = self.addAction(icon, QtWidgets.QLineEdit.TrailingPosition)
        action.triggered.connect(self.select_color)

    def select_color(self):
        curr_color = QtGui.QColor(self.text())
        new_color = QtWidgets.QColorDialog.getColor(curr_color)
        if not new_color.isValid():
            return
        self.setText(new_color.name())


class Color(SettingSelectionWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        self._colorlineedit = LineEditColorDialogWidget()
        self.layout.addWidget(self._colorlineedit)

    def setValue(self, val):
        self._colorlineedit.setText(val)

    def value(self):
        return self._colorlineedit.text()

    def reset(self):
        self.setValue(self.default)


class String(SettingSelectionWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        self._lineedit = QtWidgets.QLineEdit()
        self.layout.addWidget(self._lineedit)

    def setValue(self, val):
        self._lineedit.setText(val)

    def value(self):
        return self._lineedit.text()

    def reset(self):
        self.setValue(self.default)


class CheckBox(SettingSelectionWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        self._checkbox = QtWidgets.QCheckBox()
        self.layout.addWidget(self._checkbox)

    def setValue(self, val):
        self._checkbox.setChecked(val)

    def value(self):
        return self._checkbox.isChecked()

    def reset(self):
        self.setValue(self.default)


class SpinBox(SettingSelectionWidget):
    def __init__(self, parent, minval=0, maxval=1, stepsize=0.1, decimals=1):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)

        self._spinbox = QtWidgets.QDoubleSpinBox(self)
        self._spinbox.setRange(minval, maxval)
        self._spinbox.setSingleStep(stepsize)
        self._spinbox.setDecimals(decimals)
        self.layout.addWidget(self._spinbox)

    def setValue(self, val):
        self._spinbox.setValue(val)

    def value(self):
        return self._spinbox.value()

    def reset(self):
        if self.default is None:
            return
        self.setValue(self.default)


class FloatLineEdit(SettingSelectionWidget):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)

        self._lineedit = QtWidgets.QLineEdit(self)
        # validator = QtGui.QDoubleValidator()
        # validator.setNotation(QtGui.QDoubleValidator.Notation.StandardNotation)
        # self._lineedit.setValidator(validator)
        self.layout.addWidget(self._lineedit)

    def setValue(self, val):
        self._lineedit.setText(str(val))

    def value(self):
        return float(self._lineedit.text())

    def reset(self):
        if self.default is None:
            return
        self.setValue(self.default)


class FloatTuple(SettingSelectionWidget):
    def __init__(self, parent, nfloats=None):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)

        self._lineedits = []
        for i in range(nfloats):
            le = QtWidgets.QLineEdit(self)
            self._lineedits.append(le)
            self.layout.addWidget(le)

    def setValue(self, val):
        [self._lineedits[i].setText(str(val[i])) for i in range(len(self._lineedits))]

    def value(self):
        return [float(self._lineedits[i].text()) for i in range(len(self._lineedits))]

    def reset(self):
        if self.default is None:
            return
        self.setValue(self.default)


class SettingsTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.get_funcs = {}
        self.set_funcs = {}
        self.set_default_funcs = {}
        self.reset_funcs = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.layout.addStretch()        
        pass

    def _add_generic_setting(self, name, variable_name, setting_widget, default):
        reset_btn = QtWidgets.QPushButton(self)
        reset_btn.clicked.connect(setting_widget.reset)
        if QtWidgets.QApplication.instance().isDarkMode:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset_dark'])
        else:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset'])

        setting_widget.setDefault(default)

        self.get_funcs[variable_name] = setting_widget.value
        self.set_funcs[variable_name] = setting_widget.setValue
        self.reset_funcs[variable_name] = setting_widget.reset
        self.set_default_funcs[variable_name] = setting_widget.setDefault

        frame = QtWidgets.QFrame(self)
        layout = QtWidgets.QHBoxLayout(frame)
        layout.addWidget(QtWidgets.QLabel(name))
        layout.addWidget(setting_widget)
        layout.addWidget(reset_btn)
        layout.setStretch(0, 1)
        layout.setStretch(1, 0)
        layout.setStretch(2, 0)
        # frame.setLayout(layout)
        self.layout.addWidget(frame)

    def add_float_setting(self, 
            variable_name, 
            name,
            minval=0, 
            maxval=1, 
            stepsize=0.01,
            decimals=4,
            use_spinbox=False,
            default=0.0):
        if use_spinbox:
            setting_widget = SpinBox(self, minval, maxval, stepsize, decimals)
        else:
            setting_widget = FloatLineEdit(self)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def add_float_tuple_setting(self, 
            variable_name, 
            name,
            nfloats=None,
            default=None):
        setting_widget = FloatTuple(self, nfloats)
        if default is None:
            default = tuple([0.0] * nfloats)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def add_path_setting(self, 
            variable_name, 
            name,
            default=None):
        setting_widget = Path(self)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def add_bool_setting(self, 
            variable_name, 
            name,
            default=False):
        setting_widget = CheckBox(self)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def add_color_setting(self, 
            variable_name, 
            name,
            default='#000000'):
        setting_widget = Color(self)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def add_str_setting(self, 
            variable_name, 
            name,
            default=''):
        setting_widget = String(self)
        self._add_generic_setting(name, variable_name, setting_widget, default)

    def reset(self):
        for reset in self.reset_funcs.values():
            reset()

    def get_state(self):
        return {variable_name: data_func() for variable_name, data_func in self.get_funcs.items()}

    def set_state(self, state: dict):
        for variable_name, value in state.items():
            self.set_funcs[variable_name](value)



class SettingsSection(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.tabs = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.build()

    def add_tab(self, name):
        self.tabs[name] = SettingsTab(self)
        return self.tabs[name]

    def build(self):
        if len(self.tabs) == 1:
            self.setLayout(list(self.tabs.values())[0].layout)
            return

        self.layout = QtWidgets.QHBoxLayout(self)
        tabs = QtWidgets.QTabWidget(self)
        tabs.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout.addWidget(tabs)

        for tab_name, tab_widget in self.tabs.items():
            tabs.addTab(tab_widget, tab_name)

    def reset(self):
        for tab in self.tabs.values():
            tab.reset()

    def get_state(self):
        state = {}
        for tab_name, tab in self.tabs.items():
            state[tab_name] = tab.get_state()
        return state

    def set_state(self, state):
        for tab_name, tab_state in state.items():
            self.tabs[tab_name].set_state(tab_state)


class SettingsDialog(QtWidgets.QDialog):
    settingsChanged = QtCore.Signal()

    def __init__(self, parent, title="Select Settings"):
        super().__init__(parent)
        self.parent = parent
        self.title = title
        self.setup()
        self.build()

    @property
    def count(self):
        return len(self.setting_widgets)

    def build(self):
        def change_widget(current, previous):
            new_index = list_widget.indexFromItem(current).row()
            stack_widget.setCurrentIndex(new_index)

        self.layout = QtWidgets.QVBoxLayout(self)

        # build the title label
        title_label = QtWidgets.QLabel(self.title)
        title_label.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.layout.addWidget(title_label, stretch=0)

        # build the section tabs
        section_tabs = QtWidgets.QTabWidget(self)
        section_tabs.setTabPosition(QtWidgets.QTabWidget.West)
        section_tabs.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout.addWidget(section_tabs)
        # and build the section widgets
        items = []
        for section_name, section_widget in self.sections.items():
            section_tabs.addTab(section_widget, section_name)
        section_tabs.setCurrentIndex(0)

        # also build the buttons for the dialog
        buttons_frame = QtWidgets.QFrame()
        buttons_layout = QtWidgets.QHBoxLayout()
        buttons_frame.setLayout(buttons_layout)

        apply_button = QtWidgets.QPushButton('Apply')
        apply_button.clicked.connect(self.apply)

        save_button = QtWidgets.QPushButton('Save')
        save_button.clicked.connect(self.save)

        cancel_button = QtWidgets.QPushButton('Cancel')
        cancel_button.clicked.connect(self.reject)

        reset_button = QtWidgets.QPushButton('Reset')
        reset_button.clicked.connect(self.reset)

        buttons_layout.addWidget(apply_button)
        buttons_layout.addWidget(save_button)
        buttons_layout.addWidget(cancel_button)
        buttons_layout.addWidget(reset_button)
        self.layout.addWidget(buttons_frame)

    def add_section(self, name):
        self.sections[name] = SettingsSection(self)
        return self.sections[name]

    def setup(self):
        self.sections = {}
        self.state = {}
        self.setting_widgets = {}  # dict of setting options name: (setting-widget, default value)

        with self.add_section('Densf') as section:
            with section.add_tab('General') as tab:
                tab.add_path_setting("amsbin", "AMS Application")

        with self.add_section('Plot') as section:
            with section.add_tab('Arrows') as tab:
                tab.add_float_setting("arrow_length", 'Length', default=0.25)
                tab.add_float_setting("arrow_width", 'Width', default=0.05)
                tab.add_float_setting("arrow_head_width", 'Head Width', default=0.15)
                tab.add_float_setting("arrow_head_length", 'Head Length', default=0.4)
                tab.add_float_setting("arrow_spacing", 'Spacing', default=0.06)
                tab.add_color_setting("arrow_color", 'Color', default='#000000')

            with section.add_tab('Labels') as tab:
                tab.add_bool_setting("draw_mo_labels", 'Show MO Labels', default=False)
                tab.add_bool_setting("draw_sfo_labels", 'Show FMO Labels', default=True)
                tab.add_float_setting("orb_label_offset", 'Label Offset', default=-0.28)

            with section.add_tab('Levels') as tab:
                tab.add_float_setting("level_width", 'Width', default=0.08)
                tab.add_float_setting("level_thickness", 'Thickness', default=3.0)
                tab.add_float_setting("highlight_thickness", 'Highlight Thickness', default=2.0)
                tab.add_str_setting("mo_column_name", 'MO Column Name', default="Complex")

            with section.add_tab('Connections') as tab:
                tab.add_color_setting("OI_color", 'Orbital Interactions Color', default="#008000")
                tab.add_color_setting("PR_color", 'Pauli Repulsion Color', default="#ff0000")
                tab.add_color_setting("Sanitization_color", 'Sanitization Color', default="#bf00bf")
                tab.add_color_setting("Multiple_color", 'Multiple Color', default="#000000")
                tab.add_float_tuple_setting("alpha_range", 'Alpha Range', 2, default=(0.1, 1.0))

        with self.add_section('Algorithm') as section:
            with section.add_tab('General') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('tab General'))

        with self.add_section('PyOrbb Viewer') as section:
            with section.add_tab('Grid') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('Grid Quality'))

        self.load_defaults()
        self.reset()
        self.load_state()

    def reset(self):
        for section in self.sections.values():
            section.reset()

    def save(self):
        super().accept()
        self.settingsChanged.emit()
        self.write_state()

    def apply(self):
        self.settingsChanged.emit()

    def reject(self):
        self.set_state(self._old_state)
        self.settingsChanged.emit()
        super().reject()

    def get_state(self):
        state = {}
        for section_name, section in self.sections.items():
            state[section_name] = {}
            for tab_name, tab in section.tabs.items():
                state[section_name][tab_name] = tab.get_state()
        return state

    def get_flat_state(self):
        state = {}
        for section_name, section in self.sections.items():
            for tab_name, tab in section.tabs.items():
                state.update(tab.get_state())
        return state

    def get(self, section, tab, variable):
        return self.get_state()[section][tab][variable]

    def set(self, section, tab, variable, value):
        return self.sections[section].tabs[tab].set_funcs[variable](value)

    def set_state(self, state):
        for section_name, section_state in state.items():
            self.sections[section_name].set_state(section_state)

    def write_state(self):
        d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
        with open(os.path.join(d, f'settings_{pyfmo.__version__}.json'), 'w+') as jf:
            state = self.get_state()
            jf.write(json.dumps(state, indent=4))

    def load_state(self):
        d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
        if not os.path.exists(os.path.join(d, f'settings_{pyfmo.__version__}.json')):
            return

        with open(os.path.join(d, f'settings_{pyfmo.__version__}.json')) as jf:
            try:
                self.set_state(json.loads(jf.read()))
            except:
                raise
                self.reset()

    def load_defaults(self):
        d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
        if not os.path.exists(os.path.join(d, f'default_settings_{pyfmo.__version__}.json')):
            return

        with open(os.path.join(d, f'default_settings_{pyfmo.__version__}.json')) as jf:
            defaults = json.loads(jf.read())

        for section_name, section in self.sections.items():
            if section_name not in defaults:
                continue

            for tab_name, tab in section.tabs.items():
                if tab_name not in defaults[section_name]:
                    continue

                for variable_name, func in tab.set_default_funcs.items():
                    val = defaults[section_name][tab_name].get(variable_name, None)
                    func(val)

    def exec(self):
        self._old_state = self.get_state()
        super().exec()
