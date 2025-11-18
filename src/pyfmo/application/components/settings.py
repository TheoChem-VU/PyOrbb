from PySide6 import QtWidgets, QtCore, QtGui
import os
import platformdirs
import json
import platform



class SettingSelectionWidget(QtWidgets.QFrame):
    ...


class LineEditFileDialogWidget(QtWidgets.QLineEdit):
    def __init__(self, parent=None, filetype="filename", **filedialog_settings):
        super(LineEditFileDialogWidget, self).__init__(parent)
        self.setReadOnly(True)

        icon = QtWidgets.QApplication.style().standardIcon(QtWidgets.QStyle.SP_DirIcon)
        self.action = self.addAction(icon, QtWidgets.QLineEdit.TrailingPosition)
        self.filedialog_settings = filedialog_settings
        self.filetype = filetype
        self.action.triggered.connect(self.select_file)

    def select_file(self):
        if self.filetype == 'filename':
            path = QtWidgets.QFileDialog(self).getOpenFileName(**self.filedialog_settings)[0]
        elif self.filetype == 'existingdirectory':
            path = QtWidgets.QFileDialog(self).getExistingDirectory(**self.filedialog_settings)
        if path == "":
            return
        self.setText(path)


class Path(SettingSelectionWidget):
    def __init__(self, parent, default="", filetype=None):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        self.setLayout(self.layout)
        if platform.system() == 'Darwin':
            d = "/Applications" if os.path.exists("/Applications") else os.getcwd()
            self._filelineedit = LineEditFileDialogWidget(self, filetype="filename", caption='Select AMS application', dir=d, filter="*.app")
        else:
            self._filelineedit = LineEditFileDialogWidget(self, filetype="existingdirectory", caption='Select AMS install location', dir=os.getcwd())
        self.layout.addWidget(self._filelineedit)
        self.default = default
        self.reset()

    def setValue(self, val):
        self._filelineedit.setText(val)

    def value(self):
        return self._filelineedit.text()

    def reset(self):
        self.setValue(self.default)


class SpinBox(SettingSelectionWidget):
    def __init__(self, parent, minval=0, maxval=1, stepsize=0.1, decimals=1, default=None):
        super().__init__(parent)
        self.parent = parent
        self.layout = QtWidgets.QHBoxLayout(self)
        self.setLayout(self.layout)

        self._spinbox = QtWidgets.QDoubleSpinBox(self)
        self._spinbox.setRange(minval, maxval)
        self._spinbox.setSingleStep(stepsize)
        self._spinbox.setDecimals(decimals)
        self.layout.addWidget(self._spinbox)
        self.default = default
        self.reset()

    def setValue(self, val):
        self._spinbox.setValue(val)

    def value(self):
        return self._spinbox.value()

    def reset(self):
        if self.default is None:
            return
        self.setValue(self.default)


class SettingsTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.layout = QtWidgets.QVBoxLayout(self)
        self.get_funcs = {}
        self.set_funcs = {}
        self.reset_funcs = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.layout.addStretch()        
        pass

    def add_float_setting(self, 
            variable_name, 
            name, 
            default, 
            minval=0, 
            maxval=1, 
            stepsize=0.01,
            decimals=4):
        sb = SpinBox(self, minval, maxval, stepsize, decimals, default=default)
        reset_btn = QtWidgets.QPushButton(self)
        reset_btn.clicked.connect(sb.reset)
        if QtWidgets.QApplication.instance().isDarkMode:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset_dark'])
        else:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset'])

        self.get_funcs[variable_name] = sb.value
        self.set_funcs[variable_name] = sb.setValue
        self.reset_funcs[variable_name] = sb.reset

        layout = QtWidgets.QHBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(name))
        layout.addWidget(sb)
        layout.addWidget(reset_btn)
        layout.setStretch(0, 1)
        layout.setStretch(1, 0)
        layout.setStretch(2, 0)
        frame = QtWidgets.QFrame(self)
        frame.setLayout(layout)
        self.layout.addWidget(frame)

    def add_path_setting(self, 
            variable_name, 
            name, 
            default=""):
        p = Path(self, default=default)
        reset_btn = QtWidgets.QPushButton(self)
        reset_btn.clicked.connect(p.reset)
        if QtWidgets.QApplication.instance().isDarkMode:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset_dark'])
        else:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset'])

        self.get_funcs[variable_name] = p.value
        self.set_funcs[variable_name] = p.setValue
        self.reset_funcs[variable_name] = p.reset

        layout = QtWidgets.QHBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel(name))
        layout.addWidget(p)
        layout.addWidget(reset_btn)
        layout.setStretch(0, 1)
        layout.setStretch(1, 0)
        layout.setStretch(2, 0)
        frame = QtWidgets.QFrame(self)
        frame.setLayout(layout)
        self.layout.addWidget(frame)

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
        self.setLayout(self.layout)
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
        self.setLayout(self.layout)

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
                default = ""
                if 'AMSBIN' in os.environ:
                    if '.app' in os.environ['AMSBIN']:
                        default = os.environ['AMSBIN'].split('.app')[0] + '.app'

                tab.add_path_setting("amsbin", "AMS Application", default)

        with self.add_section('Plot') as section:
            with section.add_tab('Arrows') as tab:
                tab.add_float_setting("arrow_length", 'Length', round(.3 / 4.8280888207, 3))
                tab.add_float_setting("arrow_thickness", 'Thickness', .35)
                tab.add_float_setting("arrow_width", 'Width', .005)
                tab.add_float_setting("arrow_head_width", 'Head Width', .025)
                tab.add_float_setting("arrow_head_length", 'Head Length', round(.1 / 4.8280888207, 3))
                tab.add_float_setting("arrow_overhang", 'Overhang', .4)
                tab.add_float_setting("arrow_spacing", 'Spacing', .012)   

            with section.add_tab('Labels') as tab:
                ...

            with section.add_tab('Levels') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('tab Special'))

                # "level_width"
                # "level_thickness"
                # "highlight_thickness"

            with section.add_tab('Connections') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('tab Special'))

        with self.add_section('Algorithm') as section:
            with section.add_tab('General') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('tab General'))

        with self.add_section('PyOrbb Viewer') as section:
            with section.add_tab('Grid') as tab:
                tab.layout.addWidget(QtWidgets.QLabel('Grid Quality'))

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

    def set_state(self, state):
        for section_name, section_state in state.items():
            self.sections[section_name].set_state(section_state)

    def write_state(self):
        d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
        with open(os.path.join(d, 'settings.json'), 'w+') as jf:
            state = self.get_state()
            jf.write(json.dumps(state))

    def load_state(self):
        d = platformdirs.user_config_dir('PyOrbb', 'TheoCheM', ensure_exists=True)
        print(d)
        if not os.path.exists(os.path.join(d, 'settings.json')):
            return

        with open(os.path.join(d, 'settings.json')) as jf:
            try:
                self.set_state(json.loads(jf.read()))
            except:
                raise
                self.reset()

    def exec(self):
        self._old_state = self.get_state()
        super().exec()
