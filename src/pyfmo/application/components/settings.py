from PySide6 import QtWidgets, QtCore, QtGui



class SettingSelectionWidget(QtWidgets.QFrame):
    ...

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
        self.data_funcs = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.layout.addStretch()        
        pass

    def add_float_setting(self, 
            name, 
            variable_name, 
            default, 
            minval=0, 
            maxval=1, 
            stepsize=0.01,
            decimals=3):
        sb = SpinBox(self, minval, maxval, stepsize, decimals, default=default)
        reset_btn = QtWidgets.QPushButton(self)
        if QtWidgets.QApplication.instance().isDarkMode:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset_dark'])
        else:
            reset_btn.setIcon(QtWidgets.QApplication.instance()._ICONS['reset'])
        self.data_funcs[variable_name] = sb.value

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
        self.layout.addWidget(tabs)

        for tab_name, tab_widget in self.tabs.items():
            tabs.addTab(tab_widget, tab_name)
        self.layout.addStretch()        


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

        title_label = QtWidgets.QLabel(self.title)
        title_label.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(title_label, stretch=0)
        self.setLayout(self.layout)

        list_widget = QtWidgets.QListWidget(self)
        stack_widget = QtWidgets.QStackedWidget(self)

        splitter = QtWidgets.QSplitter(self)
        splitter.addWidget(list_widget)
        splitter.addWidget(stack_widget)
        self.layout.addWidget(splitter)

        for section_name, section_widget in self.sections.items():
            list_widget.addItem(section_name)
            stack_widget.addWidget(section_widget)

        list_widget.currentItemChanged.connect(change_widget)

    def add_section(self, name):
        self.sections[name] = SettingsSection(self)
        return self.sections[name]

    def setup(self):
        self.sections = {}
        self.state = {}
        self.setting_widgets = {}  # dict of setting options name: (setting-widget, default value)

        with self.add_section('Plot') as section:
            with section.add_tab('Arrows') as tab:

                tab.add_float_setting("arrow_length", 'Length', .3 / 4.8280888207)
                tab.add_float_setting("arrow_thickness", 'Thickness', .35)
                tab.add_float_setting("arrow_width", 'Width', .005)
                tab.add_float_setting("arrow_head_width", 'Head Width', .025)
                tab.add_float_setting("arrow_head_length", 'Head Length', .1 / 4.8280888207)
                tab.add_float_setting("arrow_overhang", 'Overhang', .4)
                tab.add_float_setting("arrow_spacing", 'Spacing', .012)    

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



"""
Settings
  |- Section1
     |- Tab1
        |- Key1: value1
        |- Key2: value2
  |- Section2
    ...
"""