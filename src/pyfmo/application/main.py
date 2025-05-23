try:
    from PySide6 import *
    has_qt = True
except ImportError:
    has_qt = False


class MainApp(QtWidgets.QApplication):
    def __post_init__(self):
        self.setStyle('Fusion')
        self.window = QtWidgets.QMainWindow()
        self.window.layout = QtWidgets.QGridLayout()
        grid_widget = QtWidgets.QWidget()
        grid_widget.setLayout(self.window.layout)
        self.window.setCentralWidget(grid_widget)
        self.window.setWindowTitle("TCViewer 2.0") 

        self.setup_tabs()

    def setup_tabs(self):
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setUsesScrollButtons(True)
        self.window.layout.addWidget(self.tabs, 0, 0, 0, 0)
        self._add_tab('test')

    def __enter__(self):
        self.__post_init__()
        return self

    def __exit__(self, *args):
        self.window.show()
        self.exec()
        self.shutdown()

    def _add_tab(self, name):
        page = QtWidgets.QWidget()
        page.layout = QtWidgets.QGridLayout()

        self.tabs.addTab(page, name)

if __name__ == '__main__':
	with MainApp() as app:
		...
