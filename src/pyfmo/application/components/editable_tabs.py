from PySide6 import QtWidgets, QtCore, QtGui
from time import perf_counter


class LineEditDialog(QtWidgets.QDialog):
    def __init__(self, parent, title=None):
        super().__init__(parent=parent)
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)

        if title is not None:
            title_label = QtWidgets.QLabel(title)
            layout.addWidget(title_label)

        self.lineedit = QtWidgets.QLineEdit()
        layout.addWidget(self.lineedit)

        buttons = QtWidgets.QFrame()
        buttons_layout = QtWidgets.QHBoxLayout()
        buttons.setLayout(buttons_layout)

        save_button = QtWidgets.QPushButton('Save')
        save_button.clicked.connect(self.accept)
        cancel_button = QtWidgets.QPushButton('Cancel')
        cancel_button.clicked.connect(self.reject)

        buttons_layout.addWidget(save_button)
        buttons_layout.addWidget(cancel_button)

        layout.addWidget(buttons)

    def open(self, default_text=""):
        self.lineedit.setText(default_text)
        self.lineedit.setFocus()
        self.lineedit.selectAll()
        code = self.exec()
        if code == QtWidgets.QDialog.Accepted:
            return self.lineedit.text()
        else:
            return None


class WindowTabs(QtWidgets.QTabWidget):
    def __init__(self, parent):
        super().__init__(parent=parent)
        self.parent = parent
        # self.setMovable(True)
        # self.setMouseTracking(True)
        self.setTabBar(WindowTabBar(self))
        # self.setTabsClosable(True)
        self.tabBarDoubleClicked.connect(self.rename_tab)
        self.rename_tab_dialog = LineEditDialog(self.parent, 'Rename tab to:')
        self.tabCloseRequested.connect(self.removeTab)

    def rename_tab(self, index):
        current_tab_name = self.tabText(index)
        new_tab_name = self.rename_tab_dialog.open(current_tab_name)
        
        if new_tab_name is not None:
            self.setTabText(index, new_tab_name)



class WindowTabBar(QtWidgets.QTabBar):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(parent)

        self.setMovable(True)
        self.setMouseTracking(True)
        self.setTabsClosable(True)

        # self.dragStartPos = QtCore.QPoint()
        # self.dragDropedPos = QtCore.QPoint()
        # self.mouseCursor = QtGui.QCursor()
        self.dragInitiated = False
        self.dragLabel = None
        self.mouseCrossedWindowTime = 0
        self.targetLabelShrinkage = None
        self.labelShrinkage = 100
        self.timerID = None
        self.mouseLeftWindow = False
        self.mouseLeftTabBar = False


    def timerEvent(self, event=None):
        super().timerEvent(event)

        index = self.parent.currentIndex()
        widg = self.parent.widget(index)
        # pixmap = widg.grab()

        time_since_crossed = perf_counter() - self.mouseCrossedWindowTime
        self.labelShrinkage = self.labelShrinkage + (self.targetLabelShrinkage - self.labelShrinkage) * time_since_crossed * 2

        # pixmap = pixmap.scaled(pixmap.width()/self.labelShrinkage, pixmap.height()/self.labelShrinkage)

        pixmap = self.orig_pixmap.scaled(self.orig_pixmap.width()/self.labelShrinkage, self.orig_pixmap.height()/self.labelShrinkage)
        rect = pixmap.rect()
        # make the pixmap transparent
        painter = QtGui.QPainter()
        painter.begin(pixmap)
        painter.setCompositionMode(QtGui.QPainter.CompositionMode_DestinationIn)
        painter.fillRect(pixmap.rect(), QtGui.QColor(0, 0, 0, 100))
        painter.end()

        # self.dragLabel.setPixmap(pixmap)
        # pixmap = self.dragLabel.pixmap()

        # self.dragLabel.resize(widg.width()/self.labelShrinkage, widg.height()/self.labelShrinkage)
        self.dragLabel.updateGeometry()
        self.dragLabel.setPixmap(pixmap)

        # this removes the window frame
        self.dragLabel.setWindowFlags(QtCore.Qt.CustomizeWindowHint)
        self.dragLabel.show()

    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)

        event.accept()
        if event.buttons() == QtCore.Qt.MouseButton.LeftButton:
            self.dragInitiated = True
            # Convert the move event into a drag
            #Create the appearance of dragging the tab content
            if self.dragLabel is None:
                self.dragLabel = QtWidgets.QLabel()
                index = self.parent.currentIndex()
                widg = self.parent.widget(index)
                pixmap = widg.grab()
                self.orig_pixmap = pixmap
                self.dragLabel.setPixmap(pixmap)
                # this makes the label transparent to mouse
                self.dragLabel.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents, True);
                self.timerID = self.startTimer(10)
        else:
            self.dragInitiated = False
            self.mouseLeftTabBar = False
            self.mouseLeftWindow = False
            if self.timerID is not None:
                self.killTimer(self.timerID)
                self.timerID = None
            # self.dragLabel = None
            return

        if self.dragInitiated:
            # if the mouse has left the tabbar:
            if not self.rect().contains(event.pos()):
                # if the mouse left the window
                if not self.parent.parent.rect().contains(event.pos()):
                    if self.targetLabelShrinkage != 6:
                        self.targetLabelShrinkage = 6
                        self.mouseCrossedWindowTime = perf_counter()
                        self.mouseLeftTabBar = True
                        self.mouseLeftWindow = True
                else:
                    if self.targetLabelShrinkage != 3:
                        self.targetLabelShrinkage = 3
                        self.mouseCrossedWindowTime = perf_counter()
                        self.mouseLeftTabBar = True
                        self.mouseLeftWindow = False
            else:
                if self.targetLabelShrinkage != 3:
                    self.mouseLeftTabBar = False
                    self.mouseLeftWindow = False
                    self.targetLabelShrinkage = 3
                    self.mouseCrossedWindowTime = perf_counter()

            self.dragLabel.move(event.globalPos())

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        
        if self.dragLabel is not None:
            self.dragLabel.hide()
        self.dragLabel = None
        self.targetLabelShrinkage = 100
        self.labelShrinkage = 100
        self.dragInitiated = False

        if self.timerID is not None:
            self.killTimer(self.timerID)
            self.timerID = None

        if self.mouseLeftWindow:
            # if the mouse is out of this window we check fi we need to make a new window
            # or add it to an existing one
            new_window = False
            window = self.get_moused_over_window(event)
            if window is None:
                new_window = True
                window = QtWidgets.QApplication.instance().add_window()

            index = self.parent.currentIndex()
            widg = self.parent.widget(index)
            new_idx = window.tabs.addTab(widg, self.parent.tabText(index))
            window.tabs.setCurrentIndex(new_idx)
            widg.setEnabled(True)
            # a new window will have by default one tab open already
            if new_window:
                window.tabs.removeTab(0)
                window.move(event.globalPos())

            if self.parent.parent.tabs.count() == 0:
                QtWidgets.QApplication.instance().remove_window(self.parent.parent)

    def get_moused_over_window(self, event):
        for window in QtWidgets.QApplication.instance().windows:
            if window is self.parent.parent:
                continue
            rect = window.rect()
            if window.rect().contains(window.mapFromGlobal(event.globalPos())):
                return window
