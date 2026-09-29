from PySide6 import QtWidgets, QtCore, QtGui
from . import rich_widgets, latex_renderer


class Spoilers(QtWidgets.QScrollArea):
    def _readBG(self):
        return self.background_color

    def _setBG(self, color):
        self.background_color = color
        self.setStyleSheet(f"QScrollArea {{ background-color: {color.name(QtGui.QColor.NameFormat.HexArgb)}; border: none; }}")

    _background_color_prop = QtCore.Property(QtGui.QColor, _readBG, _setBG)

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.setWidgetResizable(True)
        self.horizontalScrollBar().setEnabled(False)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        content = QtWidgets.QWidget(self)
        self.setWidget(content)
        self.spoilers = []

        self.background_color = QtGui.QColor('transparent')
        # vertical box layout
        self.layout = QtWidgets.QVBoxLayout(content)
        self.layout.addStretch(1)
        # self.setStyleSheet('background-color: transparent;')
        self._setBG(QtGui.QColor(7, 175, 213, 0))

    def themechange(self):
        for spoiler in self.spoilers:
            spoiler.set_pixmap()

    def addSpoiler(self, title, widget, icon=None, tooltip=None):
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(widget)
        spoiler = Spoiler(self, title, icon=icon, tooltip=tooltip)
        spoiler.setContentLayout(layout)
        self.spoilers.append(spoiler)
        self.layout.insertWidget(self.layout.count() - 1, spoiler)
        return len(self.spoilers) - 1

    def renameSpoiler(self, old_title, new_title):
        # check if there is already a spoiler with this title first:
        for spoiler in self.spoilers:
            if spoiler.title == new_title:
                raise ValueError('Cannot rename Spoiler. There is already one with the same name!')

        for spoiler in self.spoilers:
            if spoiler.title == old_title:
                spoiler.set_title(new_title)

    def empty(self):
        for i in range(self.layout.count()):
            self.remove_spoiler(0)
        self.layout.addStretch(1)

    def remove_spoiler(self, idx: int):
        widget = self.layout.takeAt(idx)

        if widget is None:
            return

        widget = widget.widget()
        if widget is None:
            return

        self.spoilers.remove(widget)
        widget.setParent(None)
        widget.destroy()

    def __len__(self):
        return len(self.spoilers)

    def emphasize(self):
        self.animation = QtCore.QPropertyAnimation(self, b"_background_color_prop")
        self.animation.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)
        self.animation.setDuration(1000)
        self.animation.setStartValue(QtGui.QColor(7, 175, 213, 100))
        self.animation.setEndValue(QtGui.QColor(7, 175, 213, 0))
        self.animation.start()



class Spoiler(QtWidgets.QWidget):
    def _readBG(self):
        return self.background_color

    def _setBG(self, color):
        self.background_color = color
        self.contentArea.setStyleSheet(f"QScrollArea {{ background-color: {color.name(QtGui.QColor.NameFormat.HexArgb)}; border: none; }}")

    _background_color_prop = QtCore.Property(QtGui.QColor, _readBG, _setBG)

    def __init__(self, parent=None, title='', animationDuration=100, icon=None, tooltip=None):
        """
        References:
            # Adapted from c++ version
            http://stackoverflow.com/questions/32476006/how-to-make-an-expandable-collapsable-section-widget-in-qt
        """
        super().__init__(parent=parent)

        self.title = title
        self.animationDuration = animationDuration
        self.toggleAnimation = QtCore.QParallelAnimationGroup()
        self.contentArea = QtWidgets.QScrollArea(self)
        if tooltip is not None:
            self.contentArea.setToolTip(tooltip)
        self.toggleButton = rich_widgets.HTMLToolButton(self)
        # self.toggleButton = QtWidgets.QToolButton(self)
        self.mainLayout = QtWidgets.QVBoxLayout()
        titleLayout = QtWidgets.QHBoxLayout()
        titleFrame = QtWidgets.QFrame()
        titleFrame.setLayout(titleLayout)

        toggleButton = self.toggleButton
        self.set_pixmap()
        # toggleButton.setStyleSheet("QToolButton { border: none; font-weight: bold; font-size: 20px; text-align: left top}")

        toggleButton.setToolButtonStyle(QtCore.Qt.ToolButtonFollowStyle)
        toggleButton.setArrowType(QtCore.Qt.NoArrow)
        # toggleButton.setText(f'{title}')

        # parent.parent.setWindowIcon(icon)
        # toggleButton.setIconSize(pixmap.size())

        toggleButton.setCheckable(True)
        toggleButton.setChecked(False)

        self.background_color = QtGui.QColor('transparent')
        self._last_toggle_forward_direction = False
        # self.contentArea.setStyleSheet("QScrollArea { background-color: transparent; border: none; }")
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
            titleLayout.addWidget(icon_lab, stretch=0)
        titleLayout.addWidget(self.toggleButton, stretch=1)
        # titleLayout.addStretch()

        mainLayout.addWidget(titleFrame)
        mainLayout.addWidget(self.contentArea)
        # mainLayout.addStretch()
        self.setLayout(self.mainLayout)

        self.toggleButton.clicked.connect(self.start_animation)

    def set_title(self, new_title: str):
        self.title = new_title
        self.set_pixmap()


    def start_animation(self, forward=True):
        if forward == self._last_toggle_forward_direction:
            return
        direction = QtCore.QAbstractAnimation.Forward if forward else QtCore.QAbstractAnimation.Backward
        self.toggleAnimation.setDirection(direction)
        self.toggleAnimation.start()
        self.toggleButton.setChecked(forward)
        self._last_toggle_forward_direction = forward

    def set_pixmap(self):
        darkmode = QtWidgets.QApplication.instance().isDarkMode
        pixmap = latex_renderer.convert_to_QPixMap("  " + self.title, darkmode=darkmode, fs=10)
        self.toggleButton.setPixmap(pixmap)

    def setContentLayout(self, contentLayout):
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

    def emphasize(self):
        self.animation = QtCore.QPropertyAnimation(self, b"_background_color_prop")
        self.animation.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)
        self.animation.setDuration(1000)
        self.animation.setStartValue(QtGui.QColor(7, 175, 213, 100))
        self.animation.setEndValue(QtGui.QColor(7, 175, 213, 0))
        self.animation.start()
