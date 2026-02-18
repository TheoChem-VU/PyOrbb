from PySide6 import QtWidgets, QtCore
from . import rich_widgets, latex_renderer

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
        self.spoilers = {}

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
        self.spoilers[title] = spoiler
        self.layout.insertWidget(self.layout.count() - 1, spoiler)

    def renameSpoiler(self, old_title, new_title):
        # self.spoilers[old_title].toggleButton.setText(new_title)
        pixmap = latex_renderer.convert_to_QPixMap("  " + new_title, darkmode=self.parent.parent.isDarkMode, fs=10)
        self.spoilers[old_title].toggleButton.setPixmap(pixmap)
        # self.spoilers[old_title].toggleButton.__lbl.setFixedHeight(pixmap.height())
        # self.spoilers[old_title].toggleButton.setIconSize(pixmap.height()//4, pixmap.width()//4)
        self.spoilers[old_title].toggleButton.update()
        self.spoilers[new_title] = self.spoilers.pop(old_title)

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

    def __len__(self):
        return len(self.spoilers)


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
        self.toggleButton = rich_widgets.HTMLToolButton(self)
        # self.toggleButton = QtWidgets.QToolButton(self)
        self.mainLayout = QtWidgets.QVBoxLayout()
        titleLayout = QtWidgets.QHBoxLayout()
        titleFrame = QtWidgets.QFrame()
        titleFrame.setLayout(titleLayout)

        toggleButton = self.toggleButton
        # toggleButton.setStyleSheet("QToolButton { border: none; font-weight: bold; font-size: 20px; text-align: left top}")

        toggleButton.setToolButtonStyle(QtCore.Qt.ToolButtonFollowStyle)
        toggleButton.setArrowType(QtCore.Qt.NoArrow)
        # toggleButton.setText(f'{title}')

        pixmap = latex_renderer.convert_to_QPixMap("  " + title, darkmode=parent.parent.parent.isDarkMode, fs=10)
        toggleButton.setPixmap(pixmap)

        # parent.parent.setWindowIcon(icon)
        # toggleButton.setIconSize(pixmap.size())

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
            titleLayout.addWidget(icon_lab, stretch=0)
        titleLayout.addWidget(self.toggleButton, stretch=1)
        # titleLayout.addStretch()

        mainLayout.addWidget(titleFrame)
        mainLayout.addWidget(self.contentArea)
        # mainLayout.addStretch()
        self.setLayout(self.mainLayout)

        def start_animation(checked):
            arrow_type = QtCore.Qt.DownArrow if checked else QtCore.Qt.RightArrow
            direction = QtCore.QAbstractAnimation.Forward if checked else QtCore.QAbstractAnimation.Backward
            # toggleButton.setArrowType(arrow_type)
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
