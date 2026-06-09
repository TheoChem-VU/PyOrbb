from PySide6 import QtWidgets, QtCore
from . import rich_widgets, latex_renderer
from functools import partial


class DrawAction(QtWidgets.QWidgetAction):
    def __init__(self, parent, text, icon, func):
        super().__init__(parent)
        self.text = text
        self.icon = icon
        self.func = func

    def createWidget(self, parent): 
        action_widget = QtWidgets.QWidget(parent)
        action_widget.setStyleSheet("margin: 5px; padding: 5px")
        action_layout = QtWidgets.QHBoxLayout(parent)
        action_layout.setContentsMargins(2, 2, 2, 2)
        action_widget.setLayout(action_layout)

        icon_pixmap = self.icon.pixmap(self.icon.actualSize(QtCore.QSize(20, 20)))
        icon_lb = QtWidgets.QLabel(parent)
        icon_lb.setPixmap(icon_pixmap)
        action_layout.addWidget(icon_lb)

        darkmode = parent.parent().parent.parent.isDarkMode
        action_btn = rich_widgets.HTMLPushButton("", parent)

        if darkmode:
            action_widget.setStyleSheet("QWidget{border: none; background-color: none;} QWidget:hover {background-color: #656565;}")
        else:
            action_widget.setStyleSheet("QWidget{border: none; background-color: none;} QWidget:hover {background-color: #f0f0f0;}")
        action_layout.addWidget(action_btn)
        action_btn.setPixmap(latex_renderer.convert_to_QPixMap(str(self.text), darkmode=darkmode))
        self.triggered.connect(self.func)
        action_btn.clicked.connect(self.func)
        return action_widget
