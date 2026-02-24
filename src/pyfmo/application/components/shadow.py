from PySide6 import QtWidgets, QtGui

_shadows = {}

def update_style():
    use_darkmode = QtWidgets.QApplication.instance().isDarkMode
    newc = QtGui.QColor(0, 0, 0, 150) if use_darkmode else QtGui.QColor(170, 170, 170, 150)
    # newc = QtGui.QColor(170, 170, 170, 150)
    for widg, shadow in _shadows.items():
        shadow.setColor(newc)
        widg.setGraphicsEffect(shadow)

def apply(widg, radius=40):
    use_darkmode = QtWidgets.QApplication.instance().isDarkMode
    newc = QtGui.QColor(0, 0, 0, 150) if use_darkmode else QtGui.QColor(170, 170, 170, 150)
    # newc = QtGui.QColor(170, 170, 170, 150)
    shadow = QtWidgets.QGraphicsDropShadowEffect(
        parent=widg, 
        blurRadius=radius, 
        xOffset=0, 
        yOffset=0, 
        color=newc)
    _shadows[widg] = shadow
    widg.setGraphicsEffect(shadow)
