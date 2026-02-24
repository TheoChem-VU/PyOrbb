from PySide6 import QtWidgets, QtGui



def apply(widg, radius=40):
    default_shadow = QtWidgets.QGraphicsDropShadowEffect(
        parent=widg, 
        blurRadius=radius, 
        xOffset=0, 
        yOffset=0, 
        color=QtGui.QColor(170, 170, 170, 200))
    widg.setGraphicsEffect(default_shadow)
