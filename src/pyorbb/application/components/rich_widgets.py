#!/usr/bin/env python
# vim: set expandtab tabstop=4 shiftwidth=4:

# PyQt Functionality Snippet by Apocalyptech
# "Licensed" in the Public Domain under CC0 1.0 Universal (CC0 1.0)
# Public Domain Dedication.  Use it however you like!
#
# https://creativecommons.org/publicdomain/zero/1.0/
# https://creativecommons.org/publicdomain/zero/1.0/legalcode

import sys
from PySide6 import QtWidgets, QtGui, QtCore

class HTMLStyle(QtWidgets.QProxyStyle):
    """
    A QProxyStyle which can be used to render HTML/Rich text inside
    QComboBoxes, QCheckBoxes and QRadioButtons.  Note that for QComboBox,
    this does NOT alter rendering of the items when you're choosing from
    the list.  For that you'll need to set an item delegate.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(QtWidgets.QStyleFactory.create('Fusion'))
        self.text_doc = QtGui.QTextDocument()

    def drawItemText(self, painter, rect, alignment, pal, enabled, text, text_role):
        """
        This is what draws the text - we use an internal QTextDocument
        to do the formatting.  The general form of this function follows the
        C++ version at https://github.com/qt/qtbase/blob/5.9/src/widgets/styles/qstyle.cpp

        Note that we completely ignore the `alignment` and `enabled` variable.
        This is always left-aligned, and does not currently support disabled
        widgets.
        """
        if not text or text == '':
            return

        # Save our current pen if we need to
        saved_pen = None
        if text_role != QtGui.QPalette.NoRole:
            saved_pen = painter.pen()
            painter.setPen(QtGui.QPen(pal.brush(text_role), saved_pen.widthF()))

        # Render the text.  There's a bit of voodoo here with the rectangles
        # and painter translation; there was various bits of finagling necessary
        # to get this to seem to work with both combo boxes and checkboxes.
        # There's probably better ways to be doing this.
        margin = 3
        painter.save()
        painter.translate(rect.left()-margin, 0)
        self.text_doc.setHtml(text)
        self.text_doc.drawContents(painter,
                QtCore.QRectF(rect.adjusted(-rect.left(), 0, -margin, 0)))
        painter.restore()

        # Restore our previous pen if we need to
        if text_role != QtGui.QPalette.NoRole:
            painter.setPen(saved_pen)

    def sizeFromContents(self, contents_type, option, size, widget=None):
        """
        For ComboBoxes, this gets called to determine the size of the list of
        options for the comboboxes.  This is too wide for our HTMLComboBox, so
        we pull in the width from there instead.
        """
        width = size.width()
        height = size.height()
        return super().sizeFromContents(contents_type,
                option,
                QtCore.QSize(width, height),
                widget)

class HTMLDelegate(QtWidgets.QStyledItemDelegate):
    """
    Class for use in a QComboBox to allow HTML Text.  I'm still a bit
    miffed that this isn't just a default part of Qt.  There's a lot of
    Google hits from people looking to do this, most suggesting
    implementing something like this, but so far I've only found this
    one actual implementation, at the end of a thread here:
    http://www.qtcentre.org/threads/62867-HTML-rich-text-delegate-and-text-centering-aligning-code-amp-pictures

    I suspect this implementation is probably heavier than we actually
    need, but it seems fairly voodooey anyway.  And keep in mind that
    after all this, you've still got to produce a completely bloody
    different solution for displaying the currently-selected item in
    the QComboBox; this is only for the list of choices.  I'm happy
    to be leaving Gtk but this kind of thing makes the move more
    bittersweet than it should be.
    """

    def __init__(self, parent=None):
        super().__init__()
        self.doc = QtGui.QTextDocument(self)

    def paint(self, painter, option, index):
        """
        Paint routine for our items in the QComboBox
        """

        # Save our painter so it can be restored later
        painter.save()

        # Copy our option var so we can make some changes without modifying
        # the underlying object
        options = QtWidgets.QStyleOptionViewItem(option)
        self.initStyleOption(options, index)

        # Add in our data to our QTextDocument
        self.doc.setHtml(options.text)

        # Acquire our style
        style = QtWidgets.QStyleFactory.create('Fusion')

        # Draw a barebones version of the control which doesn't have any
        # text specified - this is to render the background, basically, so
        # that when we're mousing over one of the items the bg changes.
        options.text = ''
        style.drawControl(QtWidgets.QStyle.CE_ItemViewItem, options, painter)

        # Grab a PaintContext and set our text color depending on if we're
        # selected or not
        ctx = QtGui.QAbstractTextDocumentLayout.PaintContext()
        if option.state & QtWidgets.QStyle.State_Selected:
            ctx.palette.setColor(QtGui.QPalette.Text, option.palette.color(
                QtGui.QPalette.Active, QtGui.QPalette.HighlightedText))
        else:
            ctx.palette.setColor(QtGui.QPalette.Text, option.palette.color(
                QtGui.QPalette.Active, QtGui.QPalette.Text))

        # Calculating some rendering geometry.
        textRect = style.subElementRect(
            QtWidgets.QStyle.SE_ItemViewItemText, options)
        textRect.adjust(3, 0, 0, 0)
        painter.translate(textRect.topLeft())
        painter.setClipRect(textRect.translated(-textRect.topLeft()))

        # Now, finally, actually render the text
        self.doc.documentLayout().draw(painter, ctx)

        # Restore our paintbrush
        painter.restore()

    def sizeHint(self, option, index):
        """
        Our size.  This actually gets called before `paint`, I think, and therefore
        is called before our text has actually been loaded into the QTextDocument,
        but apparently seems to Do The Right Thing Anyway.
        """
        return QtCore.QSize(self.doc.idealWidth(), self.doc.size().height())

# class HTMLComboBox(QtWidgets.QComboBox):
#     """
#     Custom QComboBox class to handle dealing with HTML/Rich text.  This
#     is basically just here to set a few attributes and then implement
#     a custom sizeHint.  This implementation does NOT support ComboBoxes
#     whose contents get updated - sizeHint will cache its information
#     the first time it's called and then never update it.
#     """

#     def __init__(self, parent):
#         super().__init__(parent)
#         self.setStyle(HTMLStyle())
#         self.setItemDelegate(HTMLDelegate())
#         self.stored_size = None

#         # TODO: Figure out how to actually calculate these properly
#         self.width_adjust_sizehint = 20
#         self.width_adjust_contents = -30

#     def sizeHint(self):
#         """
#         Use a QTextDocument to compute our rendered text size
#         """
#         if not self.stored_size:
#             doc = QtGui.QTextDocument()
#             model = self.model()
#             max_w = 0
#             max_h = 0
#             for rownum in range(model.rowCount()):
#                 item = model.item(rownum)
#                 doc.setHtml(item.text())
#                 size = doc.size()
#                 if size.width() > max_w:
#                     max_w = size.width()
#                 if size.height() > max_h:
#                     max_h = size.height()

#             # Need to add in a bit of padding to account for the
#             # arrow selector
#             max_w += self.width_adjust_sizehint

#             self.stored_size = QtCore.QSize(max_w, max_h)
#         return self.stored_size

#     def minimumSizeHint(self):
#         """
#         Just use the same logic as `sizeHint`
#         """
#         return self.sizeHint()

class HTMLWidgetHelper(object):
    """
    Class to enable HTML/Rich text on a "simple" Qt widget such as QCheckBox
    or QRadioButton.  The most important bit is setting the widget style to
    HTMLStyle.  The rest is all just making sure that the widget is sized
    properly; without it, the widget will be too wide.  If you don't care
    about that, you can easily just use .setStyle(HTMLStyle()) on a regular
    widget without bothering with subclassing.

    There's doubtless some corner cases we're missing here, but it works
    for my purposes.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setStyle(HTMLStyle())
        self.stored_size = None

    def sizeHint(self):
        """
        Use a QTextDocument to compute our rendered text size
        """
        if not self.stored_size:
            doc = QtGui.QTextDocument()
            doc.setHtml(self.text())
            size = doc.size()
            # Details from this derived from QCheckBox/QRadioButton sizeHint sourcecode:
            # https://github.com/qt/qtbase/blob/5.9/src/widgets/widgets/qcheckbox.cpp
            # https://github.com/qt/qtbase/blob/5.9/src/widgets/widgets/qradiobutton.cpp
            opt = QtWidgets.QStyleOptionButton()
            self.initStyleOption(opt)
            self.stored_size = QtCore.QSize(
                    size.width() + opt.iconSize.width() + 4,
                    max(size.height(), opt.iconSize.height()))
        return self.stored_size

    def minimumSizeHint(self):
        """
        Just use the same logic as `sizeHint`
        """
        return self.sizeHint()

class HTMLTableWidgetItem(HTMLWidgetHelper, QtWidgets.QTableWidgetItem):
    """
    An HTML-enabled QCheckBox.  All the actual work is done in HTMLWidgetHelper.
    We're abusing (well, using) Python's multiple inheritance since the same code
    works well for more than one widget type.
    """

class HTMLCheckBox(HTMLWidgetHelper, QtWidgets.QCheckBox):
    """
    An HTML-enabled QCheckBox.  All the actual work is done in HTMLWidgetHelper.
    We're abusing (well, using) Python's multiple inheritance since the same code
    works well for more than one widget type.
    """

class HTMLRadioButton(HTMLWidgetHelper, QtWidgets.QRadioButton):
    """
    An HTML-enabled QRadioButton.  All the actual work is done in HTMLWidgetHelper.
    We're abusing (well, using) Python's multiple inheritance since the same code
    works well for more than one widget type.
    """

class HTMLToolButton(QtWidgets.QToolButton):
    def __init__(self, parent=None, text=None):
        if parent is not None:
            super().__init__(parent)
        else:
            super().__init__()
        self.__lbl = QtWidgets.QLabel(self)
        if text is not None:
            self.__lbl.setText(text)
        self.__lyt = QtWidgets.QHBoxLayout()
        self.__lyt.setContentsMargins(0, 0, 0, 0)
        # self.__lyt.setSpacing(0)
        self.setLayout(self.__lyt)
        self.__lbl.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        self.__lbl.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Minimum,
        )
        self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.__lbl.setTextFormat(QtCore.Qt.RichText)
        self.__lyt.addWidget(self.__lbl)
        self.setStyleSheet("""
            QToolButton {
                border: none;
            }
            QToolButton::menu-button {
                subcontrol-position: bottom;
                subcontrol-origin: margin;
                background: red;
            };""")
        return

    def setText(self, text):
        self.__lbl.setText(text)
        self.updateGeometry()
        return

    def setPixmap(self, pixmap):
        self.__lbl.setPixmap(pixmap)
        self.updateGeometry()
        return

    def sizeHint(self):
        s = QtWidgets.QToolButton.sizeHint(self)
        w = self.__lbl.sizeHint()
        s.setWidth(w.width())
        s.setHeight(w.height())
        return s



class HTMLPushButton(QtWidgets.QPushButton):
    def __init__(self, text=None, parent=None):
        if parent is not None:
            super().__init__(parent)
        else:
            super().__init__()
        self.__lbl = QtWidgets.QLabel()
        if text is not None:
            self.__lbl.setText(text)
        self.__lyt = QtWidgets.QHBoxLayout()
        self.__lyt.setContentsMargins(0, 0, 0, 0)
        # self.__lyt.setSpacing(0)
        self.setLayout(self.__lyt)
        self.__lbl.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        self.__lbl.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Minimum,
        )
        self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.__lbl.setTextFormat(QtCore.Qt.RichText)
        self.__lyt.addWidget(self.__lbl)
        return

    def setText(self, text):
        self.__lbl.setText(text)
        self.updateGeometry()
        return

    def setPixmap(self, pixmap):
        self.__lbl.setPixmap(pixmap)
        self.updateGeometry()
        return

    def sizeHint(self):
        s = QtWidgets.QPushButton.sizeHint(self)
        w = self.__lbl.sizeHint()
        s.setWidth(w.width())
        s.setHeight(w.height())
        return s

# class HTMLToolButton(QtWidgets.QToolButton):
#     def __init__(self, parent=None):
#         super().__init__(parent)
#         self.setStyle(HTMLStyle())
#         self.__lbl = QtWidgets.QLabel()

#     # def setText(self, text):
#     #     self.__lbl.setText(text)
#     #     super().setText(text)

#     # def sizeHint(self):
#     #     print( self.__lbl.sizeHint(), super().sizeHint())
#     #     s = self.__lbl.sizeHint()
#     #     s.setHeight(19)
#     #     s.setWidth(s.width()*3)
#     #     return s

#     def setPixmap(self, pixmap):
#         self.__lbl.setPixmap(pixmap)
#         self.updateGeometry()
#         return