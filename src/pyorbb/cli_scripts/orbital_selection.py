import wx
import wx.grid as gridlib
import numpy as np
import ObjectListView as olv

#---------------------------------------------------------------------------

class OrbitalSelection(wx.App):
    def __init__(self, orbs):
        self.orbs = orbs
        super().__init__(None)
        self.frame = wx.Frame(None, size=(600, 500))

        self.setup_notebook()

    def start(self):
        self.frame.Show()
        # self.frame.Centre()
        self.MainLoop()

    def setup_notebook(self):
        self.notebook = wx.Notebook(self.frame, size=self.frame.GetSize())

        pages = {'complex': Page(self.notebook, self.orbs.mos.orbitals)}
        for fragment in self.orbs.fragments:
            pages[fragment] = Page(self.notebook, self.orbs.sfos.filter(fragment=fragment))

        for i, (page_name, page) in enumerate(pages.items()):
            self.notebook.InsertPage(i, page, page_name)


class FixedOLV(olv.ObjectListView):
    def _HandleSize(self, evt):
        """
        The ListView is being resized
        """
        self._PossibleFinishCellEdit()
        evt.Skip()
        self._ResizeSpaceFillingColumns()
        # Make sure our empty msg is reasonably positioned
        sz = self.GetClientSize()
        self.stEmptyListMsg.SetSize(0, sz.GetHeight()//3, sz.GetWidth(), sz.GetHeight())

    # def _SortItemsNow(self):
    #     """
    #     Sort the actual items in the list now, according to the current column and order
    #     """
    #     sortColumn = self.GetSortColumn()
    #     if not sortColumn:
    #         return

    #     secondarySortColumn = None # self.GetSecondarySortColumn()

    #     def _singleObjectComparer(col, object1, object2):
    #         value1 = col.GetValue(object1)
    #         value2 = col.GetValue(object2)

    #         try:
    #             return locale.strcoll(value1.lower(), value2.lower())
    #         except:
    #             return bool(value1 > value2)

    #     def _objectComparer(object1, object2):
    #         result = _singleObjectComparer(sortColumn, object1, object2)
    #         if secondarySortColumn and result == 0:
    #             result = _singleObjectComparer(secondarySortColumn, object1, object2)
    #         return result

    #     self.SortListItemsBy(_objectComparer)


    # def SortListItemsBy(self, cmpFunc, ascending=None):
    #     """
    #     Sort the existing list items using the given comparison function.

    #     The comparison function must accept two model objects as parameters.

    #     The primary users of this method are handlers of the SORT event that want
    #     to sort the items by their own special function.
    #     """
    #     if ascending is None:
    #         ascending = self.sortAscending

    #     def _sorter(key1, key2):
    #         cmpVal = cmpFunc(self.innerList[key1], self.innerList[key2])
    #         if ascending:
    #             return cmpVal
    #         else:
    #             return -cmpVal

    #     self.SortItems(_sorter)


    # def _HandleColumnClick(self, evt):
    #     """
    #     The user has clicked on a column title
    #     """
    #     evt.Skip()
    #     self._PossibleFinishCellEdit()

    #     # Toggle the sort column on the second click
    #     if evt.GetColumn() == self.sortColumnIndex:
    #         self.sortAscending = not self.sortAscending
    #     else:
    #         self.sortAscending = True

    #     self.SortBy(evt.GetColumn(), self.sortAscending)
    #     self._FormatAllRows()


    # def SortBy(self, newColumnIndex, ascending=True):
    #     """
    #     Sort the items by the given column
    #     """
    #     oldSortColumnIndex = self.sortColumnIndex
    #     self.sortColumnIndex = newColumnIndex
    #     self.sortAscending = ascending
    #     col = self.GetSortColumn()
    #     print(self.GetSortColumn())
    #     print(self.sortColumnIndex)
    #     values = [col.GetValue(obj) for obj in self.innerList]

    #     print(values)
    #     # value1 = col.GetValue(object1)
    #     # value2 = col.GetValue(object2)


    #     # # Let the world have a chance to sort the items
    #     # evt = OLVEvent.SortEvent(self, self.sortColumnIndex, self.sortAscending, self.IsVirtual())
    #     # self.GetEventHandler().ProcessEvent(evt)
    #     # if evt.IsVetoed():
    #     #     return

    #     # if not evt.wasHandled:
    #     #     self._SortItemsNow()
    #     order = np.argsort(values)
    #     if ascending:
    #         order = -order
    #     help(self.SortItems)

    #     self._UpdateColumnSortIndicators(self.sortColumnIndex, oldSortColumnIndex)




class Page(wx.Window):
    def __init__(self, parent, orbs):
        self.orbs = orbs
        super().__init__(parent, size=parent.GetSize())
        self.olv = FixedOLV(self, style=wx.LC_REPORT, size=parent.GetSize())
        self.write()

    def write(self):
        cols = [
            olv.ColumnDefn("Name", "left", -1, "name"),
            olv.ColumnDefn("Relative Name", "left", -1, "relative_name"),
            olv.ColumnDefn("Spin", "left", 30, "spin"),
            olv.ColumnDefn("Irrep", "left", 50, "symmetry"),
            olv.ColumnDefn("Energy (eV)", "right", -1, "energy", stringConverter=lambda o: f'{o:.2f}'),
        ]
        self.olv.SetColumns(cols)
        check_col = self.olv.CreateCheckStateColumn()
        [self.olv.Check(obj) for obj in self.orbs]

        self.olv.SetObjects(self.orbs)


class ShowGrid(gridlib.Grid):
    def __init__(self, parent, numRows, numCols):
        self.buttons = []
        gridlib.Grid.__init__(self, parent, size=parent.GetSize())
        self.DisableDragColSize()
        self.DisableDragRowSize()
        self.EnableEditing(False)
        self.CreateGrid(numRows, numCols)
        self.SetRowLabelSize(50)

        self.displayContent()
        self.handleEvents()


    def add_button(self, row, col, state=True):
        self.buttons.append((row, col, state, ButtonRenderer(state=state)))

    def displayContent(self):
        # Button coordinates
        for row, col, state, rd in self.buttons:
            attr = wx.grid.GridCellAttr()
            attr.SetRenderer(rd)
            self.SetAttr(row, col, attr)

    def handleEvents(self):
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.onCellSelected)

    def onCellSelected(self, event):
        for row, col, state, rd in self.buttons:
            if row==event.GetRow() and col==event.GetCol():
                # # Button inverted switch var
                switch_inv = {False: wx.CONTROL_FLAT, True: wx.CONTROL_CHECKED}
                # Invert button status : CONTROL_NONE (inactive) / wx.CONTROL_PRESSED (active)

                rd.buttonState = switch_inv[not rd.boolstat['boolstat']]
                # Print button status
                # self.SetCellValue(2, 2, "{0}".format({8:"inactive", 4:"active"}[self.rd.buttonState]))
                # Refresh the grid
                self.Refresh(eraseBackground=True)


class ButtonRenderer(wx.grid.PyGridCellRenderer):
    def __init__(self, state=True):
        wx.grid.GridCellRenderer.__init__(self)

        # Button neutral state
        self.boolstat = {'boolstat': state}
        # Button state inactive when initializing
        self.buttonState = wx.CONTROL_CHECKED

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        ref_rect = grid.CellToRect(row, col)
        # Define Push Button dimensions
        rect.SetSize(wx.Size(ref_rect[2], ref_rect[3]))
        # For centering the button
        rect.SetPosition(wx.Point(ref_rect[0], ref_rect[1]))
        # Draw the Push Button
        wx.RendererNative.Get().DrawCheckBox(grid, dc, rect, self.buttonState)
        # Switch button state indicator
        self.boolstat['boolstat'] = {wx.CONTROL_CHECKED: True, wx.CONTROL_FLAT: False}[self.buttonState]

#---------------------------------------------------------------------------

# if __name__ == '__main__':
#     import pyfmo


#     orbs = pyfmo.Orbitals('/Users/yumanhordijk/PhD/Programs/TheoCheM/PyFMO/calculations/PyOrb_testing_2022/MichaelAddition/substrate_cat_complex/sp.results/adf.rkf')

#     app = OrbitalSelection(orbs)
#     app.start()
#     # frame = MyFrame(orbs)
#     # frame.Show()
#     # frame.Centre()
#     # app.MainLoop()
