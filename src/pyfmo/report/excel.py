import openpyxl as xl
import pyfmo
try:
    from openpyxl.cell import get_column_letter
except ImportError:
    from openpyxl.utils import get_column_letter

from matplotlib import colormaps
from tcutility import ensure_list
import numpy as np


def to_excel(sfos1, sfos2, out_file: str = 'pyfmo.xlsx'):
    '''
    Write data about sfos1 and sfos2 to a nicely formatted excel file.
    Currently writes overlap matrices, square of overlap matrices, energy gap matrices and orbital interaction matrices.
    Both restricted and unrestricted sfos are supported.
    '''

    def make_sheet(sheet_name, sheet_title, sfos1, sfos2, values, number_format='0.00', cmap='Greens'):
        '''
        Create a new sheet and write data to it.

        Args:
            sheet_name: name to give to the sheet.
            sheet_title: title to write to cell A1 in the new sheet.
            sfos1, sfos2: lists of SFOs to be considered.
            values: matrix containing info about sfos1 and sfos2.
            number_format: format code for the values, for example, 0.00% for percentages rounded to 2 decimals.
            cmap: colormap name to color the cells by their values.
        '''
        cmap = colormaps[cmap]  # fetch the colormap from matplotlib

        sheet = wb.create_sheet(sheet_name)

        title_cell = sheet.cell(row=1, column=1, value=sheet_title)
        title_cell.font = xl.styles.Font(b=True, size=24)

        # cells that contain the names of the fragments (C2 and B3)
        f1_cell = sheet.cell(row=4, column=2, value=sfos1[0].fragment_unique_name)
        f1_cell.font = xl.styles.Font(b=True, size=16)
        f1_cell.alignment = xl.styles.Alignment(textRotation=90, horizontal="center", vertical="center")

        f2_cell = sheet.cell(row=2, column=4, value=sfos2[0].fragment_unique_name)
        f2_cell.font = xl.styles.Font(b=True, size=16)
        f2_cell.alignment = xl.styles.Alignment(horizontal="center", vertical="center")

        # the fragment names should span the whole row/column, so we merge the cells
        sheet.merge_cells(start_row=4, end_row=len(sfos1)+3, start_column=2, end_column=2)
        sheet.merge_cells(start_column=4, end_column=len(sfos2)+3, start_row=2, end_row=2)

        # we set the borders of the upper left corner to make it fit with the other borders
        sheet['C3'].border = xl.styles.Border(right=xl.styles.Side(border_style="thick"), bottom=xl.styles.Side(border_style="thick"))

        # dim_holder will be used to auto-format the columns
        dim_holder = xl.worksheet.dimensions.DimensionHolder(worksheet=sheet)
        
        for i, sfo1 in enumerate(sfos1):
            # this cell will hold the name of sfo1 in the column header
            name_cell = sheet.cell(row=i+4, column=3, value=sfo1.make_name(frag_name=False, spin=False, relative_name=False))
            name_cell.font = xl.styles.Font(b=True)
            name_cell.border = xl.styles.Border(right=xl.styles.Side(border_style="thick"))

            for j, sfo2 in enumerate(sfos2):
                # this cell will hold the name of sfo2 in the row header
                name_cell = sheet.cell(row=3, column=j+4, value=sfo2.make_name(frag_name=False, spin=False, relative_name=False))
                name_cell.font = xl.styles.Font(b=True)
                name_cell.border = xl.styles.Border(bottom=xl.styles.Side(border_style="thick"))
                name_cell.alignment = xl.styles.Alignment(horizontal="center", vertical="center")

                # normalize the value of the cell to [0, 1] to determine the color
                x = (values[i, j] - np.nanmin(values))/(np.nanmax(values) - np.nanmin(values)) * 0.8
                color = cmap(x)  # get the RGB tuple of floats from the cmap
                color = [int(x_*256) for x_ in color]  # convert the color from float to integer
                color = f'{color[0]:02x}{color[1]:02x}{color[2]:02x}'  # convert tuple of ints to a hex-code

                # set the value of the cell
                cell = sheet.cell(row=i+4, column=j+4, value=values[i, j])
                # if the value of the cell is None we should not color it (defaults to black for None-valued cells)
                if not np.isnan(values[i, j]):
                    cell.fill = xl.styles.PatternFill(start_color=color, end_color=color, fill_type="solid")

                # we draw a border if the cell is to the bottom-right of the HOMO-LUMO border
                if sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'HOMO' and sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'HOMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'), left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'HOMO':
                    cell.border = xl.styles.Border(left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'HOMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'))

                cell.number_format = number_format
                dim_holder.setdefault(get_column_letter(j+4), xl.worksheet.dimensions.ColumnDimension(sheet, min=j+4, max=j+4, bestFit=True))

        # fixing the column widths
        dim_holder['C'] = xl.worksheet.dimensions.ColumnDimension(sheet, index='C', auto_size=True)
        sheet.column_dimensions = dim_holder
        sheet.freeze_panes = sheet['D4']
        return sheet

    # open a new notebook
    wb = xl.Workbook()
    # we add a new sheet for each spin species
    spins = sorted(list(set(sfo.spin for sfo in sfos1 + sfos2)))  # becomes either ['A', 'B'] or ['AB']
    for spin in spins:
        # we sort the sfos into similar spin species and invert their order (virtual left and up, occupied right and down)
        sfos1_spin = [sfo for sfo in sfos1 if sfo.spin == spin][::-1]
        sfos2_spin = [sfo for sfo in sfos2 if sfo.spin == spin][::-1]

        # add the data we want
        name = f"Overlap {spin}" if spin != 'AB' else "Overlap"
        title = f"Overlaps (spin {spin})" if spin != 'AB' else "Overlaps"
        make_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin), number_format='0.0%')

        name = f"Overlap² {spin}" if spin != 'AB' else "Overlap²"
        title = f"Overlaps² (spin {spin})" if spin != 'AB' else "Overlaps²"
        make_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin)**2, number_format='0.0%')

        name = f"Δε {spin}" if spin != 'AB' else "Δε"
        title = f"Δε (spin {spin}) (eV)" if spin != 'AB' else "Δε (eV)"
        make_sheet(name, title, sfos1_spin, sfos2_spin, energy_gap_mat(sfos1_spin, sfos2_spin), number_format='0.00')

        name = f"Orbint {spin}" if spin != 'AB' else "Orbint"
        title = f"Orbital Interactions (spin {spin}) (1000/eV)" if spin != 'AB' else "Orbital Interactions (1000/eV)"
        oi = orbint_mat(sfos1_spin, sfos2_spin)
        oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
        make_sheet(name, title, sfos1_spin, sfos2_spin, oi, number_format='0.00')

    wb.save(out_file)



