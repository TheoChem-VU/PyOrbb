import openpyxl as xl
import pyfmo
try:
    from openpyxl.cell import get_column_letter
except ImportError:
    from openpyxl.utils import get_column_letter

from matplotlib import colormaps, colors
from tcutility import ensure_list, formula
import numpy as np
from scm import plams


CHAR_WIDTHS = {
    " ": 3,
    "!": 5,
    '"': 6,
    "#": 7,
    "$": 7,
    "%": 11,
    "&": 10,
    "'": 3,
    "(": 5,
    ")": 5,
    "*": 7,
    "+": 7,
    ",": 8,
    "-": 5,
    ".": 8,
    "/": 6,
    "0": 8,
    "1": 8,
    "2": 8,
    "3": 8,
    "4": 8,
    "5": 8,
    "6": 8,
    "7": 8,
    "8": 8,
    "9": 8,
    ":": 4,
    ";": 4,
    "<": 7,
    "=": 7,
    ">": 7,
    "?": 7,
    "@": 13,
    "A": 9,
    "B": 8,
    "C": 8,
    "D": 9,
    "E": 7,
    "F": 7,
    "G": 9,
    "H": 9,
    "I": 4,
    "J": 5,
    "K": 8,
    "L": 6,
    "M": 12,
    "N": 10,
    "O": 10,
    "P": 8,
    "Q": 10,
    "R": 8,
    "S": 7,
    "T": 7,
    "U": 9,
    "V": 9,
    "W": 13,
    "X": 8,
    "Y": 7,
    "Z": 7,
    "[": 5,
    "\\": 6,
    "]": 5,
    "^": 7,
    "_": 7,
    "`": 4,
    "a": 7,
    "b": 8,
    "c": 6,
    "d": 8,
    "e": 8,
    "f": 5,
    "g": 7,
    "h": 8,
    "i": 4,
    "j": 4,
    "k": 7,
    "l": 4,
    "m": 12,
    "n": 8,
    "o": 8,
    "p": 8,
    "q": 8,
    "r": 5,
    "s": 6,
    "t": 5,
    "u": 8,
    "v": 7,
    "w": 11,
    "x": 7,
    "y": 7,
    "z": 6,
    "{": 5,
    "|": 7,
    "}": 5,
    "~": 7,
}

def text_width(text):
    return (sum(CHAR_WIDTHS.get(char, 8) for char in str(text)) + 8.57) * 0.1318


def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


def energy_gap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1.energy - sfo2.energy))
    return np.array(ret).squeeze()


def orbint_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            if sfo1.occupation == sfo2.occupation:
                ret[-1].append(np.NaN)
            else:
                ret[-1].append((sfo1 @ sfo2)**2/abs(sfo1.energy - sfo2.energy))
    return np.array(ret).squeeze()


def coefficient_mat(sfos, mos):
    ret = []
    for sfo in ensure_list(sfos):
        ret.append([])
        for mo in ensure_list(mos):
            ret[-1].append(mo.get_coeff(sfo))
    return np.array(ret).squeeze()


def contribution_mat(orbs, sfos, mos):
    ret = []
    for sfo in ensure_list(sfos):
        ret.append([])
        for mo in ensure_list(mos):
            ret[-1].append(orbs.mulliken_contribution(sfo, mo))
    return np.array(ret).squeeze()


def get_molecules(reader):
    fragments_names = np.array(reader.read('Geometry', 'fragmenttype').split())
    coords = np.array(reader.read('Geometry', 'xyz')).reshape(-1, 3) * 0.529177249
    atoms = np.array(reader.read('Geometry', 'atomtype').split())

    order_index = np.array(ensure_list(reader.read('Geometry', 'atom order index'))[:coords.shape[0]]) - 1
    fragment_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[:coords.shape[0]]) - 1
    symbol_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[coords.shape[0]:]) - 1

    coords = coords[order_index]
    atoms = atoms[symbol_index][order_index]
    fragment = fragments_names[fragment_index][order_index]

    ret = {'complex': plams.Molecule()}
    [ret['complex'].add_atom(plams.Atom(symbol=atom, coords=coord)) for atom, coord in zip(atoms, coords)]
    for name in fragments_names:
        ret[name] = plams.Molecule()

        for atom, coord, frag in zip(atoms, coords, fragment):
            if frag != name:
                continue
            ret[name].add_atom(plams.Atom(symbol=atom, coords=coord))

    return ret


# def mo_activity_rank(orbs):
#     C = contribution_mat(orbs, tuple(orbs.sfos.sfos), tuple(orbs.mos.mos))
#     P = C * np.array([[sfo.occupation for sfo in orbs.sfos.sfos] for mo in orbs.mos.mos])
#     Pmo = np.sum(abs(P), axis=1)
#     Pmo_no_occ = abs(Pmo - np.array([mo.occupation for mo in orbs.mos.mos]))
#     mo_indices = np.argsort(Pmo_no_occ)

#     return [orbs.mos.mos[i] for i in mo_indices][::-1]


# def sfo_activity_order(orbs):
#     C = contribution_mat(orbs, orbs.sfos.sfos, orbs.mos.mos)
#     Csfo = np.sum(abs(C), axis=0)
#     sfo_indices = np.argsort(Csfo)

#     return [orbs.sfos.sfos[i] for i in sfo_indices][::-1]


# def sfo_activity_order_in_frag(orbs):
#     C = contribution_mat(orbs, orbs.sfos.sfos, orbs.mos.mos)
#     Csfo = np.sum(abs(C), axis=0)
#     sfo_indices = np.argsort(Csfo)
#     order = [orbs.sfos.sfos[i] for i in sfo_indices][::-1]
#     return {frag: [sfo for sfo in order if sfo.fragment == frag] for frag in orbs.fragments}


def to_excel(orbs, out_file: str = 'pyfmo.xlsx'):
    '''
    Write data about sfos1 and sfos2 to a nicely formatted excel file.
    Currently writes overlap matrices, square of overlap matrices, energy gap matrices and orbital interaction matrices.
    Both restricted and unrestricted sfos are supported.
    '''

    def make_matrix_sheet(sheet_name, sheet_title, sfos1, sfos2, values, number_format='0.00', cmap='Greens', use_two_scale=False, clip=None):
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
        if isinstance(cmap, str):
            cmap = colormaps[cmap]  # fetch the colormap from matplotlib

        sheet = wb.create_sheet(sheet_name)

        title_cell = sheet.cell(row=1, column=1, value=sheet_title)
        title_cell.font = xl.styles.Font(b=True, size=24)


        # cells that contain the names of the fragments (C2 and B3)
        if isinstance(sfos1[0], pyfmo.orbitals.sfo.SFO):
            frag = sfos1[0].fragment_unique_name
            f1_cell = sheet.cell(row=4, column=2, value=f'{frag} ({formula.molecule(mols[frag])})')
        else:
            f1_cell = sheet.cell(row=4, column=2, value='Complex MO')
        f1_cell.font = xl.styles.Font(b=True, size=16)
        f1_cell.alignment = xl.styles.Alignment(textRotation=90, horizontal="center", vertical="center")

        if isinstance(sfos2[0], pyfmo.orbitals.sfo.SFO):
            frag = sfos2[0].fragment_unique_name
            f2_cell = sheet.cell(row=2, column=4, value=f'{frag} ({formula.molecule(mols[frag])})')
        else:
            f2_cell = sheet.cell(row=2, column=4, value='Complex MO')
        f2_cell.font = xl.styles.Font(b=True, size=16)
        f2_cell.alignment = xl.styles.Alignment(horizontal="center", vertical="center")

        # the fragment names should span the whole row/column, so we merge the cells
        sheet.merge_cells(start_row=4, end_row=len(sfos1)+3, start_column=2, end_column=2)
        sheet.merge_cells(start_column=4, end_column=len(sfos2)+3, start_row=2, end_row=2)

        # we set the borders of the upper left corner to make it fit with the other borders
        sheet['C3'].border = xl.styles.Border(right=xl.styles.Side(border_style="thick"), bottom=xl.styles.Side(border_style="thick"))

        # normalize the data for coloring later
        clip = clip or (np.nanmin(values), np.nanmax(values))
        if use_two_scale:
            tsn = colors.TwoSlopeNorm(vcenter=0, vmin=clip[0], vmax=clip[1])
            normed_values = tsn(values)
        else:
            normed_values = (np.clip(values, *clip) - np.nanmin(np.clip(values, *clip)))/(np.nanmax(np.clip(values, *clip)) - np.nanmin(np.clip(values, *clip))) * 0.8

        column_widths = [0] * len(sfos2)
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

                column_widths[j] = max(column_widths[j], text_width(name_cell.value))
                # normalize the value of the cell to [0, 1] to determine the color

                color = cmap(normed_values[i, j])  # get the RGB tuple of floats from the cmap
                color = [int(x_*255) for x_ in color]  # convert the color from float to integer
                color = f'{color[0]:02x}{color[1]:02x}{color[2]:02x}'  # convert tuple of ints to a hex-code

                # set the value of the cell
                cell = sheet.cell(row=i+4, column=j+4, value=values[i, j])
                cell.number_format = number_format
                if not number_format.endswith('%'):
                    ndec = number_format.split('.')[1].count('0')
                    column_widths[j] = max(column_widths[j], text_width(str(round(values[i, j], ndec))))
                else:
                    ndec = number_format.split('.')[1].count('0')
                    column_widths[j] = max(column_widths[j], text_width(str(round(values[i, j]*100, ndec)) + '%'))

                # if the value of the cell is None we should not color it (defaults to black for None-valued cells)
                if not np.isnan(values[i, j]):
                    cell.fill = xl.styles.PatternFill(start_color=color, end_color=color, fill_type="solid")

                # we draw a border if the cell is to the bottom-right of the HOMO-LUMO border
                if sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO' and sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'), left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'))

        for j in range(len(sfos2)):
            sheet.column_dimensions[get_column_letter(j+4)] = xl.worksheet.dimensions.ColumnDimension(sheet, index=get_column_letter(j+4), width=column_widths[j])
        
        sheet.freeze_panes = sheet['D4']
        return sheet

    def make_key_value_table(rows, start_row, start_column, asterisks=[]):
        widths = [0] * len(rows[0])
        for i, row in enumerate(rows):
            for j, x in enumerate(row):
                try:
                    cell = sheet.cell(row=start_row + i, column=start_column + j, value=x)
                except ValueError:
                    cell = sheet.cell(row=start_row + i, column=start_column + j, value=str(x))

                if isinstance(x, float):
                    if round(x, 2) == 0.00:
                        cell.number_format = '0.00E+0'
                        widths[j] = max(widths[j], text_width('{:.2g}'.format(x)))
                    else:
                        cell.number_format = '0.00'
                        widths[j] = max(widths[j], text_width(round(x, 2)))
                elif isinstance(x, str) and '\n' in x:
                    widths[j] = max(widths[j], text_width(x.split('\n')[0]))
                else:
                    widths[j] = max(widths[j], text_width(x))

                if i == 0:
                    cell.font = xl.styles.Font(b=True)
                    cell.border = xl.styles.Border(bottom=xl.styles.Side(border_style="double"))

                elif i == (len(rows) - 1):
                    if j == 0:
                        cell.border = xl.styles.Border(left=xl.styles.Side(border_style="thin"), bottom=xl.styles.Side(border_style="thin"))
                        
                    elif j == (len(row) - 1):
                        cell.border = xl.styles.Border(right=xl.styles.Side(border_style="thin"), bottom=xl.styles.Side(border_style="thin"))

                    else:
                        cell.border = xl.styles.Border(bottom=xl.styles.Side(border_style="thin"))

                else:
                    if j == 0:
                        cell.border = xl.styles.Border(left=xl.styles.Side(border_style="thin"))

                    elif j == (len(row) - 1):
                        cell.border = xl.styles.Border(right=xl.styles.Side(border_style="thin"))

            if i > 0:
                sheet.cell(row=start_row + i, column=start_column + j + 1, value=" ")

        for j in range(len(rows[0])):
            sheet.column_dimensions[get_column_letter(start_column+j)] = xl.worksheet.dimensions.ColumnDimension(sheet, index=get_column_letter(start_column+j), width=widths[j])
        
        return len(rows) + len(asterisks) + start_row, len(rows[0]) + start_column + 1

    def make_table_sheet(sheet_name, sheet_title, rows, header):
        sheet = wb.create_sheet(sheet_name)

        title_cell = sheet.cell(row=1, column=1, value=sheet_title)
        title_cell.font = xl.styles.Font(b=True, size=24)

        widths = []
        for j, col in enumerate(header):
            cell = sheet.cell(row=3, column=j+2, value=col)
            cell.font = xl.styles.Font(b=True)
            cell.border = xl.styles.Border(bottom=xl.styles.Side(border_style="double"))
            widths.append(text_width(col))

        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                cell = sheet.cell(row=i + 4, column=j+2, value=val)
                if isinstance(val, float):
                    cell.number_format = '0.00'
                    val = round(val, 2)
                widths[j] = max(widths[j], text_width(val))

        # fixing the column widths
        dim_holder = xl.worksheet.dimensions.DimensionHolder(worksheet=sheet)
        for j in range(len(header)):
            dim_holder[get_column_letter(j+2)] = xl.worksheet.dimensions.ColumnDimension(sheet, index=get_column_letter(j+2), width=widths[j])
        sheet.column_dimensions = dim_holder
        sheet.freeze_panes = sheet['E4']

    # open a new notebook
    wb = xl.Workbook()

    # we will write some basic info about the calcualtion in the first sheet
    sheet = wb.worksheets[0]
    sheet.title = 'Info'
    # write the title cell
    title_cell = sheet.cell(row=1, column=1, value='PyFMO Analysis')
    title_cell.font = xl.styles.Font(b=True, size=24)


    # write information about the complex
    mols = get_molecules(orbs.reader)
    rows = [
        ['Complex', ''],
        ['Formula', formula.molecule(mols['complex'])],
        ['Coords', '[COPY THIS]\n     ' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols['complex']])],
        ['No. MOs', len(orbs.mos.mos)],
        ['No. occ. MOs', len([mo for mo in orbs.mos if mo.occupied])],
        ['No. virt. MOs', len([mo for mo in orbs.mos if not mo.occupied])],
        ['ΔE_int', orbs.reader.read('Energy', 'Bond Energy') * 627.503],
        ['ΔE_Pauli', orbs.reader.read('Energy', 'Pauli Total') * 627.503],
        ['ΔE_oi', orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503],
        ['ΔV_elstat', orbs.reader.read('Energy', 'elstat') * 627.503],
        ['ΔE_disp', orbs.reader.read('Energy', 'Dispersion Energy') * 627.503],
    ]
    next_row, _ = make_key_value_table(rows, 5, 2, asterisks=['EDA terms given in (kcal mol⁻¹)'])

    # write information about the fragments
    for i, fragment in enumerate(orbs.fragments):
        sfos = orbs.sfos.get_fragment_sfos(fragment)
        rows = [
            ['Fragment', fragment],
            ['Formula', formula.molecule(mols[fragment])],
            ['Coords', '[COPY THIS]\n     ' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols[fragment]])],
            ['No. SFOs', len(sfos)],
            ['No. occ. SFOs', len([sfo for sfo in sfos if sfo.occupied])],
            ['No. virt. SFOs', len([sfo for sfo in sfos if not sfo.occupied])],
        ]
        next_row, next_col = make_key_value_table(rows, next_row + 1, 2)

    # CHARGE-TRANSFER 2MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'MO1', '', 'MO2', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_two_mixing(orbs) if mix.nocc == 1]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '->', mo_names[0], '+', mo_names[1], mix.strength, mix.spin, mix.symmetry])

    next_row, next_col = make_key_value_table(rows, 5, 6)

    cell = sheet.cell(row=3, column=6, value='Two-Mixing (Charge-Transfer):')
    cell.font = xl.styles.Font(b=True, size=16)

    # PAULI 2MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'MO1', '', 'MO2', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_two_mixing(orbs) if mix.nocc == 2]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '->', mo_names[0], '+', mo_names[1], mix.strength, mix.spin, mix.symmetry])

    cell = sheet.cell(row=3, column=next_col+1, value='Two-Mixing (Pauli):')
    cell.font = xl.styles.Font(b=True, size=16)
    next_row, next_col = make_key_value_table(rows, 5, next_col+1)

    # LUMO-LUMO 2MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'MO1', '', 'MO2', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_two_mixing(orbs) if mix.nocc == 0]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '->', mo_names[0], '+', mo_names[1], mix.strength, mix.spin, mix.symmetry])

    cell = sheet.cell(row=3, column=next_col+1, value='Two-Mixing (LUMO-LUMO):')
    cell.font = xl.styles.Font(b=True, size=16)
    next_row, next_col = make_key_value_table(rows, 5, next_col+1)


    # CHARGE-TRANSFER 3MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'SFO3', '', 'MO1', '', 'MO2', '', 'MO3', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_three_mixing(orbs) if mix.nocc in [1, 2]]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '+', sfo_names[2], '->', mo_names[0], '+', mo_names[1], '+', mo_names[2], mix.strength, mix.spin, mix.symmetry])

    cell = sheet.cell(row=3, column=next_col+1, value='Three-Mixing (Charge-Transfer):')
    cell.font = xl.styles.Font(b=True, size=16)
    next_row, next_col = make_key_value_table(rows, 5, next_col+1)


    # PAULI 3MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'SFO3', '', 'MO1', '', 'MO2', '', 'MO3', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_three_mixing(orbs) if mix.nocc == 3]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '+', sfo_names[2], '->', mo_names[0], '+', mo_names[1], '+', mo_names[2], mix.strength, mix.spin, mix.symmetry])

    cell = sheet.cell(row=3, column=next_col+1, value='Three-Mixing (Pauli):')
    cell.font = xl.styles.Font(b=True, size=16)
    next_row, next_col = make_key_value_table(rows, 5, next_col+1)


    # LUMO-LUMO 3MIXING
    rows = [['Index', 'SFO1', '', 'SFO2', '', 'SFO3', '', 'MO1', '', 'MO2', '', 'MO3', 'Strength', 'Spin', 'Symmetry']]
    mixings = [mix for mix in pyfmo.analysis.closed_interactions.get_three_mixing(orbs) if mix.nocc == 0]
    for i, mix in enumerate(mixings[:25]):
        sfo_names = [sfo.make_name(frag_name=True, relative_name=False) for sfo in mix.sfos]
        mo_names = [mo.name for mo in mix.mos]
        rows.append([i, sfo_names[0], '+', sfo_names[1], '+', sfo_names[2], '->', mo_names[0], '+', mo_names[1], '+', mo_names[2], mix.strength, mix.spin, mix.symmetry])

    cell = sheet.cell(row=3, column=next_col+1, value='Three-Mixing (LUMO-LUMO):')
    cell.font = xl.styles.Font(b=True, size=16)
    next_row, next_col = make_key_value_table(rows, 5, next_col+1)

    cell = sheet.cell(row=3, column=2, value='System:')
    cell.font = xl.styles.Font(b=True, size=16)


    mo_order = pyfmo.analysis.orbital_activity.mos(orbs)
    # write a table with MO and SFO energies
    rows = []
    for mo in orbs.mos:
        rows.append([
            mo.index,
            mo.name,
            mo.relname,
            int(mo.occupation),
            mo.spin,
            mo.symmetry,
            mo.energy,
            mo_order.index(mo) + 1,
        ])

    headers = [
        'Index', 
        'Name', 
        'Relative Name', 
        'Occupation', 
        'Spin', 
        'Symmetry', 
        'Energy (eV)',
        'Activity Rank'
    ]
    make_table_sheet('MOs', 'Molecular Orbitals', rows, headers)

    sfo_order = pyfmo.analysis.orbital_activity.sfos(orbs)
    sfo_order_frag = {frag: [sfo for sfo in sfo_order if sfo.fragment == frag] for frag in orbs.fragments}
    include_site = False
    include_site_scf0 = False
    for fragment in orbs.fragments:
        rows = []
        for sfo in orbs.sfos:
            if sfo.fragment != fragment:
                continue

            rows.append([
                sfo.index,
                sfo.name,
                sfo.relname,
                int(sfo.occupation),
                sum([orbs.mulliken_contribution(mo, sfo) * mo.occupation for mo in orbs.mos if mo.occupied]),
                sfo.spin,
                sfo.symmetry,
                sfo.energy,
                sfo_order.index(sfo) + 1,
                sfo_order_frag[fragment].index(sfo) + 1,
            ])

            if hasattr(sfo, 'site_energy'):
                include_site = True

            if hasattr(sfo, 'site_energy_SCF0'):
                include_site_scf0 = True

            if include_site:
                rows[-1].append(sfo.site_energy)

            if include_site_scf0:
                rows[-1].append(sfo.site_energy_SCF0)

        headers = [
            'Index', 
            'Name', 
            'Relative Name', 
            'Occupation',
            'Gross Pop.',
            'Spin', 
            'Symmetry', 
            'Energy (eV)',
            'Activity Rank',
            'Activity Rank (in fragment)',
        ]

        if include_site:
            headers.append('Site Energy (eV)')

        if include_site_scf0:
            headers.append('Site Energy SCF0 (eV)')


        make_table_sheet(f'SFOs {fragment}', f'Fragment Orbitals for Fragment {fragment}', rows, headers)

    # we add a new sheet for each spin species
    for spin in orbs.spins:
        # we sort the sfos into similar spin species and invert their order (virtual left and up, occupied right and down)
        sfos_spin = [sfo for sfo in orbs.sfos if sfo.spin == spin]
        sfos1_spin = [sfo for sfo in sfos_spin if sfo.fragment == list(orbs.fragments)[0]]
        sfos2_spin = [sfo for sfo in sfos_spin if sfo.fragment == list(orbs.fragments)[1]]
        mos_spin = [mo for mo in orbs.mos if mo.spin == spin]

        # add the data we want
        name = f"Overlap {spin}" if spin != 'AB' else "Overlap"
        title = f"Overlaps (spin {spin})" if spin != 'AB' else "Overlaps"
        make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin), number_format='0.0%')

        name = f"Overlap² {spin}" if spin != 'AB' else "Overlap²"
        title = f"Overlaps² (spin {spin})" if spin != 'AB' else "Overlaps²"
        make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin)**2, number_format='0.0%')

        name = f"Δε {spin}" if spin != 'AB' else "Δε"
        title = f"Δε (spin {spin}) (eV)" if spin != 'AB' else "Δε (eV)"
        make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, energy_gap_mat(sfos1_spin, sfos2_spin), number_format='0.00')

        name = f"Orbint {spin}" if spin != 'AB' else "Orbint"
        title = f"Orbital Interactions (spin {spin}) (1000/eV)" if spin != 'AB' else "Orbital Interactions (1000/eV)"
        oi = orbint_mat(sfos1_spin, sfos2_spin)
        oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
        make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, oi, number_format='0.00')

        for fragment in orbs.fragments:
            cmap = colors.LinearSegmentedColormap.from_list('RdGn', ['#ad7a7fff', '#ffffffff', '#157e3bff'])
            name = f"Coefficients {fragment} {spin}" if spin != 'AB' else f"Coefficients {fragment}"
            title = f"MO Coefficients from {fragment} (spin {spin})" if spin != 'AB' else f"MO Coefficients from {fragment}"
            sfos_ = [sfo for sfo in sfos_spin if sfo.fragment == fragment]
            coeff = coefficient_mat(sfos_, mos_spin)
            make_matrix_sheet(name, title, mos_spin, sfos_, coeff.T, number_format='0.00', cmap=cmap, use_two_scale=True, clip=(-1, 1))

            name = f"Contributions {fragment} {spin}" if spin != 'AB' else f"Contributions {fragment}"
            title = f"Mulliken Contributions from {fragment} (spin {spin})" if spin != 'AB' else f"Mulliken Contributions from {fragment}"
            sfos_ = [sfo for sfo in sfos_spin if sfo.fragment == fragment]
            contribs = contribution_mat(orbs, sfos_, mos_spin)
            make_matrix_sheet(name, title, mos_spin, sfos_, contribs.T, number_format='0.00%', clip=(0, 1))

    wb.save(out_file)


if __name__ == '__main__':
    pyfmo.orbitals.Orbitals("../../test/fixtures/NH3BH3/adf.rkf").write_excel('NH3BH3.xlsx')
    # pyfmo.orbitals.Orbitals("/Users/yumanhordijk/Downloads/FeCO4CH4.adf.rkf").write_excel('FeCO4CH4.xlsx')
    # pyfmo.orbitals.Orbitals("../../test/fixtures/RadicalAddition/adf.rkf").write_excel('RadicalAddition.xlsx')
    # pyfmo.orbitals.Orbitals("../../test/fixtures/homo/FragAnal.adf.rkf").write_excel('homo.xlsx')
    # pyfmo.orbitals.Orbitals("../../test/fixtures/hetero/FragAnal.adf.rkf").write_excel('hetero.xlsx')
    # pyfmo.orbitals.Orbitals("../../test/fixtures/pentafluorophsophate/FragAnal.adf.rkf").write_excel('pentafluorophsophate.xlsx')

    # rkffile = r"D:\Users\Yuman\Desktop\PhD\PyOrb\test\fixtures\NH3BH3\adf.rkf"
    # orbs = pyfmo.Orbitals(rkffile)
    # to_excel(orbs.sfos['Acceptor'], orbs.sfos['Donor'])
