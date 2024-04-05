import openpyxl as xl
import pyfmo
try:
    from openpyxl.cell import get_column_letter
except ImportError:
    from openpyxl.utils import get_column_letter

from matplotlib import colormaps
from tcutility import ensure_list, formula
import numpy as np
from scm import plams


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


def to_excel(orbs, out_file: str = 'pyfmo.xlsx'):
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
                if sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO' and sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'), left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo2.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(left=xl.styles.Side(border_style="medium", color='808080'))
                elif sfo1.make_name(frag_name=False, spin=False, relative_name=True) == 'LUMO':
                    cell.border = xl.styles.Border(top=xl.styles.Side(border_style="medium", color='808080'))

                cell.number_format = number_format
                dim_holder.setdefault(get_column_letter(j+4), xl.worksheet.dimensions.ColumnDimension(sheet, min=j+4, max=j+4, bestFit=True))

        # fixing the column widths
        dim_holder['C'] = xl.worksheet.dimensions.ColumnDimension(sheet, index='C', auto_size=True)
        sheet.column_dimensions = dim_holder
        sheet.freeze_panes = sheet['D4']
        return sheet

    def make_key_value_table(rows, start_row, start_column, asterisks=[]):
        dim_holder = xl.worksheet.dimensions.DimensionHolder(worksheet=sheet)

        for i, (variable, value) in enumerate(rows.items()):
            cell_var = sheet.cell(row=start_row + i, column=start_column, value=variable)
            cell_val = sheet.cell(row=start_row + i, column=start_column + 1, value=value)
            if isinstance(value, float):
                    cell_val.number_format = '0.00'
            if i > 0:
                cell_pad = sheet.cell(row=start_row + i, column=start_column + 2, value=" ")

            if i == 0:
                cell_var.font = xl.styles.Font(b=True, i=True, size=16)
                cell_var.border = xl.styles.Border(bottom=xl.styles.Side(border_style="double"))

                cell_val.font = xl.styles.Font(b=True, i=True, size=16)
                cell_val.border = xl.styles.Border(bottom=xl.styles.Side(border_style="double"))

            elif i == (len(rows) - 1):
                cell_var.font = xl.styles.Font(b=True)
                cell_var.border = xl.styles.Border(left=xl.styles.Side(border_style="thin"), bottom=xl.styles.Side(border_style="thin"))
                cell_val.border = xl.styles.Border(right=xl.styles.Side(border_style="thin"), bottom=xl.styles.Side(border_style="thin"))

            else:
                cell_var.font = xl.styles.Font(b=True)
                cell_var.border = xl.styles.Border(left=xl.styles.Side(border_style="thin"))
                cell_val.border = xl.styles.Border(right=xl.styles.Side(border_style="thin"))

        dim_holder.setdefault(get_column_letter(start_column), xl.worksheet.dimensions.ColumnDimension(sheet, min=start_column, max=start_column, bestFit=True))
        dim_holder.setdefault(get_column_letter(start_column+1), xl.worksheet.dimensions.ColumnDimension(sheet, min=start_column+1, max=start_column+1, bestFit=True))

        for j, asterisk in enumerate(asterisks):
            cell = sheet.cell(row=start_row + i + j + 1, column=start_column, value=f'{"*"*(j+1)} {asterisk}')
        
        sheet.column_dimensions = dim_holder
        
        return len(rows) + len(asterisks) + start_row, start_column + 1

    def make_table_sheet(sheet_name, sheet_title, rows, header):
        sheet = wb.create_sheet(sheet_name)

        title_cell = sheet.cell(row=1, column=1, value=sheet_title)
        title_cell.font = xl.styles.Font(b=True, size=24)

        dim_holder = xl.worksheet.dimensions.DimensionHolder(worksheet=sheet)

        for j, col in enumerate(header):
            cell = sheet.cell(row=3, column=j+2, value=col)
            cell.font = xl.styles.Font(b=True)
            cell.border = xl.styles.Border(bottom=xl.styles.Side(border_style="double"))

        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                cell = sheet.cell(row=i + 4, column=j+2, value=val)
                if isinstance(val, float):
                    cell.number_format = '0.00'
                dim_holder.setdefault(get_column_letter(j+4), xl.worksheet.dimensions.ColumnDimension(sheet, min=j+4, max=j+4, bestFit=True))

        # fixing the column widths
        dim_holder['C'] = xl.worksheet.dimensions.ColumnDimension(sheet, index='C', auto_size=True)
        sheet.column_dimensions = dim_holder
        sheet.freeze_panes = sheet['D4']

    # open a new notebook
    wb = xl.Workbook()

    # we will write some basic info about the calcualtion in the first sheet
    sheet = wb.worksheets[0]
    sheet.title = 'Info'
    # write the title cell
    title_cell = sheet.cell(row=1, column=1, value='PyFMO Analysis')
    title_cell.font = xl.styles.Font(b=True, size=24)

    cell = sheet.cell(row=4, column=6, value='Here I will write the mixing situations later')

    # write information about the complex
    mols = get_molecules(orbs.reader)
    rows = {
        'Complex': '',
        'Formula': formula.molecule(mols['complex']),
        'Coords': '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols['complex']]),
        'No. MOs': len(orbs.mos.mos),
        'No. occ. MOs': len([mo for mo in orbs.mos if mo.occupied]),
        'No. virt. MOs': len([mo for mo in orbs.mos if not mo.occupied]),
        'ΔE_int': orbs.reader.read('Energy', 'Bond Energy') * 627.503,
        'ΔE_Pauli': orbs.reader.read('Energy', 'Pauli Total') * 627.503,
        'ΔE_oi': orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503,
        'ΔV_elstat': orbs.reader.read('Energy', 'elstat') * 627.503,
        'ΔE_disp': orbs.reader.read('Energy', 'Dispersion Energy') * 627.503,
    }
    next_row, _ = make_key_value_table(rows, 4, 2, asterisks=['EDA terms given in (kcal mol⁻¹)'])

    # write information about the fragments
    for i, fragment in enumerate(orbs.fragments):
        sfos = orbs.sfos.get_fragment_sfos(fragment)
        rows = {
            'Fragment': fragment,
            'Formula': formula.molecule(mols[fragment]),
            'Coords': '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols[fragment]]),
            'No. SFOs': len(sfos),
            'No. occ. SFOs': len([sfo for sfo in sfos if sfo.occupied]),
            'No. virt. SFOs': len([sfo for sfo in sfos if not sfo.occupied]),
        }
        next_row, _ = make_key_value_table(rows, next_row + 1, 2)

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
        ])

    headers = [
        'Index', 
        'Name', 
        'Relative Name', 
        'Occupation', 
        'Spin', 
        'Symmetry', 
        'Energy (eV)',
    ]
    make_table_sheet('MOs', 'Molecular Orbitals', rows, headers)

    for fragment in orbs.fragments:
        rows = []
        for sfo in orbs.sfos.get_fragment_sfos(fragment):
            rows.append([
                sfo.index,
                sfo.name,
                sfo.relname,
                int(sfo.occupation),
                sum([orbs.mulliken_contribution(mo, sfo) * mo.occupation for mo in orbs.mos if mo.occupied]),
                sfo.spin,
                sfo.symmetry,
                sfo.energy,
            ])

        headers = [
            'Index', 
            'Name', 
            'Relative Name', 
            'Occupation',
            'Gross Pop.',
            'Spin', 
            'Symmetry', 
            'Energy (eV)',
        ]
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

        for fragment in orbs.fragments:
            cmap = colors.LinearSegmentedColormap.from_list('RdGn', ['#67000dff', '#ffffffff', '#157E3AFF'])
            name = f"Coefficients {fragment} {spin}" if spin != 'AB' else f"Coefficients {fragment}"
            title = f"MO Coefficients from {fragment} (spin {spin})" if spin != 'AB' else f"MO Coefficients from {fragment}"
            sfos_ = [sfo for sfo in sfos_spin if sfo.fragment == fragment]
            coeff = coefficient_mat(sfos_, mos_spin)
            make_sheet(name, title, mos_spin, sfos_, coeff.T, number_format='0.00', cmap=cmap, use_two_scale=True)

            name = f"Contributions {fragment} {spin}" if spin != 'AB' else f"Contributions {fragment}"
            title = f"Mulliken Contributions from {fragment} (spin {spin})" if spin != 'AB' else f"Mulliken Contributions from {fragment}"
            sfos_ = [sfo for sfo in sfos_spin if sfo.fragment == fragment]
            contribs = contribution_mat(orbs, sfos_, mos_spin)
            make_sheet(name, title, mos_spin, sfos_, contribs.T, number_format='0.00%', clip=(0, 1))

    wb.save(out_file)


if __name__ == '__main__':
    import yutility

    rkffile = "../../test/fixtures/NH3BH3/adf.rkf"
    # yutility.print_kf(rkffile, True)

    orbs = pyfmo.orbitals.Orbitals(rkffile)
    orbs.write_excel()

    # rkffile = r"D:\Users\Yuman\Desktop\PhD\PyOrb\test\fixtures\NH3BH3\adf.rkf"
    # orbs = pyfmo.Orbitals(rkffile)
    # to_excel(orbs.sfos['Acceptor'], orbs.sfos['Donor'])
