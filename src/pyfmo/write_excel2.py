import xlsxwriter as xl
from tcutility import ensure_list, formula
from tcutility.report import character
from scm import plams
import numpy as np
import pyfmo  # noqa
import warnings

warnings.filterwarnings('ignore', category=UserWarning, module='xlsxwriter')



def overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


def fock_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(sfo1.fock(sfo2))
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
    sfo_idx = [sfo.index - 1 for sfo in sfos]
    mo_idx = [mo.index - 1 for mo in mos]
    return orbs.data.matrices.mulliken_contribution.total[:, sfo_idx][mo_idx, :]


def _detect_nan_rects(arr):
    rects = []
    arr = np.isnan(arr)
    for j, row in enumerate(arr):
        for i, x in enumerate(row):
            if not x:
                continue

            for rect in rects:
                is_below = (rect[0] <= j <= (rect[2] + 1)) and (rect[1] <= i <= rect[3])
                is_besides = rect[2] == j and (rect[3] + 1) == i
                if not(is_below or is_besides):
                    continue

                if is_below:
                    rect[2] = j
                    break

                elif is_besides:
                    rect[3] = i
                    break

            else:
                rects.append([j, i, j, i])

    return rects


def get_molecules(reader):
    used_regions = reader.read('Geometry', 'nr of fragments') != reader.read('Geometry', 'nr of atoms')
    fragment_indices = np.atleast_1d(reader.read('Geometry', 'fragment and atomtype index'))
    fragment_indices = fragment_indices[len(fragment_indices)//2:]
    fragment_types = np.atleast_1d(reader.read('Geometry', 'fragmenttype').split())
    if not used_regions:
        fragment_uniques = [f'{fragment_types[frag_idx-1]}:{idx+1}' for idx, frag_idx in enumerate(fragment_indices)]
        fragment_uniques = np.array(sorted(set(fragment_uniques), key=lambda fu: int(fu.split(':')[1])))
    else:
        fragment_uniques = np.array(fragment_types)

    coords = np.array(reader.read('Geometry', 'xyz')).reshape(-1, 3) * 0.529177249
    atoms = np.array(reader.read('Geometry', 'atomtype').split())

    order_index = np.array(ensure_list(reader.read('Geometry', 'atom order index'))[:coords.shape[0]]) - 1
    fragment_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[:coords.shape[0]]) - 1
    symbol_index = np.array(ensure_list(reader.read('Geometry', 'fragment and atomtype index'))[coords.shape[0]:]) - 1

    coords = coords[order_index]
    atoms = atoms[symbol_index][order_index]
    fragment = fragment_uniques[fragment_index][order_index]

    ret = {'complex': plams.Molecule()}
    [ret['complex'].add_atom(plams.Atom(symbol=atom, coords=coord)) for atom, coord in zip(atoms, coords)]
    for name in fragment_uniques:
        ret[name] = plams.Molecule()

        for atom, frag in zip(ret['complex'], fragment):
            if frag != name:
                continue
            ret[name].add_atom(plams.Atom(symbol=atom.symbol, coords=atom.coords))

    return ret



def to_excel(orbs, out_file: str = 'pyfmo.xlsx'):
    '''
    Write data about sfos1 and sfos2 to a nicely formatted excel file.
    Currently writes overlap matrices, square of overlap matrices, energy gap matrices and orbital interaction matrices.
    Both restricted and unrestricted sfos are supported.
    '''

    workbook = xl.Workbook(out_file)

    default_fmt = workbook.add_format({'num_format': '0.00'})
    title_fmt = workbook.add_format({'bold': True, 'font_size': 24})
    float_fmt = workbook.add_format({'num_format': '0.00'})
    pctg_fmt = workbook.add_format({'num_format': '0.0%'})
    bottom_border_fmt = workbook.add_format({'bottom': 6, 'bold': True, 'align': 'center'})
    bottom_right_border_fmt = workbook.add_format({'bottom': 6, 'right': 6, 'bold': True, 'align': 'center'})
    right_border_fmt = workbook.add_format({'right': 6,  'bold': True, 'align': 'right'})
    bold_centered_fmt = workbook.add_format({'bold': True, 'font_size': 16, 'align': 'center_across'})
    bold_centered_rotated_fmt = workbook.add_format({'bold': True, 'font_size': 16, 'align': 'center', 'valign': 'vcenter', 'rotation': 90})
    white_bg_all_border_fmt = workbook.add_format({'bg_color': 'white', 'border': 1})
    white_bg_top_left_border_fmt = workbook.add_format({'bg_color': 'white', 'top': 1, 'left': 1})
    table_key_fmt = workbook.add_format({'bold': True, 'left': 1})
    table_val_fmt = workbook.add_format({'bold': False, 'right': 1})
    table_val_float_fmt = workbook.add_format({'bold': False, 'right': 1, 'num_format': '0.00'})
    table_title_fmt = workbook.add_format({'bold': True, 'bottom': 6, 'font_size': 16, 'italic': True})
    table_title2_fmt = workbook.add_format({'bold': True, 'bottom': 6, 'font_size': 16})
    table_ast_fmt = workbook.add_format({'top': 1})
    table_header_fmt = workbook.add_format({'bold': True, 'bottom': 6})


    def make_matrix_sheet(sheet_name, sheet_title, orbsx, orbsy, values, number_format=default_fmt, conditional_format=None, tab_color=None):
        '''
        Create a new sheet and write data to it.
        '''
        worksheet = workbook.add_worksheet(sheet_name.replace(':', '')[:31])
        if tab_color is not None:
            worksheet.set_tab_color(tab_color)

        # write the title in the top left corner
        worksheet.write(0, 0, sheet_title, title_fmt)
        worksheet.write_blank(2, 2, '', bottom_right_border_fmt)
        column_widths = [53] * len(orbsy)

        if orbs.data.calc_info.used_regions:
            if all(isinstance(orbx, pyfmo.orbitals2.objects.SFO) for orbx in orbsx):
                frag_name = list(set(orbx.fragment for orbx in orbsx))[0]
                labelx = f'{frag_name} ({formula.molecule(mols[frag_name])})'
            else:
                labelx = 'MO'
            worksheet.merge_range(1, 3, 1, 3 + len(orbsx), labelx, bold_centered_fmt)

            if all(isinstance(orby, pyfmo.orbitals2.objects.SFO) for orby in orbsy):
                frag_name = list(set(orby.fragment for orby in orbsy))[0]
                labely = f'{frag_name} ({formula.molecule(mols[frag_name])})'
            else:
                labely = 'MO'
            worksheet.merge_range(3, 1, 2 + len(orbsy), 1, labely, bold_centered_rotated_fmt)
        else:
            worksheet.merge_range(1, 3, 1, 2 + len(orbsx), 'Orbital', bold_centered_fmt)
            worksheet.merge_range(3, 1, 2 + len(orbsy), 1, 'Orbital', bold_centered_rotated_fmt)

        # write the orbital names on the x and y axes
        for i, orbx in enumerate(orbsx):
            name = orbx.name
            if isinstance(orbx, pyfmo.orbitals2.objects.SFO) and not orbs.data.calc_info.used_regions:
                name = f'{orbx.fragment_unique}({orbx.name})'
            worksheet.write(2, 3+i, name, bottom_border_fmt)

        for i, orby in enumerate(orbsy):
            name = orby.name
            if isinstance(orby, pyfmo.orbitals2.objects.SFO) and not orbs.data.calc_info.used_regions:
                name = f'{orby.fragment_unique}({orby.name})'

            worksheet.write(3+i, 2, name, right_border_fmt)
            column_widths[i] = max(column_widths[i], character.text_width(orby.name, font_size=11))

        # write the values
        for i, row in enumerate(values.T):
            for j, x in enumerate(row):
                if np.isnan(x):
                    continue

                worksheet.write(3+i, 3+j, x, number_format)

                if number_format is float_fmt:
                    column_widths[i] = max(column_widths[i], character.text_width(f'{x: .1}', font_size=11))
                elif number_format is pctg_fmt:
                    column_widths[i] = max(column_widths[i], character.text_width(f'{x: .1%}%', font_size=11))

        # now we detect squares of NaN in the values matrix and group them together
        squares = _detect_nan_rects(values.T)
        for square in squares:
            if square[2] == len(orbsx) - 1 and square[3] == len(orbsy) - 1:
                worksheet.merge_range(square[0] + 3, square[1] + 3, square[2] + 3, square[3] + 3, '', white_bg_top_left_border_fmt)
            else:
                worksheet.merge_range(square[0] + 3, square[1] + 3, square[2] + 3, square[3] + 3, '', white_bg_all_border_fmt)

        # set the color scale so that the matrix is more readable
        if conditional_format is None:
            conditional_format = {
                'type' : '2_color_scale',
                'min_color': 'white',
                'max_color': '63be7b',
            }
        worksheet.conditional_format(3, 3, 3 + len(orbsy), 3 + len(orbsx), conditional_format)

        worksheet.freeze_panes('D4')
        worksheet.set_default_row(hide_unused_rows=True)
        # worksheet.set_default_column(hide_unused_columns=True)

        for i, width in enumerate(column_widths, start=3):
            worksheet.set_column_pixels(i, i, width)


    def make_key_value_table(rows, start_row, start_column, asterisks=[]):
        for i, (variable, value) in enumerate(rows.items()):
            if i == 0:
                sheet.write(start_row + i, start_column, variable, table_title_fmt)
                sheet.write(start_row + i, start_column + 1, value, table_title2_fmt)
                continue

            sheet.write(start_row + i, start_column, variable, table_key_fmt)

            if isinstance(value, float):
                sheet.write(start_row + i, start_column + 1, value, table_val_float_fmt)
            else:
                sheet.write(start_row + i, start_column + 1, value, table_val_fmt)

            if i > 0:
                sheet.write(start_row + i, start_column + 2, " ")

        sheet.set_column_pixels(start_column, start_column, 93)
        sheet.set_column_pixels(start_column + 1, start_column + 1, 93)

        for j, asterisk in enumerate(asterisks):
            sheet.merge_range(start_row + i + j + 1, start_column, start_row + i + j + 1, start_column + 1, f'{"*"*(j+1)} {asterisk}', table_ast_fmt)

        if len(asterisks) == 0:
            sheet.merge_range(start_row + i + 1, start_column, start_row + i + 1, start_column + 1, ' ', table_ast_fmt)
        
        return len(rows) + len(asterisks) + start_row, start_column + 1


    def make_table_sheet(sheet_name, sheet_title, rows, header, tab_color=None):
        sheet = workbook.add_worksheet(sheet_name.replace(':', '')[:31])
        if tab_color is not None:
            sheet.set_tab_color(tab_color)

        sheet.write(0, 0, sheet_title, title_fmt)

        column_widths = []
        for j, col in enumerate(header):
            sheet.write(2, j+1, col, table_header_fmt)
            column_widths.append(character.text_width(col))

        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                if isinstance(val, float):
                    sheet.write(i + 3, j + 1, val, float_fmt)
                    val = str(round(val, 2))
                else:
                    sheet.write(i + 3, j + 1, val)
                column_widths[j] = max(column_widths[j], character.text_width(val))

        for i, width in enumerate(column_widths, start=1):
            sheet.set_column_pixels(i, i, width)
        sheet.freeze_panes('A4')


    # we will write some basic info about the calcualtion in the first sheet
    sheet = workbook.add_worksheet('🛈 Info')
    sheet.set_tab_color('D6D1CD')
    # write the title cell
    sheet.write(0, 0, 'PyFMO Analysis', title_fmt)
    sheet.write(3, 5, 'Here I will write the mixing situations later')

    # write information about the complex
    mols = get_molecules(orbs.reader)
    rows = {
        'Complex': '',
        'Formula': formula.molecule(mols['complex']),
        'Coords': '[Copy This]                      \n' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols['complex']]),
        'No. MOs': len(orbs.mos),
        'No. occ. MOs': len([mo for mo in orbs.mos if mo.occupied]),
        'No. virt. MOs': len([mo for mo in orbs.mos if not mo.occupied]),
        'No. frozen cores': orbs.data.MOs.nfrozencores.total,
        'ΔE_int': orbs.reader.read('Energy', 'Bond Energy') * 627.503,
        'ΔE_Pauli': orbs.reader.read('Energy', 'Pauli Total') * 627.503,
        'ΔE_oi': orbs.reader.read('Energy', 'Orb.Int. Total') * 627.503,
        'ΔV_elstat': orbs.reader.read('Energy', 'elstat') * 627.503,
        'ΔE_disp': orbs.reader.read('Energy', 'Dispersion Energy') * 627.503,
        'Point group': orbs.reader.read('Symmetry', 'grouplabel').strip(),
    }

    next_row, _ = make_key_value_table(rows, 3, 1, asterisks=['EDA terms given in (kcal mol⁻¹)'])

    # write information about the fragments
    for i, fragment in enumerate(orbs.fragments):
        sfos = orbs.sfos.get_fragment_sfos(fragment)
        rows = {
            'Fragment': fragment,
            'Formula': formula.molecule(mols[fragment]),
            'Coords': '[Copy This]                      \n' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols[fragment]]),
            'No. SFOs': len(sfos),
            'No. occ. SFOs': len([sfo for sfo in sfos if sfo.occupied]),
            'No. virt. SFOs': len([sfo for sfo in sfos if not sfo.occupied]),
        }
        next_row, _ = make_key_value_table(rows, next_row + 1, 1)

    has_kinetic = False
    # write a table with MO and SFO energies
    rows = []
    for mo in orbs.mos:
        rows.append([
            mo.index,
            mo.name,
            mo.relative_name,
            int(mo.occupation),
            mo.spin,
            mo.symmetry,
            mo.energy,
        ])
        if mo.kinetic_energy:
            has_kinetic = True
            rows[-1].append(mo.kinetic_energy)

    headers = [
        'Index', 
        'Name', 
        'Relative Name', 
        'Occupation', 
        'Spin', 
        'Symmetry', 
        'Energy (eV)',
    ]

    if has_kinetic:
        headers.append('Kinetic Energy (eV)')

    make_table_sheet('MOs', 'Molecular Orbitals', rows, headers, tab_color='D6D1CD')

    for fragment in orbs.fragments:
        rows = []
        for sfo in orbs.sfos:
            if sfo.fragment_unique != fragment:
                continue

            rows.append([
                sfo.index,
                sfo.name,
                sfo.relative_name,
                int(sfo.occupation),
                sfo.gross_population,
                sfo.spin,
                sfo.gross_spin,
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
            'Excess Spin',
            'Symmetry', 
            'Energy (eV)',
        ]
        make_table_sheet(f'SFOs {fragment}', f'Fragment Orbitals for Fragment {fragment}', rows, headers, tab_color='D6D1CD')


    # we add a new sheet for each spin species
    for spin in orbs.sfos.spins:
        # we sort the sfos into similar spin species and invert their order (virtual left and up, occupied right and down)
        sfos_spin = [sfo for sfo in orbs.sfos if sfo.spin == spin]
        # sfos_spin = [f'{sfo.fragment_unique}({sfo})' for sfo in sfos_spin]
        sfos1_spin = [sfo for sfo in sfos_spin if sfo.fragment_unique == list(orbs.sfos.fragments)[0]]
        sfos2_spin = [sfo for sfo in sfos_spin if sfo.fragment_unique == list(orbs.sfos.fragments)[1]]
        mos_spin = [mo for mo in orbs.mos if mo.spin == spin or mo.spin == 'AB' or spin == 'AB']
        # add the data we want
        name = f"Overlap {spin}" if spin != 'AB' else "Overlap"
        title = f"Overlaps (spin {spin})" if spin != 'AB' else "Overlaps"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin, sfos_spin, overlap_mat(sfos_spin, sfos_spin), number_format=float_fmt, tab_color='FF6666')
        else:
            make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin), number_format=float_fmt, tab_color='FF6666')

        name = f"Overlap² {spin}" if spin != 'AB' else "Overlap²"
        title = f"Overlaps² (spin {spin})" if spin != 'AB' else "Overlaps²"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin, sfos_spin, overlap_mat(sfos_spin, sfos_spin)**2, number_format=float_fmt, tab_color='FF6666')
        else:
            make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, overlap_mat(sfos1_spin, sfos2_spin)**2, number_format=float_fmt, tab_color='FF6666')

        name = f"Δε {spin}" if spin != 'AB' else "Δε"
        title = f"Δε (spin {spin}) (eV)" if spin != 'AB' else "Δε (eV)"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin, sfos_spin, energy_gap_mat(sfos_spin, sfos_spin), number_format=float_fmt, tab_color='4D8B31')
        else:
            make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, energy_gap_mat(sfos1_spin, sfos2_spin), number_format=float_fmt, tab_color='4D8B31')

        name = f"Orbint {spin}" if spin != 'AB' else "Orbint"
        title = f"Orbital Interactions (spin {spin}) (1000/eV)" if spin != 'AB' else "Orbital Interactions (1000/eV)"
        if not orbs.data.calc_info.used_regions:
            oi = orbint_mat(sfos_spin, sfos_spin)
            oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
            make_matrix_sheet(name, title, sfos_spin, sfos_spin, oi, number_format=float_fmt, tab_color='4D8B31')

        else:
            oi = orbint_mat(sfos1_spin, sfos2_spin)
            oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
            make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, oi, number_format=float_fmt, tab_color='4D8B31')

        for fragment in orbs.fragments:
            fragment_name = fragment
            if ':' in fragment:
                fragment_name = f'{fragment.split(":")[0]}({fragment.split(":")[1]})'

            name = f"Coefficients {fragment_name} {spin}" if spin != 'AB' else f"Coefficients {fragment_name}"
            title = f"MO Coefficients from {fragment_name} (spin {spin})" if spin != 'AB' else f"MO Coefficients from {fragment_name}"
            sfos_ = [sfo for sfo in sfos_spin if sfo.fragment_unique == fragment]
            sfo_idx = [sfo.index - 1 for sfo in sfos_]
            mo_idx = [mo.index - 1 for mo in mos_spin]
            coeff = orbs.data.matrices.coefficients.total[:, sfo_idx][mo_idx, :]
            cnd_fmt = {
                'type': '3_color_scale',
                'min_color': 'f8696b',
                'mid_color': 'white',
                'max_color': '63be7b',
                'min_value': -1,
                'mid_value': 0,
                'max_value': 1,
                'min_type': 'num',
                'mid_type': 'num',
                'max_type': 'num',
            }
            make_matrix_sheet(name, title, sfos_, mos_spin, coeff.T, number_format=float_fmt, conditional_format=cnd_fmt, tab_color='564D80')

            name = f"Contributions {fragment_name} {spin}" if spin != 'AB' else f"Contributions {fragment_name}"
            title = f"Mulliken Contributions from {fragment_name} (spin {spin})" if spin != 'AB' else f"Mulliken Contributions from {fragment_name}"
            contribs = contribution_mat(orbs, sfos_, mos_spin)
            cnd_fmt = {
                'type': '2_color_scale',
                'min_color': 'white',
                'max_color': '63be7b',
                'min_value': 0,
                'max_value': 1,
                'min_type': 'num',
                'max_type': 'num',
            }
            make_matrix_sheet(name, title, sfos_, mos_spin, contribs.T, number_format=pctg_fmt, conditional_format=cnd_fmt, tab_color='058ED9')

    for spin in orbs.sfos.spins:
        if not has_kinetic:
            break

        name = f"Fock {spin}" if spin != 'AB' else "Fock"
        title = f"Fock (spin {spin})" if spin != 'AB' else "Fock"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin, sfos_spin, fock_mat(sfos_spin, sfos_spin), number_format=float_fmt, tab_color='FF6666')
        else:
            make_matrix_sheet(name, title, sfos1_spin, sfos2_spin, fock_mat(sfos1_spin, sfos2_spin), number_format=float_fmt, tab_color='FF6666')

    workbook.close()



if __name__ == '__main__':
    import pyfmo  # noqa
    from tcutility import timer
    from time import perf_counter
    import matplotlib.pyplot as plt

    norbs = []
    load_time_new = []
    load_time_old = []
    # for alkyl in ['C1', 'C2', 'C3', 'C4', 'C5']:
    for alkyl in ['C1']:
        load_time_new.append([])
        load_time_old.append([])
        for _ in range(1):
            with timer.timer('new_orbitals.load_rkf'):
                start = perf_counter()
                orbs = pyfmo.orbitals2.objects.Orbitals(f'../../calculations/PyOrb_testing_2022/Alkyl/{alkyl}/EDA.results/adf.rkf')
                orbs.write_excel2()
                load_time_new[-1].append(perf_counter() - start)
                norbs.append(len(orbs.mos))
            # with timer.timer('new_orbitals.write_excel'):
            #     orbs.write_excel2()

            # with timer.timer('old_orbitals.load_rkf'):
            #     start = perf_counter()
            #     orbs_old = pyfmo.orbitals.Orbitals(f'../../calculations/PyOrb_testing_2022/Alkyl/{alkyl}/EDA.results/adf.rkf')
            #     contribution_mat(orbs_old, orbs_old.sfos.sfos, orbs_old.mos.mos)
            #     orbs_old.write_excel()
            #     load_time_old[-1].append(perf_counter() - start)
            # # with timer.timer('old_orbitals.mulliken_analysis'):
            # #     contribution_mat(orbs_old, orbs_old.sfos.sfos, orbs_old.mos.mos)
            # # with timer.timer('old_orbitals.write_excel'):
            #     orbs_old.write_excel()

    print(norbs, load_time_new)
    plt.scatter(norbs, np.array(load_time_new))
    plt.scatter(norbs, np.array(load_time_old))

    plt.title('Comparison Old and New method (incl. writing excel file)')
    plt.plot(np.unique(norbs), [np.mean(times) for times in load_time_new], label='New')
    plt.plot(np.unique(norbs), [np.mean(times) for times in load_time_old], label='Old')
    plt.yscale('log')
    plt.xlabel('No. MOs')
    plt.ylabel('Loading time (s)')

    plt.legend()
    plt.show()
