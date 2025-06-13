import xlsxwriter as xl
from tcutility import ensure_list, formula
from tcutility.report import character
from scm import plams
import numpy as np
import pyfmo  # noqa
import warnings

warnings.filterwarnings('ignore', category=UserWarning, module='xlsxwriter')



def _overlap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append((sfo1 @ sfo2))
    return np.atleast_2d(np.array(ret).squeeze())


def _fock_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(sfo1.fock(sfo2))
    return np.atleast_2d(np.array(ret).squeeze())


def _energy_gap_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1.energy - sfo2.energy))
    return np.atleast_2d(np.array(ret).squeeze())


def _orbint_mat(sfos1, sfos2):
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            if sfo1.occupation == sfo2.occupation:
                ret[-1].append(np.NaN)
            else:
                ret[-1].append((sfo1 @ sfo2)**2/abs(sfo1.energy - sfo2.energy))
    return np.atleast_2d(np.array(ret).squeeze())


def _contribution_mat(orbs, sfos, mos):
    idx_offset = 0 if mos[0].spin in ['AB', 'A'] else len(mos)
    sfo_idx = [sfo.index - 1 + idx_offset for sfo in sfos]
    mo_idx = [mo.index - 1 + idx_offset for mo in mos]
    C = orbs.data.matrices.mulliken_contribution.total[:, sfo_idx][mo_idx, :]
    return C


def _get_molecules(reader):
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


def to_excel(orbs: pyfmo.Orbitals, out_file: str = 'pyfmo.xlsx', sfo_energy_type: str = 'energy'):
    '''
    Write data about orbitals and general information about the calculation to a nicely formatted excel file.

    Args:
        orbs: the orbitals object to write the Excel file for.
        out_file: the path to write the Excel file to.
        sfo_energy_type: the energy type for the |SFO| objects to use for analysing orbital interactions.
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
    table_key_fmt = workbook.add_format({'bold': True})
    table_val_fmt = workbook.add_format({'bold': False})
    table_val_float_fmt = workbook.add_format({'bold': False, 'num_format': '0.00'})
    table_val_pctg_fmt = workbook.add_format({'bold': False, 'num_format': '0.0%'})
    table_title_fmt = workbook.add_format({'bold': True, 'bottom': 6, 'font_size': 16, 'italic': True})
    table_title2_fmt = workbook.add_format({'bold': True, 'bottom': 6, 'font_size': 16})
    table_ast_fmt = workbook.add_format({'top': 1})
    table_header_fmt = workbook.add_format({'bold': True, 'bottom': 6})
    table_left = workbook.add_format({'left': 1})
    table_right = workbook.add_format({'right': 1})


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
            if all(isinstance(orbx, pyfmo.orbitals.objects.SFO) for orbx in orbsx):
                frag_name = list(set(orbx.fragment_unique for orbx in orbsx))[0]
                labelx = f'{frag_name} ({formula.molecule(mols[frag_name])})'
            else:
                labelx = 'MO'
            worksheet.merge_range(1, 3, 1, 3 + len(orbsx), labelx, bold_centered_fmt)

            if all(isinstance(orby, pyfmo.orbitals.objects.SFO) for orby in orbsy):
                frag_name = list(set(orby.fragment_unique for orby in orbsy))[0]
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
            if isinstance(orbx, pyfmo.orbitals.objects.SFO) and not orbs.data.calc_info.used_regions:
                name = f'{orbx.fragment_unique}({orbx.name})'
            worksheet.write(2, 3+i, name, bottom_border_fmt)

        for i, orby in enumerate(orbsy):
            name = orby.name
            if isinstance(orby, pyfmo.orbitals.objects.SFO) and not orbs.data.calc_info.used_regions:
                name = f'{orby.fragment_unique}({orby.name})'

            worksheet.write(3+i, 2, name, right_border_fmt)
            column_widths[i] = max(column_widths[i], character.text_width(orby.name, font_size=11))

        # write the values
        for i, row in enumerate(values.T):
            for j, x in enumerate(row):
                if np.isnan(x):
                    worksheet.write(3+i, 3+j, 0)
                    continue

                worksheet.write(3+i, 3+j, x, number_format)

                if number_format is float_fmt:
                    column_widths[i] = max(column_widths[i], character.text_width(f'{x: .1}', font_size=11))
                elif number_format is pctg_fmt:
                    column_widths[i] = max(column_widths[i], character.text_width(f'{x: .1%}%', font_size=11))

        worksheet.hide_zero()
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

            sheet.write(start_row + i, start_column - 1, ' ', table_right)
            sheet.write(start_row + i, start_column + 2, ' ', table_left)

            sheet.write(start_row + i, start_column, variable, table_key_fmt)

            if isinstance(value, float):
                sheet.write(start_row + i, start_column + 1, value, table_val_float_fmt)
            else:
                sheet.write(start_row + i, start_column + 1, value, table_val_fmt)

        sheet.set_column_pixels(start_column, start_column, 93)
        sheet.set_column_pixels(start_column + 1, start_column + 1, 93)

        for j, asterisk in enumerate(asterisks):
            fmt = None
            if j == 0:
                fmt = table_ast_fmt
            sheet.merge_range(start_row + i + j + 1, start_column, start_row + i + j + 1, start_column + 1, f'{"*"*(j+1)} {asterisk}', fmt)

        if len(asterisks) == 0:
            sheet.merge_range(start_row + i + 1, start_column, start_row + i + 1, start_column + 1, ' ', table_ast_fmt)
        
        return len(rows) + len(asterisks) + start_row, start_column + 1


    def make_array_table(rows, title, start_row, start_column, asterisks=[], col_fmts={}, header=[]):
        # sheet.write(start_row, start_column, title, table_title_fmt)
        sheet.merge_range(start_row, start_column, start_row, start_column + len(rows[0]) - 1, title, table_title_fmt)

        for i, head in enumerate(header):
            sheet.write(start_row + 1, start_column + i, head, table_key_fmt)

        sheet.write(start_row + 1, start_column - 1, ' ', table_right)
        sheet.write(start_row + 1, start_column + len(rows[0]), ' ', table_left)
        for i, row in enumerate(rows):
            sheet.write(start_row + i + 2, start_column - 1, ' ', table_right)
            sheet.write(start_row + i + 2, start_column + len(rows[0]), ' ', table_left)

            for j, part in enumerate(row):

                if j in col_fmts:
                    fmt = col_fmts[j]
                elif isinstance(part, float):
                    fmt = table_val_float_fmt
                else:
                    fmt = table_val_fmt

                sheet.write(start_row + i + 2, start_column + j, part, fmt)
                sheet.set_column_pixels(start_column + j, start_column + j, 93)

        for j, asterisk in enumerate(asterisks):
            fmt = None
            if j == 0:
                fmt = table_ast_fmt
            sheet.merge_range(start_row + i + j + 3, start_column, start_row + i + j + 3, start_column + len(rows[0]) - 1, asterisk, fmt)

        return len(rows) + len(asterisks) + start_row + 1, start_column + len(rows[0]) - 1


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
            sheet.set_column_pixels(i, i, width + 20)

        sheet.freeze_panes('A4')
        sheet.autofilter(2, 1, 1+len(rows), len(header))


    # we will write some basic info about the calcualtion in the first sheet
    sheet = workbook.add_worksheet('🛈 Info')
    sheet.set_tab_color('058ED9')
    # write the title cell
    sheet.write(0, 0, 'PyFMO Analysis', title_fmt)

    # write information about the complex
    mols = _get_molecules(orbs.reader)
    rows = {
        'Complex': '',
        'Formula': formula.molecule(mols['complex']),
        'Coords': '[Copy This]                      \n' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols['complex']]),
        'Nº MOs': len(orbs.mos),
        'Nº occ. MOs': len([mo for mo in orbs.mos if mo.occupied]),
        'Nº virt. MOs': len([mo for mo in orbs.mos if not mo.occupied]),
        'Nº frozen cores': orbs.data.MOs.nfrozencores.total,
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
        sfos = orbs.sfos.filter(fragment=fragment)
        rows = {
            'Fragment': fragment,
            'Formula': formula.molecule(mols[fragment]),
            'Coords': '[Copy This]                      \n' + '\n'.join([f'{atom.symbol}\t{atom.x}\t{atom.y}\t{atom.z}' for atom in mols[fragment]]),
            'Nº SFOs': len(sfos),
            'Nº occ. SFOs': len([sfo for sfo in sfos if sfo.occupied]),
            'Nº virt. SFOs': len([sfo for sfo in sfos if not sfo.occupied]),
        }
        next_row, next_col = make_key_value_table(rows, next_row + 2, 1)


    mixer = pyfmo.analysis.mixing.Mixer(orbs, energy_type=sfo_energy_type)
    mixes = mixer.orbital_interactions(N=100)
    # write information about the mixing
    rows = [(str(mix.sfos[0]),
             str(mix.sfos[1]),
             str(mix.mos[0]),
             str(mix.mos[1]),
             mix.strength, 
             mix.fraction, 
             abs(mix.sfos[0] @ mix.sfos[1]), 
             abs(getattr(mix.sfos[0], mix.energy_type) - getattr(mix.sfos[1], mix.energy_type)),
             abs(mix.sfos[0].mulliken_contribution(mix.mos[0])),
             abs(mix.sfos[1].mulliken_contribution(mix.mos[0])),
             abs(mix.sfos[0].mulliken_contribution(mix.mos[1])),
             abs(mix.sfos[1].mulliken_contribution(mix.mos[1])),
             all(
                 [abs(mix.sfos[0].mulliken_contribution(mix.mos[0])) >= 0.02,
                 abs(mix.sfos[1].mulliken_contribution(mix.mos[0])) >= 0.02,
                 abs(mix.sfos[0].mulliken_contribution(mix.mos[1])) >= 0.02,
                 abs(mix.sfos[1].mulliken_contribution(mix.mos[1])) >= 0.02])) for mix in mixes]
    _, next_col = make_array_table(rows, 
            'Orbital Interactions',
            3,
            next_col + 3,
            col_fmts={5: table_val_pctg_fmt},
            header=['SFO1', 'SFO2', 'MO1', 'MO2', 'ΔE⁽²⁾ (kcal mol⁻¹)', 'Frac.*', 'S', 'Δε (eV)**', 'SFO1->MO1', 'SFO2->MO1', 'SFO1->MO2', 'SFO2->MO2', 'Check'],
            asterisks=[
                '* Frac. represents the relative amount of orbital interaction explained by this interaction',
                f'** SFO energy type: {sfo_energy_type}'])

    # write information about the mixing
    mixes = mixer.pauli_repulsions(N=100)
    rows = [(str(mix.sfos[0]),
             str(mix.sfos[1]),
             str(mix.mos[0]),
             str(mix.mos[1]),
             mix.strength, 
             mix.fraction,
             abs(mix.sfos[0] @ mix.sfos[1]),
             abs(mix.sfos[0].mulliken_contribution(mix.mos[0])),
             abs(mix.sfos[1].mulliken_contribution(mix.mos[0])),
             abs(mix.sfos[0].mulliken_contribution(mix.mos[1])),
             abs(mix.sfos[1].mulliken_contribution(mix.mos[1])),
             all(
                 [abs(mix.sfos[0].mulliken_contribution(mix.mos[0])) >= 0.02,
                 abs(mix.sfos[1].mulliken_contribution(mix.mos[0])) >= 0.02,
                 abs(mix.sfos[0].mulliken_contribution(mix.mos[1])) >= 0.02,
                 abs(mix.sfos[1].mulliken_contribution(mix.mos[1])) >= 0.02])) for mix in mixes]
    _, next_col = make_array_table(rows, 
            'Pauli Repulsive Interactions',
            3,
            next_col + 3,
            col_fmts={5: table_val_pctg_fmt},
            header=['SFO1', 'SFO2', 'MO1', 'MO2', 'ΔE⁽⁴⁾ (kcal mol⁻¹)', 'Frac.*', 'S', 'SFO1->MO1', 'SFO2->MO1', 'SFO1->MO2', 'SFO2->MO2', 'Check'],
            asterisks=['* Frac. represents the relative amount of Pauli repulsion explained by this interaction'])

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

        for fragment in orbs.fragments:
            rows[-1].append(mo.fragment_character(fragment))

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

    for fragment in orbs.fragments:
        headers.append(f'{fragment} Character')

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

    sfos_spin = {spin: [sfo for sfo in orbs.sfos if sfo.spin == spin] for spin in orbs.sfos.spins}
    sfos1_spin = {spin: [sfo for sfo in sfos_spin[spin] if sfo.fragment_unique == list(orbs.sfos.fragments)[0]] for spin in orbs.sfos.spins}
    sfos2_spin = {spin: [sfo for sfo in sfos_spin[spin] if sfo.fragment_unique == list(orbs.sfos.fragments)[1]] for spin in orbs.sfos.spins}
    mos_spin = {spin: [mo for mo in orbs.mos if mo.spin == spin or mo.spin == 'AB' or spin == 'AB'] for spin in orbs.sfos.spins}
    spin_names = {'A': '𝛼', 'B': '𝛽'}
    # we add a new sheet for each spin species
    for spin in orbs.sfos.spins:
        cnd_fmt = {
            'type': '3_color_scale',
            'min_color': '63be7b',
            'mid_color': 'white',
            'max_color': '63be7b',
            'min_value': -1,
            'mid_value': 0,
            'max_value': 1,
            'min_type': 'num',
            'mid_type': 'num',
            'max_type': 'num',
        }
        # add the data we want
        name = f"S {spin_names[spin]}" if spin != 'AB' else "S"
        title = f"Overlaps (spin {spin_names[spin]})" if spin != 'AB' else "Overlaps"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin[spin], sfos_spin[spin], _overlap_mat(sfos_spin[spin], sfos_spin[spin]), number_format=float_fmt, tab_color='FF6666', conditional_format=cnd_fmt)
        else:
            make_matrix_sheet(name, title, sfos1_spin[spin], sfos2_spin[spin], _overlap_mat(sfos1_spin[spin], sfos2_spin[spin]), number_format=float_fmt, tab_color='FF6666', conditional_format=cnd_fmt)

    for spin in orbs.sfos.spins:
        name = f"S² {spin_names[spin]}" if spin != 'AB' else "S²"
        title = f"Overlaps² (spin {spin_names[spin]})" if spin != 'AB' else "Overlaps²"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin[spin], sfos_spin[spin], _overlap_mat(sfos_spin[spin], sfos_spin[spin])**2, number_format=float_fmt, tab_color='FF6666')
        else:
            make_matrix_sheet(name, title, sfos1_spin[spin], sfos2_spin[spin], _overlap_mat(sfos1_spin[spin], sfos2_spin[spin])**2, number_format=float_fmt, tab_color='FF6666')

    for spin in orbs.sfos.spins:
        name = f"S²_occ {spin_names[spin]}" if spin != 'AB' else "S²_occ"
        title = f"Pauli Overlaps² (spin {spin_names[spin]})" if spin != 'AB' else "Pauli Overlaps²"
        if not orbs.data.calc_info.used_regions:
            _sfos = [sfo for sfo in sfos_spin[spin] if sfo.occupation > 0]
            make_matrix_sheet(name, title, _sfos, _sfos, _overlap_mat(_sfos, _sfos)**2, number_format=float_fmt, tab_color='FF6666')
        else:
            _sfos1 = [sfo for sfo in sfos1_spin[spin] if sfo.occupation > 0]
            _sfos2 = [sfo for sfo in sfos2_spin[spin] if sfo.occupation > 0]
            if len(_sfos1) == 0 or len(_sfos2) == 0:
                continue
            make_matrix_sheet(name, title, _sfos1, _sfos2, _overlap_mat(_sfos1, _sfos2)**2, number_format=float_fmt, tab_color='FF6666')

    for spin in orbs.sfos.spins:
        name = f"Δε {spin_names[spin]}" if spin != 'AB' else "Δε"
        title = f"Orbital Energy Gap (spin {spin_names[spin]}) (eV)" if spin != 'AB' else "Orbital Energy Gap (eV)"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin[spin], sfos_spin[spin], _energy_gap_mat(sfos_spin[spin], sfos_spin[spin]), number_format=float_fmt, tab_color='4D8B31')
        else:
            make_matrix_sheet(name, title, sfos1_spin[spin], sfos2_spin[spin], _energy_gap_mat(sfos1_spin[spin], sfos2_spin[spin]), number_format=float_fmt, tab_color='4D8B31')

    for spin in orbs.sfos.spins:
        name = f"OI {spin_names[spin]}" if spin != 'AB' else "OI"
        title = f"Orbital Interactions (spin {spin_names[spin]}) (1000/eV)" if spin != 'AB' else "Orbital Interactions (1000/eV)"
        if not orbs.data.calc_info.used_regions:
            oi = _orbint_mat(sfos_spin[spin], sfos_spin[spin])
            oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
            make_matrix_sheet(name, title, sfos_spin[spin], sfos_spin[spin], oi, number_format=float_fmt, tab_color='4D8B31')
        else:
            oi = _orbint_mat(sfos1_spin[spin], sfos2_spin[spin])
            oi[~np.isnan(oi)] *= 1000  # in the case of orbital interactions, there is a mask applied to the matrix and we want to multiply each value with 1000 for easier reading
            make_matrix_sheet(name, title, sfos1_spin[spin], sfos2_spin[spin], oi, number_format=float_fmt, tab_color='4D8B31')

    for spin in orbs.sfos.spins:
        for fragment in orbs.fragments:
            fragment_name = fragment
            if ':' in fragment:
                fragment_name = f'{fragment.split(":")[0]}({fragment.split(":")[1]})'

            name = f"Coeff {fragment_name} {spin_names[spin]}" if spin != 'AB' else f"Coeff {fragment_name}"
            title = f"MO Coefficients from {fragment_name} (spin {spin_names[spin]})" if spin != 'AB' else f"MO Coefficients from {fragment_name}"
            sfos_ = [sfo for sfo in sfos_spin[spin] if sfo.fragment_unique == fragment]
            idx_offset = 0 if spin in ['A', 'AB'] else len(mos_spin[spin])
            sfo_idx = [sfo.index - 1 + idx_offset for sfo in sfos_]
            mo_idx = [mo.index - 1 + idx_offset for mo in mos_spin[spin]]
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
            make_matrix_sheet(name, title, sfos_, mos_spin[spin], coeff.T, number_format=float_fmt, conditional_format=cnd_fmt, tab_color='564D80')

    for spin in orbs.sfos.spins:
        for fragment in orbs.fragments:
            fragment_name = fragment
            if ':' in fragment:
                fragment_name = f'{fragment.split(":")[0]}({fragment.split(":")[1]})'

            name = f"Contr {fragment_name} {spin_names[spin]}" if spin != 'AB' else f"Contr {fragment_name}"
            title = f"Mulliken Contributions from {fragment_name} (spin {spin_names[spin]})" if spin != 'AB' else f"Mulliken Contributions from {fragment_name}"
            sfos_ = [sfo for sfo in sfos_spin[spin] if sfo.fragment_unique == fragment]
            contribs = _contribution_mat(orbs, sfos_, mos_spin[spin])
            cnd_fmt = {
                'type': '2_color_scale',
                'min_color': 'white',
                'max_color': '63be7b',
                'min_value': 0,
                'max_value': 1,
                'min_type': 'num',
                'max_type': 'num',
            }
            make_matrix_sheet(name, title, sfos_, mos_spin[spin], contribs.T, number_format=pctg_fmt, conditional_format=cnd_fmt, tab_color='058ED9')

    for spin in orbs.sfos.spins:
        if not has_kinetic:
            break

        name = f"F {spin_names[spin]}" if spin != 'AB' else "F"
        title = f"Fock (spin {spin_names[spin]})" if spin != 'AB' else "Fock"
        if not orbs.data.calc_info.used_regions:
            make_matrix_sheet(name, title, sfos_spin[spin], sfos_spin[spin], _fock_mat(sfos_spin[spin], sfos_spin[spin]), number_format=float_fmt, tab_color='FF6666')
        else:
            make_matrix_sheet(name, title, sfos1_spin[spin], sfos2_spin[spin], _fock_mat(sfos1_spin[spin], sfos2_spin[spin]), number_format=float_fmt, tab_color='FF6666')

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
                orbs = pyfmo.orbitals.objects.Orbitals(f'../../calculations/PyOrb_testing_2022/Alkyl/{alkyl}/EDA.results/adf.rkf')
                orbs.write_excel2()
                load_time_new[-1].append(perf_counter() - start)
                norbs.append(len(orbs.mos))
            # with timer.timer('new_orbitals.write_excel'):
            #     orbs.write_excel2()

            # with timer.timer('old_orbitals.load_rkf'):
            #     start = perf_counter()
            #     orbs_old = pyfmo.orbitals.Orbitals(f'../../calculations/PyOrb_testing_2022/Alkyl/{alkyl}/EDA.results/adf.rkf')
            #     _contribution_mat(orbs_old, orbs_old.sfos.sfos, orbs_old.mos.mos)
            #     orbs_old.write_excel()
            #     load_time_old[-1].append(perf_counter() - start)
            # # with timer.timer('old_orbitals.mulliken_analysis'):
            # #     _contribution_mat(orbs_old, orbs_old.sfos.sfos, orbs_old.mos.mos)
            # # with timer.timer('old_orbitals.write_excel'):
            #     orbs_old.write_excel()

    print(norbs, load_time_new)
    plt.scatter(norbs, np.array(load_time_new))
    plt.scatter(norbs, np.array(load_time_old))

    plt.title('Comparison Old and New method (incl. writing excel file)')
    plt.plot(np.unique(norbs), [np.mean(times) for times in load_time_new], label='New')
    plt.plot(np.unique(norbs), [np.mean(times) for times in load_time_old], label='Old')
    plt.yscale('log')
    plt.xlabel('Nº MOs')
    plt.ylabel('Loading time (s)')

    plt.legend()
    plt.show()
