from . import orbitals  # noqa
Orbitals = orbitals.objects.Orbitals
from . import plotting  # noqa
from . import analysis  # noqa
from . import application  # noqa
import re

IRREP_TRANSLATION_LATEX = {
    "E1:1": r"E$_\mathrm{1}^\mathrm{1}$",
    "E1:2": r"E$_\mathrm{1}^\mathrm{2}$",
    "EE1:1": r"E$^{\mathrm{1}\prime}_\mathrm{1}$",
    "EE1:2": r"E$^{\mathrm{2}\prime}_\mathrm{1}$",
    "EEE1:1": r"E$^{\mathrm{1}\prime\prime}_\mathrm{1}$",
    "EEE1:2": r"E$^{\mathrm{2}\prime\prime}_\mathrm{1}$",
    "AA": r"A$^{\prime}$",
    "AAA": r"A$^{\prime\prime}$",
    "AA1": r"A$^{\prime}_\mathrm{1}$",
    "AA2": r"A$^{\prime}_\mathrm{2}$",
    "AAA1": r"A$^{\prime\prime}_\mathrm{1}$",
    "AAA2": r"A$^{\prime\prime}_\mathrm{2}$",
    "A": r"A",
    "A1": r"A$_\mathrm{1}$",
    "A2": r"A$_\mathrm{2}$",
    "B1": r"B$_\mathrm{1}$",
    "B2": r"B$_\mathrm{2}$",
    "T1": r"T$_\mathrm{1}$",
    "T2": r"T$_\mathrm{2}$",
    "E": r"E",
    
    "SIGMA": r"$\Sigma$",
    "PI": r"$\Pi$",
    "DELTA": r"$\Delta$",
    "PHI": r"$\Phi$",

    "S": r"$s$",
    "P": r"$p$",
    "D": r"$d$",
    "F": r"$f$",
}


IRREP_TRANSLATION_HTML = {
    "E1": "E<sub>1</sub>",
    "E1:1": "E<sub>1</sub><sup>1</sup>",
    "E1:2": "E<sub>1</sub><sup>2</sup>",
    "EE1": "E<sup>1</sup>′",
    "EE1:1": "E<sup>1</sup>′<sub>1</sub>",
    "EE1:2": "E<sup>2</sup>′<sub>1</sub>",
    "EEE1": "E<sup>1</sup>′′",
    "EEE1:1": "E<sup>1</sup>′′<sub>1</sub>",
    "EEE1:2": "E<sup>2</sup>′′<sub>1</sub>",
    "AA": "A′",
    "AAA": "A″",
    "AA1": "A′<sub>1</sub>",
    "AA2": "A′<sub>2</sub>",
    "AAA1": "A′′<sub>1</sub>",
    "AAA2": "A′′<sub>2</sub>",
    "A": "A",
    "A1": "A<sub>1</sub>",
    "A2": "A<sub>2</sub>",
    "B1": "B<sub>1</sub>",
    "B2": "B<sub>2</sub>",
    "T1": "T<sub>1</sub>",
    "T2": "T<sub>2</sub>",
    "E": "E",
    
    "SIGMA": "<i>σ</i>",
    "PI": "<i>π</i>",
    "DELTA": "<i>δ</i>",
    "PHI": "<i>φ</i>",

    "S": "<i>s</i>",
    "P": "<i>p</i>",
    "D": "<i>d</i>",
    "F": "<i>f</i>",
}

def translate_irrep_label(symm_label: str, mode='latex', use_formatting=True) -> str:
    if mode == 'latex':
        if symm_label in IRREP_TRANSLATION_LATEX:
            return IRREP_TRANSLATION_LATEX[symm_label]

        if ':' in symm_label:
            subspecies, dimension = symm_label.split(':')
        else:
            subspecies = symm_label
            dimension = None
        
        if '.' in subspecies:
            subspecies, parity = subspecies.split('.')
        else:
            parity = None

        s = IRREP_TRANSLATION_LATEX.get(subspecies, subspecies)
        if parity is not None:
            s += rf'$_\mathrm{{{parity}}}$'
            
        if dimension is not None:

            if re.search(r"\d+", dimension): 
                dimension = re.sub(r"([xyz])(\d+)", r"\1^{\2}", dimension)
            s += rf'$_{{{dimension}}}$'

        return s

    if mode == 'html':
        if symm_label in IRREP_TRANSLATION_HTML:
            s = IRREP_TRANSLATION_HTML[symm_label]
        else:
            if ':' in symm_label:
                subspecies, dimension = symm_label.split(':')
            else:
                subspecies = symm_label
                dimension = None
            
            if '.' in subspecies:
                subspecies, parity = subspecies.split('.')
            else:
                parity = None

            s = IRREP_TRANSLATION_HTML.get(subspecies, subspecies)
            if parity is not None:
                s += f'<sub>{parity}</sub>'
            if dimension is not None:

                if re.search(r"\d+", dimension): 
                    dimension = re.sub(r"([xyz])(\d+)", r"\1<sup>\2</sup>", dimension)
                s += f'<sub><i>{dimension}</i></sub>'

        if not use_formatting:
            for format_tag in ['<i>', '</i>', '<sub>', '</sub>', '<sup>', '</sup>']:
                s = s.replace(format_tag, '')

        return s

    return symm_label


def generate_label(orb, mode='latex', use_formatting=True):
    if mode == 'latex':
        spin_part = {
            'A': r'$\alpha$',
            'B': r'$\beta$'
        }.get(orb.spin, '')
    elif mode == 'html':
        spin_part = {
            'A': '<i>α</i>',
            'B': '<i>β</i>'
        }.get(orb.spin, '')

    if isinstance(orb, orbitals.objects.MO):
        orb_name = f'{orb.name}{spin_part}'
        orb_name = orb_name.replace(orb.symmetry, translate_irrep_label(orb.symmetry, mode=mode))
    else:
        if orb.spin == 'AB':
            orb_name = orb.name
        else:
            orb_name = f'{orb.name}{spin_part}'

        if orb.subspecies.startswith('P:'):
            principal_qn = orb_name.split(':')[0][:-1]
            orb_name = orb_name.replace(principal_qn, str(int(principal_qn) + 1), 1)
        if orb.subspecies.startswith('D:'):
            principal_qn = orb_name.split(':')[0][:-1]
            orb_name = orb_name.replace(principal_qn, str(int(principal_qn) + 2), 1)
        if orb.subspecies.startswith('F:'):
            principal_qn = orb_name.split(':')[0][:-1]
            orb_name = orb_name.replace(principal_qn, str(int(principal_qn) + 3), 1)

        orb_name = orb_name.replace(orb.subspecies, translate_irrep_label(orb.subspecies, mode=mode))

    if not use_formatting:
        for format_tag in ['<i>', '</i>', '<sub>', '</sub>', '<sup>', '</sup>']:
            orb_name = orb_name.replace(format_tag, '')

    return orb_name

