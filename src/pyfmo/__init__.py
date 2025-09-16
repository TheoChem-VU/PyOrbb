from . import orbitals  # noqa
Orbitals = orbitals.objects.Orbitals
from . import plotting  # noqa
from . import analysis  # noqa
from . import application  # noqa
import re

IRREP_TRANSLATION_LATEX = {
    "E1:1": r"E$_\mathrm{1}^\mathrm{1}$",
    "E1:2": r"E$_\mathrm{1}^\mathrm{2}$",
    "AA": r"A$^{\prime}$",
    "AAA": r"A$^{\prime\prime}$",
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
    "E1:1": "E<sub>1</sub><sup>1<sup>",
    "E1:2": "E<sub>1</sub><sup>2<sup>",
    "AA": "A′",
    "AAA": "A″",
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

def translate_irrep_label(symm_label: str, mode='latex') -> str:
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
            return IRREP_TRANSLATION_HTML[symm_label]

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

        return s

    return symm_label
