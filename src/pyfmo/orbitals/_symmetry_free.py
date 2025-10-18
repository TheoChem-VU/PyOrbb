from typing import List


def get_irreps(reader: str) -> List[str]:
	symm_symbol = reader.read('Geometry', 'grouplabel').strip()

	if symm_symbol == 'NOSYM':
		return ['A']
	elif symm_symbol == 'ATOM':
		symbol_starts = ['S', 'P', 'D', 'F']
	elif symm_symbol == 'D(LIN)':
		symbol_starts = ['SIGMA.g', 'SIGMA.u', 'PI.g', 'PI.u', 'DELTA.g', 'DELTA.u', 'PHI.g', 'PHI.u']
	else:
		return

	irreps = []
	reader_sections = reader._sections.keys()
	for start in symbol_starts:
		irreps.extend([sec for sec in reader_sections if sec.startswith(start)])

	return irreps



def get_ncbs(reader: str) -> List[str]:
	irreps = get_irreps(reader)
	try:
		ncbs = [reader.read(irrep, 'ncbas') for irrep in irreps]
	except KeyError:
		ncbs = [0 for _ in irreps]
		
	return ncbs