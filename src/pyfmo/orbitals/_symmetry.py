from typing import List


def get_irreps(reader: str) -> List[str]:
	symm_symbol = reader.read('Geometry', 'grouplabel').strip()

	match symm_symbol:
		case 'NOSYM':
			return ['A']
		case 'ATOM':
			symbol_starts = ['S', 'P', 'D', 'F']
		case 'D(LIN)':
			symbol_starts = ['SIGMA.g', 'SIGMA.u', 'PI.g', 'PI.u', 'DELTA.g', 'DELTA.u', 'PHI.g', 'PHI.u']
		case _:
			return

	irreps = []
	reader_sections = reader._sections.keys()
	for start in symbol_starts:
		irreps.extend([sec for sec in reader_sections if sec.startswith(start)])

	return irreps
