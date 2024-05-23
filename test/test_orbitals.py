import pyfmo
import os 

j = os.path.join

## orbital load test -- radicals ##

def test_load_orbitals():
	pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))


def test_energy1():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['1A_A'].energy, 2) == -381.81


def test_energy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert round(orbs.mos['SOMO_A'].energy, 2) == -6.23


def test_energy3():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO_A'] == orbs.mos['13A_A']


def test_energy4():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SOMO'] == orbs.mos['13A_A']

def test_energy5():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	assert orbs.mos['SUMO'] == orbs.mos['13A_B']

## symmetry test for E' symmetries ##

def test_symm_orbitals():
	pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))

## test to check degeneracy and number of degenerate orbitals

def test_degeneracy():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.mos['5EE1:1'].degenerate is True

def test_degeneracy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.mos['1AA2'].degenerate is False

def test_n_degeneracy():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.mos['5EE1:1'].n_degenerate == 2

def test_n_degeneracy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.mos['1AA2'].n_degenerate == 1

#test for energies with symmetry labels and relative names in the HOMO-LUMO nomenclature
def test_symm_energy1():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert round(orbs.mos['2EEE1:2'].energy, 3) == -10.420

def test_symm_energy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert round(orbs.mos['LUMO+4'].energy, 3) == 0.465
## test in de FMOs 
def test_LUMO_lowered():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert round(orbs.sfos['2(LUMO)'].energy, 3) == -11.189

def test_FMO_degeneracy():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.sfos['2(HOMO-1)'].degenerate is True 

def test_FMO_degeneracy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.sfos['1(3S)'].degenerate is False

def test_FMO_n_degeneracy():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.sfos['2(HOMO-1)'].n_degenerate == 2

def test_FMO_n_degeneracy2():
	orbs = pyfmo.orbitals.Orbitals(j(os.path.split(__file__)[0], 'fixtures','pentafluorophsophate','FragAnal.adf.rkf'))
	assert orbs.sfos['1(2P:x)'].n_degenerate == 3



def test_SCF_site_energy_1():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert round(orbs.sfos['Na(12A)'].site_energy_SCF0, 4) == 5.1043

def test_SCF_site_energy_2():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert round(orbs.sfos['Cl(9A)'].site_energy, 4) == -9.4931

def test_SCF_site_energy_3():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].site_energy and orbs.sfos['Na(1A)'].site_energy != orbs.sfos['Na(1A)'].site_energy_SCF0
	# assert orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].site_energy_SCF0 and orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].energy_approximation

def test_Fock_matrix_diagonal():
	rkf = j(os.path.split(__file__)[0], 'fixtures','FMAT_SFO', 'nh4.adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','FMAT_SFO', 'nh4.adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert orbs.sfos['nh3(3A1)'].energy != orbs.sfos['nh3(3A1)'].site_energy_SCF0
	assert orbs.sfos['nh3(3A1)'].site_energy == orbs.sfos['nh3(3A1)'].site_energy_SCF0


def test_fragments_list():
	rkf = j(os.path.split(__file__)[0], 'fixtures','NH3BH3','adf.rkf')
	sfos = pyfmo.orbitals.sfo.SFOs(kfpath=rkf)
	unsorted_fragments = [sfo.fragment_unique_name for sfo in sfos]
	assert len(unsorted_fragments) != len(sfos.fragments)


if __name__ == '__main__':
	import pytest
	print(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	pytest.main()

