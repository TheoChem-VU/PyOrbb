import pyfmo
import os 

j = os.path.join


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



def test_SCF_stage_energy_1():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert round(orbs.sfos['Na(12A)'].site_energy_SCF0, 4) == 3.4891

def test_SCF_stage_energy_2():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert round(orbs.sfos['Cl(9A)'].site_energy, 4) == -9.4931

def test_SCF_stage_energy_3():
	rkf = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex', 'adf.rkf')
	path_SCF0 = j(os.path.split(__file__)[0], 'fixtures','SAOP_DZP', 'complex_SCF0', 'adf.rkf')
	orbs = pyfmo.orbitals.Orbitals(rkf, path_SCF0=path_SCF0)
	assert orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].site_energy and orbs.sfos['Na(1A)'].site_energy != orbs.sfos['Na(1A)'].site_energy_SCF0
	# assert orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].site_energy_SCF0 and orbs.sfos['Na(1A)'].energy != orbs.sfos['Na(1A)'].energy_approximation



if __name__ == '__main__':
	import pytest
	print(j(os.path.split(__file__)[0], 'fixtures','RadicalAddition','adf.rkf'))
	pytest.main()

