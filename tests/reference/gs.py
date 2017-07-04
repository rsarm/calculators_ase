
import sys
from ase.io              import read
from ase.calculators.obc import OBC
from ase                 import units


mol=read(sys.argv[1])


calc=OBC() #['gaff','ghemical','mmff94','mmff94s','uff']

calc.parameters['ff']  = 'mmff94'
calc.parameters['ref'] = 'AMPTRB10.mol2'
calc.setup_ff()

mol.set_calculator(calc)


# Relaxing geometry with autoopt()
calc.autoopt(mol)

print '\nEnergy by terms (kcal/mol):-----------------'
e_terms_dict=calc.get_potential_energy_terms(mol,energy_unit=units.kcal/units.mol)

print ''
