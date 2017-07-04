
import sys
from ase.io              import read
from ase.calculators.obc import OBC
from ase                 import units


mol=read(sys.argv[1] + '.mol')


calc=OBC()

calc.parameters['ref'] = sys.argv[1] + '.mol' #mol
calc.parameters['ff']  = 'mmff94'
calc.setup_ff()

mol.set_calculator(calc)

print '%20s     ' % (sys.argv[1]),
print '   %12.6f' % (mol.get_potential_energy() / units.kcal*units.mol),

# Optimize
calc.autoopt(mol)

print '%12.6f'    % (mol.get_potential_energy() / units.kcal*units.mol)
