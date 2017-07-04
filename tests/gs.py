
import sys
from ase.io              import read
from ase.calculators.obc import OBC
from ase                 import units


mol=read(sys.argv[1])


calc=OBC() #['gaff','ghemical','mmff94','mmff94s','uff']

calc.parameters['ref'] = mol
calc.parameters['ff']    = 'mmff94'
calc.setup_ff()

mol.set_calculator(calc)

print 'Initial Energy    = %12.6f eV'   % mol.get_potential_energy()
print 'Initial Max Force = %12.6f eV/A' % mol.get_forces().max()


print '\nRelaxing geometry with autoopt()\n'
calc.autoopt(mol)

print 'Final Energy      = %12.6f eV'   % mol.get_potential_energy()
print 'Final Max Force   = %12.6f eV/A' % mol.get_forces().max()


print '\nEnergy by terms (eV):-----------------------'
e_terms_dict=calc.get_potential_energy_terms(mol,energy_unit=units.eV)


print '\nEnergy by terms (kcal/mol):-----------------'
e_terms_dict=calc.get_potential_energy_terms(mol,energy_unit=units.kcal/units.mol)
