
import sys
from ase.io              import read
from ase.calculators.obc import OBC



mol=read(sys.argv[1])


calc=OBC(atoms=mol,ff='uff') #['gaff','ghemical','mmff94','mmff94s','uff']


#mol.set_calculator(calc)

print mol.get_potential_energy()
print mol.get_forces()


#Optimize geometry with autoopt()
calc.autoopt(mol)

print mol.get_potential_energy()
print mol.get_forces()
