

from ase.io              import read
from ase.calculators.obc import OBC



mol=read('ccc.xyz')


calc=OBC()
calc.parameters['ff']='uff' #['gaff','ghemical','mmff94','mmff94s','uff']


mol.set_calculator(calc)

print mol.get_potential_energy()
print mol.get_forces()

