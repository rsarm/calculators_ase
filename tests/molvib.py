


import sys
import numpy as np

from ase.calculators.obc import OBC
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.atoms import Atoms
from ase.io    import read



# Run without command-line arguments to get the H2 molecule
# or pass a xyz file as command-line argument.
if len(sys.argv)==1:
   mol=Atoms('HH')
   mol.positions=np.array([[0.70,0.00,0.00],[0.00,0.00,0.00]])
else:
   mol=read(sys.argv[1])

calc=OBC() #['gaff','ghemical','mmff94','mmff94s','uff']

calc.parameters['atoms'] = mol
calc.parameters['ff']    = 'uff'
calc.setup_ff()

mol.set_calculator(calc)


# Relax
#BFGS(mol).run(fmax=0.005)
calc.autoopt(mol)

vib = Vibrations(mol)

vib.run()


vib.summary()

vib.clean()
