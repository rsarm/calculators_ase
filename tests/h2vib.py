


import sys
import numpy as np

from ase.calculators.obc import OBC
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.atoms import Atoms



mol=Atoms('HH')
mol.positions=np.array([[0.70,0.00,0.00],[0.00,0.00,0.00]])


calc=OBC()
calc.parameters['ff']='uff'

mol.set_calculator(calc)


# Relax
BFGS(mol).run(fmax=0.005)


vib = Vibrations(mol)

vib.run()


vib.summary()

vib.clean()
