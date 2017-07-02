
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.io import read
from ase import units
from ase.calculators.obc import OBC
from ase.optimize             import QuasiNewton
#from ase.constraints import FixAtoms


import sys
import numpy as np






atoms=read(sys.argv[1])


calc=OBC(atoms=atoms,ff='uff')
#calc.parameters['ff']='uff'

#atoms.set_calculator(calc)

#Relaxing
#QuasiNewton(atoms).run(fmax=0.005);
#or with autoopt
calc.autoopt(atoms)



# Fixing the atoms
#n_fixed=int(sys.argv[2])
#c = FixAtoms(indices=range(n_fixed,atoms.get_number_of_atoms()))
#atoms.set_constraint(c)

####################################################################
############ Only for the Formating of the output ##################
####################################################################
def printenergy(a,f):
    """Function to print the potential and kinetic energy."""
    epot = a.get_potential_energy()
    ekin = a.get_kinetic_energy()
    print >>f, 'Mol', epot, ekin

def write_formated(numpy_array,_format=None):
  """Write numpy arrays in a nice format and
     without brackets.

     Usage:
     >>> x=np.random.rand(5)
     >>> print write_formated(x)
     >>> print read_xyz.write_formated(x, "%12.3f")
  """

  if _format==None:
    return ' '.join("%12.8f" % (k) for k in numpy_array)
  else:
    return ' '.join( _format % (k) for k in numpy_array)
####################################################################
####################################################################




f=open('traj_' + sys.argv[1],'w')

print >>f, atoms.get_number_of_atoms()
printenergy(atoms,f)
for s,atm in enumerate(atoms.positions):
  print >>f,atoms.get_chemical_symbols()[s],\
            write_formated(atm),\
            write_formated(atoms.get_forces()[s])


MaxwellBoltzmannDistribution(atoms, 9000 * units.kB)
dyn = VelocityVerlet(atoms, .08 * units.fs)   # 10.2 * units.f

for i in range(1000):

    dyn.run(1)

    print >>f, atoms.get_number_of_atoms()
    printenergy(atoms,f)
    for s,atm in enumerate(atoms.positions):
      print >>f,atoms.get_chemical_symbols()[s],\
                write_formated(atm),\
                write_formated(atoms.get_forces()[s])


f.close()
