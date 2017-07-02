"""This module defines an ASE calculator 'obc'.

This calculator uses the force fields implemented in OpenBabel
to compute energies and numerical forces from them.

"""

import numpy as np

from ase.calculators.calculator import Calculator
from ase                        import units

import pybel








def ase2xyz(atoms):
    """Return a string containing the xyz file.

    input:

    * m :: ase Atoms object.

    """

    N=atoms.get_number_of_atoms()

    xyz=str(N)+'\nMol  0.0  0.0\n'

    # Loop first over the first N-1 atoms.
    # The last line of the xyz file is added after the loop
    # to avoid the an empty line at the end of the file.
    for i in range(N-1):
        xyz+=atoms.get_chemical_symbols()[i]+ ' '.join("%16.8f" % n for n in atoms.positions[i])+'\n'

    return xyz + atoms.get_chemical_symbols()[N-1]+ ' '.join("%16.8f" % n for n in atoms.positions[N-1])





#In [120]: for i,a in enumerate(pybel.ob.OBMolAtomIter(m.OBMol)):
#     ...:     a.SetVector(i*3*0.2,0.1,0.0)
#
#In [121]: ff.SetCoordinates(m.OBMol)
#Out[121]: True
#
#In [122]: ff.Energy()
#Out[122]: 868.463711181507




def ob_forcefield(atoms,ffname):
    """Sets up the OB force field."""

    #Setting up openbabel
    obmol = pybel.ob.OBMol()
    obconversion = pybel.ob.OBConversion()
    obconversion.SetInAndOutFormats("xyz", "mol")

    # Openbabel reads the mol as .xyz string
    obconversion.ReadString(obmol,ase2xyz(atoms))


    # Get the FF
    ff=pybel.ob.OBForceField.FindForceField(ffname)
    #ff.SetLogLevel(pybel.ob.OBFF_LOGLVL_HIGH)
    #ff.SetLogToStdErr()
    ff.Setup(obmol)

    return (ff,obmol)






def e_ff(atoms, obmol, obff, ffunitfactor):
    """Potential energy computed with a force field of openbabel.

    * atoms :: ase Atoms object
    """

    for i,a in enumerate(pybel.ob.OBMolAtomIter(obmol)):
        a.SetVector(atoms.positions[i,0],
                    atoms.positions[i,1],
                    atoms.positions[i,2])

    obff.SetCoordinates(obmol)

    # The factor 4.1572299 is R/2 (half the gas constant)
    return obff.Energy() * units.kcal / units.mol * ffunitfactor # Result in eV




def _e_ff(atoms,ffname,ffunitfactor):
    """Potential energy computed with a force field of openbabel.

    * atoms :: ase Atoms object
    """

    obff=ob_forcefield(atoms,ffname)[0]

    # The factor 4.1572299 is R/2 (half the gas constant)
    return obff.Energy() * units.kcal / units.mol * ffunitfactor # Result in eV






def _numeric_force(atoms, a, i, ffname,ffunitfactor, d=0.0001):
    """Evaluate force along i'th axis on a'th atom using finite difference.

    This will trigger two calls to get_potential_energy(), with atom a moved
    plus/minus d in the i'th axial direction, respectively.

    Based in ase.calculators.test
    """
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus  = e_ff(atoms,ffname,ffunitfactor)
    atoms.positions[a, i] -= 2 * d
    eminus = e_ff(atoms,ffname,ffunitfactor)
    atoms.positions[a, i] = p0


    return (eminus - eplus) / (2 * d)


def numeric_force(atoms, a, i, obmol, obff, ffunitfactor, d=0.0001):
    """Evaluate force along i'th axis on a'th atom using finite difference.

    This will trigger two calls to get_potential_energy(), with atom a moved
    plus/minus d in the i'th axial direction, respectively.

    Based in ase.calculators.test
    """
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus  = e_ff(atoms, obmol, obff, ffunitfactor)
    atoms.positions[a, i] -= 2 * d
    eminus = e_ff(atoms, obmol, obff, ffunitfactor)
    atoms.positions[a, i] = p0


    return (eminus - eplus) / (2 * d)




def _f_ff(atoms,ffname,ffunitfactor):
    """Evaluate numeric forces.

    Args:
    * atoms     :: ASE Atoms object
    """

    return np.array([[numeric_force(atoms, a, i, ffname,ffunitfactor)
            for i in range(3)] for a in range(len(atoms))])
#




def f_ff(atoms, obmol, obff, ffunitfactor):
    """Evaluate numeric forces.

    Args:
    * atoms     :: ASE Atoms object
    """

    return np.array([[numeric_force(atoms, a, i, obmol, obff, ffunitfactor, d=0.0001)
            for i in range(3)] for a in range(len(atoms))])
#




class OBC(Calculator):
    """OBC.

    Example:
    -------

    from ase.io              import read
    from ase.calculators.obc import OBC

    mol=read('mol.xyz')

    calc=OBC()
    calc.parameters['ff']='uff' #['gaff','ghemical','mmff94','mmff94s','uff']


    mol.set_calculator(calc)

    print mol.get_potential_energy()
    print mol.get_forces()

    #Optimize geometry with autoopt()
    calc.autoopt(mol)
    print mol.get_potential_energy()
    print mol.get_forces()
    """



    implemented_properties = ['energy', 'forces']

    default_parameters = {'ff'  : 'uff'}

    nolabel = True




    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

        self.obff, self.obmol = ob_forcefield(kwargs['atoms'],kwargs['ff'])






    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        ffname  = self.parameters.ff

        # Fixing the units for each ForceField
        # Results are given in eV.
        # The factor 4.1572299 is R/2 (half the gas constant)
        if ffname == 'uff'    or ffname=='gaff':
            ffunitfactor = 0.24054479161712947364 #=  1.0 / 4.15722990
        if ffname == 'ghemical':
            ffunitfactor = 0.12027239580856473682 #=  0.5 / 4.15722990
        if ffname == 'mmff94' or ffname == 'mmff94s' :
            ffunitfactor = 1.
        #alphas  = self.parameters.alphas
        #gamma   = self.parameters.gamma

        #self.results['energy'] = e_ff(self.atoms,ffname,ffunitfactor)
        #self.results['forces'] = f_ff(self.atoms,ffname,ffunitfactor)

        self.results['energy'] = e_ff(self.atoms, self.obmol, self.obff, ffunitfactor)
        self.results['forces'] = f_ff(self.atoms, self.obmol, self.obff, ffunitfactor)





    def autoopt(self,atoms,nsteps=500,loglevel=0,algo='sd'):
        """Geometry optimization using the optimizer implemented
        in OpenBabel.
        """

        obff,obmol=ob_forcefield(atoms,self.parameters.ff)

        # The highest is 3 (pybel.ob.OBFF_LOGLVL_HIGH)
        obff.SetLogLevel(loglevel)
        obff.SetLogToStdErr()

        # Optimization algorithms
        if algo.lower()=='cg' or algo.lower()=='conjugategradients':
            ff.ConjugateGradients(nsteps)
        if algo.lower()=='sd' or algo.lower()=='steepestdescent':
            obff.SteepestDescent(nsteps)
        #
        #obff.GetCoordinates(obmol)
        obff.UpdateCoordinates(obmol)

        opt_positions=[[a.x(),a.y(),a.z()] for a in pybel.ob.OBMolAtomIter(obmol)]

        atoms.set_positions(np.array(opt_positions))
