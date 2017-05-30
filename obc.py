"""This module defines an ASE calculator 'obc'.

"""

import numpy as np

from ase.calculators.calculator import Calculator
from ase                        import units

import pybel






def ase2xyz(m):
    """Return a string containing the xyz file.

    input:

    * m :: ase Atoms object.

    """

    N=m.get_number_of_atoms()

    xyz=str(N)+'\nMol  0.0  0.0\n'

    for i in range(N-1):
        xyz+=m.get_chemical_symbols()[i]+ ' '.join("%16.8f" % n for n in m.positions[i])+'\n'

    return xyz + m.get_chemical_symbols()[N-1]+ ' '.join("%16.8f" % n for n in m.positions[N-1])



def e_ff(mol,ffname,ffunitfactor):
    """Potential energy computed with a force field of openbabel.

    * mol :: ase Atoms object
    """

    #Setting up openbabel
    obmol = pybel.ob.OBMol()
    obconversion = pybel.ob.OBConversion()
    obconversion.SetInAndOutFormats("xyz", "mol")

    # Openbabel reads the mol as .xyz string
    obconversion.ReadString(obmol,ase2xyz(mol))


    # Get the ff 
    ff=pybel.ob.OBForceField.FindForceField(ffname)
    #ff.SetLogLevel(pybel.ob.OBFF_LOGLVL_HIGH) 
    #ff.SetLogToStdErr() 
    ff.Setup(obmol)

    # The factor 4.1572299 is R/2 (half the gas constant)
    return ff.Energy() * units.kcal / units.mol * ffunitfactor # Result in eV





def numeric_force(atoms, a, i, ffname,ffunitfactor, d=0.0001):
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




def f_ff(atoms,ffname,ffunitfactor):
    """Evaluate numeric forces.

    Args:
    * atoms     :: ASE Atoms object
    """

    return np.array([[numeric_force(atoms, a, i, ffname,ffunitfactor)
            for i in range(3)] for a in range(len(atoms))])







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
    """



    implemented_properties = ['energy', 'forces']

    default_parameters = {'ff'  : 'uff'}

    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges', 'magmoms']):
        Calculator.calculate(self, atoms, properties, system_changes)

        ffname  = self.parameters.ff

        # Fixing the units for each ForceField
        # Results are given in eV.
        # The factor 4.1572299 is R/2 (half the gas constant)
        if ffname == 'uff' or ffname=='gaff':
            ffunitfactor = 0.24054479161712947364 #=  1.0 / 4.15722990
        if ffname == 'ghemical':
            ffunitfactor = 0.12027239580856473682 #=  0.5 / 4.15722990
        if ffname == 'mmff94' or ffname == 'mmff94s' :
            ffunitfactor = 1.
        #alphas  = self.parameters.alphas
        #gamma   = self.parameters.gamma

        self.results['energy'] = e_ff(self.atoms,ffname,ffunitfactor)
        self.results['forces'] = f_ff(self.atoms,ffname,ffunitfactor)

