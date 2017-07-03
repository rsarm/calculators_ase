"""This module defines an ASE calculator 'OBC' (OpenBabel Calculator).

This calculator uses the force fields implemented in OpenBabel
to compute energies and numerical forces.

The implemented force fields are:
gaff
ghemical
mmff94
mmff94s
uff

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
    # to avoid an empty line at the end of the file.
    for i in range(N-1):
        xyz+=atoms.get_chemical_symbols()[i]+\
        ' '.join("%16.8f" % n for n in atoms.positions[i])+'\n'

    return xyz+atoms.get_chemical_symbols()[N-1]+\
        ' '.join("%16.8f" % n for n in atoms.positions[N-1])
#




def extref2str(_file):
    if _file==None:
        return None

    f=open(_file,'r')
    lines=f.readlines()
    f.close()

    extref_str=''
    for s in lines:
        extref_str+=s

    fileformat = _file.split('.')[-1]

    return {fileformat: extref_str}





def ob_forcefield(atoms,ffname,extref=None):
    """Sets up the OB force field."""

    # Setting up openbabel
    obmol = pybel.ob.OBMol()
    obconversion = pybel.ob.OBConversion()

    if extref==None:
        obconversion.SetInAndOutFormats("xyz", "mol2")
        # Openbabel reads the mol as .xyz string
        obconversion.ReadString(obmol,ase2xyz(atoms))
    else:
        fileformat=extref.keys()[0]
        obconversion.SetInAndOutFormats(fileformat, "mol2")
        obconversion.ReadString(obmol,extref[fileformat])


    # Get the FF
    ff=pybel.ob.OBForceField.FindForceField(ffname)

    ff.SetLogLevel(0) # pybel.ob.OBFF_LOGLVL_HIGH
    ff.SetLogToStdErr()

    ff.Setup(obmol)

    return (ff,obmol)






def e_ff(atoms, obmol, obff, ffunitfactor):
    """Potential energy computed with a force field of openbabel.

    * atoms :: ase Atoms object
    """

    # Update self.obmol coordinates.
    for i,a in enumerate(pybel.ob.OBMolAtomIter(obmol)):
        a.SetVector(atoms.positions[i,0],
                    atoms.positions[i,1],
                    atoms.positions[i,2])

    # Update coordinates in the forcefield.
    obff.SetCoordinates(obmol)

    # The factor 4.1572299 is R/2 (half the gas constant)
    return obff.Energy() * ffunitfactor # Result in eV







def numeric_force(atoms, a, i, obmol, obff, ffunitfactor, d):
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







def f_ff(atoms, obmol, obff, ffunitfactor, nfd):
    """Evaluate numeric forces.

    Args:
    * atoms     :: ASE Atoms object
    """

    return np.array([[numeric_force(atoms, a, i, obmol, obff, ffunitfactor, nfd)
            for i in range(3)] for a in range(len(atoms))])
#




class OBC(Calculator):
    """OBC.

    Example:
    -------

    from ase.io              import read
    from ase.calculators.obc import OBC
    from ase                 import units

    mol=read('mol.xyz')

    calc=OBC()

    calc.parameters['atoms'] = mol
    calc.parameters['ff']    = 'mmff94' # 'gaff', 'ghemical', 'mmff94s', 'uff'

    # Set up the calculator:
    calc.setup_ff()

    mol.set_calculator(calc)

    print mol.get_potential_energy()
    print mol.get_forces()

    # Optimize geometry with autoopt()
    calc.autoopt(mol)
    print mol.get_potential_energy()
    print mol.get_forces()
    """


    implemented_properties = ['energy', 'forces']

    default_parameters = {'nfd'    : 0.0001, # Numeric force displacement lenght
                          'ff'     : 'uff' , # FF name
                          'atoms'  : None  , # reference molecule as ase.Atoms object
                          'extref' : None    # reference molecule as file (any format)
                         }

    nolabel = True




    def __init__(self,**kwargs):
        Calculator.__init__(self, **kwargs)





    def setup_ff(self):
        """Unitializes the FF.

        Note: the parameter 'extref' overrides 'atom'.

        """

        if self.parameters.atoms==None and self.parameters.extref==None:
            raise ValueError(
            """Either 'atoms' or 'extref' need to be specified:
            calc.parameters['atoms']  = atom  # ase.Atoms object
            calc.parameters['extref'] = 'my_ref_molecule.mol2' # Reference molecule
            """)

        # Fixing the units for each ForceField
        # Results are given in eV.
        # The factor 4.1572299 is R/2 (half the gas constant)
        if self.parameters.ff == 'uff'    or self.parameters.ff=='gaff':
            self.ffunitfactor = 0.24054479161712947364 # = 1.0 / 4.15722990
        #
        if self.parameters.ff == 'ghemical':
            self.ffunitfactor = 0.12027239580856473682 # = 0.5 / 4.15722990
        #
        if self.parameters.ff == 'mmff94' or self.parameters.ff == 'mmff94s' :
            self.ffunitfactor = 1.
        #
        self.ffunitfactor = self.ffunitfactor * units.kcal/units.mol

        # Unitializes the FF
        self.obff, self.obmol = ob_forcefield(self.parameters.atoms,
                                              self.parameters.ff,
                                              extref2str(self.parameters.extref)
                                             )





    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        self.results['energy'] = e_ff(self.atoms, self.obmol, self.obff, self.ffunitfactor)
        self.results['forces'] = f_ff(self.atoms, self.obmol, self.obff, self.ffunitfactor, self.parameters.nfd)





    def autoopt(self,atoms,nsteps=500,loglevel=0,algo='sd'):
        """Geometry optimization using the optimizers implemented
        in OpenBabel.
        """

        # The highest is 3 (pybel.ob.OBFF_LOGLVL_HIGH)
        self.obff.SetLogLevel(loglevel)
        self.obff.SetLogToStdErr()

        # Optimization algorithms
        if algo.lower()=='cg' or algo.lower()=='conjugategradients':
            self.obff.ConjugateGradients(nsteps)
        if algo.lower()=='sd' or algo.lower()=='steepestdescent':
            self.obff.SteepestDescent(nsteps)
        #

        #obff.GetCoordinates(obmol)
        self.obff.UpdateCoordinates(self.obmol)

        opt_positions=[[a.x(),a.y(),a.z()] for a in pybel.ob.OBMolAtomIter(self.obmol)]

        atoms.set_positions(np.array(opt_positions))

        self.obff.SetLogLevel(0) # setting loglevel back to 0 again




    def get_potential_energy_terms(self,atoms,prt=True,energy_unit=units.eV):
        """Energy splitted by terms.

        This is a redefinition of get_potential_energies(),
        with a pourpose different to what it's originally intended.

        """

        terms_dict={
        'Bond Stretching'      : self.obff.E_Bond()         *self.ffunitfactor/energy_unit,
        'Angle Bending'        : self.obff.E_Angle()        *self.ffunitfactor/energy_unit,
        'Stretch-Bending'      : self.obff.E_StrBnd()       *self.ffunitfactor/energy_unit,
        'Out-Of-Plane Bending' : self.obff.E_OOP()          *self.ffunitfactor/energy_unit,
        'Torsional'            : self.obff.E_Torsion()      *self.ffunitfactor/energy_unit,
        'Van der Waals'        : self.obff.E_VDW()          *self.ffunitfactor/energy_unit,
        'Electrostatic'        : self.obff.E_Electrostatic()*self.ffunitfactor/energy_unit
        }

        if prt==True:
            for k in terms_dict.keys():
                print '%30s %12.6f' % (k,terms_dict[k])


        etot = sum([terms_dict[k] for k in terms_dict.keys()])

        print '--------------------------------------------'
        print '%30s %12.6f'  % ('E Total',etot)

        return terms_dict















############
# Old code #
############

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




def _f_ff(atoms,ffname,ffunitfactor):
    """Evaluate numeric forces.

    Args:
    * atoms     :: ASE Atoms object
    """

    return np.array([[numeric_force(atoms, a, i, ffname,ffunitfactor)
            for i in range(3)] for a in range(len(atoms))])
#
