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








def ase_to_xyz(atoms):
    """Return a string containing the xyz file.

    input:

    * m :: ase.Atoms object.

    """

    N=atoms.get_number_of_atoms()

    xyz=str(N)+'\nMol  0.0  0.0\n'

    for i in range(N):
        xyz+=atoms.get_chemical_symbols()[i]+\
        ' '.join("%16.8f" % n for n in atoms.positions[i])+'\n'

    return xyz[:-1] # To avoid an empy line at the end.






def set_reference(ref):
    """Returns a dictionary with the format and
    reference molecule as string.

    input:

    * ref :: ase.Atoms object or path to a molecule file.

    """

    # If the reference is an ase.Atoms object:
    if hasattr(ref,'get_chemical_symbols'):
        return {'format':'xyz',
                'molstr': ase_to_xyz(ref)}
    #
    # If reference is a molecule file:
    elif type(ref)==str:
        f=open(ref,'r')
        lines=f.readlines()
        f.close()

        ref_str=''
        for s in lines:
            ref_str+=s

        return {'format': ref.split('.')[-1],
                'molstr': ref_str}






def ob_forcefield(reference_dict, ffname):
    """Sets up the OB force field.

    input:

    * reference_dict :: python dictionary, like
                        {'format':'xyz',
                         'molstr':'2\nMol\nH 0.0 0.0 0.0\nH 0.7 0.0 0.0'}.
    * ffname         :: string with ff name, like 'uff'.

    """

    # Setting up openbabel
    obmol = pybel.ob.OBMol()
    obconversion = pybel.ob.OBConversion()

    obconversion.SetInAndOutFormats(reference_dict['format'], "mol2")
    obconversion.ReadString(obmol,  reference_dict['molstr'])


    # Get the FF
    ff=pybel.ob.OBForceField.FindForceField(ffname)

    ff.SetLogLevel(0) # pybel.ob.OBFF_LOGLVL_HIGH
    ff.SetLogToStdErr()

    ff.Setup(obmol)

    return (ff,obmol)






def e_ff(atoms, obmol, obff):
    """Potential energy computed with a force field of openbabel.

    input:

    * atoms        :: ase.Atoms object.
    * obmol        :: OpenBabel OBMol object.
    * obff         :: Initialized OpenBabel forcefield object.

    """

    # Update self.obmol coordinates.
    for i,a in enumerate(pybel.ob.OBMolAtomIter(obmol)):
        a.SetVector(atoms.positions[i,0],
                    atoms.positions[i,1],
                    atoms.positions[i,2])

    # Update coordinates in the forcefield.
    obff.SetCoordinates(obmol)

    # The factor 4.1572299 is R/2 (half the gas constant)
    return obff.Energy() * units.kcal / units.mol # Result in eV







def numeric_force(atoms, a, i, obmol, obff, d):
    """Evaluate force along i'th axis on a'th atom using finite difference.

    This will trigger two calls to get_potential_energy(), with atom a moved
    plus/minus d in the i'th axial direction, respectively.

    Based in ase.calculators.test

    input:

    * atoms        :: ase.Atoms object.
    * a            :: Atom index (int).
    * i            :: Coordinate index (0,1 or 2).
    * obmol        :: OpenBabel OBMol object.
    * obff         :: Initialized OpenBabel forcefield object.
    * d            :: Numeric force displacement lenght (float).

    """

    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus  = e_ff(atoms, obmol, obff)
    atoms.positions[a, i] -= 2 * d
    eminus = e_ff(atoms, obmol, obff)
    atoms.positions[a, i] = p0


    return (eminus - eplus) / (2 * d)







def f_ff(atoms, obmol, obff, d):
    """Evaluate numeric forces.

    input:

    * atoms        :: ase.Atoms object.
    * obmol        :: OpenBabel OBMol object.
    * obff         :: Initialized OpenBabel forcefield object.
    * d            :: Numeric force displacement lenght (float).

    """

    return np.array([[numeric_force(atoms, a, i, obmol, obff, d)
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

    default_parameters = {'nfd' : 0.0001, # Numeric force displacement lenght
                          'ff'  : 'uff' , # FF name
                          'ref' : None    # reference molecule as file (any format) 
                                          # or ase.Atoms object
                         }

    nolabel = True




    def __init__(self,**kwargs):
        Calculator.__init__(self, **kwargs)





    def setup_ff(self):
        """Initializes the OpenBabel FF."""

        if self.parameters.ref==None:
            raise ValueError(
            """Parameter 'ref' needs to be specified:
            Example:
            calc.parameters['ref']  = atom                  # ase.Atoms object.
            calc.parameters['ref'] = 'my_ref_molecule.mol2' # Reference file.
            """)


        # Initializing the OpenBabel FF
        self.obff, self.obmol = ob_forcefield(set_reference(self.parameters.ref),
                                              self.parameters.ff)





    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        self.results['energy'] = e_ff(self.atoms, self.obmol, self.obff)
        self.results['forces'] = f_ff(self.atoms, self.obmol, self.obff, self.parameters.nfd)





    def autoopt(self,atoms,nsteps=500,loglevel=0,algo='sd'):
        """Geometry optimization using the optimizers implemented
        in OpenBabel.

        input:

        * atoms    :: ase.Atoms object.
        * nsteps   :: Number of optimization steps (int).
        * loglevel :: loglevel for OpenBabel FF (0,1,2 or 3).
        * algo     :: Algorithms for minimization
                      'conjugategradients' (short 'cg') or
                      'steepestdescent' (short 'sd')

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

        input:

        * atoms       :: ase.Atoms object.
        * prt         :: Print energy terms or not (bool)
        * energy_unit :: Output and print units (values defined in ase.units)

        """

        terms_dict={
        'Bond Stretching'      : self.obff.E_Bond()         * units.kcal / units.mol /energy_unit,
        'Angle Bending'        : self.obff.E_Angle()        * units.kcal / units.mol /energy_unit,
        'Stretch-Bending'      : self.obff.E_StrBnd()       * units.kcal / units.mol /energy_unit,
        'Out-Of-Plane Bending' : self.obff.E_OOP()          * units.kcal / units.mol /energy_unit,
        'Torsional'            : self.obff.E_Torsion()      * units.kcal / units.mol /energy_unit,
        'Van der Waals'        : self.obff.E_VDW()          * units.kcal / units.mol /energy_unit,
        'Electrostatic'        : self.obff.E_Electrostatic()* units.kcal / units.mol /energy_unit
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


        ## Fixing the units for each ForceField
        ## Results are given in eV.
        ## The factor 4.1572299 is R/2 (half the gas constant)
        #if self.parameters.ff == 'uff'    or self.parameters.ff=='gaff':
            #self.ffunitfactor = 1. #0.24054479161712947364 # = 1.0 / 4.15722990
        ##
        #if self.parameters.ff == 'ghemical':
            #self.ffunitfactor = 1 #0.12027239580856473682 # = 0.5 / 4.15722990
        ##
        #if self.parameters.ff == 'mmff94' or self.parameters.ff == 'mmff94s' :
            #self.ffunitfactor = 1.
        ##
        #self.ffunitfactor = self.ffunitfactor * units.kcal/units.mol


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
