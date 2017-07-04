
import sys

def write_xyz(lines,nat,label):
    """xxx."""

    fout=open(label+'.xyz','w')

    print >>fout, nat,'\nMol 0.0 0.0'
    for lm in lines[l+6:l+6+nat]:
        lms=lm.split()
        print >>fout, lms[1][0],lms[2],lms[3],lms[4]

    fout.close()


def write_mol2(lines,nat,label):
    """xxx."""

    fout=open(label+'.xyz','w')

    print >>fout, nat,'\nMol 0.0 0.0'
    for lm in lines[l+6:l+6+nat]:
        lms=lm.split()
        print >>fout, lms[1][0],lms[2],lms[3],lms[4]

    fout.close()



fin=open(sys.argv[1],'r')
lines=fin.readlines()
fin.close()


l=0
while(l<len(lines)):
    if lines[l].split()[0]=='@<TRIPOS>MOLECULE':
        label=lines[l+1].split()[0]
        nat=int(lines[l+2].split()[0])

        write_xyz(lines,nat,label)

        l=l+nat

    l=l+1
