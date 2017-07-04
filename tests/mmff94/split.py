
import sys
import pybel

def write(moliter,form):
    for m in moliter:
        m.write(form, m.title + '.' + form)




#moliter = pybel.readfile('mol2',sys.argv[1])
#write(moliter,'xyz')

moliter = pybel.readfile('mol2',sys.argv[1])
write(moliter,'mol')








