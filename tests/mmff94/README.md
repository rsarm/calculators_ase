Tests from

http://www.ccl.net/cca/data/MMFF94/MMFF94_dative.mol2

and

http://www.ccl.net/cca/data/MMFF94/MMFF94.energies


To run the tests, download both files (MMFF94_dative.mol2 and MMFF94.energies)
and then run **split.py MMFF94_dative.mol2** to get the xyz files from the
list *MMFF94_dative.mol2*. (the header of MMFF94.energies should be removed)

Then run **run.sh** to get a table with the OBC energies compared to those
on the file *MMFF94.energies*
