

echo  "            Name                 E1_ref       E2_ref     Name           E1_calc         Er_calc"

for i in `less energies.dat | awk '{print $1}'` ; do
    calc=`python gs.py ${i}.xyz`
    echo -n "$calc   "

    grep $i energies.dat
done
