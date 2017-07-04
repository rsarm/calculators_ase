

echo  "            Name                 E1_calc      E2_calc    Name           E1_ref          Er_ref"

for i in `less energies.dat | awk '{print $1}'` ; do
    calc=`python gs.py ${i}.xyz`
    echo -n "$calc   "

    grep $i energies.dat
done
