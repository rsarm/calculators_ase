
#      AGLYSL01       26.87431        26.87431               AGLYSL01            26.874216    26.874213
echo  "            Name                 E1_ref       E2_ref     Name           E1_calc         Er_calc"
#            AGLYSL01            26.874216    26.874213   AGLYSL01       26.87431        26.87431

for i in `less energies.dat | awk '{print $1}'` ; do
    calc=`python gs.py ${i}.xyz`
    echo -n "$calc   "

    grep $i energies.dat
done
