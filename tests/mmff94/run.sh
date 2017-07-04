

echo  "            Name                 OBC         OBC (rel)     Name           OPTIMOL         BatchMin"
echo  "---------------------------------------------------------------------------------------------------"

for i in `less MMFF94.energies | awk '{print $1}'` ; do
    calc=`python gs.py $i`
    echo -n "$calc   | "

    grep $i MMFF94.energies 
done
