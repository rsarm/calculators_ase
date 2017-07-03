
for i in `less energies.dat | awk '{print $1}'` ; do
    python gs.py ${i}.xyz;
done
