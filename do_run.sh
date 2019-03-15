todo_list=`awk '{print $0}' delve_paths `


for full_path in $todo_list
do

req=`echo $full_path |awk -F "/" '{print $6}' |awk -F "r" '{print $2}'`
att=`echo $full_path |awk -F "/" '{print $NF}' |awk -F "0" '{print $2}'`
exp=`echo $full_path |awk -F "/" '{print $(NF-1)}' |awk -F "D00" '{print $2}'`
python make_red_catlist.py $full_path ./
python DECADE-expCalib_Y3apass.py --expnum $exp --reqnum $req --attnum $att 

mv *ZP*csv delve_ZP/
mv STD* delve_ZP/
mv Obj* delve_ZP/
mv Merg* delve_ZP/

rm Zero* delve_ZP/ 
rm *Obj*
rm *.fits
rm *match*
rm *std*

done
