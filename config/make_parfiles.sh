#! /bin/sh

while read run vdrift tb gain
do
   echo "Run : $run  $vdrift  $tb $gain"
    


sed "s/0.020/$vdrift/g" parTex.txt > parTex_$run.par
sed -i "s/410/$tb/g"  parTex_$run.par  
sed -i "s/1.00/$gain/g"  parTex_$run.par  


done < list_vdrift.txt


