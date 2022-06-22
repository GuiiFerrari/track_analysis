#! /bin/sh

while read name
do
    #echo "filename : $name"
    



./pATTPCexe $name
#echo "./pATTPCexe $name"

#sleep 5 &


done < files_17F.txt

#exit
