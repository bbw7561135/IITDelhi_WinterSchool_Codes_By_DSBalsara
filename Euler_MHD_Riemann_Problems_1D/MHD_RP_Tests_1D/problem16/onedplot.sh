#!/bin/bash

nx=400

xmin=-0.5
xmax=0.5

nimages=1

xcoord='x'

datafilename='rh01 pr01 vx01 vy01 vz01 by01 bz01'
fortran_compiler='gfortran'
gnuplot='gnuplot'




#-------------------------------------------------------------------
#-------------------------------------------------------------------

rm -f *.jpeg

file='readdata_1d.F'
rm -f readdata_temp.F
sed -e "s;nxnum;${nx};g" \
    -e "s;nslice;${nimages};g" \
    -e "s;"xmin1";${xmin};g" \
    -e "s;"xmax1";${xmax};g" \
       "$file" > readdata_temp.F

${fortran_compiler} readdata_temp.F -o readdata


for datafilename1 in $datafilename
do

rm -f *_rewrite
rm -f min max
rm -f filelist

echo ${datafilename1}
./readdata ${datafilename1}

read mindata < min
read maxdata < max

rm -r min max 

echo "datamin = "$mindata, "datamax = "$maxdata


while read datafilename1;
do
  echo $datafilename1

file='plotscript_1d'

rm -f plotscript_temp
sed -e "s;"outputfile";"${datafilename1}.jpeg";g" \
    -e "s;"datafile";"${datafilename1}_rewrite";g" \
    -e "s;"datamin1";${mindata};g" \
    -e "s;"datamax1";${maxdata};g" \
    -e "s;"xmin1";${xmin};g" \
    -e "s;"xmax1";${xmax};g" \
    -e "s;"xcoordinate";${xcoord};g" \
    -e "s;"variablename";${datafilename1};g" \
       "$file" > plotscript_temp

${gnuplot} plotscript_temp
rm -f plotscript_temp

done < filelist
rm -f filelist
rm -f *_rewrite
rm -f readdata_temp.F

done
rm -f readdata

exit 0
