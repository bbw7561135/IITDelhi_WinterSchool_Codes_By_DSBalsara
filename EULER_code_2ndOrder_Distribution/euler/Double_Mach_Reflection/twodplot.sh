#!/bin/bash

nx=96
ny=96

xmin=-0.5
xmax=0.5
ymin=-0.5
ymax=0.5

nimages=15

xcoord='x'
ycoord='y'

datafilename='rhoa prsa vlxa vlya vlza'
fortran_compiler='gfortran'
gnuplot='gnuplot'




#-------------------------------------------------------------------
#-------------------------------------------------------------------

rm -f *.jpeg

cp readdata_2d_src readdata_2d.F

file='readdata_2d.F'
rm -f readdata_temp.F
sed -e "s;nxnum;${nx};g" \
    -e "s;nynum;${ny};g" \
    -e "s;nslice;${nimages};g" \
    -e "s;"xmin1";${xmin};g" \
    -e "s;"xmax1";${xmax};g" \
    -e "s;"ymin1";${ymin};g" \
    -e "s;"ymax1";${ymax};g" \
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

file='plotscript_2d'

rm -f plotscript_temp
sed -e "s;"outputfile";"${datafilename1}.jpeg";g" \
    -e "s;"datafile";"${datafilename1}_rewrite";g" \
    -e "s;"datamin1";${mindata};g" \
    -e "s;"datamax1";${maxdata};g" \
    -e "s;"xmin1";${xmin};g" \
    -e "s;"xmax1";${xmax};g" \
    -e "s;"ymin1";${ymin};g" \
    -e "s;"ymax1";${ymax};g" \
    -e "s;"xcoordinate";${xcoord};g" \
    -e "s;"ycoordinate";${ycoord};g" \
    -e "s;"variablename";${datafilename1};g" \
       "$file" > plotscript_temp

${gnuplot} plotscript_temp
rm -f plotscript_temp
#./fixbb ${datafilename1}.eps
#convert -density 1000 -flatten ${datafilename1}.eps ${datafilename1}.jpeg

done < filelist
rm -f filelist
rm -f *_rewrite
rm -f readdata_temp.F

done
rm -f readdata

exit 0
