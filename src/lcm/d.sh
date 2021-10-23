#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as d.sh
# chmod +x d.sh
# ./d.sh
# checked in https://www.shellcheck.net/




printf "make pgm files \n"
gcc d.c -lm -Wall -march=native -fopenmp

if [ $? -ne 0 ]
then
    echo ERROR: compilation failed !!!!!!
    exit 1
fi


export  OMP_DISPLAY_ENV="TRUE"
printf "display OMP info \n"

printf "run the compiled program\n"
time ./a.out > a.txt

export  OMP_DISPLAY_ENV="FALSE"

printf "change Image Magic settings\n"
export MAGICK_WIDTH_LIMIT=100MP
export MAGICK_HEIGHT_LIMIT=100MP

printf "convert all pgm files to png using Image Magic v 6 convert \n"
# for all pgm files in this directory
for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename "$file" .pgm)
  # convert  using ImageMagic
  convert "${b}".pgm -resize 2000x2000 "${b}".png
  echo "$file"
done


printf "delete all pgm files \n"
rm ./*.pgm

 
echo OK

printf "info about software \n"
bash --version
make -v
gcc --version
convert -version
convert -list resource
# end
