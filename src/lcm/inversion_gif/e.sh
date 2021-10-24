#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as e.sh
# chmod +x e.sh
# ./e.sh
# checked in https://www.shellcheck.net/




printf "make pgm files \n"
gcc e.c -lm -Wall -march=native -fopenmp

if [ $? -ne 0 ]
then
    echo ERROR: compilation failed !!!!!!
    exit 1
fi


export  OMP_DISPLAY_ENV="TRUE"
printf "display OMP info \n"

printf "run the compiled program\n"
time ./a.out > e.txt

export  OMP_DISPLAY_ENV="FALSE"

printf "change Image Magic settings\n"
export MAGICK_WIDTH_LIMIT=100MP
export MAGICK_HEIGHT_LIMIT=100MP

printf "convert all pgm files to png using Image Magic v 6 convert \n"






for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename -- "$file" .pgm)
  # convert from pgm to gif and add text ( level ) using ImageMagic
  # https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html
  convert "$file" -pointsize 50 -annotate +10+100 "${b:0:4}" "${b}".gif
  echo "$file"
done
 
# convert gif files to animated gif
# https://stackoverflow.com/questions/69691394/how-to-convert-images-with-negative-number-as-a-name-to-animation-video
readarray -t files < <(printf '%s\n' *.gif | LC_ALL=C sort -n)
convert "${files[@]}" -resize 600x600 a600.gif




printf "delete all pgm files \n"
rm ./*.pgm

 
echo OK

printf "info about software \n"
echo "$SHELL"
bash --version
make -v
gcc --version
convert -version
convert -list resource
# end



