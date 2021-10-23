#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as d.sh
# chmod +x d.sh
# ./d.sh
# checked in https://www.shellcheck.net/




printf "resize all p files to png using Image Magic v 6 convert \n"
# for all pgm files in this directory
for file in *.png ; do
  # b is name of file without extension
  b=$(basename "$file" .png)
  # convert  using ImageMagic
  convert "${b}".png -resize 600x600 "${b}_600".png
  echo "$file"
done


