#!/bin/bash

if [ $# = 0 ] ; then
   echo "  USAGE : $(basename $0 ) HARM-file units"
   echo " Plot the amplitude of S2 S2 and M2 constituents from the HARM-file."
   echo " units is one of m dm dm or mm "
   exit 0
fi

file=$1
unit=$2

export BIMG_PALDIR=$DEVGIT/CHART_COUPE/PALDIR


cat << eof > markampli
0
5
10
15
20
25
30
eof

for wave in S2 S1 M2 ; do
  case $unit in 
  ( 'm'    ) scal=1    ;;
  ( 'dm'   ) scal=10   ;;
  ( 'cm'   ) scal=100  ;;
  ( 'mm'   ) scal=1000 ;;
  esac

chart -clrdata $file  -clrmark markampli -pixel -clrvar ${wave}_A -clrscale $scal -p nrl.pal \
    -clrxypal 0.1 0.95  0.1 0.2 -format PALETTE i2 -title "NATL60-MJM165 (1 month)  $wave" -o $wave.cgm \
    -string 0.5 0.12 1.1 0 " (${unit}) "
ctrans -d sun -res 800x800  $wave.cgm > toto.sun
convert -quality 100 -density 300 toto.sun ${wave}_${file%.nc}.png

done
