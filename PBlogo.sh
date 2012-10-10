#! /bin/bash

# 2012 - P. Poulain, A. G. de Brevern

usage="Usage: PBlogo.sh file.transfac [start end]"

# check number of arguments
if [ $# -eq 0 ]
then
    echo "Not enough argument"
    echo $usage
    exit
fi

# check first argument is an existing regular file
if [ ! -f $1 ]
then
    echo "$1 is not a regular file"
    echo $usage
    exit
fi

echo "calling WebLogo with $1 to build PBlogo"

# data
outfile=${1%transfac}logo.pdf
title=${1%.PB.transfac}

# default weblogo
weblogo \
--format pdf \
--errorbars NO \
--fineprint "PBlogo" \
--title $title \
--color "#1240AB" d      "strand main" \
--color "#1240AB" abcdef "strand others" \
--color "#0BD500" ghij "coil" \
--color "#FD0006" m     "helix main" \
--color "#FD0006" klnop "helix others" \
--composition none \
--datatype transfac \
-s large -f $1 -o $outfile

# weblogo with lower and upper residue index
if [ $# -eq 3 -a $2 -lt $3 ]
then
outfile=${1%transfac}logo.$2-$3.pdf
weblogo \
--format pdf \
--errorbars NO \
--fineprint "PBlogo" \
--title $title \
--color "#1240AB" d      "strand main" \
--color "#1240AB" abcdef "strand others" \
--color "#0BD500" ghij "coil" \
--color "#FD0006" m     "helix main" \
--color "#FD0006" klnop "helix others" \
--lower $2 \
--upper $3 \
--composition none \
--datatype transfac \
-s large -f $1 -o $outfile
fi

echo "wrote $outfile"

#--composition "{'a':3.89,'b':4.41,'c':8.12,'d':18.85,'e':2.45,'f':6.68,'g':1.15,'h':2.40,'i':1.86,'j':0.83,'k':5.45,'l':5.46,'m':30.22,'n':1.99,'o':2.77,'p':3.47}" \
