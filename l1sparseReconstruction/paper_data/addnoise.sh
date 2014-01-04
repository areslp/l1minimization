#!/bin/sh
for file in `ls *.obj| grep -v noise`
do
    echo $file
    ./addnoise.exe $file 0.1
    ./addnoise.exe $file 0.2
    ./addnoise.exe $file 0.3
    ./addnoise.exe $file 0.4
    ./addnoise.exe $file 0.5
done
