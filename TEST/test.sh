#!/bin/bash

FILES=` ls *.f90 `

for mat in $FILES
do 
echo ${mat}
diff ${mat} /home/lsg/HSSPACK_Serial/SRC/${mat}
done 
