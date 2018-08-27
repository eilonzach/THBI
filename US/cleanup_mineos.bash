#!/bin/bash

shopt -s nullglob

files=( *.asc* )
if (( ${#files[@]} )); then
    
for i in *_*.asc
 do 
 j=`echo $i | sed s/.asc//g`
 echo $j
 rm $j*
done

for i in *_*.q
 do 
 j=`echo $i | sed s/.q//g`
 echo $j
 rm $j*
done

fi
