#!/bin/bash

for i in *.asc1
 do 
 j=`echo $i | sed s/.asc1//g`
 echo $j
done
