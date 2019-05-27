#!/bin/bash

########################
# 80c = 710nc  ( 294K) #
# 81c = 711nc  ( 600K) #
# 82c = 712nc  ( 900K) #
# 83c = 713nc  (1200K) #
# 84c = 714nc  (2500K) #
# 85c = 715nc  (   0K) #
########################

LIB=$1

# find the current temperature
for ((TEMP=710; TEMP<=715; TEMP++))
do
	if grep -q $TEMP ./inputfile/inventory.inp
	then
		THIS=$TEMP
		break
	fi
done

# change the temperature
if [ ! -z $THIS ]; then
TTOT='1,$s/'$THIS'/'$LIB'/g'
sed -i $TTOT ./inputfile/CE_mat.inp
sed -i $TTOT ./inputfile/inventory.inp
echo LIBRARY FROM $THIS TO $LIB
fi


	






