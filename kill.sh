#!/bin/bash


# A way to increase number by 1
# 	- var=$((var+1))
# 	- ((var=var+1))
# 	- ((var+=1))
# 	- ((var++))

# ${arr[*]} 	: all of the items in the array
# ${!arr[*]}	: all of the indice in the array
# ${#arr[*]}	: number of iems in the array
# ${#arr[0]}	: length of item zero


###############################################################################
# input variables
N0=$1
N1=$2

if [ -z $N0 ]
then
	echo " Specify a node number"
	exit 1
fi

if [ -z $N1 ]
then
	N1=$N0
fi


###############################################################################
# reading a mpi log file
ps -ef | grep 'mpirun' > mpilog.txt
I=0
while read line
do
	script[I]=$line
	((I++))
done < mpilog.txt
rm -f mpilog.txt


###############################################################################
# klling the job
nl=${#script[*]}
for (( I=0; I<nl; I++ )); do
	str=($(echo "${script[$I]}" | tr ' ' '\n'))
	for  (( J=$N0; J<=$N1; J++)); do
		if [ "${str[9]}" == "spacen0"${J} ]
		then
			kill -9 ${str[1]}
			echo "[${str[9]}] job ${str[1]} has been killed"
		fi
	done
done


