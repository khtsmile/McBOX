#!/bin/bash

if [ -z $1 ]; then
	N0=01
elif [ $1 == rm ]; then
	rm -f log*
	exit 1
else
	N0=0$1
fi

if [ -z $2 ]; then
	NN=$N0
else
	NN=${N0}0$2
fi

tail -f log${NN}.out
