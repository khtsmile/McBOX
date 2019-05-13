#! /bin/bash

tot_batch=50
mkdir -p output

rm -f ./output/*

for ((i=1; i<=$tot_batch; i++)); do 
	echo "=============================================="
	echo "               $i th BATCH RUN                "
	echo "=============================================="
	fname=batch_$i
	file=./output/$fname
	touch $file
	./run_McBOX
	cat spectral_data >> $file
	
done 
