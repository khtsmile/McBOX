for I in {1..5}
do
	mpirun -machinefile spacen01 -np 1 ./McCrit_beta.out 28 > log${I}.out &
	wait
done
