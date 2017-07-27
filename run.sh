#!/bin/bash
Replicates=10
WorldId=0

for i in {0..32}
do
	while read line
	do
		for ((i=0; i<Replicates; i++))
		do
			nohup time ./model $WorldId $line &
			((WorldId++))
		done
	done < config_${i}.txt

	sleep 2.5h
done