Replicates=2
WorldId=0
while read line
do
	for ((i=0; i<Replicates; i++))
	do
		nohup ./model $WorldId $line &
		((WorldId++))
	done
done < config.txt
