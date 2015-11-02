#!/bin/bash

make all

graph="rome99.gr"

for i in `seq 1 11`;
do
	read S T
	./asplot "plotfiles/asline"$i $S $T < $graph
	./bdvplot "plotfiles/bdvline"$i $S $T < $graph
	./blplot "plotfiles/blline"$i $S $T < $graph
done
