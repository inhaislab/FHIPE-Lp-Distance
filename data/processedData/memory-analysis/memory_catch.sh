#!/bin/bash

N="
8
16
32
64
128
"
M="10"
P="
2
4
6
8
10
"

Mode=("encrypt" "keygen")
for mode in ${Mode[@]} ; do

	for n in $N ; do

		for p in $P ; do
	
			echo "processing N=$n, M=$M, P=$p..."

			filename="memory-$n-$M-$p.log"

			valgrind --tool=massif --massif-out-file=$mode/$filename ./p-norm-ipe-test $mode $n $M $p  
			
			if [ -f "$mode/$filename" ]; then
				ms_print $mode/memory-$n-$M-$p.log > $mode/memory-$n-$M-$p-visual.txt 
			else
				echo "valgrind not export file name : $filename"
			fi

		done
	done
done

echo "done"
