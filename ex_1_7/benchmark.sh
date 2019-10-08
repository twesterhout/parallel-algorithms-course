#!/bin/bash

PREFIX="build"

main()
{
	# declare -a BOUNDS=(
	# 	1000 10000 100000 1000000 2000000
	# 	3000000 4000000 5000000 6000000 7000000
	# 	8000000 9000000 10000000 11000000 12000000
	# 	13000000 14000000 15000000 16000000 17000000
	# 	18000000 19000000 20000000
	# )
	declare -a BOUNDS=(
		100000000 200000000
		300000000 400000000
		500000000 600000000
		700000000 800000000
		900000000 1000000000
	)
	declare -a PROCS=(1 2 3 4)
	declare -r block_size=5000000
	declare -r PARALLEL=1
	for upper_bound in ${BOUNDS[@]}; do
		if [ $PARALLEL -eq 1 ]; then
			printf '%lu\t' $upper_bound
			for num_procs in ${PROCS[@]}; do
				{ time -p mpirun -n $num_procs \
					"$PREFIX/generate" -b $block_size -q $upper_bound; } 2>&1 \
						| grep real | cut -d' ' -f 2
			done | tr '\n' '\t'
			printf '\n'
		else
			printf '%lu\t' $upper_bound
			{ time -p "$PREFIX/generate" -q $upper_bound; } 2>&1 \
					| grep real | cut -d' ' -f 2
		fi
	done
}

main
