for i in `seq 0 40`;
do
	omega=$(python -c "print(1.92+$i/1000.0)")

	./poisson -w $omega -c 100 -r 100 -t 100 --SOR -n 0
	sleep 1
done
