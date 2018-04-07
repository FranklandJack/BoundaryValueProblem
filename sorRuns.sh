for i in `seq 0 20`;
do
	omega=$(python -c "print(1+$i/20.0)")

	./poisson -w $omega -c 100 -r 100 -t 100 --SOR -n 0
	sleep 1
done
