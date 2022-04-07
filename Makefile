CC = sw9gcc

build/ssptrsv : build/master.o build/slave.o build/slave_cross.o build/slave_2skew.o
	$(CC) -mhybrid $^ -o $@

build/master.o : src/master.c
	$(CC) -O3 -mhost -o $@ -c $^
build/slave_cross.o : src/slave_cross.c
	$(CC) -O3 -mslave -msimd -mftz -o $@ -c $^
build/slave.o : src/slave.c
	$(CC) -O3 -mslave -msimd -mftz -o $@ -c $^
build/slave_2skew.o : src/slave_2skew.c
	$(CC) -O3 -mslave -msimd -mftz -o $@ -c $^

run: build/ssptrsv
	bsub -b -I -q q_test_ts -n 1 -cgsp 64 -share_size 12000 -host_stack 1024 ./build/ssptrsv

clean :
	rm build/*.o build/ssptrsv
