# Please change the "-q" option for different queues
bsub -b -I -q q_test_ts -n 1 -cgsp 64 -share_size 12000 -host_stack 1024 ./build/ssptrsv
