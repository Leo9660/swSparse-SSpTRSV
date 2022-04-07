# swSparse-SSpTRSV
SpTRSV for structured-grid on SW26010-Pro

# Compile and run our code
sh ./build.sh

# Notification
The main purpose of this code repository is to show the code. You can only compile and run this code on SW26010-Pro. Some header files we use here is unique on the platform.

# About our code
`master.c` is for threads spawn and correctness test.
`slave.c`, `slave_2skew.c` and `slave_cross.c` are code for CPEs, written with athread interfaces (like POSIX), in three different cases.
