all: test

lut_transform_tb.o:lut_transform_tb.cpp
	g++ -O2 -Wall -c $<

ref.o:ref.c
	gcc -O2 -Wall -c $<

uut.o:uut.c
	gcc -O1 -Wall -c $<

checksum.o:checksum.c
	gcc -O1 -Wall -c $<

test: lut_transform_tb.o ref.o uut.o checksum.o
	g++ -s  $^ -o $@