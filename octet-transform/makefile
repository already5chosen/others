all: test test2 test3

lut_transform_tb.o:lut_transform_tb.cpp
	g++ -O2 -Wall -c $<

ref.o:ref.c
	gcc -O2 -Wall -c $<

uut.o:uut.c
	gcc -O1 -Wall -c $<

uut2.o:uut2.c
	gcc -O1 -Wall -c $<

uut3.o:uut3.c
	gcc -O1 -Wall -c $<

checksum.o:checksum.c
	gcc -O1 -Wall -c $<

test: lut_transform_tb.o ref.o uut.o checksum.o
	g++ -s  $^ -o $@

test2: lut_transform_tb.o ref.o uut2.o checksum.o
	g++ -s  $^ -o $@
	
test3: lut_transform_tb.o ref.o uut3.o checksum.o
	g++ -s  $^ -o $@