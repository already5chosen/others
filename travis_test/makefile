cc = clang
CC = clang++
COPT = -Wall -O2 -mavx

all: travis_test.exe  travis_test_wlp.exe

travis_test.o: travis_test.cpp
	${CC} -c ${COPT} $< -o $@ 

travis_test_wlp.o: travis_test.cpp
	${CC} -c ${COPT} -DWINDOWS_LARGE_PAGES $< -o $@ 

rigth_rite.o: rigth_rite.c
	${cc} -c ${COPT} $< -o $@ 

memset32_c.o: memset32_c.c
	${cc} -c ${COPT} $< -o $@ 

travis_test.exe: travis_test.o memset32_c.o
	${CC} -s ${COPT} $^ -o $@ 

travis_test_wlp.exe: travis_test_wlp.o memset32_c.o rigth_rite.o
	${CC} -s ${COPT} $^ -o $@ 
