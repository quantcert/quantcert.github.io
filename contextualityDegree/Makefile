.PHONY: tests

dependencies:
	make -C lib/BCpackage-0.40 bc2cnf > /dev/null 
	cd ./lib/kissat_gb && \
	if [ ! -d ./build ]; then \
		./configure; \
	fi
	make -C ./lib/kissat_gb > /dev/null

main.out: ./lib/*.c
	gcc -std=c99 -Wall -g -Wextra -Ofast -fopenmp lib/qontextium.c -lm -o main.out

all: dependencies main.out

tests: ./lib/*.c dependencies
	gcc -std=c99 -Wall -g -Wextra -Ofast -fopenmp tests/test.c -lm -o tests.out
	./tests.out 2> /dev/null

clean: 
	-make -C ./lib/kissat_gb/ clean
	-make -C ./lib/BCpackage-0.40/ clean 
	-rm -rf ./lib/kissat_gb/build
	-rm -f main.out
	-rm -f tests.out
