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

clean: 
	make -C ./lib/kissat_gb/ clean
	make -C ./lib/BCpackage-0.40/ clean 
	rm -f main.out
