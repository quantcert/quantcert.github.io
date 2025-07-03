.PHONY: tests

# Compiler and flags
CC = gcc
CFLAGS = -Wall -g -Wextra -O3 -fopenmp -I./include -std=c99

# Source files
SRC_FILES = src/constants.c src/bit_vector.c src/hashset.c src/bv.c src/complex_int.c src/quantum_assignment.c src/config_checker.c src/contextuality_degree.c src/quadrics.c src/hypergram.c src/cayley_hexagon.c
OBJ_FILES = $(SRC_FILES:src/%.c=build/%.o)

# Create build directory if it doesn't exist
$(shell mkdir -p build)

dependencies:
	make -C external/BCpackage-0.40 bc2cnf > /dev/null 
	cd ./external/kissat_gb && \
	if [ ! -d ./build ]; then \
		./configure; \
	fi
	make -C ./external/kissat_gb > /dev/null

# Target executable
main.out: $(OBJ_FILES) build/qontextium.o
	$(CC) $(CFLAGS) $(OBJ_FILES) build/qontextium.o -lm -o main.out

# Compile source files to object files
build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

all: dependencies main.out

tests: $(OBJ_FILES) dependencies
	$(CC) $(CFLAGS) $(OBJ_FILES) tests/test.c -lm -o tests.out
	./tests.out 2> /dev/null

clean: 
	rm -f build/*.o
	-make -C ./external/kissat_gb/ clean
	-make -C ./external/BCpackage-0.40/ clean 
	-rm -rf ./external/kissat_gb/build
	-rm -f main.out
	-rm -f tests.out

example.out: $(OBJ_FILES) build/example.o
	$(CC) $(CFLAGS) $(OBJ_FILES) build/example.o -lm -o example.out
