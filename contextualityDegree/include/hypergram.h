/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file hypergram.c
 * @brief This file contains the code used to manipulate the hypergrams
 * described in the paper
 * An abstract structure determines the contextuality degree of observable-based Kochen-Specker proofs
 * by Axel Muller and Alain Giorgetti
 */

#ifndef _SUB_GEO
#define _SUB_GEO

#include "contextuality_degree.h"
#include "bit_vector.h"

/**
 * @brief Context-commutation structure
 * Represents a geometry from which the contextuality degree can be computed
 */
typedef struct hypergram {
    size_t** geometries;/*array of geometries 0 is the end of an array*/
    size_t cpt_geometries;
    size_t cpt_points;
    size_t max_points_per_geometry;
    bit_matrix commutation_matrix;
    bv* assignment;
    size_t n_qubits;
} hypergram;

/**
 * @return true if the gram matrix is valid
 * i.e. it is square and symmetric
 */
bool is_gram_matrix_valid(bit_matrix bm);

/**
 * @brief Returns the incidence matrix of the hypergraph of the contexts of a hypergram
 *
 * @param obss Set of observables, will be sorted as an output
 * @param size Number of observables
 * @param basis (output) the basis
 *
 * @return int Cardinal of the basis
 */
bit_matrix incidence_matrix(hypergram ccs);

/**
 * @brief Checks if a hypergram is valid
 * 
 * @return true If the hypergram is valid
 * i.e. the commutation matrix is valid and the geometries are valid
 */
bool is_hypergram_valid(hypergram ccs);

/**
 * @brief Parses a gram matrix from a file (csv style)
 * 
 * ex:
 * 0,1,1
 * 1,0,1
 * 1,1,0
 * 
 * or
 * 
 * 011
 * 101
 * 110
 * 
 * @param f 
 * @return bit_matrix 
 */
bit_matrix parse_gram_matrix(FILE *f);

/**
 * @brief Parses a hypergraph from a file
 * 
 * ex:
 * 1,2,3
 * 4,5,6
 * 7,8,9
 * 
 * no zero allowed
 * 
 * @param f 
 * @return CCS 
 */
size_t **parse_geometries(FILE *f, size_t *cpt_geometries, size_t *points_per_geometry, size_t *cpt_points);

/**
 * @brief Frees the memory allocated for a hypergram
 *
 * @param ccs Hypergram to free
 */
void hypergram_free(hypergram ccs);

/**
 * @brief Creates a hypergram from a two files
 * The first file contains the geometries
 * ex:
 * 1,2,3
 * 4,5,6
 * 7,8,9
 *
 * The second file contains the gram matrix
 * ex:
 * 011
 * 101
 * 110
 *
 * @param geometries_file
 * @param gram_file
 * @return hypergram or (hypergram){0} if the hypergram is not valid
 */
hypergram hypergram_create_from_file(FILE *geometries_file, FILE *gram_file);

/**
 * @brief Creates a hypergram from a gram matrix file
 * 
 * @param gram_file 
 * @return hypergram 
 */
hypergram hypergram_create_from_gram_file(FILE *gram_file);

/**
 * @brief Prints a hypergram
 * 
 * @param ccs Hypergram to print
 */
void hypergram_print(hypergram ccs);

/**
 * @brief Creates a list of all the observables present in a quantum assignment object
 * 
 * @param qa 
 * @param res 
 * @return size_t 
 */
size_t quantum_assignment_compute_assignment(quantum_assignment qa, bv **res);


/**
 * @brief Transforms a quantum assignment into a hypergram
 */
hypergram quantum_assignment_to_hypergram(quantum_assignment qa);



/**
 * @brief Finds a set bit in a gram matrix
 * 
 * @param bm Bit matrix
 * @param pivot_i (output) index of the row of the first set bit
 * @param pivot_j (output) index of the column of the first set bit
 * 
 * @return true if a set bit was found, false otherwise
 */
bool gram_matrix_find_set_bit(bit_matrix bm, bv *assignment, int n_qubits, size_t *pivot_i, size_t *pivot_j);


/**
 * @brief Finds the assignment of a CCS from its independent gram matrix
 * by following the algorithm described in the paper
 * 
 * Contrary to the paper, instead of using a mutable matrix alone, we use a static matrix 
 * and the assignment to keep track of the anticommutations
 * 
 * @param gram_matrix Independent gram matrix
 * @param assignment (output) Assignment
 * @return int Number of qubits
 */
int pauli_assignment_from_anticommutations(bit_matrix gram_matrix, bv *assignment);

/**
 * @brief Computes the assignment of a hypergram
 * 
 * @param ccs Hypergram
 */
void hypergram_compute_assignment(hypergram *ccs);

/**
 * @brief Transforms a hypergram into a quantum assignment
 */
quantum_assignment hypergram_to_quantum_assignment(hypergram ccs);


/**
 * @brief Converts a gram matrix to a quantum assignment by performing an exhaustive search
 * of hyperedges compatible with the gram matrix
 * 
 * @param gram Bit matrix
 * @return quantum_assignment 
 */
quantum_assignment gram_matrix_to_quantum_assignment(bit_matrix gram);


#endif // _SUB_GEO