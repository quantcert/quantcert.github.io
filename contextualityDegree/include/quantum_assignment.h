/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file quantum_assignment.c
 * @brief Data structure and functions to handle quantum assignments
 * 
 */


#ifndef QUANTUM_ASSIGNMENT
#define QUANTUM_ASSIGNMENT

#include "bv.h"

/**
 * @brief quantum assignment used to check contextuality
 *
 * @param geometry_indices lists all the geometries of the geometry array to take into account to compute the contextuality degree
 * @param geometries list of all (not necessary) the geometries
 * @param cpt_geometries number of geometries checked
 * @param points_per_geometry number of observables in each geometry
 * @param n_qubits number of qubits per observable
 *
 */
typedef struct
{
    size_t *geometry_indices;
    bv **geometries;
    size_t cpt_geometries;
    size_t points_per_geometry;
    int n_qubits;

    bool *lines_negativity;
} quantum_assignment;

bool quantum_assignment_autofill_indices(quantum_assignment* qa);

/**
 * @brief returns true if the product of all the observables is minus the identity,
 * and false if it is the identity matrix
 * If the product is not the identity modulo some phase, an error is printed and
 * the result can be anything
 *
 * @param geometry array of observables
 * @param size size of the observable array
 * @param n_qubits number of qubits we work on
 * @param verbose if true prints each step of the operation
 */
bool is_negative_custom(bv geometry[], int size, int n_qubits, bool verbose, FILE *output);
bool is_negative(bv geometry[], int size, int n_qubits);

/**
 * @brief Prints a quantum assignment
 * 
 * @param qa 
 */
void print_quantum_assignment(quantum_assignment* qa);

/**
 * @brief Prints a quantum assignment to a file in CSV format
 * 
 * @param qa 
 * @param output 
 */
void quantum_assignment_to_CSV(quantum_assignment qa, FILE *output);

/**
 * @brief writes clauses to a file to be read by ortools
 *
 * @param quantum_assignment
 *
 * @param points_per_geometry
 */
void quantum_assignment_to_ortools(quantum_assignment qa, FILE *output);


/**
 * @brief Prints a quantum assignment as a readable configuration to a file
 *
 * @param qa
 * @param output
 */
void quantum_assignment_print_to_file(quantum_assignment *qa, FILE *output);

/**
 * @brief Computes the negativite contexts of a quantum assignment
 * 
 * @param qa 
 */
bool quantum_assignment_compute_negativity(quantum_assignment* qa);

/**
 * @brief Conts the number of negative lines in a quantum assignment
 * 
 * @param qa quantum assignment
 * @return int 
 */
int negative_lines_count(quantum_assignment* qa);

/**
 * @brief Returns a quantum assignment where the contexts are those that are invalid
 * given a classical assignment
 * 
 * @param qa 
 * @param bool_sol classical assignment
 * @param validity if true the contexts are those that are valid
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_from_invalid_contexts(quantum_assignment qa,bool* bool_sol,bool validity);


// Function to determine the maximum number of columns
void parse_matrix_dimensions(FILE *file, size_t *rows, size_t *cols);

/**
 * @brief parses a file containing a quantum assignment
 * 
 * @param file 
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_parse(FILE *file);

/**
 * @brief merges two quantum assignments
 * 
 * @param qa1 
 * @param qa2
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_merge(quantum_assignment qa1,quantum_assignment qa2);

/**
 * @brief frees a quantum assignment
 * 
 * @param qa 
 */
void free_quantum_assignment(quantum_assignment* qa);

#endif //QUANTUM_ASSIGNMENT