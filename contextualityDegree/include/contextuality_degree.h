/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file contextuality_degree.h
 * @brief contains methods to compute the contextuality degree of a quantum assignment
 */
#ifndef CDEGREE
#define CDEGREE 1

#include "bv.h"
#include "quantum_assignment.h"

#define DISABLED_PARAMETER -1.0f

typedef enum
{
    SAT_SOLVER,
    RETRIEVE_SOLUTION,
    INVALID_LINES_HEURISTIC_SOLVER,
} solver_mode;

extern solver_mode global_solver_mode;
extern size_t global_heuristic_iterations; //maximum number of iterations for the heuristic method
extern float global_heuristic_flip_probability;  // probability of choosing a random assignment in the heuristic method
extern float global_heuristic_threshold;    // threshold for the heuristic method


/**
 * @brief return true with probability p
 * 
 * @param p 
 * @return true 
 * @return false 
 */
bool rand_float(float p);

/**
 * @brief prints a bool array grouping every 5 bits as a character
 * 
 * @param arr 
 * @param size 
 */
void print_bool(bool* arr,size_t size);

/**
 * @brief parses a string representing a bool array grouping every 5 bits as a character
 * 
 * @param arr 
 * @param size 
 */
void parse_bool(bool* arr,size_t size);


/**
 * @brief returns the hamming distance between a quantum assignment and a given
 * classical solution in a boolean form
 * 
 * @param qa 
 * @param bool_sol boolean array containing a solution
 * @param output file descriptor where the computation is printed (NULL if no output)
 * @return int hamming distance between the quantum assignment and the solution
 */
int check_contextuality_solution(quantum_assignment* qa,bool* bool_sol,FILE* output);

/**
 * @brief returns the hamming distance between a quantum assignment and a given 
 * classical solution using a SAT solver
 * 
 * @param qa 
 * @param bc2cnf_file name of the file containing the BC clauses
 * @param sat_file name of the file containing the SAT clauses
 * @param output file descriptor where the computation is printed (NULL if no output)
 * @param ret_sol boolean array containing the solution
 * @return int hamming distance between the quantum assignment and the solution
 */
int compute_contextuality_solution(quantum_assignment* qa,char* bc2cnf_file,char* sat_file,FILE* output,bool* ret_sol);


/**
 * @brief Writes a geometry into a file to be read by the bc2cnf program.
 * Computes the negativity of the geometry (the sum of all the observable must be the identity modulo some phase)
 * 
 * @param geometry array of geometries
 * @param l index of the geometry we want to print
*/
void write_line(FILE *f,bv* geometry,bool negative,int points_per_geometry);

/**
 * @brief Computes the maximum of lines a point contains
 * 
 * @param qa 
 * @return int 
 */
int max_line_per_point(quantum_assignment* qa);

/**
 * @brief Computes for each observable the contexts it is present in
 * 
 * @param qa 
 * @return int** A matrix containing for each observable the contexts it is present in (end of a context is marked by -1) 
 */
int **compute_contexts_per_obs(quantum_assignment *qa, bool print_solution);

/**
 * @brief method of finding a minimal hamming distance that successively flips the 
 * value of values of classical assignments contained in the most invalid lines
 * 
 * (No context can contain more than once the same observable)
 * 
 * @param qa input quantum assignment
 * @param ret_sol best solution found
 * @return int minimal hamming distance found
 */
int geometry_contextuality_degree_max_invalid_heuristics(quantum_assignment* qa,bool print_solution,bool* ret_sol);


/**
 * @brief Returns the contextuality degree of a list of geometries
 * 
 * @param geometry_indices lists all the geometries of the geometry array to take into account to compute the contextuality degree
 * @param geometries list of all (not necessary) the geometries
 * @param cpt_geometries number of geometries checked
 * @param points_per_geometry number of observables in each geometry
 * @param n_qubits number of qubits per observable
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
*/
int geometry_SAT_contextuality_degree(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,bool* ret_sol);


/**
 * @brief Returns the contextuality degree of a list of geometries with a given method
 * 
 * @param qa 
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual(only for the SAT solver)
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
 * @param mode method used to compute the contextuality degree
 * @param bool_sol if not NULL, the solution is stored in this array
 * @return int minimal hamming distance found
 */
int geometry_contextuality_degree_custom(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,solver_mode mode,bool* bool_sol);
/**
 * @brief Returns the contextuality degree of a quantum assignment
 * 
 * @param qa 
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
 * @param bool_sol if not NULL, the solution is stored in this array
 * @return int minimal hamming distance found
 */
int geometry_contextuality_degree(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,bool* bool_sol);
#endif //CDEGREE