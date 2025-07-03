/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file quadrics.c
 * @brief C functions to generate various geometries, among which quadrics.
 */

#ifndef QUADRICS
#define QUADRICS

#include "bv.h"
#include "contextuality_degree.h"




/**All these formula are described in Contextuality degree of quadrics
in multi-qubit symplectic polar spaces*/
#define NB_HYPERBOLICS(N) ((pow4(N)+pow2(N))/2)
#define NB_ELLIPTICS(N) ((pow4(N)-pow2(N))/2)
#define NB_QUADRICS(N) (NB_HYPERBOLICS(N)+NB_ELLIPTICS(N))
#define NB_PERPSETS(N) (pow4(N)-1)
int NB_GENERATORS(int N);

/**
 * @brief Counts the number of totally isotropic subspaces
 * up to 8 qubits, before the numbers get too large 
 * (see QPL2023)
 * 
 * @param N number of qubits
 * @param K dimension of the subspace
 * @return size_t 
 */
size_t NB_SUBSPACES(int N,int K);

#define NB_OBS_PER_HYPERBOLIC(N) (NB_HYPERBOLICS(N)-1)
#define NB_OBS_PER_ELLIPTIC(N) (NB_ELLIPTICS(N)-1)
#define NB_OBS_PER_QUADRIC(N) (NB_OBS_PER_HYPERBOLIC(N)+NB_OBS_PER_ELLIPTIC(N))
#define NB_OBS_PER_PERPSET(N) ((pow4(N)/2)-1)
#define NB_OBS_PER_GENERATOR(N) (pow2(N)-1)

#define NB_LINES_PER_HYPERBOLIC(N) ((NB_OBS_PER_HYPERBOLIC(N))*(NB_OBS_PER_HYPERBOLIC(N-1))/3)
#define NB_LINES_PER_ELLIPTIC(N) ((NB_OBS_PER_ELLIPTIC(N)*(NB_OBS_PER_ELLIPTIC(N-1)))/3)
#define NB_LINES_PER_QUADRIC(N) (NB_LINES_PER_HYPERBOLIC(N)+NB_LINES_PER_ELLIPTIC(N))
#define NB_LINES_PER_PERPSET(N) (pow4((N)-1)-1)

#define NO_LINE (0)

/**all 6 2-spreads as lists of indexes of the W2 doily lines*/
extern size_t index_two_spreads[NB_TWO_SPREADS_PER_DOILY][NB_LINES_TWO_SPREAD];


/**
 * @brief inserts a line index in the incidence map of observables
 * 
 * @param obs 
 * @param index 
 */
void insert_line_index(bv obs, size_t index,size_t** lines_indices,int n_qubits);
/**
 * @brief generates all the lines for a given number of qubits
 * 
 * @param lines resulting array of lines of observables(The index 0 is NOT a line !)
 * @param lines_indices indicates all the lines each observable belongs to (all lists of lines are sorted)
 * @param n_qubits number of qubits of the geometry
*/
quantum_assignment generate_total_lines(size_t*** lines_indices,int n_qubits);
/**
 * @brief returns the index of the leftmost set bit, or -1 if it doesn't exist
 * 
 * @param obs 
 * @return int 
 */
int bv_left_most(bv obs);
/**
 * @brief Comparison function to sort a bv array in DESCENDING order
 */
int bv_anti_cmp(const void * first, const void * second );
/**
 * @brief computes the basis of a set of observables
 * using Gaussian elimination
 * 
 * @param obss set of observables, will be sorted as an output
 * @param size number of observables
 * @param basis (output) the basis
 * 
 * @return int cardinal of the basis
 */
int find_basis(bv *obss,int size,bv** basis);
/**
 * @brief Get the basis combination of an observable given that 
 * the basis is sorted in descending order
 * 
 * @param obs 
 * @param basis 
 * @param size 
 * @return uint64_t 
 */
uint64_t get_basis_combination(bv obs,bv *basis,int size);

/**
 * @brief Callback function used to get all generators from the commuting function
*/
bool generate_subspace(bv bv1[],int size);

/**
 * @brief internal function used to generate all mutually commuting elements
 * 
 * @param bv1 list of observable to generate
 * @param index current index to generate observables on
 * @param size size of the observable list
 * @param n_qubits number of qubits per observable
 * 
 * @return true iff the iteration has to stop
 */
bool commuting_rec(bv bv1[],uint32_t index,uint32_t size,int n_qubits,quantum_assignment* res);

/**
 * @brief Generates all possible sets of mutually commuting observables
 * 
 * @param size number of mutually commuting observables
 * @param n_qubits number of qubits per observable
 */
quantum_assignment commuting(uint32_t size,int n_qubits);
/**
 * @brief Tests the conjecture stating that any subspace with k >= 4 has for each 
 * of its points an even number of negative lines passing through it
 * 
 * @param n_qubits 
 * @param k 
 * @param bv1 
 */
void subspace_neg_lines_count(int n_qubits,int k,bv bv1[]);


void subspaces_rec(int n_qubits,int k,int l,bv bv1[]);
/**
 * @brief Generates all subspaces of a given dimension
 * 
 * /!\ The geometries need to be freed after use
 * 
 * @param n_qubits number of qubits
 * @param k dimension of the subspace
 * 
 * @return quantum_assignment 
 */
quantum_assignment subspaces(int n_qubits, int k);
quantum_assignment affine_planes(quantum_assignment planes);
/**
 * @brief function generating the lines following a given form from a source observable and a given form
 * (even though for perpsets an additionnal condition is that the source observable must belong to each line)
 * 
 * WARNING : the array geometries_indices is allocated on the heap and must be freed after use
 * 
 * @param obs source observable of the geometry
 * @param form function determining if a point belongs to a geometry
 * @param lines_indices array specifying for each point all the lines it belongs to
 * @param n_qubits number of qubits of the geometry
 * @param lines array of all the lines of n qubits
 * 
 * @param lines_res array of all the lines of the resulting geometry
*/
quantum_assignment zero_locus(bv obs, unsigned int (*form)(bv, bv, int),size_t** lines_indices,int n_qubits,bv** lines,size_t* lines_res,bool complement);
quantum_assignment perpset(bv obs,size_t** lines_indices,int n_qubits,bv** lines,bool complement);
quantum_assignment quadric(bv obs,size_t** lines_indices,int n_qubits,bv** lines,bool complement);

/**
 * @brief Generates a 2-spread by removing the given spread from the doily
 * 
 * @param w2_doily_lines 
 * @param lines_indices 
 */
void generate_two_spread(size_t* lines_indices);
/**
 * @brief Builds recursively a spread by checking every possible combination of 5 lines such that they 
 * don't contain a common point (and since 5*3 = 15 the 5 lines contain all the points of the doily)
 * 
 * @param w2_doily_lines 
 * @param incidence_list 
 * @param lines_indices 
 * @param depth 
 */
void doily_spreads_rec(bv w2_doily_lines[][NB_POINTS_PER_LINE],uint16_t incidence_list,size_t* lines_indices,int depth);
/**
 * @brief generates the 6 spreads of the W2 doily in order to build the corresponding 2-spreads
 * 
 * 
 * 
 */
void doily_spreads(bv W2_doily_lines[][NB_POINTS_PER_LINE]);

#endif //QUADRICS
