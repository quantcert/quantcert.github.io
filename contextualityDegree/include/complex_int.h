/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/

/**
 * @file complex_int.h
 * @brief used to perform operations on complex numbers
 * and 2*2 matrices of complex numbers
 */

#ifndef MY_COMPLEX
#define MY_COMPLEX 1

#include <stdio.h>

#include "bv.h"

/**
 * @brief Discrete complex number
 * stored as two signed integers
 */
typedef struct{
    signed int real;
    signed int imag;
}W2_complex;

/**
 * @brief 2*2 matrix
 * 
 * used for pauli matrix multiplication
 */
typedef struct{
    W2_complex a;
    W2_complex b;
    W2_complex c;
    W2_complex d;
}pauli_matrix;

extern pauli_matrix PAULI_MATS[N_PAULI];

/*every complex number used by the pauli matrices*/

void init_complex();

/**
 * Returns the Pauli matrix corresponding to a given observable in bv form
 * @param o one-qubit observable representing a pauli gate
*/
pauli_matrix get_matrix(bv o);

/**
 * Returns true iff the 2 numbers in parameter are equal
*/
bool eq(W2_complex a,W2_complex b);
/**
 * Returns true iff the number given in parameter is 0
*/
bool is_zero(W2_complex a);

/**
 * @brief adds 2 complex numbers 
 */
W2_complex c_add(W2_complex a,W2_complex b);

/**
 * @brief multiplies 2 complex numbers 
 */
W2_complex c_mul(W2_complex a,W2_complex b);

/**
 * @brief multiplies two 2*2 complex matrices
 * @return pauli_matrix 
 */
pauli_matrix matrix_mult(pauli_matrix m1,pauli_matrix m2);



/**
 * @brief prints a complex number
 * 
 * @param n 
 */
void print_complex_to_file(W2_complex n,FILE* output);
void print_complex(W2_complex n);

/**
 * @brief prints a 2*2 matrix
 * 
 * @param m matrix to print
 */
void print_mat(pauli_matrix m);

/** 
 * @brief Returns the phase of an identity matrix
 * Prints an error if the given matrix is not the identity modulo some phase
 * 
*/
W2_complex get_id_matrix_phase(pauli_matrix mat);

/**
 * @brief gives the multiple of the identity matrix
 * given by the product of three pauli
 * matrices given by their index (I,X,Y or Z)
 * ONLY when the product is a multiple of the identity matrix
 * 
 * otherwise the result is not defined and an error is printed
 * 
 * @param x1 
 * @param x2 
 * @param x3 
 * @return W2_complex 
 */
W2_complex phase(size_t x1,size_t x2,size_t x3);

#endif //MY_COMPLEX
