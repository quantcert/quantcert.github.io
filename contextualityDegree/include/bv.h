/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file bv.c
 * @brief Contains functions about the bitvector used to represent observables
 *
 * Some functions are based on [DHGMS22] (https://doi.org/10.1088/1751-8121/aca36f)
 * and [HDS22] (https://doi.org/10.1038/s41598-022-13079-3)
 *
 */

#ifndef MY_BV
#define MY_BV 1/*safe guard*/

#include <time.h>

#include "constants.h"

/*bit vector representation of an observable*/
typedef uint32_t bv;
/*even or odd part of an observable*/
typedef uint32_t word;

/* Binary representation of each 1-qubit observable */
#define I 0b00
#define X 0b01
#define Y 0b11
#define Z 0b10

/* In the bitvector 0bRSTU, RS is the Z part and TU is the X part. */
#define II 0b0000//0
#define IX 0b0001//1
#define XI 0b0010//2
#define XX 0b0011//3
#define IZ 0b0100//4
#define IY 0b0101//5
#define XZ 0b0110//6
#define XY 0b0111//7
#define ZI 0b1000//8
#define ZX 0b1001//9
#define YI 0b1010//10
#define YX 0b1011//11
#define ZZ 0b1100//12
#define ZY 0b1101//13
#define YZ 0b1110//14
#define YY 0b1111//15

#define BV_LIMIT_CUSTOM(N) ((bv)pow4(N))//represents the limit value (not included) for an observable bit vector
#define BV_LIMIT BV_LIMIT_CUSTOM(N_QUBITS)

//#define mult(a,b) ((a) & (b))//multiplies an observable by another
//#define add(a,b)  ((a) ^ (b))//applies an observable to the other

#define Qplus ^ //applies an observable to the other
#define Qtimes &

extern char *HEXAD[][DOILY_SIZE];

/**
 * @brief returns the number of bits set to one on a bit vector
 * 
 * @param bv1 
 * @return unsigned int 
 */
unsigned int popcnt(bv bv1);

/**
 * @brief returns the factorial of a given number
 * 
 * @param a 
 * @return int 
 */
int factorial(int a);
/**
 * @brief Get the even part of the observable
 * 
 * @param bv1 
 * @param n_qubits 
 * @return word 
 */
word get_Z(bv bv1,int n_qubits);
/**
 * @brief Get the odd part of the observable
 * 
 * @param bv1 
 * @param n_qubits 
 * @return word 
 */
word get_X(bv bv1,int n_qubits);

/**
 * @brief merge the x and z parts into an bv
 * 
 * @param bv1 
 * @param n_qubits number of qubits in the bit vector
 * @return bv
 */
bv to_index_custom(word z,word x,int n_qubits);

/**
 * @brief returns 1 iff the number of bits in the word
 * is odd
 * (compilation optimization removes automatically unnecessary ifs)
 * bit parity trick found in http://graphics.stanford.edu/~seander/bithacks.html#ParityLookupTable
 * @param n 
 * @param n_qubits 
 * @return int 
 */
int bit_parity(word n,const int n_qubits);

word innerProductVector(bv i1,bv i2,int n_qubits);

/**
 * @brief symplectic product of 2 observables
 * 
 * bitvector parity in O(log2(N))
*/
unsigned int innerProduct_custom(bv i1,bv i2,int n_qubits);
unsigned int innerProduct(bv i1,bv i2);

/**
 * @brief Base quadratic form of an observable
 *
 * Q_0(bv1) = x1x2 + x3x4 + ... + x_{N-1}x_N
 *
 * See [DHGMS22] Section 5.1
 *
 */
unsigned int baseQuadraticFormVector(bv i, int n_qubits);

/**
 * @brief Quadratic of bv1 given a base
 * 
 * Q_base(bv1) = Q_0(bv1) + <bv1,base>
 * 
 * See [DHGMS22] Section 5.1
 * 
 */
unsigned int quadraticForm_custom(bv base, bv bv1, int n_qubits);

unsigned int quadraticForm(bv base, bv bv1);

/**
 * @brief same as inner product but with 2 qubits (faster)
 * 
 * @param i1 
 * @param i2 
 * @return int 
 */
unsigned int innerProduct_W2(bv i1,bv i2);

/**
 * @brief Performs the transvection of the point q from the base p
 * 
 * T_p(q) = q + <p,q>p
 * 
 * See [HDS22] Section 1
 */
bv transvection(bv p,bv q,int n_qubits);

/**
 * @brief Returns true if the observable is symmetric i.e. if its number of Y's is even
 * 
 * @param bv1 
 * @return int 
 */
bool is_symmetric(bv bv1,int n_qubits);

/**
 * @brief Prints every bit of a bitvector
 * 
 * @param b 
 */
void print_bits_custom(bv b,int n_bits);

/**
 * @brief gives the number of Is in an observable of size N_QUBITS
 * 
 * @param bv1 
 * @return int 
 */
int n_I_custom(bv bv1,int n_qubits);
int n_I(bv bv1);

/**
 * @brief returns a 2bit long integer corresponding to the selected gate
 *                     
 * example : get_gate(IXZZ,1) = X
 *                     ^
*/
unsigned int get_gate(bv i,unsigned int gate,int n_qubits);

bv set_gate_custom(bv bv1,bv obs,int index,int n_qubits);

/**
 * @brief Sets the gate "obs" to the quantum gate "index" of the observable "bv1" 
 * 
 * @param bv1 
 * @param obs 
 * @param gate 
 * @return bv 
 */
bv set_gate(bv bv1,bv obs,int index);

/**
 * @brief Prints a bitvector in the standard output (used to write files)
 * 
 * @param bv1 
 * @param N_QUBITS 
 */
void print_BV_to_file(bv bv1,int n_qubits,FILE* output);

/**
 * @brief prints a bit vector with its pauli gates in the console
 * 
 * @param bv1 bit vector to print
 */
void print_BV_custom(bv bv1,int n_qubits);
void print_BV_W2(bv bv1);

/**
 * @brief gives for each character its corresponding qubit
 * 
 * @param gate 
 * @return unsigned int 
 */
unsigned int char_to_gate(char gate);

/**
 * @brief Transforms a string into a qubit bit vector
 * 
 * @param str observable in string form
 * @param n_qubits number of qubits in the observable
 * @return bv 
 */
bv str_to_bv_custom(char *str,int n_qubits);

/**
 * @brief tTransforms a string into a N_QUBITS qubit bit vector
 * 
 * @param str observable in string form
 * @return bv 
 */
bv str_to_bv(char *str);

bv extend_bv(bv bv1,size_t from_qubits,size_t to_qubits);

#endif //MY_BV
