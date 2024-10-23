/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/

/**
 * @file complex.c
 * @brief used to perform operations on complex numbers
 * and 2*2 matrices of complex numbers
 */

#ifndef MY_COMPLEX
#define MY_COMPLEX 1

#include <stdio.h>

#include "constants.c"
#include "bv.c"

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

/*every complex number used by the pauli matrices*/


/*
 * every pauli matrix
 */
pauli_matrix PAULI_MATS[N_PAULI];

void init_complex(){
    W2_complex w_0 = ((W2_complex){0, 0});//0
    W2_complex wP1 = ((W2_complex){1, 0});//1
    W2_complex wM1 = ((W2_complex){-1,0});//-1
    W2_complex wPI = ((W2_complex){0, 1});//i
    W2_complex wMI = ((W2_complex){0,-1});//-i

    PAULI_MATS[I] = (pauli_matrix){//I
    wP1,w_0,
    w_0,wP1
    };
    PAULI_MATS[X] = (pauli_matrix){//Z
    wP1,w_0,
    w_0 ,wM1
    };
    PAULI_MATS[Y] = (pauli_matrix){//X
    w_0,wP1,
    wP1,w_0
    };
    PAULI_MATS[Z] = (pauli_matrix){//Y
    w_0,wMI,
    wPI,w_0
    };
}

/**
 * Returns the Pauli matrix corresponding to a given observable in bv form
 * @param o one-qubit observable representing a pauli gate
*/
pauli_matrix get_matrix(bv o){
    if(o > Y)print("This matrix does not exist!");
    return PAULI_MATS[o];
}

/**
 * Returns true iff the 2 numbers in parameter are equal
*/
bool eq(W2_complex a,W2_complex b){
    return a.real == b.real && a.imag == b.imag;
}

/**
 * Returns true iff the number given in parameter is 0
*/
bool is_zero(W2_complex a){
    return a.real == 0 && a.imag == 0;
}

/**
 * @brief adds 2 complex numbers 
 */
W2_complex c_add(W2_complex a,W2_complex b){
    return (W2_complex){a.real+b.real,a.imag+b.imag};
}

/**
 * @brief multiplies 2 complex numbers 
 */
W2_complex c_mul(W2_complex a,W2_complex b){
    int rr = a.real*b.real;
    int ri = a.real*b.imag;
    int ir = a.imag*b.real;
    int ii = a.imag*b.imag;
    //printf(":%d %d %d %d@ %d %d:",rr,ri,ir,ii,rr-ii,ri+ir);
    return (W2_complex){rr-ii,ri+ir};
}

/**
 * @brief multiplies two 2*2 complex matrices
 * @return pauli_matrix 
 */
pauli_matrix matrix_mult(pauli_matrix m1,pauli_matrix m2){
    W2_complex a = c_add(c_mul(m1.a,m2.a),c_mul(m1.c,m2.b));
    W2_complex b = c_add(c_mul(m1.b,m2.a),c_mul(m1.d,m2.b));
    W2_complex c = c_add(c_mul(m1.a,m2.c),c_mul(m1.c,m2.d));
    W2_complex d = c_add(c_mul(m1.b,m2.c),c_mul(m1.d,m2.d));
    
    return (pauli_matrix){
        a,b,
        c,d
    };
}



/**
 * @brief prints a complex number
 * 
 * @param n 
 */
void print_complex_to_file(W2_complex n,FILE* output){

    if(n.real != 0)fprintf(output,"%+d",n.real);
    else if(n.imag == 1){fprintf(output,"+i");return;}
    else if(n.imag == -1){fprintf(output,"-i");return;}
    else if(n.imag != 0)fprintf(output,"%+di",n.imag);
    else fprintf(output," 0");
}
void print_complex(W2_complex n){print_complex_to_file(n,stderr);}

/**
 * @brief prints a 2*2 matrix
 * 
 * @param m matrix to print
 */
void print_mat(pauli_matrix m){
    print("[");
    print_complex(m.a);
    print(",");
    print_complex(m.b);
    print("\n ");
    print_complex(m.c);
    print(",");
    print_complex(m.d);
    print("]\n\n");
}

/** 
 * @brief Returns the phase of an identity matrix
 * Prints an error if the given matrix is not the identity modulo some phase
 * 
*/
W2_complex get_id_matrix_phase(pauli_matrix mat){
    if(!eq(mat.a,mat.d) ||  !is_zero(mat.b) || !is_zero(mat.c) || is_zero(mat.a)){
        print("incorrect value!");
        is_done = true;
    }
    /*
    * One property of any matrix which is x times I is that
    * there is only one eigenvalue which is x
    */    
    return mat.a;
}

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
W2_complex phase(size_t x1,size_t x2,size_t x3){

    pauli_matrix res =  matrix_mult(PAULI_MATS[x1],
                        matrix_mult(PAULI_MATS[x2],
                                    PAULI_MATS[x3]));
    

   return get_id_matrix_phase(res);
}

#endif //MY_COMPLEX
