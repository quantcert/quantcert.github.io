/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/

/**
 * @file complex_int.c
 * @brief used to perform operations on complex numbers
 * and 2*2 matrices of complex numbers
 */
#include "complex_int.h"

/*
 * every pauli matrix
 */
#define W2_COMPLEX(r, i) ((W2_complex){r, i})
#define PAULI_MAT(a, b, c, d) ((pauli_matrix){a, b, c, d})

pauli_matrix PAULI_MATS[N_PAULI] = {
    // I matrix: [[1,0],[0,1]]
    [I] = {{1, 0}, {0, 0}, {0, 0}, {1, 0}},
    // X matrix: [[0,1],[1,0]]
    [X] = {{0, 0}, {1, 0}, {1, 0}, {0, 0}},
    // Y matrix: [[0,-i],[i,0]]
    [Y] = {{0, 0}, {0, -1}, {0, 1}, {0, 0}},
    // Z matrix: [[1,0],[0,-1]]
    [Z] = {{1, 0}, {0, 0}, {0, 0}, {-1, 0}}
};

// void init_complex(){
//     W2_complex w_0 = ((W2_complex){0, 0}); // 0
//     W2_complex wP1 = ((W2_complex){1, 0});//1
//     W2_complex wM1 = ((W2_complex){-1,0});//-1
//     W2_complex wPI = ((W2_complex){0, 1});//i
//     W2_complex wMI = ((W2_complex){0,-1});//-i

//     PAULI_MATS[I] = (pauli_matrix){//I
//     wP1,w_0,
//     w_0,wP1
//     };
//     PAULI_MATS[X] = (pauli_matrix){//Z
//     wP1,w_0,
//     w_0 ,wM1
//     };
//     PAULI_MATS[Y] = (pauli_matrix){//X
//     w_0,wP1,
//     wP1,w_0
//     };
//     PAULI_MATS[Z] = (pauli_matrix){//Y
//     w_0,wMI,
//     wPI,w_0
//     };
// }

pauli_matrix get_matrix(bv o){
    if(o > Y)print("This matrix does not exist!");
    return PAULI_MATS[o];
}

bool eq(W2_complex a,W2_complex b){
    return a.real == b.real && a.imag == b.imag;
}

bool is_zero(W2_complex a){
    return a.real == 0 && a.imag == 0;
}

W2_complex c_add(W2_complex a,W2_complex b){
    return (W2_complex){a.real+b.real,a.imag+b.imag};
}

W2_complex c_mul(W2_complex a,W2_complex b){
    int rr = a.real*b.real;
    int ri = a.real*b.imag;
    int ir = a.imag*b.real;
    int ii = a.imag*b.imag;
    
    return (W2_complex){rr-ii,ri+ir};
}

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

void print_complex_to_file(W2_complex n,FILE* output){

    if(n.real != 0)fprintf(output,"%+d",n.real);
    else if(n.imag == 1){fprintf(output,"+i");return;}
    else if(n.imag == -1){fprintf(output,"-i");return;}
    else if(n.imag != 0)fprintf(output,"%+di",n.imag);
    else fprintf(output," 0");
}
void print_complex(W2_complex n){print_complex_to_file(n,stderr);}

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

W2_complex phase(size_t x1,size_t x2,size_t x3){

    pauli_matrix res =  matrix_mult(PAULI_MATS[x1],
                        matrix_mult(PAULI_MATS[x2],
                                    PAULI_MATS[x3]));
    

   return get_id_matrix_phase(res);
}
