/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file constants.c
 * @brief contains all the constants used in the program
 */

#ifndef MY_CONSTS
#define MY_CONSTS 1/*safe guard*/

#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <omp.h>

#ifndef N_QUBITS
    #define N_QUBITS 4//number of qubits in an observable
#endif
#ifndef SYMMETRIES
    #define SYMMETRIES false //if true only looks for doilies modulo known symmetries

#endif
#ifndef SKIP_OVOID_SYMMETRIES
    #define SKIP_OVOID_SYMMETRIES false//skips known symmetries that generate the same signatures
// but gives every sub_signature but with wrong numbers(much faster)
#endif
#ifndef HASHABLE_SRC
    #define HASHABLE_SRC "signature.c"
#endif
#ifndef N_CONFIGS
    #define N_CONFIGS 12//number of possible configurations per signature
#endif

#define DOILY_PROCESSING (N_CONFIGS == 12)//true iff doilies are processed

/////////////////////EDITABLE VARIABLES/////////////////////////////////////

#define IMPLEMENT_SIGINT true //if true, generates a results table when receiving a SIGINT signal
#define PRINT_PROGRESSION (N_QUBITS>=5)//if true prints the progression in % of the program
#define PRINT_SIGNATURE_COUNT false//if true always prints the number of sub-signatures found

#define MULTI_THREAD true //enables multi-thread(if you don't want debug prints issues)
#define SUB_NEGATIVE_TYPES DOILY_PROCESSING//if true differenciates types 7A/7B and 8A/8B
#define AUTOMATIC_PROCEDURE false//if true uses an automatic function to trace the path of the procedure

#define HEURISTIC_NUM_THREADS (MULTI_THREAD?200:1)//number of threads to use for the heuristic method
#define HEURISTIC_MAX_ITERATIONS 100000//maximum number of iterations for the heuristic method

#define DOUBLE_CHECK false//if true, checks the validity of every doily generated (has always been valid anyway)
#define DANGEROUS_OPTIMIZATIONS false//unproved optimizations
#define SEE_GRAPH true//if true, allows the user to go through the graph of the contextuality degree of a geometry
#define DETERMINISTIC false//if true, randomness will always use the same seed (must be used with MULTI_THREAD turned off)

//the 3 following options work properly when SKIP_OVOID_SYMMETRIES is DISABLED (since it is used to find the correct numbers)
#define BUFFERIZED false //bufferizes prints to file to accelerate the program(default:false)
#define BURNSIDER false//if true uses the burnside lemma, otherwise the faster lexicographic order(default:false)
#define XZ_SYMMETRY SYMMETRIES //checks the symmetry between X and Z gates 
// /!\ burnside doesn't work with ZX symmetry

#define W_PERM      SYMMETRIES //removes qubits permutations for doilies (IXXI <=> XIXI)

enum program_code{
    CODE_DOILY,
    CODE_DOILY_DFA,
    CODE_HEXAD,
    CODE_TWO_SPREAD,
    CODE_FANO
};
#define PROGRAM_MODE CODE_DOILY_DFA//what program to use


///////////////////NON-EDITABLE VARIABLES//////////////////////////////////////////////////////

#define pow2(a)   ( 1 << (a))//2^a (useful since it is done ahead of compilation)
#define pow4(a)   (pow2(2*(a)))//4^a
#define mask(a)   (pow2(a)-1)//applies a bit mask of size a

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define print(args...) fprintf(stderr, ##args)/*prints to the error output*/

/* MATH */

#define NOILY_SIZE(N) (pow4((N))-1) //4^N pauli observables minus the identity
#define DOILY_SIZE  (NOILY_SIZE(2)+1) //number of observables in the W(3,2) doily (we don't use the index 0)
#define TROILY_SIZE NOILY_SIZE(3) //size of a W5 (troily)

#define NB_LINES_CUSTOM(N) (((pow4(N)-1)*(pow4(N-1)-1))/3) /* = (4^n-1)(4^(n-1)-1)/3, see [dHG+21] */
#define NB_OBSERVABLES_CUSTOM(N) (pow4(N)-1)/*number of observables for N qubits*/
#define NB_LINES_PER_POINT_CUSTOM(N) ((NB_LINES_CUSTOM(N)*NB_POINTS_PER_LINE)/NB_OBSERVABLES_CUSTOM(N))

#define NB_LINES            NB_LINES_CUSTOM(N_QUBITS)
#define NB_OBSERVABLES      NB_OBSERVABLES_CUSTOM(N_QUBITS)
#define NB_LINES_PER_POINT  NB_LINES_PER_POINT_CUSTOM(N_QUBITS)

#define NB_TWO_SPREADS_PER_DOILY OVOID_INDEXES_COUNT /*dual of the ovoid*/
#define SPREAD_SIZE OVOID_COUNT /*dual of the ovoid*/
#define NB_LINES_TWO_SPREAD (NB_LINES_DOILY-SPREAD_SIZE)/*a 2-spread is a doily without a spread*/

/* GF(2) constants */

#define OVOID_COUNT 5//number of observables in an ovoid
#define N_COMMUTING 6//number of observables commuting with one another in the W(3,2) doily
#define INDEX_DOILY_SIZE 2//number of qubits in the doily used for indexes

#define NB_POINTS_PER_LINE 3//number of points in a line of a geometry
#define NB_POINTS_PER_AFFINE 4//number of points in an affine space of a geometry
#define N_PAULI 4//number of pauli gates
#define TRIAD_COUNT 3
#define NB_LINES_DOILY  NB_LINES_CUSTOM(INDEX_DOILY_SIZE)//(DOILY_SIZE - 1)
#define OVOID_INDEXES_COUNT 6
#define RES_TAB_SIZE pow2(3*MIN(N_QUBITS+1,5))//size of the results table : either 2*(8^N) or 2*(8^6) as an upper bound since it is above 2*(15+N choose N)
#define PAIRS_SIZE 1000
#define MAX_THREADS 64
#define FIND_MODE false
#define FIND_PAIRS (false && (N_QUBITS == 4))
#define ONE_DOILY_PER_CONFIGURATION false//stops the production of doilies when the 12 configurations have been found
#define DOILY_CONTEXTUALITY false//checks the contextuality degree of doilies and its 6 2-spreads
#define MERMIN_PERES_SQUARE false//generates mermin-peres squares from doilies /!\ only works with MULTI_THREAD = false
#define CHECK_SUBSPACES_LINES_EVEN false//if true prints the number of negative lines per point in each subspace

#define DOILY_MAX_TWO_POW 4//4 bits(16 values) are necessary to represent the number of observables per type in a doily

#define HEPTAD_SIZE 7//size of a Conwell heptad
#define TROILY_MAX_TWO_POW 6//6 bits(64 values) are necessary to represent the number of observables per type in a doily
#define TROILIY_COPIES 288//number of copies of the same troily generated
#define INDEX_TROILY_SIZE 3//number of qubits per observable in a troily
#define FANO_PLANE_SIZE 7//number of observables in a fano plane
#define N_TROILY_FANO_PLANES ((pow2(1)+1)*(pow2(2)+1)*(pow2(3)+1))//See contextuality degree of quadrics

#define BGET(b, i) (((b)>>(i))&1 )//gets a boolean value from a bit vector integer
#define BSET(b, i) ( (b) |= (1<<(i)))//sets a boolean value on a bit vector integer

/*default functions when openmp is not set*/
#ifndef _OPENMP
    void omp_set_num_threads(int n){}
    int omp_get_num_threads(){return 1;}
    int omp_get_thread_num(){return 1;}
#endif
/*outer function headers to avoid warnings*/
void program_end();

__ssize_t getline(char ** restrict,  size_t * restrict,  FILE * restrict);

struct sigaction;

volatile sig_atomic_t is_done = false;//to be set to true when SIGINT is triggered

/**
 * @brief triggered when the sigint signal is activated
 * sets the ending flag
 */
void sigint() {
  print("\ninterruption : closing file\n");
  //fclose(stdout);
  is_done = true;
  //end_permut();
}

void implement_sigint(){
    /*generates a result array before interruption
    we redirect the SIGINT signal to a program terminating function*/
    struct sigaction action;
    memset(&action, 0, sizeof(action));
    action.sa_handler = sigint;
    sigaction(SIGINT, &action, NULL);
}

void** init_matrix(size_t dim1,size_t dim2,size_t size){
    void** res = calloc(dim1,sizeof(void*));
    void* sub_array = calloc(dim1*dim2,size);
    for (size_t i = 0; i < dim1; i++){
        res[i] = (sub_array+(i*dim2*size));/**size since it cannot know the type*/
    }
    if(res == NULL || sub_array == NULL)print("matrix allocation error : %p %p\n",res,sub_array);
    return res;
}

void free_matrix(void* mat){//for compilation purposes
    if(mat == NULL)return;
    free(((void**)mat)[0]);
    free(mat);
}

void main_header(){
    if (!MULTI_THREAD)omp_set_num_threads(1);
    if (BUFFERIZED)setvbuf(stdout, NULL, _IOFBF, 0);      //bufferizes the output to gain speed
    if(IMPLEMENT_SIGINT)implement_sigint();
}



#endif //MY_CONSTS
