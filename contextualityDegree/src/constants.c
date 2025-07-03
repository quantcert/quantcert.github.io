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
#include "constants.h"

volatile sig_atomic_t is_done = false; // to be set to true when SIGINT is triggered
bool global_interact_with_user = true; // if true, the program will interact with the user

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
    if(dim1 == 0 || dim2 == 0 || size == 0)return NULL;
    void** res = dim1>0?calloc(dim1,sizeof(void*)):NULL;
    void* sub_array = calloc(dim1*dim2,size);
    for (size_t i = 0; i < dim1; i++){
        res[i] = (sub_array+(i*dim2*size));/**size since it cannot know the type*/
    }
    if(res == NULL || sub_array == NULL)print("matrix allocation error : %ld %ld\n",dim1,dim2);
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
    if (IMPLEMENT_SIGINT)implement_sigint();
}

static unsigned int x = 0;
#pragma omp threadprivate(x)

unsigned int fast_random(){
    if(x == 0)x = omp_get_thread_num() + (DETERMINISTIC ? (0) : (time(NULL)));
    
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}
