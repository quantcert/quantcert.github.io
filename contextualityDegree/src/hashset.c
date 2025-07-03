/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
/**********************************************************************************/
/**
 * @file hashset.c
 * @brief Hashset of fixed-size bitsets
 */

#include <stdlib.h>
#include <stdint.h>

#include "hashset.h"

#include "constants.h"

#define BUCKET_SIZE 9973//2339/*prime number, good ratio size/speed*/
//#define LIST_SIZE 1//1000
#define N_BIT_INT 64
#define PERMUT_N_BITS ELT_SIZE*N_BIT_INT
//16*64 = 4^5
//16*4  = 4^4

//size_t max = 0;/*unused, maximum size of a list*/

hash_set hash_set_init(){
    hash_set list;
    list.size = 1000;
    list.max = 0;
    list.list = calloc(BUCKET_SIZE*list.size,sizeof(hash_set_bitset));
    return list;
}

bool hash_set_bitset_is_null(hash_set_bitset permut)
{
    for (size_t i = 0; i < ELT_SIZE; i++)if(permut.v[i] != 0)return false;
    return true;
}

bool hash_set_bitset_equals(hash_set_bitset p1, hash_set_bitset p2){
    for (size_t i = 0; i < ELT_SIZE; i++)if(p1.v[i] != p2.v[i])return false;
    return true;
}

bool hash_set_bitset_sup(hash_set_bitset p1,hash_set_bitset p2){
    for (int i = ELT_SIZE-1; i >= 0; i--)if(p1.v[i] != p2.v[i])return (p1.v[i] > p2.v[i]);
    return false;
}

int hash_set_bitset_bitcount(hash_set_bitset p1)
{
    int res = 0;
    for (size_t i = 0; i < ELT_SIZE; i++)res += __builtin_popcount(p1.v[i])+__builtin_popcount(p1.v[i]>>32);
    return res;
}

hash_set_bitset hash_set_bitset_init(){
    return (hash_set_bitset){v:{0}};
}

void hash_set_bitset_add(hash_set_bitset *permut,unsigned int number){
    permut->v[number/N_BIT_INT] |= 1L<<(number%N_BIT_INT);
}

bool hash_set_bitset_get(hash_set_bitset permut,int number){
    return BGET(permut.v[number/N_BIT_INT],number%N_BIT_INT);
}

bool hash_set_exists(hash_set* list,hash_set_bitset permut){

    hash_set_store_type hash = 1;/*the hash of the permutation*/
    for (size_t i = 0; i < ELT_SIZE; i++){
        hash = (hash*((permut.v[i]%BUCKET_SIZE)+1))%BUCKET_SIZE;/*+1 because we dont want to multiply by 0*/
    }
    size_t size = list->size;
    int index = size*hash;/*each hash has a list of 2000 capacity*/

    for (size_t i = 0; i < size; i++)
    {
        hash_set_bitset test_permut = list->list[index+i];

        if(hash_set_bitset_is_null(test_permut)){/*if there is place, we store it*/
            list->list[index+i] = permut;
            return false;/*it was not found, we had to store it*/
        }
        if (hash_set_bitset_equals(permut, test_permut))
            return true;                  /*we found it*/
        if(i > list->max){list->max = i;} /*checks the max size of a list*/
    }
    print(" couac:%ld ",hash);/* if we exceed the capacity of a bucket list*/
    return false;
}

void hash_set_free(hash_set* list){
    if(list->max >= 2)print("\nmaximum bucket list size:%lu; ",list->max);/*prints the max size of a bucket list*/
    free(list->list);
}
