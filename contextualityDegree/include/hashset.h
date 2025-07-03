/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
/**********************************************************************************/
/**
 * @file hashset.h
 * @brief Hashset of fixed-size bitsets
 */
#ifndef MY_DHSET
#define MY_DHSET

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define ELT_SIZE 5 //(1L<<2*(MAX(N_QUBITS-3,0)))//16 if 5qubits//4 if 4qubits ...

typedef uint64_t hash_set_store_type;

/**
 * @brief a 4-qubit doily contains 15 out of the (4^4)-1 = 255 possible observables
 * so since we look for permutations of such observables, we can store each 
 * observable of the permutation on one bit
 * thus the object used for this purpose is composed of 4 64-bit integers
 */
typedef struct {hash_set_store_type v[ELT_SIZE];} hash_set_bitset;
typedef struct {
    hash_set_bitset *list;
    size_t size;
    size_t max;
}hash_set;

/**
 * @brief initializes the hashset of permutations
 * 
 */
hash_set hash_set_init();

/**
 * @brief returns true iff a permutation is empty
 * 
 * @param permut 
 * @return true 
 * @return false 
 */
bool hash_set_bitset_is_null(hash_set_bitset permut);

/**
 * @brief returns true iff the 2 permutations are equal
 * 
 * @param p1 
 * @param p2 
 * @return true 
 * @return false 
 */
bool hash_set_bitset_equals(hash_set_bitset p1, hash_set_bitset p2);

bool hash_set_bitset_sup(hash_set_bitset p1,hash_set_bitset p2);

int hash_set_bitset_bitcount(hash_set_bitset p1);

/**
 * @brief creates an empty doily permutation
 * 
 * @return hash_set_bitset 
 */
hash_set_bitset hash_set_bitset_init();

/**
 * @brief adds an observable index in the permutation by toggleing the
 * right corresponding bit on one of the integers
 * 
 * @param permut 
 * @param number 
 */
void hash_set_bitset_add(hash_set_bitset *permut,unsigned int number);

/**
 * @return true iff the number is in the permut
 */
bool hash_set_bitset_get(hash_set_bitset permut,int number);

/**
 * @brief if the permutation doesn't exist, it is added to the hashset 
 * and returns true
 * else it returns false
 * 
 * @param permut 
 * @return true 
 * @return false 
 */
bool hash_set_exists(hash_set* list,hash_set_bitset permut);

void hash_set_free(hash_set* list);

#endif //MY_DHSET