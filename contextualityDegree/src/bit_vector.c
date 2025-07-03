/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file bitset.c
 * @brief A set of bits stored in an array of integers, and a matrix of bit sets
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "bit_vector.h"

#include "constants.h"

bit_vector bit_set_create(size_t size,bit_set_type* bits){
    bit_vector bs;
    bs.size = size;

    if(bits != NULL) {
        bs.bits = bits;
        return bs;
    }
    bs.bits = (bit_set_type*)calloc(BIT_SET_ARR_SIZE(bs), sizeof(bit_set_type));
    if (bs.bits == NULL) {
        fprintf(stderr, "Failed to allocate memory for bits [%ld]\n",size);
        return (bit_vector){0};
    }
    return bs;
}

void bit_set_free(bit_vector bs) {
    free(bs.bits);
}

void bit_set_set_bit(bit_vector bs, size_t index, bool value) {
    if (index >= bs.size) {
        print("Invalid bit_vector or index (set) [%ld,%ld]\n",index,bs.size);
        return;
    }

    size_t byteIndex = index / (sizeof(bit_set_type)*8);
    size_t bitIndex = index % (sizeof(bit_set_type)*8);

    if (value) {
        bs.bits[byteIndex] |= (1ULL << bitIndex);
    } else {
        bs.bits[byteIndex] &= ~(1ULL << bitIndex);
    }
}

bool bit_set_get_bit(bit_vector bs, size_t index) {
    if (index >= bs.size) {
        print("Invalid bit_vector or index (get) [%ld,%ld] ",index,bs.size);
        return false;
    }

    size_t byteIndex = index / (sizeof(bit_set_type)*8);
    size_t bitIndex = index % (sizeof(bit_set_type)*8);

    return (bs.bits[byteIndex] & (1ULL << bitIndex)) != 0;
}

size_t bit_set_cardinality(bit_vector bs) {
    size_t count = 0;
    for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs); i++) {
        bit_set_type bits = bs.bits[i];
        for (size_t j = 0; j < sizeof(bit_set_type)*8; j++) {
            count += BGET(bits, j);
        }
    }
    return count;
}

bool bit_set_is_empty(bit_vector bs) {
    for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs); i++)if (bs.bits[i] != 0)return false;
    return true;
}

bit_vector bit_set_op(bit_vector bs1, bit_vector bs2,char op,bit_set_type* bits) {
    if (bs1.size != bs2.size) {
        fprintf(stderr, "BitSets must be of the same size\n");
        return (bit_vector){0};
    }

    bit_vector result = bit_set_create(bs1.size, bits);
    
    switch (op) {
        case '&':
            for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs1); i++)result.bits[i] = bs1.bits[i] & bs2.bits[i];
            break;
        case '|':
            for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs1); i++)result.bits[i] = bs1.bits[i] | bs2.bits[i];
            break;
        case '^':
            for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs1); i++)result.bits[i] = bs1.bits[i] ^ bs2.bits[i];
            break;
        default:
            print("Invalid operator\n");
            bit_set_free(result);
            return (bit_vector){0};
    }

    return result;
}

bit_vector bit_set_apply_op(bit_vector* bs1, bit_vector bs2,char op) {
    return bit_set_op(*bs1, bs2, op, bs1->bits);
}

bit_vector copy_bitset(bit_vector bs,bit_set_type* bits) {
    bit_vector copy = bit_set_create(bs.size, bits);
    memcpy(copy.bits, bs.bits,BIT_SET_ARR_SIZE(bs)*sizeof(bit_set_type));
    return copy;
}

void print_bitset(bit_vector bs) {
    for (size_t i = 0; i < bs.size; i++) {
        print("%d", bit_set_get_bit(bs, i));
    }
    // print("(");
    // for (size_t i = 0; i < 64; i++) {
    //     print("%d", bit_set_get_bit(bs, i));
    // }
    // print(")");
}

// Provided bit_vector definition and utility functions here...

int bit_set_left_most(bit_vector bs) {
    for (long i = bs.size-1; i >= 0; i--) {
        if (bit_set_get_bit(bs, i)) {
            return i;
        }
    }
    return -1; // No set bit found
}

bit_matrix bit_matrix_create(size_t size, size_t col_size) {

    bit_matrix bm = {
        .size = size,
        .bit_sets = (bit_vector *)calloc(size, sizeof(bit_vector)),
        .bits = (bit_set_type **)init_matrix(size, (size_t)BIT_SIZE(col_size), sizeof(bit_set_type))};

    for (size_t i = 0; i < size; i++)bm.bit_sets[i] = bit_set_create(col_size, bm.bits[i]);

    return bm;
}

void bit_matrix_set_bit(bit_matrix bm, size_t row, size_t col, bool value) {
    bit_set_set_bit(bm.bit_sets[row], col, value);
}

bool bit_matrix_get_bit(bit_matrix bm, size_t row, size_t col) {
    return bit_set_get_bit(bm.bit_sets[row], col);
}

bit_matrix copy_bit_matrix(bit_matrix bm) {
    bit_matrix copy = bit_matrix_create(bm.size, bm.bit_sets[0].size);
    for (size_t i = 0; i < bm.size; i++)copy.bit_sets[i] = copy_bitset(bm.bit_sets[i], copy.bits[i]);
    return copy;
}

void bit_matrix_print(bit_matrix bm) {
    for (size_t i = 0; i < bm.size; i++) {
        print_bitset(bm.bit_sets[i]);
        //print(" %d",bit_set_left_most(bm.bit_sets[i]));
        print("\n");
    }
}

bool bit_matrix_equals(bit_matrix bm1, bit_matrix bm2) {
    if (bm1.size != bm2.size){
        print("mismatch row size : %ld != %ld\n",bm1.size,bm2.size);
        return false;
    }
    
    for (size_t i = 0; i < bm1.size; i++) {
        if (bm1.bit_sets[i].size != bm2.bit_sets[i].size){
            print("mismatch col size at %ld : %ld != %ld\n",i,bm1.bit_sets[i].size,bm2.bit_sets[i].size);
            return false;
        }
        
        for (size_t j = 0; j < bm1.bit_sets[i].size; j++) {
            if (bit_set_get_bit(bm1.bit_sets[i], j) != bit_set_get_bit(bm2.bit_sets[i], j)){
                print("mismatch at %ld,%ld : %d != %d\n",i,j,bit_set_get_bit(bm1.bit_sets[i], j),bit_set_get_bit(bm2.bit_sets[i], j));
                return false;
            }
        }
    }
    return true;
}

bit_matrix bit_matrix_product(bit_matrix bm1, bit_matrix bm2) {
    if (bm1.bit_sets[0].size != bm2.size) {
        print("BitMatrix column size must match row size of second matrix: %ld != %ld\n", bm1.bit_sets[0].size, bm2.size);
        return (bit_matrix){0};
    }

    bit_matrix result = bit_matrix_create(bm1.size, bm2.bit_sets[0].size);

    for (size_t i = 0; i < bm1.size; i++) {
        for (size_t j = 0; j < bm2.bit_sets[0].size; j++) {
            bool value = false;
            for (size_t k = 0; k < bm1.bit_sets[0].size; k++) {
                value ^= bit_matrix_get_bit(bm1, i, k) & bit_matrix_get_bit(bm2, k, j);
            }
            bit_matrix_set_bit(result, i, j, value);
        }
    }

    return result;
}

bool bit_matrix_is_empty(bit_matrix bm) {
    for (size_t i = 0; i < bm.size; i++) {
        if (!bit_set_is_empty(bm.bit_sets[i])) {
            return false;
        }
    }
    return true;
}

void bit_matrix_free(bit_matrix bm) {

    /*No bitset_free because the data is in bits*/
    free(bm.bit_sets);
    free_matrix(bm.bits);
}