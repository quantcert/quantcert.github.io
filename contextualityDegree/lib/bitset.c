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

#ifndef BITSET_H
#define BITSET_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#define BIT_SIZE(b) (((b) / (sizeof(bit_set_type)*8))+1)
#define BIT_SET_ARR_SIZE(bs) (BIT_SIZE(bs.size))

/** type of the bit set */
typedef uint64_t bit_set_type;

/**
 * @brief a set of bits stored in an array of integers
 * 
 * @param size Size of the bit set
 * @param bits memory adress to use for the bit set
 * @return bit_set New bit set
 */
typedef struct {
    bit_set_type* bits;
    size_t size;
} bit_set;

/**
 * @brief Creates a new bit set
 * 
 * @param size Size of the bit set
 * @param bits memory adress to use for the bit set
 * @return bit_set New bit set
 */
bit_set bit_set_create(size_t size,bit_set_type* bits){
    bit_set bs;
    bs.size = size;

    if(bits != NULL) {
        bs.bits = bits;
        return bs;
    }
    bs.bits = (bit_set_type*)calloc(BIT_SET_ARR_SIZE(bs), sizeof(bit_set_type));
    if (bs.bits == NULL) {
        fprintf(stderr, "Failed to allocate memory for bits [%ld]\n",size);
        return (bit_set){0};
    }
    return bs;
}

/**
 * @brief Frees the memory allocated for a bit set
 * 
 * @param bs Bit set to free
 */
void bit_set_free(bit_set bs) {
    free(bs.bits);
}

/**
 * @brief Sets a bit in a bit set
 * 
 * @param bs Bit set to modify
 * @param index Index of the bit to modify
 * @param value Value to set the bit to
 */
void bit_set_set_bit(bit_set bs, size_t index, bool value) {
    if (index >= bs.size) {
        print("Invalid bs or index (set) [%ld,%ld]\n",index,bs.size);
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

/**
 * @brief Returns the value of a bit in a bit set
 * 
 * @param bs Bit set to examine
 * @param index Index of the bit to examine
 * @return true If the bit is set
 * @return false If the bit is not set
 */
bool bit_set_get_bit(bit_set bs, size_t index) {
    if (index >= bs.size) {
        print("Invalid bs or index (get) [%ld,%ld] ",index,bs.size);
        return false;
    }

    size_t byteIndex = index / (sizeof(bit_set_type)*8);
    size_t bitIndex = index % (sizeof(bit_set_type)*8);

    return (bs.bits[byteIndex] & (1ULL << bitIndex)) != 0;
}

/**
 * @brief Returns the number of set bits in a bit set
 *
 * @param bs Bit set to examine
 * @return size_t number of set bits in a bit set
 */
size_t bit_set_cardinality(bit_set bs) {
    size_t count = 0;
    for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs); i++) {
        bit_set_type bits = bs.bits[i];
        for (size_t j = 0; j < sizeof(bit_set_type)*8; j++) {
            count += BGET(bits, j);
        }
    }
    return count;
}

/**
 * @brief Checks if a bit set is empty
 * 
 * @param bs Bit set to check
 * @return true If the bit set is empty
 * @return false If the bit set is not empty
 */
bool bit_set_is_empty(bit_set bs) {
    for (size_t i = 0; i < BIT_SET_ARR_SIZE(bs); i++)if (bs.bits[i] != 0)return false;
    return true;
}

/**
 * @brief Performs a bitwise operation on two bit sets
 * 
 * @param bs1 First bit set
 * @param bs2 Second bit set
 * @param op Operator ('&', '|', '^')
 * @return bit_set Result of the operation
 */
bit_set bit_set_op(bit_set bs1, bit_set bs2,char op,bit_set_type* bits) {
    if (bs1.size != bs2.size) {
        fprintf(stderr, "BitSets must be of the same size\n");
        return (bit_set){0};
    }

    bit_set result = bit_set_create(bs1.size, bits);
    
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
            return (bit_set){0};
    }

    return result;
}

/**
 * @brief Applies a bitwise operation of two bit sets to the first bit set
 * 
 * @param bs1 First bit set where the result is stored
 * @param bs2 Second bit set
 * @param op Operator ('&', '|', '^')
 * @return bit_set Result of the operation
 */
bit_set bit_set_apply_op(bit_set* bs1, bit_set bs2,char op) {
    return bit_set_op(*bs1, bs2, op, bs1->bits);
}

/**
 * @brief Copies a bit set
 * 
 * @param bs Bit set to copy
 * @return bit_set Copy of the bit set
 */
bit_set copy_bitset(bit_set bs,bit_set_type* bits) {
    bit_set copy = bit_set_create(bs.size, bits);
    memcpy(copy.bits, bs.bits,BIT_SET_ARR_SIZE(bs)*sizeof(bit_set_type));
    return copy;
}

/**
 * @brief Prints a bit set
 * 
 * @param bs Bit set to print
 */
void print_bitset(bit_set bs) {
    for (size_t i = 0; i < bs.size; i++) {
        print("%d", bit_set_get_bit(bs, i));
    }
    // print("(");
    // for (size_t i = 0; i < 64; i++) {
    //     print("%d", bit_set_get_bit(bs, i));
    // }
    // print(")");
}

// Provided bit_set definition and utility functions here...

/**
 * @brief Returns the index of the leftmost set bit, or -1 if it doesn't exist
 * 
 * @param bs bit_set to examine
 * @return int Index of the leftmost set bit, or -1 if none
 */
int bit_set_left_most(bit_set bs) {
    for (long i = bs.size-1; i >= 0; i--) {
        if (bit_set_get_bit(bs, i)) {
            return i;
        }
    }
    return -1; // No set bit found
}

/**
 * @brief A matrix of bit sets
 * @param bit_sets Array of bit sets
 * @param size Number of bit sets (number of rows)
 * @param bits Array of bit set data
 */
typedef struct {
    bit_set* bit_sets;
    size_t size;
    bit_set_type** bits;
} bit_matrix;

/**
 * @brief Creates a new bit matrix
 * 
 * @param size Number of rows
 * @param bit_set_size number of columns
 * @return bit_matrix New bit matrix
 */
bit_matrix bit_matrix_create(size_t size, size_t col_size) {

    bit_matrix bm = {
        .size = size,
        .bit_sets = (bit_set *)calloc(size, sizeof(bit_set)),
        .bits = (bit_set_type **)init_matrix(size, (size_t)BIT_SIZE(col_size), sizeof(bit_set_type))};

    for (size_t i = 0; i < size; i++)bm.bit_sets[i] = bit_set_create(col_size, bm.bits[i]);

    return bm;
}

/**
 * @brief Sets a bit in a bit matrix
 * 
 * @param bm Bit matrix to modify
 * @param row Row of the bit to modify
 * @param col Column of the bit to modify
 * @param value Value to set the bit to
 */
void bit_matrix_set_bit(bit_matrix bm, size_t row, size_t col, bool value) {
    bit_set_set_bit(bm.bit_sets[row], col, value);
}

/**
 * @brief Returns the value of a bit in a bit matrix
 * 
 * @param bm Bit matrix to examine
 * @param row Row of the bit to examine
 * @param col Column of the bit to examine
 * @return true If the bit is set
 * @return false If the bit is not set
 */
bool bit_matrix_get_bit(bit_matrix bm, size_t row, size_t col) {
    return bit_set_get_bit(bm.bit_sets[row], col);
}

bit_matrix copy_bit_matrix(bit_matrix bm) {
    bit_matrix copy = bit_matrix_create(bm.size, bm.bit_sets[0].size);
    for (size_t i = 0; i < bm.size; i++)copy.bit_sets[i] = copy_bitset(bm.bit_sets[i], copy.bits[i]);
    return copy;
}

/**
 * @brief Prints a bit matrix
 * 
 * @param bm Bit matrix to print
 */
void bit_matrix_print(bit_matrix bm) {
    for (size_t i = 0; i < bm.size; i++) {
        print_bitset(bm.bit_sets[i]);
        //print(" %d",bit_set_left_most(bm.bit_sets[i]));
        print("\n");
    }
}

/**
 * @brief Checks if two bit matrices are equal
 * 
 * @param bm1 First bit matrix
 * @param bm2 Second bit matrix
 * @return true If the bit matrices are equal
 * @return false If the bit matrices are not equal
 */
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
        fprintf(stderr, "BitMatrix column size must match row size of second matrix\n");
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

/**
 * @brief Frees the memory allocated for a bit matrix
 * 
 * @param bm Bit matrix to free
 */
void bit_matrix_free(bit_matrix bm) {

    /*No bitset_free because the data is in bits*/
    free(bm.bit_sets);
    free_matrix(bm.bits);
}

#endif //BITSET_H