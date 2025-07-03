/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file bit_vector.h
 * @brief A set of bits stored in an array of integers, and a matrix of bit sets
 */

#ifndef BITSET_H
#define BITSET_H

#define BIT_SIZE(b) (((b) / (sizeof(bit_set_type)*8))+1)
#define BIT_SET_ARR_SIZE(bs) (BIT_SIZE(bs.size))

/** type of the bit set */
typedef uint64_t bit_set_type;

/**
 * @brief a set of bits stored in an array of integers
 * 
 * @param size Size of the bit set
 * @param bits memory adress to use for the bit set
 * @return bit_vector New bit set
 */
typedef struct {
    bit_set_type* bits;
    size_t size;
} bit_vector;

/**
 * @brief A matrix of bit sets
 * @param bit_sets Array of bit sets
 * @param size Number of bit sets (number of rows)
 * @param bits Array of bit set data
 */
typedef struct
{
    bit_vector *bit_sets;
    size_t size;
    bit_set_type **bits;
} bit_matrix;

/**
 * @brief Creates a new bit set
 * 
 * @param size Size of the bit set
 * @param bits memory adress to use for the bit set
 * @return bit_vector New bit set
 */
bit_vector bit_set_create(size_t size,bit_set_type* bits);

/**
 * @brief Frees the memory allocated for a bit set
 * 
 * @param bs Bit set to free
 */
void bit_set_free(bit_vector bs);

/**
 * @brief Sets a bit in a bit set
 * 
 * @param bs Bit set to modify
 * @param index Index of the bit to modify
 * @param value Value to set the bit to
 */
void bit_set_set_bit(bit_vector bs, size_t index, bool value);

/**
 * @brief Returns the value of a bit in a bit set
 * 
 * @param bs Bit set to examine
 * @param index Index of the bit to examine
 * @return true If the bit is set
 * @return false If the bit is not set
 */
bool bit_set_get_bit(bit_vector bs, size_t index);

/**
 * @brief Returns the number of set bits in a bit set
 *
 * @param bs Bit set to examine
 * @return size_t number of set bits in a bit set
 */
size_t bit_set_cardinality(bit_vector bs);

/**
 * @brief Checks if a bit set is empty
 * 
 * @param bs Bit set to check
 * @return true If the bit set is empty
 * @return false If the bit set is not empty
 */
bool bit_set_is_empty(bit_vector bs);

/**
 * @brief Performs a bitwise operation on two bit sets
 * 
 * @param bs1 First bit set
 * @param bs2 Second bit set
 * @param op Operator ('&', '|', '^')
 * @return bit_vector Result of the operation
 */
bit_vector bit_set_op(bit_vector bs1, bit_vector bs2,char op,bit_set_type* bits);

/**
 * @brief Applies a bitwise operation of two bit sets to the first bit set
 * 
 * @param bs1 First bit set where the result is stored
 * @param bs2 Second bit set
 * @param op Operator ('&', '|', '^')
 * @return bit_vector Result of the operation
 */
bit_vector bit_set_apply_op(bit_vector* bs1, bit_vector bs2,char op);
/**
 * @brief Copies a bit set
 * 
 * @param bs Bit set to copy
 * @return bit_vector Copy of the bit set
 */
bit_vector copy_bitset(bit_vector bs,bit_set_type* bits);

/**
 * @brief Prints a bit set
 * 
 * @param bs Bit set to print
 */
void print_bitset(bit_vector bs);

// Provided bit_vector definition and utility functions here...

/**
 * @brief Returns the index of the leftmost set bit, or -1 if it doesn't exist
 * 
 * @param bs bit_vector to examine
 * @return int Index of the leftmost set bit, or -1 if none
 */
int bit_set_left_most(bit_vector bs);

/**
 * @brief Creates a new bit matrix
 * 
 * @param size Number of rows
 * @param bit_set_size number of columns
 * @return bit_matrix New bit matrix
 */
bit_matrix bit_matrix_create(size_t size, size_t col_size);

/**
 * @brief Sets a bit in a bit matrix
 * 
 * @param bm Bit matrix to modify
 * @param row Row of the bit to modify
 * @param col Column of the bit to modify
 * @param value Value to set the bit to
 */
void bit_matrix_set_bit(bit_matrix bm, size_t row, size_t col, bool value);

/**
 * @brief Returns the value of a bit in a bit matrix
 * 
 * @param bm Bit matrix to examine
 * @param row Row of the bit to examine
 * @param col Column of the bit to examine
 * @return true If the bit is set
 * @return false If the bit is not set
 */
bool bit_matrix_get_bit(bit_matrix bm, size_t row, size_t col);

bit_matrix copy_bit_matrix(bit_matrix bm);

/**
 * @brief Prints a bit matrix
 * 
 * @param bm Bit matrix to print
 */
void bit_matrix_print(bit_matrix bm);

/**
 * @brief Checks if two bit matrices are equal
 * 
 * @param bm1 First bit matrix
 * @param bm2 Second bit matrix
 * @return true If the bit matrices are equal
 * @return false If the bit matrices are not equal
 */
bool bit_matrix_equals(bit_matrix bm1, bit_matrix bm2);

bit_matrix bit_matrix_product(bit_matrix bm1, bit_matrix bm2);

bool bit_matrix_is_empty(bit_matrix bm);

/**
 * @brief Frees the memory allocated for a bit matrix
 * 
 * @param bm Bit matrix to free
 */
void bit_matrix_free(bit_matrix bm);

#endif //BITSET_H