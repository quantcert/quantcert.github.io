/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file config_checker.h
 * @brief Contains functions to check the structure of the invalid lines of a quantum assignment
 * There are multiple filters you can use to visualize the contextual graphs :

NOTHING
- This filter always returns false, meaning no lines will pass through this filter. It effectively filters out all lines.
ALL
- This filter always returns true, meaning all lines will pass through this filter.
POINT DEGREE
- This filter checks if any of the points in the geometry (for the given line i) match a specific point degree specified by the user.
OBSERVABLE VALUE
- This filter checks if any of the points in the geometry for the line i match a specific observable value given by the user
LINE DEGREE
- This filter checks if the sorted degrees of the points in the geometry for each line exactly match the degrees specified by the user
SYMMETRIC
- This filter checks if all points in the geometry for line i are symmetric. An observable is symmetric if the number of Y's is even.
 *
 */
#ifndef CONFIG_CHECKER
#define CONFIG_CHECKER

#include "quantum_assignment.h"

#define CONFIG_CHECKER_STR_BUFFER_SIZE 4096

typedef enum
{
    MODE_NOTHING,
    MODE_ALL,
    MODE_POINT_DEGREE,
    MODE_OBSERVABLE_VALUE,
    MODE_LINE_DEGREE,
    MODE_SYMMETRIC
} mode;

bool mode_nothing_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);
bool mode_all_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);
bool mode_point_degree_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);
bool mode_observable_value_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);
bool mode_line_degree_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);
bool mode_symmetric_filter(quantum_assignment qa,int i,int* param,int* sorted_tab);

extern bool (*mode_filters[])(quantum_assignment,int,int*,int*);

void swap(int *a, int *b);

void compute_invalid_lines(quantum_assignment qa, bool *bool_sol, bool **invalid_lines, int **number_of_invalid_lines);

/**
 * @brief internal function to find the index of a line in a table
 *
 * @param table
 * @param line
 * @param capacity
 * @return int
 */
int get_table_index(int **table, int *line, size_t nb_points_per_line);

int compare_integers_increasing(const void *a, const void *b);

/**
 * @brief Counts the complexity degree of a set of
 * categories based on the provided category_list and category_count, iterating
 * through configurations up to CONFIG_MAX_DEG. If to_print is true, it outputs
 * the category details and their counts, using data from the quantum_assignment
 * structure to determine the number of points per geometry.
 */
void print_category_list(int **category_list, int *category_count, quantum_assignment qa);

mode ask_mode(quantum_assignment qa, int *param);

void set_order(quantum_assignment qa, int *obs_specific_type_degree,int index,int *order);

/**
 *  processes geometric data and filters invalid lines based on the specified mode, 
 *  It can also optionally prints detailed information about the filtered 
 * lines and their degrees relative to the selected filter and the overall invalid 
 * configuration.
 * 
 * @return int the complexity degree of the filtered lines(number of different types of invalid lines
 * in terms of their degrees)
 */
int filter_lines(bool to_print,quantum_assignment qa, bool *invalid_lines, int *number_of_invalid_lines, mode m, int *param, bool *line_filter, int **obs_specific_type_degree,int** category_count, int ***category_list);

/**
 * @brief Checks the geometric structure of the invalid lines of a quantum assignment
 *
 * @param qa
 * @param bool_sol given solution
 * @param to_print if true prints the results
 * @param line_filter if not NULL, only checks the lines that are true in the filter
 * @return int the number of different types of invalid lines
 */
int check_structure(quantum_assignment* qa, bool *bool_sol, bool to_print, bool *line_filter);

#endif //CONFIG_CHECKER