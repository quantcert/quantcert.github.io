/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file config_checker.c
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

#include "config_checker.h"

#include "quantum_assignment.h"

#define CONFIG_MAX_DEG 20000000

bool (*mode_filters[])(quantum_assignment, int, int *, int *) = {
    mode_nothing_filter,
    mode_all_filter,
    mode_point_degree_filter,
    mode_observable_value_filter,
    mode_line_degree_filter,
    mode_symmetric_filter};

bool mode_nothing_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    //avoid warnings on all variables
    (void)qa;(void)i;(void)param;(void)sorted_tab;
    return false;
}
bool mode_all_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    (void)qa;(void)i;(void)param;(void)sorted_tab;
    return true;
}
bool mode_point_degree_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    (void)i;
    for (size_t j = 0; j < qa.points_per_geometry; j++)
        if (sorted_tab[j] == param[0])
            return true;
    return false;
}
bool mode_observable_value_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    (void)sorted_tab;
    for (size_t j = 0; j < qa.points_per_geometry; j++)
        if (qa.geometries[qa.geometry_indices[i]][j] == (bv)param[0])
            return true;
    return false;
}
bool mode_line_degree_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    (void)i;
    for (size_t j = 0; j < qa.points_per_geometry; j++)
        if (sorted_tab[j] != param[j])
            return false;
    return true;
}
bool mode_symmetric_filter(quantum_assignment qa,int i,int* param,int* sorted_tab){
    (void)param;(void)sorted_tab;
    for (size_t j = 0; j < qa.points_per_geometry; j++)
        if (!is_symmetric(qa.geometries[qa.geometry_indices[i]][j], qa.n_qubits))
            return false;
    return true;
}

void swap(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void compute_invalid_lines(quantum_assignment qa, bool *bool_sol, bool **invalid_lines, int **number_of_invalid_lines)
{
    /*true for each invalid line*/
    (*invalid_lines) = (bool *)calloc(qa.cpt_geometries, sizeof(bool));
    if (!(*invalid_lines)){
        print("ERROR: invalid lines allocation\n");
        return;
    }
    /*number of invalid lines per point*/
    (*number_of_invalid_lines) = (int *)calloc(BV_LIMIT_CUSTOM(qa.n_qubits), sizeof(int));
    if (!(*number_of_invalid_lines)){
        print("ERROR: invalid lines allocation\n");
        free(*invalid_lines);
        return;
    }

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        bool test = false;
        for (size_t j = 0; j < qa.points_per_geometry; j++)
        {
            test ^= bool_sol[qa.geometries[qa.geometry_indices[i]][j]];
        }
        (*invalid_lines)[i] = (test != qa.lines_negativity[i]);
        if (!(*invalid_lines)[i])
            continue;

        for (size_t j = 0; j < qa.points_per_geometry; j++)
        {
            (*number_of_invalid_lines)[qa.geometries[qa.geometry_indices[i]][j]]++;
        }
    }
}

int get_table_index(int **table, int *line, size_t nb_points_per_line)
{
    for (size_t i = 0; i < CONFIG_MAX_DEG; i++)
    {
        bool valid = true;
        bool all_zeroes = true;
        for (size_t j = 0; j < nb_points_per_line; j++)
        {
            if (table[i][j] != 0)
                all_zeroes = false; // we check at the same time if the entry is empty
            if (table[i][j] != line[j])
            {
                valid = false;
                if (!all_zeroes)
                    break;
            }
        }
        if (valid)
            return i; // if found return index
        if (all_zeroes)
        { // else we create that entry
            for (size_t j = 0; j < nb_points_per_line; j++)
            {
                table[i][j] = line[j];
            }
            return i;
        }
    }
    print("\n\nERROR: table overflow\n");

    return -1;
}

int compare_integers_increasing(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

void print_category_list(int **category_list, int *category_count, quantum_assignment qa)
{
    print("\n\nLine degree types:\n");
    
    for (int i = 0; i < CONFIG_MAX_DEG; i++)
    {
        bool all_zeroes = true;
        
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            if (category_list[i][j] != 0)all_zeroes = false;
        }
        if (all_zeroes)break;         

        print("[");
        for (size_t j = 0; j < qa.points_per_geometry; j++)
        {
            print("%d", category_list[i][j]);
            if (j < qa.points_per_geometry - 1)
                print(",");
        }
        print("] (x%d)\n", category_count[i]);
    }
}

mode ask_mode(quantum_assignment qa, int *param)
{
    int mode = MODE_NOTHING;
    int nul = 0;
    print("\n\nType search mode (0->stop,1->all,2->point degree,3->observable value,4->line degree,5->symmetric): ");
    nul += scanf("%d", &mode);
    if (mode == MODE_POINT_DEGREE)
    {
        print("\npoint degree: ");
        nul += scanf("%d", &param[0]);
    }
    else if (mode == MODE_OBSERVABLE_VALUE)
    {
        char str[CONFIG_CHECKER_STR_BUFFER_SIZE];
        print("\nobservable value: ");
        nul += scanf("%s", str);
        param[0] = str_to_bv_custom(str, qa.n_qubits);
    }
    else if (mode == MODE_LINE_DEGREE)
    {
        for (size_t i = 0; i < qa.points_per_geometry; i++)
        {
            print("point %ld: ", i + 1);
            nul += scanf("%d", &param[i]);
        }
    }
    if (!nul)
        print("\noops, wrong text!\n");
    return mode;
}

void set_order(quantum_assignment qa, int *obs_specific_type_degree,int index,int *order)
{
    for (size_t i = 0; i < qa.points_per_geometry; i++)
        order[i] = i;

    /*ordering the points in each context (bubble sort)*/
    for (size_t j = 0; j < qa.points_per_geometry - 1; j++)
    {
        for (size_t k = 0; k < qa.points_per_geometry - j - 1; k++)
        {
            if (obs_specific_type_degree[qa.geometries[qa.geometry_indices[index]][order[k + 1]]] <
                obs_specific_type_degree[qa.geometries[qa.geometry_indices[index]][order[k]]])
            {
                swap(&order[k + 1], &order[k]);
            }
        }
    }
}

int filter_lines(bool to_print,quantum_assignment qa, bool *invalid_lines, int *number_of_invalid_lines, mode m, int *param, bool *line_filter, int **obs_specific_type_degree,int** category_count, int ***category_list)
{
    bool *is_line_in_specific_type = calloc(qa.cpt_geometries, sizeof(bool));

    *obs_specific_type_degree = calloc(BV_LIMIT_CUSTOM(qa.n_qubits), sizeof(int));
    if (!(*obs_specific_type_degree))
    {
        print("ERROR: obs_specific_type_degree allocation\n");
        return -1;
    }

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        if (!invalid_lines[i] || (line_filter != NULL && !line_filter[qa.geometry_indices[i]]))
            continue;

        int sorted_tab[qa.points_per_geometry];
        for (size_t j = 0; j < qa.points_per_geometry; j++)
        {
            sorted_tab[j] = number_of_invalid_lines[qa.geometries[qa.geometry_indices[i]][j]];
            if (sorted_tab[j] >= CONFIG_MAX_DEG)
                print("\noverflow error!(%d)\n", sorted_tab[j]);
        }
        qsort(sorted_tab, qa.points_per_geometry, sizeof(int), compare_integers_increasing);

        if (!mode_filters[m](qa,i,param,sorted_tab))continue;

        is_line_in_specific_type[i] = true;
        for (size_t j = 0; j < qa.points_per_geometry; j++)
        {
            bv bv1 = qa.geometries[qa.geometry_indices[i]][j];
            //if(m == MODE_POINT_DEGREE && number_of_invalid_lines[bv1] != param[0])continue;
            (*obs_specific_type_degree)[bv1]++;
        }
    }

    if (to_print){
        print("\nOn the left: degree of the points relative to the selected filter");
        print("\nOn the right: degree of the points relative to the whole invalid configuration\n");
        print("\n\nspecific / obs / general");
    }
    /*number of lines with that category*/
    *category_count = calloc(CONFIG_MAX_DEG, sizeof(int));
    /*for each type of invalid line, determined by the degrees of its vertices, we count the
    number of lines with that type*/
    *category_list = (int **)init_matrix(CONFIG_MAX_DEG, qa.points_per_geometry, sizeof(int));

    
    

    int max_index = 0;

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        int order[qa.points_per_geometry];
        set_order(qa, *obs_specific_type_degree,i, order);
        // for (size_t i = 0; i < qa.points_per_geometry; i++)
        //     order[i] = i;

        // /*ordering the points in each context (bubble sort)*/
        // for (size_t j = 0; j < qa.points_per_geometry - 1; j++){
        //     for (size_t k = 0; k < qa.points_per_geometry - j - 1; k++){
        //         if ((*obs_specific_type_degree)[qa.geometries[qa.geometry_indices[i]][order[k + 1]]] <
        //             (*obs_specific_type_degree)[qa.geometries[qa.geometry_indices[i]][order[k]]]){
        //             swap(&order[k + 1], &order[k]);
        //         }
        //     }
        // }

        if (!is_line_in_specific_type[i])continue;

        int line[qa.points_per_geometry];
        for (size_t j = 0; j < qa.points_per_geometry; j++)
            line[j] = (*obs_specific_type_degree)[qa.geometries[qa.geometry_indices[i]][order[j]]];
        int index = get_table_index(*category_list, line, qa.points_per_geometry);
        (*category_count)[index]++;
        if (index > max_index)max_index = index;

        if (!to_print)continue;

        print("\n");
        for (size_t j = 0; j < qa.points_per_geometry; j++)
            print("%d ", (*obs_specific_type_degree)[qa.geometries[qa.geometry_indices[i]][order[j]]]);
        for (size_t j = 0; j < qa.points_per_geometry; j++)
            print_BV_custom(qa.geometries[qa.geometry_indices[i]][order[j]], qa.n_qubits);
        for (size_t j = 0; j < qa.points_per_geometry; j++)
            print(" %d", number_of_invalid_lines[qa.geometries[qa.geometry_indices[i]][order[j]]]);
    }
    return max_index+1;
}

int check_structure(quantum_assignment* qa, bool *bool_sol, bool to_print, bool *line_filter)
{

    quantum_assignment_compute_negativity(qa);

    
    
    int ccpt = 0;
    if (to_print && line_filter != NULL)
        for (size_t i = 0; i < (size_t)NB_LINES_CUSTOM(qa->n_qubits); i++)
            ccpt++;

    

    bool *invalid_lines = NULL;
    int *number_of_invalid_lines = NULL;
    compute_invalid_lines(*qa, bool_sol, &invalid_lines, &number_of_invalid_lines);

    mode m = 1;
    int param[qa->points_per_geometry];
    if (to_print && global_interact_with_user)m = ask_mode(*qa, param);

    int *obs_specific_type_degree = NULL;
    int *category_count = NULL;
    int **category_list = NULL;

    int complexity_degree = filter_lines(to_print, *qa, invalid_lines, number_of_invalid_lines, m, param, line_filter, &obs_specific_type_degree, &category_count, &category_list);
    //second_filter(to_print, obs_specific_type_degree, is_line_in_specific_type, number_of_invalid_lines, qa, );
    //print("[[%d]]\n", category_list[0][0]);

    /*for each number of invalid lines, we count the number of points with that number of invalid lines*/
    int *number_of_obs_with_degree = calloc(qa->cpt_geometries, sizeof(int));

    for (bv i = 1; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)
        number_of_obs_with_degree[number_of_invalid_lines[i]]++;

    if (to_print)print("\nPoint degree types:\n");

    for (size_t i = 1; i < qa->cpt_geometries; i++)
    {
        if (to_print && number_of_obs_with_degree[i] != 0)
            print("degree %ld (x%d);", i, number_of_obs_with_degree[i]);
    }

    if (to_print)print("\n\nspecific degrees : ");

    /*for each present point degree, we count the number of points with that degree*/
    int *cpt_specific_degrees = calloc(qa->cpt_geometries, sizeof(int));

    for (bv j = 0; j < BV_LIMIT_CUSTOM(qa->n_qubits); j++)
        cpt_specific_degrees[obs_specific_type_degree[j]]++;

    if (to_print)
        for (size_t i = 1; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)
            if (cpt_specific_degrees[i] != 0)
                print("%ld(x%d)", i, cpt_specific_degrees[i]);

    int vertices_count = 0;
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)
        if (number_of_invalid_lines[i] > 0)
            vertices_count++;

    if (to_print)
        print("\n\nNumber of vertices in the invalid configuration: %d ; \n", vertices_count);

    if(to_print)print_category_list(category_list, category_count, *qa);

    if (to_print
    )print("Number of different line types: %d\n", complexity_degree);

    free_matrix(category_list);
    free(category_count);
    free(number_of_invalid_lines);
    free(invalid_lines);
    free(obs_specific_type_degree);
    free(number_of_obs_with_degree);
    free(cpt_specific_degrees);

    return complexity_degree;
}