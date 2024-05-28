/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2.1                                             */
/**********************************************************************************/
/**
 * @file qontextium.c
 * @brief Main file of the qontextium program
 */

#ifndef QONTEXTIUM
#define QONTEXTIUM

#include "quadrics.c"

#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool
    SET_HYPERBOLICS = false,
    SET_ELLIPTICS = false,
    SET_PERPSETS = false,
    SET_SUBSPACES = false,
    SET_AFFINE = false,
    SET_IMPORT = false;


int main(int argc, char **argv)
{

    main_header();

    print("You can stop the program at any moment by pressing Ctrl+C\n");

    init_complex();

    int VARQ = N_QUBITS;

    // print("%d %d %d\n",NB_HYPERBOLICS(VARQ),NB_OBS_PER_HYPERBOLIC(VARQ),NB_LINES_PER_HYPERBOLIC(VARQ));
    // print("%d %d %d\n",NB_ELLIPTICS(VARQ),NB_OBS_PER_ELLIPTIC(VARQ),NB_LINES_PER_ELLIPTIC(VARQ));

    print("generating lines...");
    size_t **lines_indices = NULL;
    quantum_assignment lines_qa = generate_total_lines(&lines_indices, VARQ);
    quantum_assignment import_qa = (quantum_assignment){0};
    print(" done !\n    ");

    int k_subspaces = -1;
    bool complement = false;

    for (int i = 1; i < argc; i++)
    {
        bool is_zerolocus = false;
        if (strcmp(argv[i], "--hyperbolics") == 0)
        {
            SET_HYPERBOLICS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--elliptics") == 0)
        {
            SET_ELLIPTICS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--perpsets") == 0)
        {
            SET_PERPSETS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--affine") == 0)
        {
            SET_AFFINE = true;
        }
        else if (strcmp(argv[i], "--subspaces") == 0)
        {
            SET_SUBSPACES = true;
            if (i < argc - 1)
                k_subspaces = strtod(argv[i + 1], NULL);
            if (k_subspaces == -1)
                k_subspaces = VARQ - 1;
            i++;
        }
        if (is_zerolocus)
        {
            if (i < argc - 1 && strcmp(argv[i + 1], "--complement") == 0)
            {
                complement = true;
                i++;
            }
        }

        if (strcmp(argv[i], "--solver") == 0)
        {
            print("solver choosed : ");
            i++;
            if (i < argc && strcmp(argv[i], "sat") == 0)
            {
                global_solver_mode = SAT_SOLVER;
                print("sat\n");
            }
            else if (i < argc && strcmp(argv[i], "retrieve") == 0)
            {
                global_solver_mode = RETRIEVE_SOLUTION;
                print("retrieve\n");
            }
            else
            {
                print("none\n");
            }
        }

        if (strcmp(argv[i], "--import") == 0)
        {
            i++;
            if (i < argc)
            {
                FILE *f = fopen(argv[i], "r");
                if (f == NULL)
                {
                    print("file not found\n");
                    return 0;
                }
                import_qa = quantum_assignment_parse(f);
                fclose(f);
                SET_IMPORT = true;
            }
        }
    }
    print("h:%d e:%d p:%d s:%d k:%d a:%d i:%d\n", SET_HYPERBOLICS, SET_ELLIPTICS, SET_PERPSETS, SET_SUBSPACES, k_subspaces, SET_AFFINE, SET_IMPORT);

    
    
    
    
    
    
    

    bool *my_bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

    // #pragma omp parallel for shared(is_done)
    for (bv i = I; i < (bv)BV_LIMIT_CUSTOM(VARQ); i++)
    {
        if (is_done)
            continue;

        bool *bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

        quantum_assignment qa = (quantum_assignment){0};

        if (SET_PERPSETS && i != I)
        {
            qa = perpset(i, lines_indices, VARQ, lines_qa.geometries, complement);
            print("per:%d;", geometry_contextuality_degree(qa, false, true, false, NULL));
        }
        if (SET_HYPERBOLICS || SET_ELLIPTICS)
        {
            qa = quadric(i, lines_indices, VARQ, lines_qa.geometries, complement);
            if ((qa.cpt_geometries == (size_t)NB_LINES_PER_HYPERBOLIC(qa.n_qubits) && SET_HYPERBOLICS) ||
                (qa.cpt_geometries == (size_t)NB_LINES_PER_ELLIPTIC(qa.n_qubits) && SET_ELLIPTICS))
            {
                // print_quantum_assignment(qa);

                // printf("n_neg:%d\n",negative_lines_count(qa));
                print("test_sol:%d\n", geometry_contextuality_degree(qa, false, true, false, my_bool_sol));
                

                
                
                
                

                
                
                
            }

            // break;
        }
        free_quantum_assignment(&qa);
        free(bool_sol);
    }

    // test_CK_noIs(lines_qa, my_bool_sol, 0);free(my_bool_sol);return 0;/*takes a YYYY quadric*/

    if (SET_SUBSPACES)
    {

        bool *bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

        int c_deg_generators = (k_subspaces == 1) ? geometry_contextuality_degree(lines_qa, false, true, false, bool_sol) : subspaces_contextuality_degree(VARQ, k_subspaces);
        print("\nsubspace(n=%d,k=%d): contextuality degree : %d\n", VARQ, k_subspaces, c_deg_generators);

        free(bool_sol);
    }
    if (SET_AFFINE)
    {
        quantum_assignment planes = subspaces(VARQ, 2);
        quantum_assignment affine = affine_planes(planes);
        print("n_neg:%d/%ld\n", negative_lines_count(affine), affine.cpt_geometries);
        print("affine:%d\n", geometry_contextuality_degree(affine, false, true, false, NULL));
    }
    if (SET_IMPORT)
    {
        print("imported configuration:\n");
        print_quantum_assignment(import_qa);
        print("contextuality degree :%d\n", geometry_contextuality_degree(import_qa, false, true, false, NULL));
    }

    free_quantum_assignment(&lines_qa);
    free_matrix(lines_indices);
    free(lines_qa.geometry_indices);

    return 0;
}

#endif //QONTEXTIUM