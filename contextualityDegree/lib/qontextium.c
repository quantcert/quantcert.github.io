/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
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

    SET_ALL_QUADRICS = false,

    SET_SUBSPACES = false,
    SET_AFFINE = false,
    SET_IMPORT = false;


int main(int argc, char **argv)
{

    print("   ____              __            __  _               \n"
          "  / __ \\____  ____  / /____  _  __/ /_(_)_  ______ ___ \n"
          " / / / / __ \\/ __ \\/ __/ _ \\| |/_/ __/ / / / / __ `__ \\\n"
          "/ /_/ / /_/ / / / / /_/  __/>  </ /_/ / /_/ / / / / / /\n"
          "\\___\\_\\____/_/ /_/\\__/\\___/_/|_|\\__/_/\\__,_/_/ /_/ /_/ \n"

          "Welcome to Qontextium!\n\n"
          "Qontextium is a program for estimating the contextuality degree \n"
          "of quantum configurations. Developed by Axel Muller and Alain Giorgetti\n"
          "at Université de Franche-Comté, CNRS, institut FEMTO-ST\n"
          "Happy analyzing quantum configurations with Qontextium!\n\n\n");

    main_header();

    print("You can stop the program at any moment by pressing Ctrl+C\n\n");

    init_complex();

    int VARQ = atoi(argv[1]);


    // print("%d %d %d\n",NB_HYPERBOLICS(VARQ),NB_OBS_PER_HYPERBOLIC(VARQ),NB_LINES_PER_HYPERBOLIC(VARQ));
    // print("%d %d %d\n",NB_ELLIPTICS(VARQ),NB_OBS_PER_ELLIPTIC(VARQ),NB_LINES_PER_ELLIPTIC(VARQ));

    print("generating lines...");
    size_t **lines_indices = NULL;
    quantum_assignment lines_qa = generate_total_lines(&lines_indices, VARQ);
    quantum_assignment import_qa = (quantum_assignment){0};
    print(" done !\n    ");

    int k_subspaces = -1;
    bool complement = false;

    for (int i = 2; i < argc; i++)
    {
        bool is_zerolocus = false;
        if (strcmp(argv[i], "--hyperbolic") == 0)
        {
            print("configuration :hyperbolic quadric\n");
            SET_HYPERBOLICS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--elliptic") == 0)
        {
            print("configuration :elliptic quadric\n");
            SET_ELLIPTICS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--perpset") == 0)
        {
            print("configuration :perpset\n");
            SET_PERPSETS = true;
            is_zerolocus = true;
        }
        else if (strcmp(argv[i], "--affine") == 0)
        {
            print("configuration :affine planes\n");
            SET_AFFINE = true;
        }
        else if (strcmp(argv[i], "--subspaces") == 0)
        {
            print("configuration :subspaces\n");
            SET_SUBSPACES = true;
            if (i < argc - 1)
                k_subspaces = strtod(argv[i + 1], NULL);
            if (k_subspaces <= 0 || k_subspaces >= VARQ){
                print("\nInvalid number of dimensions: %d\n", k_subspaces);
                return 0;
            }
            i++;
        }
        if (is_zerolocus)
        {
            if (i < argc - 1 && strcmp(argv[i + 1], "--complement") == 0)
            {
                print("\tcomplement\n");
                complement = true;
                i++;
            }
        }

        if (strcmp(argv[i], "--solver") == 0)
        {
            print("selected solver :");
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
            else if (i < argc && strcmp(argv[i], "heuristic") == 0)
            {
                global_solver_mode = INVALID_LINES_HEURISTIC_SOLVER;
                print("heuristic\n");
            }
            else
            {
                print("none\n");
            }
        }

        else if (strcmp(argv[i], "--import") == 0)
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
        else if (strcmp(argv[i], "--all") == 0)
        {
            SET_ALL_QUADRICS = true;
        }
    }
    

    
    
    
    
    
    
    

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
            print("best Hamming distance found:%d;", geometry_contextuality_degree(&qa, false, true, false, NULL));
            if (!SET_ALL_QUADRICS)break;
        }
        if (SET_HYPERBOLICS || SET_ELLIPTICS)
        {
            qa = quadric(i, lines_indices, VARQ, lines_qa.geometries, complement);
            if ((qa.cpt_geometries == (size_t)NB_LINES_PER_HYPERBOLIC(qa.n_qubits) && SET_HYPERBOLICS) ||
                (qa.cpt_geometries == (size_t)NB_LINES_PER_ELLIPTIC(qa.n_qubits) && SET_ELLIPTICS))
            {
                // print_quantum_assignment(qa);

                // printf("n_neg:%d\n",negative_lines_count(qa));
                print("best Hamming distance found: %d\n", geometry_contextuality_degree(&qa, false, true, false, my_bool_sol));
                

                
                
                
                

                if(!SET_ALL_QUADRICS)break;
                
                
            }

            // break;
        }
        free_quantum_assignment(&qa);
        free(bool_sol);
    }

    

    if (SET_SUBSPACES)
    {

        bool *bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

        int c_deg;

        if (k_subspaces == 1){
            c_deg = geometry_contextuality_degree(&lines_qa, false, true, false, bool_sol);
        }else{
            print("\nGenerating subspaces...(expected: at most %ld * %d)\n\n", NB_SUBSPACES(VARQ, k_subspaces), NB_OBS_PER_GENERATOR(k_subspaces + 1));

            quantum_assignment qa = subspaces(VARQ, k_subspaces);
            print("\n\nGeometries generated((%ld contexts,%d obs. per cont.,%d qubits,%ld neg.lines))\n\nnow checking contextuality...\n\n", _subspaces_current_index, NB_OBS_PER_GENERATOR(k_subspaces + 1), VARQ, _subspaces_n_negative);
            
            bool sol[_subspaces_current_index];
            c_deg = geometry_contextuality_degree(&qa, false, true, true, sol);
            free_matrix(qa.geometries);
        }
        print("\nsubspace(%d qubits,%d dimensions):Hamming distance found: %d\n", VARQ, k_subspaces, c_deg);
        free(bool_sol);
    }
    if (SET_AFFINE){
        quantum_assignment planes = subspaces(VARQ, 2);
        quantum_assignment affine = affine_planes(planes);
        print("n_neg:%d/%ld\n", negative_lines_count(&affine), affine.cpt_geometries);
        print("affine:%d\n", geometry_contextuality_degree(&affine, false, true, false, NULL));

        
        
    }
    if (SET_IMPORT){
        print("imported configuration:\n");
        print_quantum_assignment(&import_qa);
        print("Hamming distance found :%d\n", geometry_contextuality_degree(&import_qa, false, true, false, NULL));
        free_quantum_assignment(&import_qa);
        free_matrix(import_qa.geometries);
    }

    free_quantum_assignment(&lines_qa);
    free_matrix(lines_indices);
    free_matrix(lines_qa.geometries);
    free(my_bool_sol);

    return 0;
}

#endif //QONTEXTIUM