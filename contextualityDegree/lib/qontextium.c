/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file qontextium.c
 * @brief Main file of the qontextium program
 */

#ifndef QONTEXTIUM
#define QONTEXTIUM

#include "quadrics.c"
#include "hypergram.c"
#include "cayley_hexagon.c"

#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum
{
    EXPORT_ALL,
    EXPORT_VALID,
    EXPORT_INVALID
} EXPORT_TYPE;

bool
    SET_HYPERBOLICS = false,
    SET_ELLIPTICS = false,
    SET_PERPSETS = false,
    SET_HEXAGONS = false,
    SET_SKEW_HEXAGONS = false,

    SET_ALL_QUADRICS = false,

    SET_SUBSPACES = false,
    SET_AFFINE = false,
    SET_IMPORT_ASSIGNMENT = false,
    SET_IMPORT_HYPERGRAM = false,
    
    SET_EXPORT = false;

EXPORT_TYPE export_type = EXPORT_ALL;




/**
 * @brief wrapper for the geometry_contextuality_degree function
 * that prints the contexts depending on the user choice
 * 
 * @param qa 
 * @param contextuality_only 
 * @param print_solution 
 * @param optimistic 
 * @param bool_sol 
 * @return int 
 */
int geometry_contextuality_degree_and_print(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,bool* bool_sol){
    int deg = geometry_contextuality_degree(qa, contextuality_only, print_solution, optimistic, bool_sol);
    
    if(SET_EXPORT){
        if(export_type == EXPORT_ALL){
            quantum_assignment_to_CSV(*qa, stdout);
        }else{
            quantum_assignment line_invalids = quantum_assignment_from_invalid_contexts(*qa, bool_sol, export_type == EXPORT_VALID);
            quantum_assignment_to_CSV(line_invalids, stdout);
            free_quantum_assignment(&line_invalids);
        }
    }

    print("best Hamming distance found: %d\n", deg);

    return deg;
}


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

    int VARQ = -1;

    // print("%d %d %d\n",NB_HYPERBOLICS(VARQ),NB_OBS_PER_HYPERBOLIC(VARQ),NB_LINES_PER_HYPERBOLIC(VARQ));
    // print("%d %d %d\n",NB_ELLIPTICS(VARQ),NB_OBS_PER_ELLIPTIC(VARQ),NB_LINES_PER_ELLIPTIC(VARQ));

    size_t **lines_indices = NULL;
    quantum_assignment lines_qa = {0};//generate_total_lines(&lines_indices, VARQ);
    quantum_assignment import_qa = (quantum_assignment){0};
    print(" done !\n    ");

    int k_subspaces = -1;
    bool complement = false;

    for (int i = 1; i < argc; i++)
    {
        bool is_zerolocus = false;
        bool generate_lines = false;
        if (strcmp(argv[i], "--hyperbolic") == 0)
        {
            print("configuration :hyperbolic quadric\n");
            SET_HYPERBOLICS = true;
            is_zerolocus = true;
            generate_lines = true;
        }
        else if (strcmp(argv[i], "--elliptic") == 0)
        {
            print("configuration :elliptic quadric\n");
            SET_ELLIPTICS = true;
            is_zerolocus = true;
            generate_lines = true;
        }
        else if (strcmp(argv[i], "--perpset") == 0)
        {
            print("configuration :perpset\n");
            SET_PERPSETS = true;
            is_zerolocus = true;
            generate_lines = true;
        }
        else if (strcmp(argv[i], "--affine") == 0)
        {
            print("configuration :affine planes\n");
            SET_AFFINE = true;
            generate_lines = true;
        }
        else if (strcmp(argv[i], "--hexagon") == 0)
        {
            print("configuration :split cayley hexagons\nembedding: ");
            SET_HEXAGONS = true;
            VARQ = 3;
            i++;
            if (i < argc && strcmp(argv[i], "skew") == 0)
            {
                SET_SKEW_HEXAGONS = true;
                print("skew\n");
            }else{
                print("classical\n");
            }
        }
        else if (strcmp(argv[i], "--subspaces") == 0)
        {
            print("configuration :subspaces\n");
            SET_SUBSPACES = true;
            generate_lines = true;
            if (i < argc - 1)
                k_subspaces = strtod(argv[i + 1], NULL);
            if (k_subspaces <= 0){
                print("\nInvalid number of dimensions: %d\n", k_subspaces);
                return 0;
            }
            i++;
        }
        if(generate_lines){
            if (i < argc)
            {
                VARQ = atoi(argv[i + 1]);
                
                if(VARQ <= 0){
                    print("\nInvalid number of qubits: %d\n", VARQ);
                    return 0;
                }
                if (SET_SUBSPACES && k_subspaces >= VARQ)
                {
                    print("\nInvalid number of dimensions: %d >= %d\n", k_subspaces, VARQ);
                    return 0;
                }
                lines_qa = generate_total_lines(&lines_indices, VARQ);
                i++;
                print(" VARQ:%d\n", VARQ);
            }
        }
        if (is_zerolocus || SET_HEXAGONS)
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
                if(strcmp(argv[i],"assignment") == 0){
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
                        SET_IMPORT_ASSIGNMENT = true;
                    }else{
                        print("no file specified\n");
                    }
                }else if(strcmp(argv[i],"hypergram") == 0){
                    //two files needed
                    if (i + 2 < argc)
                    {
                        FILE *f_hypergrpaph = fopen(argv[i + 1], "r");
                        FILE *f_gram = fopen(argv[i + 2], "r");
                        if (f_hypergrpaph == NULL){
                            print("hypergrpaph file not found\n");
                            return 0;
                        }
                        if (f_gram == NULL){
                            print("gram file not found\n");
                            return 0;
                        }
                        hypergram import_hg = hypergram_create_from_file(f_hypergrpaph, f_gram);
                        if (import_hg.geometries == NULL)
                        {
                            print("error while importing hypergram\n");
                            return 0;
                        }
                        hypergram_compute_assignment(&import_hg);
                        print("\nnumber of qubits:%ld\n",import_hg.n_qubits);
                        import_qa = hypergram_to_quantum_assignment(import_hg);
                        hypergram_free(import_hg);
                        fclose(f_hypergrpaph);
                        fclose(f_gram);
                        SET_IMPORT_HYPERGRAM = true;
                    }
                    else
                    {
                        print("two files needed\n");
                    }
                }
                    
            }
        }
        else if (strcmp(argv[i], "--all") == 0)
        {
            SET_ALL_QUADRICS = true;
        }
        else if (strcmp(argv[i], "--export") == 0)
        {
            SET_EXPORT = true;
            i++;

            if (i < argc)
            {
                if (strcmp(argv[i], "all") == 0){
                    export_type = EXPORT_ALL;
                }
                else if (strcmp(argv[i], "valid") == 0){
                    export_type = EXPORT_VALID;
                }
                else if (strcmp(argv[i], "invalid") == 0)
                {
                    export_type = EXPORT_INVALID;
                }
            }
        }
        else if (strcmp(argv[i], "--no-interaction") == 0){
            global_interact_with_user = false;
        }
        else if (strcmp(argv[i], "--heuristic-iter") == 0){
            i++;
            if (i < argc){
                global_heuristic_iterations = atoi(argv[i]);
            }
        }
    }


    bool *my_bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

    // #pragma omp parallel for shared(is_done)
    for (bv i = I; i < (bv)BV_LIMIT_CUSTOM(VARQ); i++)
    {
        if (is_done)continue;

        bool *bool_sol = calloc(BV_LIMIT_CUSTOM(VARQ), sizeof(bool));

        quantum_assignment qa = (quantum_assignment){0};

        if (SET_PERPSETS && i != I)
        {
            qa = perpset(i, lines_indices, VARQ, lines_qa.geometries, complement);
            geometry_contextuality_degree_and_print(&qa, false, true, false, NULL);
            print_quantum_assignment(&qa);
            if (!SET_ALL_QUADRICS)break;
        }
        if (SET_HYPERBOLICS || SET_ELLIPTICS)
        {
            qa = quadric(i, lines_indices, VARQ, lines_qa.geometries, complement);
            size_t expected_size = complement?lines_qa.cpt_geometries - qa.cpt_geometries:qa.cpt_geometries;
            if ((expected_size == (size_t)NB_LINES_PER_HYPERBOLIC(qa.n_qubits) && SET_HYPERBOLICS) ||
                (expected_size == (size_t)NB_LINES_PER_ELLIPTIC(qa.n_qubits) && SET_ELLIPTICS))
            {
                //print_quantum_assignment(&qa);

                // printf("n_neg:%d\n",negative_lines_count(qa));
                geometry_contextuality_degree_and_print(&qa, false, true, false, my_bool_sol);
                
                
                
                
                
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

        print("\nsubspace(%d qubits,%d dimensions)\n", VARQ, k_subspaces);
        if (k_subspaces == 1){
            geometry_contextuality_degree_and_print(&lines_qa, false, true, false, bool_sol);
        }else{
            print("\nGenerating subspaces...(expected: at most %ld * %d)\n\n", NB_SUBSPACES(VARQ, k_subspaces), NB_OBS_PER_GENERATOR(k_subspaces + 1));

            quantum_assignment qa = subspaces(VARQ, k_subspaces);
            print("\n\nGeometries generated((%ld contexts,%d obs. per cont.,%d qubits,%ld neg.lines))\n\nnow checking contextuality...\n\n", _subspaces_current_index, NB_OBS_PER_GENERATOR(k_subspaces + 1), VARQ, _subspaces_n_negative);
            
            bool sol[_subspaces_current_index];
            geometry_contextuality_degree_and_print(&qa, false, true, true, sol);
            free_matrix(qa.geometries);
        }
        free(bool_sol);
        //print_quantum_assignment(&lines_qa);
    }
    if (SET_AFFINE){
        quantum_assignment planes = subspaces(VARQ, 2);
        quantum_assignment affine = affine_planes(planes);
        print("negative lines:%d/%ld\n", negative_lines_count(&affine), affine.cpt_geometries);
        geometry_contextuality_degree_and_print(&affine, false, true, false, NULL);

        
        
    }
    if(SET_HEXAGONS){
        VARQ = 3;
        lines_qa = generate_total_lines(&lines_indices, VARQ);
        quantum_assignment qa = {0}, complement_qa = {0};

        if(!SET_SKEW_HEXAGONS){
            while (next_classical_cayley_hexagon(lines_qa, lines_indices, &qa, &complement_qa) && !is_done){

                geometry_contextuality_degree_and_print(complement ? (&complement_qa) : (&qa), false, true, false, my_bool_sol);
                if(!SET_ALL_QUADRICS)break;
            }
        }else{
            while (next_skew_cayley_hexagon(lines_qa, lines_indices, &qa, &complement_qa) && !is_done){

                geometry_contextuality_degree_and_print(complement ? (&complement_qa) : (&qa), false, true, false, my_bool_sol);
                if (!SET_ALL_QUADRICS)break;
            }
        }
        free_quantum_assignment(&qa);
        free_quantum_assignment(&complement_qa);
        free_cayley_hexagons();
    }
    if (SET_IMPORT_ASSIGNMENT || SET_IMPORT_HYPERGRAM){
        print("imported configuration:\n");
        print_quantum_assignment(&import_qa);
        geometry_contextuality_degree_and_print(&import_qa, false, true, false, NULL);
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