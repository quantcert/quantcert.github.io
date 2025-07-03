/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
/**********************************************************************************/
/**
 * @file test.c
 * @author Axel Muller
*/
#ifndef TEST_C
#define TEST_C

#include "constants.h"
#include "quadrics.h"
#include "hypergram.h"
#include "cayley_hexagon.h"
#include "complex_int.h"

size_t n_passed = 0; //number of passed tests
size_t n_failed = 0; //number of failed tests

bool assert_true(bool condition, const char *message){
    if(condition){
        printf("\u2705 - pass:\t\t %s\n", message); 
        n_passed++; 
    }else{
        printf("\u274C - fail:\t\t %s\n", message);
        n_failed++;
    }
    return condition;
}

bool assert_equal(int a, int b, const char *message){
    if(a == b){
        printf("\u2705 - pass: %d == %d\t %s\n", a, b, message);
        n_passed++;
    }else{
        printf("\u274C - fail: %d != %d\t %s\n", a, b, message);
        n_failed++;
    }
    return a == b;
}

void print_summary(){
    printf("\n\nTests summary:\n");
    printf("Passed\u2705: %zu\n", n_passed);
    printf("Failed\u274C: %zu\n", n_failed);
    if(n_failed > 0){
        printf("Some tests failed, please check the output.\n");
    }else{
        printf("All tests passed successfully!\n");
    }
}

int main(){

    global_interact_with_user = false;
    global_heuristic_iterations = 1000;

    main_header();

    int VARQ = 3;

    size_t **lines_indices = NULL;
    quantum_assignment lines_qa_three = generate_total_lines(&lines_indices, VARQ);

    /////////////////////////////

    assert_equal(lines_qa_three.cpt_geometries,315, 
    "315 lines for 3 qubits");

    /////////////////////////////

    quantum_assignment doily = subspaces(2,1);
    int doily_deg = geometry_contextuality_degree_custom(&doily, false, false, false, SAT_SOLVER, NULL);
    assert_equal(doily_deg,3, 
    "Doily has contextuality degree 3");
    
    /////////////////////////////

    int heuristic_doily_deg = geometry_contextuality_degree_custom(&doily, false, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);
    assert_equal(heuristic_doily_deg,3,
    "Doily has contextuality degree 3 with heuristic");

    free_quantum_assignment(&doily);

    /////////////////////////////

    FILE *grid_file = fopen("./misc/qa_grid.txt", "r");
    assert_true(grid_file != NULL, 
    "grid file opened");

    /////////////////////////////

    quantum_assignment import_qa = quantum_assignment_parse(grid_file);
    assert_equal(import_qa.cpt_geometries,6, 
    "imported grid has 6 geometries");
    
    fclose(grid_file);

    /////////////////////////////

    bool* bool_sol = calloc(BV_LIMIT_CUSTOM(import_qa.n_qubits), sizeof(bool));
    int import_deg = geometry_contextuality_degree_custom(&import_qa, false, false, false, SAT_SOLVER, bool_sol);
    assert_equal(import_deg,1, 
    "imported grid has contextuality degree 1");

    /////////////////////////////

    int heuristic_deg = geometry_contextuality_degree_custom(&import_qa, false, false, true, INVALID_LINES_HEURISTIC_SOLVER, bool_sol);
    assert_equal(heuristic_deg,1,
    "imported grid has contextuality degree 1 with heuristic");
    
    /////////////////////////////

    assert_equal(1, check_contextuality_solution(&import_qa, bool_sol, stderr), "checker is ok for grid");

    /////////////////////////////

    assert_equal(1, check_structure(&import_qa,bool_sol,false, NULL), 
    "structure is ok for grid");

    /////////////////////////////

    quantum_assignment troily = subspaces(3,1);
    int troily_heuristic_deg = geometry_contextuality_degree_custom(&troily, false, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);
    assert_equal(troily_heuristic_deg,63,
    "Troily has contextuality degree 63 with heuristic");

    free_quantum_assignment(&troily);

    /////////////////////////////

    FILE *mer_hypergraph = fopen("./misc/grid.hypergraph.txt", "r");
    FILE *mer_gram = fopen("./misc/grid.gram.txt", "r");

    if(mer_hypergraph == NULL || mer_gram == NULL){
        printf("Error: could not open files\n");
        return 1;
    }
    hypergram mermin_hg = hypergram_create_from_file(mer_hypergraph, mer_gram);

    assert_equal(mermin_hg.cpt_geometries, 6, "The imported grid  + hyperedges has 6 geometries");
    hypergram_print(mermin_hg);

    /////////////////////////////

    //rewind the  gram file to read the assignment
    rewind(mer_gram);

    hypergram mermin_hg_r = hypergram_create_from_gram_file(mer_gram);
    assert_equal(mermin_hg_r.cpt_geometries, 6, "The imported grid has 6 geometries from gram file");

    /////////////////////////////

    quantum_assignment mermin_qa = hypergram_to_quantum_assignment(mermin_hg);

    assert_equal(mermin_qa.n_qubits, 2, "The imported grid has 2 qubits");

    /////////////////////////////

    int mermin_deg = geometry_contextuality_degree_custom(&mermin_qa, false, false, false, INVALID_LINES_HEURISTIC_SOLVER, NULL);

    assert_equal(mermin_deg,1,"Mermin square has contextuality degree 1");

    fclose(mer_hypergraph);
    fclose(mer_gram);
    hypergram_free(mermin_hg);

    /////////////////////////////

    quantum_assignment hexagon = {0},complement_hexagon = {0};
    
    int cpt = 0;

    while (cpt < 1 && next_classical_cayley_hexagon(lines_qa_three, lines_indices, &hexagon, &complement_hexagon))
    {
        int c_degree = geometry_contextuality_degree_custom(&hexagon, true, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);
        int c_degree_comp = geometry_contextuality_degree_custom(&complement_hexagon, true, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);

        assert_equal(c_degree,0,"Hexagon has contextuality degree 0");
        assert_equal(c_degree_comp,0,"Complement hexagon has contextuality degree 0");
        cpt++;
    }

    cpt = 0;

    while (cpt < 2 && next_skew_cayley_hexagon(lines_qa_three, lines_indices, &hexagon, &complement_hexagon))
    {
        int c_degree = geometry_contextuality_degree_custom(&hexagon, true, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);
        int c_degree_comp = geometry_contextuality_degree_custom(&complement_hexagon, true, false, true, INVALID_LINES_HEURISTIC_SOLVER, NULL);

        assert_equal(c_degree,0,"Skew Hexagon has estimated c. degree 0");
        assert_equal(c_degree_comp,24,"Skew Complement hexagon has estimated c. degree 24");

        cpt++;
    }

    free_quantum_assignment(&hexagon);
    free_quantum_assignment(&complement_hexagon);

    free_cayley_hexagons();

    /////////////////////////////

    free_quantum_assignment(&import_qa);
    free_matrix(import_qa.geometries);
    free(bool_sol);

    /////////////////////////////

    free_quantum_assignment(&lines_qa_three);
    free_matrix(lines_indices);
    free_matrix(lines_qa_three.geometries);

    print_summary();

    return 0;
}

#endif // TEST_C