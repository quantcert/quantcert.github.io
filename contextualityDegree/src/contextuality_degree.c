/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file contextuality_degree.c
 * @brief contains methods to compute the contextuality degree of a quantum assignment
 */

#include "contextuality_degree.h"

#include <regex.h>
#include <errno.h>
#include <math.h>
#include "bv.h"
#include "quantum_assignment.h"
#include "config_checker.h"

#define HEURISTIC_NUM_THREADS (MULTI_THREAD ? (9) : (1)) // number of threads to use for the heuristic method

solver_mode global_solver_mode = SAT_SOLVER;
size_t global_heuristic_iterations = 10000; //maximum number of iterations for the heuristic method
float global_heuristic_flip_probability = 0.95;   //probability of choosing a random assignment in the heuristic method
float global_heuristic_threshold = DISABLED_PARAMETER;          // threshold for the heuristic method


bool rand_float(float p) {
    bool r = ((float)fast_random() / (float)UINT32_MAX) <= p;
    return r;
}


void print_bool(bool* arr,size_t size){
    for (size_t i = 0; i < size/4; i++)
    {
        char c = 'a';
        for (size_t j = 0; j < 4; j++)
        {
            size_t index = i*4+j;
            if(index >= size)return;
            c += (arr[index])<<j;
            //print("%d %d,",j,arr[index]);
        }
        print("%c",c);
    }
    print("\n");
}

void parse_bool(bool* arr,size_t size){
    print("enter a solution : ");
    char str[/* (size/4)+1+1 */10000];
    int ret = scanf("%s",str);
    if(ret == 0)return;
    for (size_t i = 0; i < size/4; i++)
    {
        char c = str[i];
        if(c < 'a' || c > 'z')return;
        c -= 'a';
        for (size_t j = 0; j < 4; j++)
        {
            size_t index = i*4+j;
            if(index >= size)return;
            arr[index] = BGET(c,j);
        }
    }
}


int check_contextuality_solution(quantum_assignment* qa,bool* bool_sol,FILE* output){

    quantum_assignment_compute_negativity(qa);

    int c_deg = 0,neg = 0;

    bool to_print = (output != NULL);

    if(to_print){
        fprintf(output,"\n___________________________\n");
        fprintf(output,"context => quantum / classical\n");
    }
    for (size_t i = 0; i < qa->cpt_geometries; i++)
    {
        if(to_print)
            for (size_t j = 0; j < qa->points_per_geometry; j++){
                bv bv1 = qa->geometries[qa->geometry_indices[i]][j];
                if (bv1 == I)break;
                print_BV_to_file(bv1, qa->n_qubits, output);
                fprintf(output, "(%s)", bool_sol[bv1] ? "-1" : "+1");
            }
        if(to_print)fprintf(output,"  =>  ");

        bool is_neg = qa->lines_negativity[i];
        if(is_neg)neg++;

        

        bool test = false;
        for (size_t j = 0; j < qa->points_per_geometry; j++){
            bv bv1 = qa->geometries[qa->geometry_indices[i]][j];
            if(bv1 == I)break;
            test ^= bool_sol[bv1];
            
        }
        // Add debug prints
        int increment = (test != is_neg) ? 1 : 0;
        c_deg += increment;

        if(to_print){
            fprintf(output," %s / %s : %s", is_neg ?"-1":"+1", test?"-1":"+1",(test == is_neg)?"valid    ":"invalid");
            if(test != is_neg)fprintf(output," %d",c_deg);
            
            fprintf(output,"\n");
        }

    }
    if(to_print){
        fprintf(output, "\nnegative lines :%d\nHamming distance :%d\n", neg, c_deg);
        
    }

    
    return c_deg;
}

int compute_contextuality_solution(quantum_assignment* qa,char* bc2cnf_file,char* sat_file,FILE* output,bool* ret_sol){

    FILE *fp;
    char line[CONFIG_CHECKER_STR_BUFFER_SIZE];
    regex_t sat_re;
    regmatch_t rm[3];

    char *tofind_sat = "(-?[0-9]+)";

    if (regcomp(&sat_re, tofind_sat, REG_EXTENDED) != 0){
        print("Failed to compile regex '%s'\n", tofind_sat);return -1;
    }
    //print("Number of captured expressions: %zu\n", sat_re.re_nsub);
    bool *sat_neg = calloc(100 * qa->cpt_geometries, sizeof(bool));

    fp = fopen(sat_file,"r");

    if(fp == 0){
        print("no solution found !\n");
        return -1;
    }
    while ((fgets(line, CONFIG_CHECKER_STR_BUFFER_SIZE, fp)) != NULL)
    {
        line[strcspn(line, "\n")] = '\0';

        char *sub_line = line;
        while(regexec(&sat_re, sub_line, 2, rm, 0) == 0){

            if(rm[1].rm_eo == 0)break;
            /*rm[0] is the line*/
            char int_str[100];
            sprintf(int_str,"%.*s", (int)(rm[1].rm_eo - rm[1].rm_so), sub_line+ rm[1].rm_so);
            long res = strtol(int_str,NULL,10);

            sat_neg[(res<0)?(-res):res] = !(res<0);

            sub_line += rm[1].rm_eo;
        }
    } 

    fclose(fp);

    bool *bool_sol = calloc(BV_LIMIT_CUSTOM(qa->n_qubits), sizeof(bool));

    regfree(&sat_re);
    regex_t re;

    char *tofind = "^c v([0-9]*) <-> ([0-9]*)$";

    if (regcomp(&re, tofind, REG_EXTENDED) != 0){
        print("Failed to compile regex '%s'\n", tofind);return -1;
    }
    //print("Regex: %s\n", tofind);
    //print("Number of captured expressions: %zu\n", re.re_nsub);
    fp = fopen(bc2cnf_file, "r");
    if (fp == 0){
        print( "Failed to open file %s (%d: %s)\n", bc2cnf_file, errno, strerror(errno));return -1;
    }

    while ((fgets(line, CONFIG_CHECKER_STR_BUFFER_SIZE, fp)) != NULL)
    {

        line[strcspn(line, "\n")] = '\0';

        if (regexec(&re, line, 3, rm, 0) == 0){
            /*rm[0] is the line*/
            char from_str[100],to_str[100];

            sprintf(from_str,"%.*s", (int)(rm[1].rm_eo - rm[1].rm_so), line+ rm[1].rm_so);
            sprintf(to_str,"%.*s", (int)(rm[2].rm_eo - rm[2].rm_so), line+ rm[2].rm_so);

            long from = strtol(from_str,NULL,10);
            long to = strtol(to_str,NULL,10);
            //print("((%ld,%ld))",from,to);
            bool_sol[from] = sat_neg[to];
        }
    }
    fclose(fp);
    /////////////////////printing the solution///////////////////////////////////////

    int res = check_contextuality_solution(qa,bool_sol,output);

    if(ret_sol != NULL)
        for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)ret_sol[i] = bool_sol[i];

    free(sat_neg);
    free(bool_sol);
    regfree(&re);

    return res;
}

void write_line(FILE *f,bv* geometry,bool negative,int points_per_geometry){
    /*edge case where the context is empty*/
    if(geometry[0] == I){
        fprintf(f,"T\n");
        return;
    }
    /*We compute the "negativeness" of a geometry before printing it as -1^x to
    solve it as a linear problem, which is why the expected sum is 1(odd) for negative geometries
    and 0(even) for the positive ones*/
    if(negative)
         fprintf(f,"ODD(");
    else fprintf(f,"EVEN(");

    /*We print all the observables of the geometry*/
    bool begin = true;
    for (int i = 0; i < points_per_geometry; i++){
        if(geometry[i] == I)break;//print("ERROR WRITE LINE;");
        if(begin){begin = false;}else{fprintf(f,",");}
        fprintf(f,"v%d",geometry[i]);
    }
    fprintf(f,")\n");
}

int max_line_per_point(quantum_assignment* qa){

    int cpt[BV_LIMIT_CUSTOM(qa->n_qubits)];
    int max = 0;
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)cpt[i] = 0;
    
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        for (size_t j = 0; j < qa->points_per_geometry; j++){
            bv bv1 = qa->geometries[qa->geometry_indices[i]][j];
            if(bv1 == I)break;
            int* test = &(cpt[bv1]);
            (*test)++;
            if(*test > max){
                max = *test;
                //print("max_i=%d\n",bv1);
            }
        }
    }
    return max;
}

int **compute_contexts_per_obs(quantum_assignment *qa, bool print_solution){

    const size_t max_line_per_obs = max_line_per_point(qa);
    // print("avg_line_per_obs:%ld\n",avg_line_per_obs);
    if (print_solution)print(".");
    

    /*stores for every observable the contexts it is present in*/
    int **line_per_obs = (int **)init_matrix(BV_LIMIT_CUSTOM(qa->n_qubits), max_line_per_obs, sizeof(int));
    if (print_solution)print(".");
    /*initializes the matrix with empty lines*/
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)
        for (size_t j = 0; j < max_line_per_obs; j++)
            line_per_obs[i][j] = -1;

    if (print_solution)print("l(x%ld)", qa->cpt_geometries);

    for (size_t i = 0; i < qa->cpt_geometries; i++)
    {
        if (i % 1000000 == 0 && i != 0 && print_solution)print("%ld,", i);

        for (size_t j = 0; j < qa->points_per_geometry; j++)
        {
            bv bv1 = qa->geometries[qa->geometry_indices[i]][j];
            if (bv1 == I)continue;/*An identity gate is the end of a context*/

            for (size_t k = 0; k < max_line_per_obs; k++){/* we look for the first empty line to place the context*/
                int *line = &(line_per_obs[bv1][k]);
                if (*line == -1){
                    *line = i;
                    break;
                }
                if (k == max_line_per_obs - 1)print("linecpt overflow !\n");
            }
        }
    }
    return line_per_obs;
}

int geometry_contextuality_degree_max_invalid_heuristics(quantum_assignment* qa,bool print_solution,bool* ret_sol){

    //initializes timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    if(print_solution)print("\n.");
    

    bool min_sol[BV_LIMIT_CUSTOM(qa->n_qubits)];
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)min_sol[i] = false;

    

    size_t max_line_per_obs = max_line_per_point(qa);
    int **line_per_obs = compute_contexts_per_obs(qa,print_solution);

    if (print_solution)print(".\n");

    int global_min = qa->cpt_geometries;//negative_lines_count(qa);
    int test_th_global_min = global_min;

    int nth_sol_global = 0;
    

    bool auto_threshold = global_heuristic_threshold == DISABLED_PARAMETER;

    float max_threshold = 1.0f;
    float optimal_threshold = auto_threshold?0.85f:global_heuristic_threshold;
    float min_threshold = 0.0f;
    float range_radius = 1.0f;


    #pragma omp parallel num_threads(HEURISTIC_NUM_THREADS)
    {
        int n_invalid[BV_LIMIT_CUSTOM(qa->n_qubits)];
        bool bool_sol[BV_LIMIT_CUSTOM(qa->n_qubits)];
        

        for (size_t j = 0; j < BV_LIMIT_CUSTOM(qa->n_qubits); j++){
            n_invalid[j] = 0;
            bool_sol[j] = false;//fast_random() % 2;
        }

        int hamming_test = check_contextuality_solution(qa,bool_sol,NULL);

        for (size_t i = 0; i < qa->cpt_geometries; i++){/*At the beginning, the number of invalid contexts for each observable is the number of negative contexts*/
            bool sign = false;
            for (size_t j = 0; j < qa->points_per_geometry; j++){
                sign ^= bool_sol[qa->geometries[qa->geometry_indices[i]][j]];
            }
            
            if(qa->lines_negativity[i] != sign){  
            
                for (size_t j = 0; j < qa->points_per_geometry; j++){
                    bv bv1 = qa->geometries[qa->geometry_indices[i]][j];
                    if(bv1 == I)break;
                    n_invalid[bv1]++;
                }
            }
        }
        int current_max = 0;

        int n_neg = 0;

        for (size_t cpt = 0; cpt < global_heuristic_iterations && !is_done; cpt++){
            
            
            
            int old_current_max = current_max;
            current_max = 0;/*maximum number of invalid contexts found for a single observable*/
            
            /*if the last thread ran 20 iterations, we shrink the range of possible 
            threshold(theta) values (if theta is not already specified)*/
            if(cpt%20 == 20-1 && auto_threshold && omp_get_thread_num() == omp_get_num_threads()-1){
                #pragma omp critical
                {
                    range_radius *= 0.5f;
                    max_threshold = MIN(optimal_threshold + range_radius,1.0f);
                    min_threshold = MAX(optimal_threshold - range_radius,0.0f);
                    
                    if (range_radius < 0.01f /*0.0f*/)
                    {/*we eventually reset the range*/
                        range_radius = 1.0f;
                        optimal_threshold = 0.5f;
                        test_th_global_min = qa->cpt_geometries;
                    }
                }
            }
            float threshold_select;
            float rand_select = global_heuristic_flip_probability;
            
            if (auto_threshold){/*every thread gets its different threshold*/
                float th_ratio = (float)omp_get_thread_num() / omp_get_num_threads();
                threshold_select = min_threshold + th_ratio * (max_threshold - min_threshold);
                
                if (omp_get_num_threads() == 1) threshold_select = optimal_threshold;
            }
            else{
                threshold_select = global_heuristic_threshold;
            }
            
            #pragma omp critical
            {
                //print("%d,",hamming_test);
                if(hamming_test < test_th_global_min){
                    test_th_global_min = hamming_test;
                    //print("test_th_global_min : %d\n",test_th_global_min);
                    optimal_threshold = threshold_select;
                }
                if(hamming_test <= global_min){/*if a thread found a lower bound that the current best one*/
                    
                    if (hamming_test < global_min && print_solution){
                        clock_gettime(CLOCK_MONOTONIC, &end);
                        double time_taken = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec)) * 1e-9;
                        print("current Hamming distance : %d, %.2fs\n", hamming_test, time_taken);
                        //print("(%.02f,%d[%ld]) ", time_taken, hamming_test,cpt);   
                        
                        check_structure(qa, bool_sol, false, NULL);
                    }    
                    global_min = hamming_test;
                    //print("(%.2f)",optimal_threshold);
                    
                    
                    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++){
                        
                        min_sol[i] = bool_sol[i];
                        
                    }
                    nth_sol_global++;
                    
                }
            }
            //if (hamming_test <= global_min && hamming_test < 100000 /* 134700 */)check_structure(qa, bool_sol, false, NULL);
            if(global_min == 0)break;/*if a fully valid solution has been found*/

            for (size_t i = I + 1; i < BV_LIMIT_CUSTOM(qa->n_qubits); i++)
            { /*for each observable*/
                /*selects the assignments that won't be flipped*/
                if (!(n_invalid[i] > old_current_max * threshold_select && rand_float(rand_select) /* && n_I_custom(i, qa->n_qubits) %2 == 0 */))
                    continue; 
                bool_sol[i] ^= true;
                n_neg += (bool_sol[i] ? +1 : -1);

                for (size_t j = 0; j < max_line_per_obs && line_per_obs[i][j] != -1; j++) /*for each context in an obs*/
                {
                    size_t line_index = line_per_obs[i][j];
                    size_t geo_line = qa->geometry_indices[line_index];
                    bv *line = qa->geometries[geo_line];
                    bool test = false;
                    for (size_t k = 0; k < qa->points_per_geometry; k++) /*We compute the new sign of the classical assignment*/
                    {
                        if ((line[k] == I))break;
                        test ^= bool_sol[line[k]];
                    }
                    /*we add -1 if the line is negative and +1 if it is positive*/
                    int adder = 1 - 2 * (test == qa->lines_negativity[line_index]); 

                    for (size_t k = 0; k < qa->points_per_geometry; k++) /*We update the number of invalid contexts for each observable*/
                    {
                        if ((line[k] == I))break;
                        int *inv = &(n_invalid[line[k]]);
                        *inv += adder;
                        
                        if (*inv > current_max)current_max = *inv;
                    }
                    hamming_test += adder; /*We update the hamming distance dynamically*/
                    
                }
            }
        }
        if (print_solution)print(".");
    }
    if(is_done)is_done = false;
    free_matrix(line_per_obs);
    /*if wanted, the solution is copied to the given array*/
    if(ret_sol != NULL)for(size_t i = 0; i < (size_t)BV_LIMIT_CUSTOM(qa->n_qubits); i++)ret_sol[i] = min_sol[i];

    if (print_solution)print("\nHamming distance found: %d\nepsilon (if minimal): %.3f\n", global_min, 2.0f * (float)global_min / (float)qa->cpt_geometries);

    return global_min;
}


int geometry_SAT_contextuality_degree(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,bool* ret_sol){

    if(qa->cpt_geometries == 0)return -1;

    // initializes timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    if(print_solution)print("\nlines : %ld ; negative lines : %d\n",qa->cpt_geometries,negative_lines_count(qa));
    /*heuristics*/
    

    bool no_ret_sol = (ret_sol == NULL);

    int th_num = omp_get_thread_num();

    char bc2cnf_file[128] = {0},sat_log[128] = {0},bc2cnf_cp[128] = {0},sat_cp[128] = {0},sat_command[4096],cp_command[4096];
    char* heuristic = (optimistic?"--sat":"--unsat");
    sprintf(bc2cnf_file,"tmp%d.txt",th_num);
    sprintf(sat_log,"tmp%d.log",th_num);
    sprintf(bc2cnf_cp,"bak%d.txt",th_num);
    sprintf(sat_cp,"bak%d.log",th_num);
    
    sprintf(sat_command, "./external/BCpackage-0.40/bc2cnf -nosimplify > %s && ./external/kissat_gb/build/kissat %s -q %s > %s",
            bc2cnf_file,
            heuristic,
            bc2cnf_file,
            sat_log);
    sprintf(cp_command,"cp %s %s && cp %s %s",bc2cnf_file,bc2cnf_cp,sat_log,sat_cp);

    /*We check the contextuality degree */
    /*First we test the maximal contextuality degree possible,check its satisfiabiliy, then 
    decrement it until it is no longer satisfiable*/
    int neg_lines = negative_lines_count(qa);
    if(neg_lines == 0 && false){
        if(print_solution)printf("positive config. ");
        return 0;
    }
    if (no_ret_sol)ret_sol = calloc(BV_LIMIT_CUSTOM(qa->n_qubits), sizeof(bool));

    int start_degree = contextuality_only?0:/* cpt_geometries - 1 ;*/neg_lines;
    
    int c_degree_test = start_degree;
    int hamming_distance = c_degree_test;

    if(print_solution)print("\nStarting SAT computation...\n");
    /*While we haven't tested all possible degrees*/
    while(c_degree_test >= 0 && !is_done){
        /*pipeline process : we use bc2cnf to transform the problem into a DIMACS CNF form*/
        //
        FILE* sat_fp = popen(sat_command,"w");

        /*We print each line of the problem into the pipeline*/

        /*For G geometries and a potential degree D, is it possible to satisfy at least G-D 
        equations and at most G ones.*/
        fprintf(sat_fp,"BC1.1\nASSIGN ");
        if(!contextuality_only)fprintf(sat_fp, "[%ld,%ld](", qa->cpt_geometries - c_degree_test, qa->cpt_geometries);

        for (size_t i = 0; i < qa->cpt_geometries; i++){
            if(i != 0)fprintf(sat_fp,",");
            write_line(sat_fp,qa->geometries[qa->geometry_indices[i/* (i==cpt_geometries-1)?i:(i+1-2*(i%2)) */]],qa->lines_negativity[i],qa->points_per_geometry);
        }

        if(!contextuality_only)fprintf(sat_fp,")");
        fprintf(sat_fp,";");

        /*the exit status tells us if the problem is SAT or not*/
        int close_status = pclose(sat_fp);

        int status = WEXITSTATUS(close_status);
        bool is_sat_test = status == 10;//else 20//is_sat(f);
        bool expected_output = status == 10 || status == 20;
        // if(print_solution)print("%d,",c_degree_test);
        /*if the exit status is something else than what is expected, we print it*/
        
        if (!expected_output){
            print("status error: {%d}", status);
            
            return hamming_distance;
        }
        if (expected_output  && !is_sat_test){
            if(print_solution)print("\nContextuality degree found: %d\nepsilon: %.3f", hamming_distance, 2.0f * (float)hamming_distance / (float)qa->cpt_geometries);
            break;
        }
        if(system(cp_command) != EXIT_SUCCESS)print("copy error!");

        if(!contextuality_only && expected_output){
            hamming_distance = compute_contextuality_solution(qa,bc2cnf_cp,sat_cp,NULL,ret_sol);
            if(hamming_distance == -1)print("\nunexpected hamming distance error\n");
            if(print_solution){
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time_taken = (end.tv_sec - start.tv_sec) + ( (end.tv_nsec - start.tv_nsec)) * 1e-9;
                print("current Hamming distance: %d, %.2fs\n", hamming_distance, time_taken);
            }

            c_degree_test/* --;// */ = hamming_distance-1;
        }else{
            c_degree_test--;
        }
    };

    /*if the tested degree fails, then we go back to the one above*/

    /*if a solution is found*/
    #pragma omp critical
    {
        //compute_contextuality_solution(qa,bc2cnf_cp,sat_cp,print_solution?stdout:NULL,ret_sol);
        check_contextuality_solution(qa,ret_sol,NULL);
    }

    remove(bc2cnf_file);
    remove(sat_log);
    remove(bc2cnf_cp);
    remove(sat_cp);
    if(no_ret_sol)free(ret_sol);

    return hamming_distance;
}

int geometry_contextuality_degree_custom(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,solver_mode mode,bool* bool_sol){

    quantum_assignment_autofill_indices(qa);
    quantum_assignment_compute_negativity(qa);

    bool has_bool_sol = bool_sol != NULL;

    int c_degree = -1;


    if(!has_bool_sol)bool_sol = calloc(BV_LIMIT_CUSTOM(qa->n_qubits),sizeof(bool));

    
    

    switch (mode){
    case SAT_SOLVER:c_degree = geometry_SAT_contextuality_degree(qa,contextuality_only,print_solution,optimistic,bool_sol);break;
    case RETRIEVE_SOLUTION:
        parse_bool(bool_sol,BV_LIMIT_CUSTOM(qa->n_qubits));
        c_degree = check_contextuality_solution(qa,bool_sol,NULL);break;
    case INVALID_LINES_HEURISTIC_SOLVER:c_degree = geometry_contextuality_degree_max_invalid_heuristics(qa, print_solution, bool_sol);break;
    default:break;
    }



    if(!has_bool_sol)free(bool_sol);
    return c_degree;
}

int geometry_contextuality_degree(quantum_assignment* qa,bool contextuality_only,bool print_solution,bool optimistic,bool* bool_sol){
    return geometry_contextuality_degree_custom(qa,contextuality_only,print_solution,optimistic,global_solver_mode,bool_sol);
}    