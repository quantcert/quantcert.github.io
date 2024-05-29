/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
/**********************************************************************************/
/**
 *
 * @file contextuality_degree.c
 * @author Axel Muller
 * @brief
 * @version 0.1
 * @date 2023-11-17
 *
 */
#ifndef CDEGREE
#define CDEGREE 1

#include "bv.c"
#include "complex.c"
#include <regex.h>
#include <errno.h>
#include <math.h>

#define STR_BUFFER_SIZE 4096

/**
 * @brief quantum assignment used to check contextuality
 * 
 * @param geometry_indices lists all the geometries of the geometry array to take into account to compute the contextuality degree
 * @param geometries list of all (not necessary) the geometries
 * @param cpt_geometries number of geometries checked
 * @param points_per_geometry number of observables in each geometry
 * @param n_qubits number of qubits per observable
 * 
 */
typedef struct{
    size_t* geometry_indices;
    bv** geometries;
    size_t cpt_geometries;
    size_t points_per_geometry;
    int n_qubits;

    bool* lines_negativity;
} quantum_assignment;

typedef enum
{
    SAT_SOLVER,
    RETRIEVE_SOLUTION,
} solver_mode;

int global_solver_mode = SAT_SOLVER;


/**
 * @brief return true with probability p
 * 
 * @param p 
 * @return true 
 * @return false 
 */
bool rand_float(float p) {
    return ((float)rand() / RAND_MAX) <= p;
}

/**
 * @brief prints a bool array grouping every 5 bits as a character
 * 
 * @param arr 
 * @param size 
 */
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

/**
 * @brief parses a string representing a bool array grouping every 5 bits as a character
 * 
 * @param arr 
 * @param size 
 */
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

/**
 * @brief returns true if the product of all the observables is minus the identity,
 * and false if it is the identity matrix 
 * If the product is not the identity modulo some phase, an error is printed and
 * the result can be anything
 * 
 * @param geometry array of observables
 * @param size size of the observable array
 * @param n_qubits number of qubits we work on
 * @param verbose if true prints each step of the operation
*/
bool is_negative_custom(bv geometry[],int size,int n_qubits,bool verbose,FILE* output){
    pauli_matrix mat = get_matrix(I);

    /*for each qubit*/
    for (int i = 0; i < n_qubits; i++){
        /*for each observable*/

        pauli_matrix printed = get_matrix(I);

        for (int j = size-1; j >= 0; j--){/*multiplication is performed from right to left*/
            /*we compute the product of the nth qubits of each observable
            (it should be the identity modulo some phase)
            and then we compute the product of all these matrices for all qubits*/
            bv gate = get_gate(geometry[j],i,n_qubits);
            mat = matrix_mult(mat,get_matrix(gate));
            if(verbose || true)printed = matrix_mult(printed,get_matrix(gate));
        }  
        
        W2_complex phase = get_id_matrix_phase(printed);
        if(verbose){
            print_complex_to_file(phase,output);
            if(i != n_qubits-1)fprintf(output," * ");
        }
        //if(phase.imag != 0 || phase.real != 1)print("ERROR : QUBIT PRODUCT IS NOT +I : %d+%di; ",phase.real,phase.imag);//uncomment to check for total positivity
    }
    
    
    /*we check that the phase is correct*/
    W2_complex phase = get_id_matrix_phase(mat);

    if(verbose){
        fprintf(output," = ");
        print_complex_to_file(phase,output);
    }

    if(phase.imag != 0 || phase.real == 0){
        print("geometry phase error");
    }
    return phase.real == -1;
}
bool is_negative(bv geometry[],int size,int n_qubits){return is_negative_custom(geometry,size,n_qubits,false,stderr);}

bool quantum_assignment_autofill_indices(quantum_assignment* qa){
    bool has_no_indices = qa->geometry_indices == NULL;
    if(has_no_indices){
        qa->geometry_indices = calloc(qa->cpt_geometries,sizeof(size_t));
        for (size_t i = 0; i < qa->cpt_geometries; i++)qa->geometry_indices[i] = i;
    }
    return has_no_indices;
}

/**
 * @brief Prints a quantum assignment
 * 
 * @param qa 
 */
void print_quantum_assignment(quantum_assignment qa){
    bool no_indices = quantum_assignment_autofill_indices(&qa);
    
    print("geometries : \n");
    for (size_t i = 0; i < qa.cpt_geometries; i++){
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            print_BV_custom(qa.geometries[qa.geometry_indices[i]][j],qa.n_qubits);
            //print("%ld ",qa.geometries[qa.geometry_indices[i]][j]);
        }
        print("%c",is_negative_custom(qa.geometries[qa.geometry_indices[i]],qa.points_per_geometry,qa.n_qubits,false,NULL)?'-':'+');
        print("\n");
    }
    print("\n");

    if(no_indices)free(qa.geometry_indices);
}

/**
 * @brief Computes the negativite contexts of a quantum assignment
 * 
 * @param qa 
 */
bool compute_negativity(quantum_assignment* qa){

    if(qa->lines_negativity != NULL)return false;

    /* bool alloc_indices =  */quantum_assignment_autofill_indices(qa);

    qa->lines_negativity = calloc(qa->cpt_geometries,sizeof(bool));
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        qa->lines_negativity[i] = is_negative_custom(qa->geometries[qa->geometry_indices[i]],qa->points_per_geometry,qa->n_qubits,false,NULL);
    }

    // if(alloc_indices){
    //     free(qa->geometry_indices);
    //     qa->geometry_indices = NULL;
    // }

    return true;
}

/**
 * @brief Conts the number of negative lines in a quantum assignment
 * 
 * @param qa quantum assignment
 * @return int 
 */
int negative_lines_count(quantum_assignment qa){
    bool no_indices = quantum_assignment_autofill_indices(&qa);
    bool no_negativity = compute_negativity(&qa);
    
    int res = 0;
    for (size_t i = 0; i < qa.cpt_geometries; i++){
        if(qa.lines_negativity[i])res++;
    }
    if(no_indices)free(qa.geometry_indices);
    if(no_negativity)free(qa.lines_negativity);

    return res;
}

/**
 * @brief frees a quantum assignment
 * 
 * @param qa 
 */
void free_quantum_assignment(quantum_assignment* qa){
    free(qa->geometry_indices);
    free(qa->lines_negativity);
    qa->geometry_indices = NULL;
    qa->lines_negativity = NULL;
}


void swap(int* a,int* b){
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void compute_invalid_lines(quantum_assignment qa,bool* bool_sol,bool** invalid_lines,int** number_of_invalid_lines){
    /*true for each invalid line*/
    (*invalid_lines) = (bool *)calloc(qa.cpt_geometries, sizeof(bool));
    /*number of invalid lines per point*/
    (*number_of_invalid_lines) = (int *)calloc(BV_LIMIT_CUSTOM(qa.n_qubits), sizeof(int));

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

/**
 * @brief internal function to find the index of a line in a table
 * 
 * @param table 
 * @param line 
 * @param capacity 
 * @return int 
 */
int get_table_index(int** table,int* line,size_t capacity,size_t nb_points_per_line){
    for (size_t i = 0; i < capacity; i++)
    {
        bool valid = true;
        bool all_zeroes = true;
        for (size_t j = 0; j < nb_points_per_line; j++)
        {
            if(table[i][j] != 0)all_zeroes = false;//we check at the same time if the entry is empty
            if(table[i][j] != line[j]){
                valid = false;
                if(!all_zeroes)break;
            }
        }
        if(valid)return i;//if found return index
        if(all_zeroes){//else we create that entry
            for (size_t j = 0; j < nb_points_per_line; j++)
            {
                table[i][j] = line[j];
            }
            return i;
        }
    }
    
    return -1;
}

int compare_integers_increasing(const void *a, const void *b){
    return (*(int *)a - *(int *)b);
}

int compute_category_list(bool to_print,int** category_list,int* category_count,int MAX_DEG,quantum_assignment qa){
    int complexity_degree = 0;

    if (to_print)
        print("\n\nLine degree types:\n");
    for (int i = 0; i < MAX_DEG; i++)
    {
        bool all_zeroes = true;
        for (size_t j = 0; j < qa.points_per_geometry; j++)
            if (category_list[i][j] != 0)
                all_zeroes = false;
        if (all_zeroes)
        {
            complexity_degree = i;
            break;
        }
        if (to_print)
        {
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
    return complexity_degree;
}

/**
 * @brief Checks the geometric structure of the invalid lines of a quantum assignment
 * 
 * @param qa 
 * @param bool_sol given solution
 * @param to_print if true prints the results
 * @param line_filter if not NULL, only checks the lines that are true in the filter
 * @return int the number of different types of invalid lines
 */
int check_structure(quantum_assignment qa,bool* bool_sol,bool to_print,bool *line_filter){
    
    compute_negativity(&qa);
    
    
    
    int ccpt = 0;
    if(to_print && line_filter != NULL)for (size_t i = 0; i < (size_t)NB_LINES_CUSTOM(qa.n_qubits); i++)ccpt++;
    
    

    bool* invalid_lines = NULL;
    int* number_of_invalid_lines = NULL;
    compute_invalid_lines(qa,bool_sol,&invalid_lines,&number_of_invalid_lines);

    /*for each number of invalid lines, we count the number of points with that number of invalid lines*/
    int number_of_obs_with_degree[BV_LIMIT_CUSTOM(qa.n_qubits)];
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)number_of_obs_with_degree[i] = 0;
    
    for (bv i = 1; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++){
        number_of_obs_with_degree[number_of_invalid_lines[i]]++;
    }
    if(to_print)print("\nPoint degree types:\n");
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++){
        if (number_of_obs_with_degree[i] != 0){
            if (to_print)
                print("degree %ld (x%d);", i, number_of_obs_with_degree[i]);
        }
    }

    int MAX_DEG = 20000000;

    /*for each type of invalid line, determined by the degrees of its vertices, we count the 
    number of lines with that type*/
    int** category_list = (int**)init_matrix(MAX_DEG,qa.points_per_geometry,sizeof(int));
    /*number of lines with that category*/
    int* category_count = calloc(MAX_DEG,sizeof(int));

    
    

    for (int i = 0; i < MAX_DEG; i++){
        category_count[i] = 0;
        for (size_t j = 0; j < qa.points_per_geometry; j++)category_list[i][j] = 0;
    }
    /*in case we want to check the presence of a specific type of line*/
    bool is_line_in_specific_type[qa.cpt_geometries];
    for (size_t i = 0; i < qa.cpt_geometries; i++)is_line_in_specific_type[i] = false;
    
    int obs_specific_type_degree[BV_LIMIT_CUSTOM(qa.n_qubits)];
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)obs_specific_type_degree[i] = 0;

    int mode = 1;

    int param[qa.points_per_geometry];

    if(to_print && SEE_GRAPH){
        int nul = 0;
        print("\n\nType search mode (0->stop,1->all,2->point degree,3->observable value,4->line degree,5->symmetric): ");
        nul += scanf("%d",&mode);
        if(mode == 2 || mode == 3){
            print((mode==2)?"point degree : ":"observable value: ");
            nul += scanf("%d",&param[0]);
        }else if(mode == 4){
            print("point 1 : ");
            nul += scanf("%d",&param[0]);
            print("point 2 : ");
            nul += scanf("%d",&param[1]);
            print("point 3 : ");
            nul += scanf("%d",&param[2]);
        }
        if(!nul)print("\noops, wrong text!\n");
    }

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        if(!invalid_lines[i] || (line_filter != NULL && !line_filter[qa.geometry_indices[i]]))continue;

        int sorted_tab[qa.points_per_geometry];
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            sorted_tab[j] = number_of_invalid_lines[qa.geometries[qa.geometry_indices[i]][j]];
            if(sorted_tab[j] >= MAX_DEG)print("\n!!%d!!\n",sorted_tab[j]);
        }
        qsort(sorted_tab, qa.points_per_geometry, sizeof(int), compare_integers_increasing);

        if(to_print || true){

            bool pass = true;
            if (mode == 0)pass = false;
            if(mode == 2){
                pass = false;
                for (size_t j = 0; j < qa.points_per_geometry; j++)if(sorted_tab[j] != param[0])pass = true;
            }
            if(mode == 3){
                pass = false;
                for (size_t j = 0; j < qa.points_per_geometry; j++)if(qa.geometries[qa.geometry_indices[i]][j] == (bv)param[0])pass = true;
            }
            if(mode == 4){
                for (size_t j = 0; j < qa.points_per_geometry; j++)if(sorted_tab[j] != param[j])pass = false;
            }
            if(mode == 5){
                for (size_t j = 0; j < qa.points_per_geometry; j++)if(!is_symmetric(qa.geometries[qa.geometry_indices[i]][j],qa.n_qubits))pass = false;
            }
            if(!pass)continue;

            is_line_in_specific_type[i] = true;
            for (size_t j = 0; j < qa.points_per_geometry; j++)
            {
                obs_specific_type_degree[qa.geometries[qa.geometry_indices[i]][j]]++;
            }
        }
    }

    if(to_print)print("\n\nspecific degrees : ");

    /*for each present point degree, we count the number of points with that degree*/
    int cpt_specific_degrees[BV_LIMIT_CUSTOM(qa.n_qubits)];
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)cpt_specific_degrees[i] = 0;
    
    for (bv j = 0; j < BV_LIMIT_CUSTOM(qa.n_qubits); j++)cpt_specific_degrees[obs_specific_type_degree[j]]++;
    
    if(to_print)for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)if(cpt_specific_degrees[i] != 0)print("%ld(x%d)",i,cpt_specific_degrees[i]);

    if(to_print)print("\n\nspecific / obs / general");

    int* tab_to_sort = obs_specific_type_degree;//number_of_invalid_lines;

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        int order[qa.points_per_geometry];
        for (size_t i = 0; i < qa.points_per_geometry; i++)order[i] = i;

        /*ordering the points in each context (bubble sort)*/
        for (size_t j = 0; j < qa.points_per_geometry - 1; j++) {
            for (size_t k = 0; k < qa.points_per_geometry - j - 1; k++) {
                if (tab_to_sort[qa.geometries[qa.geometry_indices[i]][order[k + 1]]] <
                    tab_to_sort[qa.geometries[qa.geometry_indices[i]][order[k]]]) {
                    swap(&order[k + 1], &order[k]);
                }
            }
        }

        if(!is_line_in_specific_type[i])continue;

        int line[qa.points_per_geometry];
        for (size_t j = 0; j < qa.points_per_geometry; j++)line[j] = tab_to_sort[qa.geometries[qa.geometry_indices[i]][order[j]]];
        int index = get_table_index(category_list,line,MAX_DEG,qa.points_per_geometry);
        category_count[index]++;

        if(!to_print)continue;

        print("\n");
        for (size_t j = 0; j < qa.points_per_geometry; j++)print("%d ",obs_specific_type_degree[qa.geometries[qa.geometry_indices[i]][order[j]]]);
        for (size_t j = 0; j < qa.points_per_geometry; j++)print_BV(qa.geometries[qa.geometry_indices[i]][order[j]]);
        for (size_t j = 0; j < qa.points_per_geometry; j++)print(" %d",number_of_invalid_lines[qa.geometries[qa.geometry_indices[i]][order[j]]]);

    }
    

    int vertices_count = 0;
    for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)
        if (number_of_invalid_lines[i] > 0)
            vertices_count++;

    if (to_print)
        print("\n\nNumber of vertices in the invalid configuration: %d ; \n", vertices_count);

    int complexity_degree = compute_category_list(to_print,category_list,category_count,MAX_DEG,qa);

    if(to_print
    )print("\nNumber of different line types: %d\n",complexity_degree);

    free_matrix(category_list);
    free(category_count);
    free(number_of_invalid_lines);
    free(invalid_lines);

    return complexity_degree;
}

/**
 * @brief returns the hamming distance between a quantum assignment and a given
 * classical solution in a boolean form
 * 
 * @param qa 
 * @param bool_sol boolean array containing a solution
 * @param output file descriptor where the computation is printed (NULL if no output)
 * @return int hamming distance between the quantum assignment and the solution
 */
int check_contextuality_solution(quantum_assignment qa,bool* bool_sol,FILE* output){

    bool alloc_neg = compute_negativity(&qa);

    int c_deg = 0,neg = 0;

    bool to_print = (output != NULL);

    if(to_print){
        fprintf(output,"\n___________________________\n");
        fprintf(output,"context => quantum / classical\n");
    }
    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        if(to_print)for (size_t j = 0; j < qa.points_per_geometry; j++){
            bv bv1 = qa.geometries[qa.geometry_indices[i]][j];
            if(bv1 == I)break;
            print_BV_to_file(bv1,qa.n_qubits,output);
            fprintf(output,"(%s)",bool_sol[bv1]?"-1":"+1");
        }
        if(to_print)fprintf(output,"  =>  ");
        
        bool is_neg = qa.lines_negativity[i];
        if(is_neg)neg++;

        

        bool test = false;
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            bv bv1 = qa.geometries[qa.geometry_indices[i]][j];
            if(bv1 == I)break;
            test ^= bool_sol[bv1];
            
        }    
        c_deg += (test != is_neg);

        if(to_print){
            fprintf(output," %s / %s : %s", is_neg ?"-1":"+1", test?"-1":"+1",(test == is_neg)?"valid    ":"invalid");
            if(test != is_neg)fprintf(output," %d",c_deg);
            
            fprintf(output,"\n");
        }

    }
    if(to_print){
        fprintf(output, "\nnegative lines : %d\nHamming distance : %d\n", neg, c_deg);
        
    }

    

    if(alloc_neg)free(qa.lines_negativity);

    return c_deg;

}

// Function to determine the maximum number of columns
void parse_matrix_dimensions(FILE *file, size_t *rows, size_t *cols)
{
    char line[STR_BUFFER_SIZE];
    *rows = 0;

    while (fgets(line, sizeof(line), file))
    {
        size_t count = 0;
        char *token = strtok(line, ",");
        while (token){
            count++;
            token = strtok(NULL, ",");
        }
        if (count > *cols)*cols = count;
        if (count > 1)(*rows)++;
    }
}

/**
 * @brief parses a file containing a quantum assignment
 * 
 * @param file 
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_parse(FILE *file){

    quantum_assignment qa = {0};

    print("Parsing file\n");

    parse_matrix_dimensions(file,&qa.cpt_geometries,&qa.points_per_geometry);

    if ((int)qa.points_per_geometry == -1){
        fclose(file);return qa;
    }
    print("size : %ldx%ld\n",qa.cpt_geometries,qa.points_per_geometry);

    rewind(file); // Reset the file pointer to the beginning

    // Count the number of lines (rows)
    qa.geometries = (bv**)init_matrix(qa.cpt_geometries, qa.points_per_geometry , sizeof(bv));

    char line[STR_BUFFER_SIZE];
    // Read the file again and fill the matrix
    size_t row = 0;
    while (fgets(line, sizeof(line), file) && row < qa.cpt_geometries)
    {
        int col = 0;
        char *token = strtok(line, ",");
        while (token){
            if (qa.n_qubits == 0 && strchr(token, ' ') == NULL && strlen(token) > 1)qa.n_qubits = strlen(token);
            qa.geometries[row][col++] = str_to_bv_custom(token, qa.n_qubits);
            token = strtok(NULL, ",");
        }
        if(col > 1)row++;
    }

    quantum_assignment_autofill_indices(&qa);
    compute_negativity(&qa);

    return qa;
}

/**
 * @brief returns the hamming distance between a quantum assignment and a given 
 * classical solution using a SAT solver
 * 
 * @param qa 
 * @param bc2cnf_file name of the file containing the BC clauses
 * @param sat_file name of the file containing the SAT clauses
 * @param output file descriptor where the computation is printed (NULL if no output)
 * @param ret_sol boolean array containing the solution
 * @return int hamming distance between the quantum assignment and the solution
 */
int compute_contextuality_solution(quantum_assignment qa,char* bc2cnf_file,char* sat_file,FILE* output,bool* ret_sol){

    FILE *fp;
    char line[STR_BUFFER_SIZE];
    regex_t sat_re;
    regmatch_t rm[3];

    char *tofind_sat = "(-?[0-9]+)";

    if (regcomp(&sat_re, tofind_sat, REG_EXTENDED) != 0){
        print("Failed to compile regex '%s'\n", tofind_sat);return -1;
    }
    //print("Number of captured expressions: %zu\n", sat_re.re_nsub);
    bool* sat_neg = calloc(100*qa.cpt_geometries,sizeof(bool));
    
    fp = fopen(sat_file,"r");

    if(fp == 0){
        print("no solution found !\n");
        return -1;
    }
    while ((fgets(line, STR_BUFFER_SIZE, fp)) != NULL)
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

    bool* bool_sol = calloc(BV_LIMIT_CUSTOM(qa.n_qubits),sizeof(bool));
    
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

    while ((fgets(line, STR_BUFFER_SIZE, fp)) != NULL){

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

    if(ret_sol != NULL)for (size_t i = 0; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++)ret_sol[i] = bool_sol[i];

    free(sat_neg);
    free(bool_sol);
    regfree(&re);

    return res;
}


/**
 * @brief Writes a geometry into a file to be read by the bc2cnf program.
 * Computes the negativity of the geometry (the sum of all the observable must be the identity modulo some phase)
 * 
 * @param geometry array of geometries
 * @param l index of the geometry we want to print
*/
void write_line(FILE *f,bv* geometry,bool negative,int points_per_geometry){

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



/**
 * @brief Returns the contextuality degree of a list of geometries
 * 
 * @param geometry_indices lists all the geometries of the geometry array to take into account to compute the contextuality degree
 * @param geometries list of all (not necessary) the geometries
 * @param cpt_geometries number of geometries checked
 * @param points_per_geometry number of observables in each geometry
 * @param n_qubits number of qubits per observable
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
*/
int geometry_SAT_contextuality_degree(quantum_assignment qa,bool contextuality_only,bool print_solution,bool optimistic,bool* ret_sol){

    if(qa.cpt_geometries == 0)return -1;

    /*spasm use example*/
    
    
    
    
    if(print_solution)print("\nlines : %ld ; negative lines : %d\n",qa.cpt_geometries,negative_lines_count(qa));
    /*heuristics*/
    

    bool no_ret_sol = (ret_sol == NULL);

    int th_num = omp_get_thread_num();

    char bc2cnf_file[128] = {0},sat_log[128] = {0},bc2cnf_cp[128] = {0},sat_cp[128] = {0},sat_command[4096],cp_command[4096];
    char* heuristic = (optimistic?"--sat":"--unsat");
    sprintf(bc2cnf_file,".tmp%d.txt",th_num);
    sprintf(sat_log,".tmp%d.log",th_num);
    sprintf(bc2cnf_cp,".bak%d.txt",th_num);
    sprintf(sat_cp,".bak%d.log",th_num);
    
    sprintf(sat_command,"./lib/BCpackage-0.40/bc2cnf -nosimplify > %s && ./lib/kissat_gb/build/kissat %s -q %s > %s",
        bc2cnf_file,
        heuristic,
        bc2cnf_file,
        sat_log);
    sprintf(cp_command,"cp %s %s && cp %s %s",bc2cnf_file,bc2cnf_cp,sat_log,sat_cp);

    /*We check the contextuality degree */
    /*First we test the maximal contextuality degree possible,check its satisfiabiliy, then 
    decremeent it until it is no longer satisfiable*/
    int neg_lines = negative_lines_count(qa);
    if(neg_lines == 0 && false){
        if(print_solution)printf("positive config. ");
        return 0;
    }
    if (no_ret_sol)ret_sol = calloc(BV_LIMIT_CUSTOM(qa.n_qubits), sizeof(bool));

    int start_degree = contextuality_only?0:/* cpt_geometries - 1 ;*/neg_lines;
    
    int c_degree_test = start_degree;
    int hamming_distance = c_degree_test;

    print("\nStarting SAT computation...\n");

    /*While we haven't tested all possible degrees*/
    while(c_degree_test >= 0 && !is_done){
        /*pipeline process : we use bc2cnf to transform the problem into a DIMACS CNF form*/
        //
        FILE* sat_fp = popen(sat_command,"w");

        /*We print each line of the problem into the pipeline*/

        /*For G geometries and a potential degree D, is it possible to satisfy at least G-D 
        equations and at most G ones.*/
        fprintf(sat_fp,"BC1.1\nASSIGN ");
        if(!contextuality_only)fprintf(sat_fp,"[%ld,%ld](",qa.cpt_geometries-c_degree_test,qa.cpt_geometries);

        for (size_t i = 0; i < qa.cpt_geometries; i++){
            if(i != 0)fprintf(sat_fp,",");
            
            write_line(sat_fp,qa.geometries[qa.geometry_indices[i/* (i==cpt_geometries-1)?i:(i+1-2*(i%2)) */]],qa.lines_negativity[i],qa.points_per_geometry);
        }

        if(!contextuality_only)fprintf(sat_fp,")");
        fprintf(sat_fp,";");

        /*the exit status tells us if the problem is SAT or not*/
        int status = WEXITSTATUS(pclose(sat_fp));
        bool is_sat_test = status == 10;//else 20//is_sat(f);
        bool expected_output = status == 10 || status == 20;
        // if(print_solution)print("%d,",c_degree_test);
        /*if the exit status is something else than what is expected, we print it*/
        if (!expected_output) print("status error: {%d}", status);
        if (expected_output  && !is_sat_test){
            if(print_solution)print("\nContextuality degree found : %d\n",hamming_distance);
            break;
        }
        if(system(cp_command) != EXIT_SUCCESS)print("copy error!");

        if(!contextuality_only && expected_output){
            hamming_distance = compute_contextuality_solution(qa,bc2cnf_cp,sat_cp,NULL,ret_sol);
            if(hamming_distance == -1)print("\nunexpected hamming distance error\n");
            if(print_solution)print("current Hamming distance: %d\n",hamming_distance);
            c_degree_test/* --;// */ = hamming_distance-1;
        }else{
            c_degree_test--;
        }
    };
    print("Computation done\n");
    /*if the tested degree fails, then we go back to the one above*/

    /*if a solution is found*/
    #pragma omp critical
    {
        //compute_contextuality_solution(qa,bc2cnf_cp,sat_cp,print_solution?stdout:NULL,ret_sol);
        check_contextuality_solution(qa,ret_sol,print_solution?stdout:NULL);
    }


    remove(bc2cnf_file);
    remove(sat_log);
    remove(bc2cnf_cp);
    remove(sat_cp);
    if(no_ret_sol)free(ret_sol);

    return hamming_distance;
}

/**
 * @brief Returns a quantum assignment where the contexts are those that are invalid
 * given a classical assignment
 * 
 * @param qa 
 * @param bool_sol classical assignment
 * @param validity if true the contexts are those that are valid
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_from_invalid_contexts(quantum_assignment qa,bool* bool_sol,bool validity){
    quantum_assignment res = qa;
    res.geometry_indices = calloc(qa.cpt_geometries,sizeof(size_t));
    res.cpt_geometries = 0;
    res.lines_negativity = NULL;

    for (size_t i = 0; i < qa.cpt_geometries; i++)
    {
        bool classical_negativity = false;
        for (size_t j = 0; j < qa.points_per_geometry && qa.geometries[qa.geometry_indices[i]][j] != I; j++){
            classical_negativity ^= bool_sol[qa.geometries[qa.geometry_indices[i]][j]];
        }
        
        if((classical_negativity ^ !(validity)) == is_negative_custom(qa.geometries[qa.geometry_indices[i]],qa.points_per_geometry,qa.n_qubits,false,NULL)){
            res.geometry_indices[res.cpt_geometries] = qa.geometry_indices[i];
            res.cpt_geometries++;
        }
    }

    compute_negativity(&res);

    return res;
}

/**
 * @brief Returns the contextuality degree of a list of geometries with a given method
 * 
 * @param qa 
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual(only for the SAT solver)
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
 * @param mode method used to compute the contextuality degree
 * @param bool_sol if not NULL, the solution is stored in this array
 * @return int minimal hamming distance found
 */
int geometry_contextuality_degree_custom(quantum_assignment qa,bool contextuality_only,bool print_solution,bool optimistic,solver_mode mode,bool* bool_sol){

    bool has_indices = !quantum_assignment_autofill_indices(&qa);
    bool has_negative = !compute_negativity(&qa);

    bool has_bool_sol = bool_sol != NULL;

    int c_degree = -1;

    if(!DETERMINISTIC)srand(time(NULL));


    if(!has_bool_sol)bool_sol = calloc(BV_LIMIT_CUSTOM(qa.n_qubits),sizeof(bool));

    
    
    
    
    
        switch (mode){
        case SAT_SOLVER:c_degree = geometry_SAT_contextuality_degree(qa,contextuality_only,print_solution,optimistic,bool_sol);break;
        case RETRIEVE_SOLUTION:
            parse_bool(bool_sol,BV_LIMIT_CUSTOM(qa.n_qubits));
            c_degree = check_contextuality_solution(qa,bool_sol,print_solution?stdout:NULL);break;
        default:break;
        }

    

    
    
    

    

    

    
    
    
    
    
    

    
    
    
    
    
    int continued = true;

    while (true)
    {
        if(SEE_GRAPH){
            print("\nDo you want to see the contextual configuration ? (0:no, 1:yes): ");
            
            if(scanf("%d",&continued) == 0)print("oops wrong text");
        }else{
            continued = false;
        }
        if(!continued)break;
        check_structure(qa, bool_sol, stdout, NULL); 
        
    }

    if(print_solution){
        print("\nsolution code : ");
        print_bool(bool_sol, BV_LIMIT_CUSTOM(qa.n_qubits));
        print("\n");
    }
    //check_contextuality_solution(qa,bool_sol,stderr);

    if(!has_indices)free(qa.geometry_indices);
    if(!has_negative)free(qa.lines_negativity);
    if(!has_bool_sol)free(bool_sol);
    return c_degree;
}

/**
 * @brief Returns the contextuality degree of a quantum assignment
 * 
 * @param qa 
 * @param contextuality_only if true doesn't compute the degree but only wether or not the geometry is contextual
 * @param print_solution if true prints the solution if one is found
 * @param optimistic enables a sat solver heuristic making it faster to find a solution IFF there is one
 * @param bool_sol if not NULL, the solution is stored in this array
 * @return int minimal hamming distance found
 */
int geometry_contextuality_degree(quantum_assignment qa,bool contextuality_only,bool print_solution,bool optimistic,bool* bool_sol){
    return geometry_contextuality_degree_custom(qa,contextuality_only,print_solution,optimistic,global_solver_mode,bool_sol);
}    

#endif //CDEGREE