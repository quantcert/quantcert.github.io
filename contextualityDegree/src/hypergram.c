/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file hypergram.c
 * @brief This file contains the code used to manipulate the hypergrams
 * described in the paper
 * An abstract structure determines the contextuality degree of observable-based Kochen-Specker proofs
 * by Axel Muller and Alain Giorgetti
 */


#include "hypergram.h"

#include <unistd.h>
#include "constants.h"
#include "complex_int.h"
#include "quadrics.h"

bool is_gram_matrix_valid(bit_matrix bm)
{
    if (bm.size != bm.bit_sets[0].size)
    {
        print("error : the matrix is not square : %ldx%ld\n", bm.size, bm.bit_sets[0].size);
        return false;
    }
    for (size_t i = 0; i < bm.size; i++)
    {
        if(bit_matrix_get_bit(bm,i,i) != 0){
            print("error : the matrix is not valid : %ldx%ld\n",i,i);
            return false;
        }
        for (size_t j = 0; j < bm.size; j++)
        {
            if (bit_matrix_get_bit(bm,i, j) != bit_matrix_get_bit(bm,j, i))
            {
                print("error : the matrix is not symmetric : %ldx%ld\n", i, j);
                return false;
            }
        }
    }
    
    return true;
}

bit_matrix incidence_matrix(hypergram ccs)
{

    size_t height = ccs.commutation_matrix.size;
    

    bit_matrix bm = bit_matrix_create(ccs.cpt_geometries, height/*cpt_points+1*/);

    for (size_t i = 0; i < ccs.cpt_geometries; i++)
    {
        for (size_t j = 0; j < ccs.max_points_per_geometry; j++)
        {
            size_t index = ccs.geometries[i][j];
            if (index == 0)break;
            bit_matrix_set_bit(bm, i, index, true);
        }
    }
    return bm;
}

bool is_hypergram_valid(hypergram ccs)
{
    print("number of points = %ld\nnumber of contexts = %ld\n", ccs.cpt_points, ccs.cpt_geometries);

    if (!is_gram_matrix_valid(ccs.commutation_matrix)){
        print("error : the gram matrix is not valid\n");
        return false;
    }

    bit_matrix incidence = incidence_matrix(ccs);

    bit_matrix product = bit_matrix_product(incidence,ccs.commutation_matrix);

    //The hypergram is valid iff H.G = 0

    if(!bit_matrix_is_empty(product)){
        print("error : the hypergram is not asssignable\n");
        return false;
    }
    
    return true;
}

bit_matrix parse_gram_matrix(FILE* f){
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    size_t cpt = 0;
    bit_matrix bm = (bit_matrix){0};

    while ((read = getline(&line, &len, f)) != -1) {
        if(cpt == 0){
            //here we count the number of lines to allocate the matrix
            size_t size = 0;
            for (ssize_t i = 0; i < read; i++)
                if (line[i] == '0' || line[i] == '1')size++;
            bm = bit_matrix_create(size+1,size+1);
            print("matrix size : %ldx%ld\n",bm.size,bm.bit_sets[0].size);
        }
        //now read every number between commas (either 0 or 1)
        size_t index = 0;
        for (ssize_t i = 0; i < read; i++)
        {
            char c = line[i];
            if(c != '0' && c != '1')continue;
            bit_matrix_set_bit(bm,cpt+1,index+1,(c-'0'));
            index++;
        }
        cpt++;
        if(cpt == bm.size)break;
    }
    free(line);
    return bm;
}

size_t** parse_geometries(FILE* f,size_t* cpt_geometries,size_t* points_per_geometry,size_t* cpt_points){
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    size_t cpt = 0;

    *cpt_points = 0;
    *cpt_geometries = 0;
    *points_per_geometry = 0;
    /*finding the line with the most numbers (counting the commas)*/
    while ((read = getline(&line, &len, f)) != -1) {
        size_t cpt_commas = 1;
        for (ssize_t i = 0; i < read; i++)if(line[i] == ',')cpt_commas++;
        if(cpt_commas > *points_per_geometry)*points_per_geometry = cpt_commas;
        (*cpt_geometries)++;
        
    }
    print("cpt_geometries = %ld ; points_per_geometry = %ld\n",*cpt_geometries,*points_per_geometry);
    size_t** geometries = (size_t**)init_matrix(*cpt_geometries,*points_per_geometry +1,sizeof(size_t));

    /*now we can fill the geometries by using strtok*/
    fseek(f,0,SEEK_SET);

    while ((read = getline(&line, &len, f)) != -1) {
        size_t index = 0;
        char* token = strtok(line,",");
        while (token != NULL)
        {
            int n = atoi(token);
            geometries[cpt][index++] = n;
            if(n > (int)*cpt_points)*cpt_points = n;
            token = strtok(NULL,",");
        }
        cpt++;
        
    }

    free(line);
    return geometries;
}

void hypergram_free(hypergram ccs)
{
    free_matrix(ccs.geometries);
    bit_matrix_free(ccs.commutation_matrix);
    free(ccs.assignment);
}

hypergram hypergram_create_from_file(FILE *geometries_file, FILE *gram_file){
    hypergram res = {
        .cpt_points = 0,
        .n_qubits = 0
    };

    res.geometries = parse_geometries(geometries_file, &res.cpt_geometries, &res.max_points_per_geometry, &res.cpt_points);
    res.commutation_matrix = parse_gram_matrix(gram_file);

    if (!is_hypergram_valid(res))
    {
        print("error : the hypergram is not valid\n");
        hypergram_free(res);
        return (hypergram){0};
    }

    return res;
}

/**
 * @brief creates a hypergram from a gram file, WITH an assignment
 * 
 * @param gram_file 
 * @return hypergram 
 */
hypergram hypergram_create_from_gram_file(FILE *gram_file){

    bit_matrix bm = parse_gram_matrix(gram_file);

    quantum_assignment qa = gram_matrix_to_quantum_assignment(bm);

    hypergram res = quantum_assignment_to_hypergram(qa);
    

    if (!is_hypergram_valid(res))
    {
        print("error : the hypergram is not valid\n");

        hypergram_free(res);
        return (hypergram){0};
    }

    print_quantum_assignment(&qa);

    return res;
}

void hypergram_print(hypergram ccs) {
    print("geometries : \n");
    for (size_t i = 0; i < ccs.cpt_geometries; i++){
        for (size_t j = 0; j < ccs.max_points_per_geometry && ccs.geometries[i][j] != 0; j++){
            print("%ld ",ccs.geometries[i][j]);
        }
        print("\n");
    }
    print("commutation matrix : \n");
    bit_matrix_print(ccs.commutation_matrix);
    if(ccs.assignment == NULL)return;
    print("assignment : \n");
    for (bv i = I; i <= ccs.cpt_points; i++)print_BV_custom(ccs.assignment[i],ccs.n_qubits);
}

size_t quantum_assignment_compute_assignment(quantum_assignment qa,bv** res)
{
    (*res) = calloc(BV_LIMIT_CUSTOM(qa.n_qubits),sizeof(bv));

    //we will use this bit_vector to keep track of the observables we have already seen
    bit_vector bs = bit_set_create(BV_LIMIT_CUSTOM(qa.n_qubits),NULL);

    for (size_t i = 0; i < qa.cpt_geometries; i++){
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            bit_set_set_bit(bs,qa.geometries[qa.geometry_indices[i]][j],true);
        }
    }
    //we will now fill the res array with the observables we have seen (set bits in the bit_vector)
    size_t index = 1;
    for (bv i = I; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++){
        if (bit_set_get_bit(bs, i))(*res)[index++] = i;
    }

    bit_set_free(bs);

    return index;
}


hypergram quantum_assignment_to_hypergram(quantum_assignment qa) {
    hypergram res = {
        .geometries = (size_t**)init_matrix(qa.cpt_geometries,qa.points_per_geometry,sizeof(size_t)),
        .cpt_geometries = qa.cpt_geometries,
        .max_points_per_geometry = qa.points_per_geometry,
        .n_qubits = qa.n_qubits
    };

    
    res.cpt_points = quantum_assignment_compute_assignment(qa, &res.assignment);
    size_t* reversed_assignment = calloc(BV_LIMIT_CUSTOM(qa.n_qubits),sizeof(size_t));
    for (size_t i = 0; i < res.cpt_points; i++)reversed_assignment[res.assignment[i]] = i;

    for (size_t i = 0; i < qa.cpt_geometries; i++){
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            res.geometries[i][j] = reversed_assignment[qa.geometries[qa.geometry_indices[i]][j]];
        }
    }

    //computing the commutation relations of the commutation matrix from the assignment
    res.commutation_matrix = bit_matrix_create(res.cpt_points,res.cpt_points);
    
    print("number of points = %ld\n",res.cpt_points);
    print("matrix size : %ldx%ld\n",res.commutation_matrix.size,res.commutation_matrix.bit_sets[0].size);

    for (bv i = 0; i < res.cpt_points; i++){
        for (bv j = 0; j < res.cpt_points; j++){
            bit_matrix_set_bit(res.commutation_matrix, i, j, innerProduct_custom(res.assignment[i], res.assignment[j], qa.n_qubits));
        }
    }

    free(reversed_assignment);
    return res;
}



bool gram_matrix_find_set_bit(bit_matrix bm,bv* assignment, int n_qubits, size_t *pivot_i, size_t *pivot_j){
    for (size_t i = 0; i < bm.size; i++){
        for (size_t j = 0; j < i; j++){
            if (bit_matrix_get_bit(bm, i, j) ^ innerProduct_custom(assignment[i], assignment[j], n_qubits)){
                *pivot_i = i;
                *pivot_j = j;
                return true;
            }
        }
    }
    return false;
}


int pauli_assignment_from_anticommutations(bit_matrix gram_matrix, bv *assignment)
{
    //clean the assignment
    for (size_t i = 0; i < gram_matrix.size; i++)assignment[i] = I;
    size_t n_qubits_max = MIN((gram_matrix.size / 2) + 1,15);

    size_t ith_qubit=0;
    size_t pivot_i = 0, pivot_j = 0;
    while (gram_matrix_find_set_bit(gram_matrix, assignment, n_qubits_max, &pivot_i, &pivot_j))
    {
        
        for (size_t j = 0; j < gram_matrix.size; j++)
        { /* pivot_i = Z,pivot_j = X*/
            bv qubit = I;

            bv bv_pivot_i = assignment[pivot_i];
            bv bv_pivot_j = assignment[pivot_j];
            bv bv_j = assignment[j];

            if (bit_matrix_get_bit(gram_matrix, j, pivot_i) ^ innerProduct_custom(bv_j, bv_pivot_i, n_qubits_max)){
                qubit = qubit Qplus X;
            }
            if (bit_matrix_get_bit(gram_matrix, j, pivot_j) ^ innerProduct_custom(bv_j, bv_pivot_j, n_qubits_max)){
                qubit = qubit Qplus Z;
            }
            
            
            assignment[j] = set_gate_custom(assignment[j], qubit, n_qubits_max - ith_qubit - 1, n_qubits_max);
            
        }
        
        ith_qubit++;
        if(ith_qubit >= n_qubits_max){
            print("Error: too many qubits needed for the assignment\n");exit(0);
        }
    }

    for (size_t i = 0; i < gram_matrix.size; i++){
        assignment[i] = extend_bv(assignment[i], n_qubits_max, ith_qubit);
        
    }
    
    return ith_qubit;
}

void hypergram_compute_assignment(hypergram *ccs){
    ccs->assignment = calloc(ccs->cpt_points+1,sizeof(bv));
    ccs->n_qubits = pauli_assignment_from_anticommutations(ccs->commutation_matrix, ccs->assignment);
}

quantum_assignment hypergram_to_quantum_assignment(hypergram ccs)
{

    quantum_assignment res = {
        .geometries = (bv **)init_matrix(ccs.cpt_geometries, ccs.max_points_per_geometry, sizeof(bv)),
        .cpt_geometries = ccs.cpt_geometries,
        .points_per_geometry = ccs.max_points_per_geometry};
    if (ccs.assignment == NULL)hypergram_compute_assignment(&ccs);
    res.n_qubits = ccs.n_qubits;
    
    for (size_t i = 0; i < ccs.cpt_geometries; i++){
        for (size_t j = 0; j < ccs.max_points_per_geometry && ccs.geometries[i][j] != 0; j++){
            res.geometries[i][j] = ccs.assignment[ccs.geometries[i][j]];
        }
    }

    quantum_assignment_autofill_indices(&res);
    quantum_assignment_compute_negativity(&res);
    
    return res;
}


void iterate_contexts_rec(hypergram* ccs, size_t ctx_size, size_t point_index, size_t* capacity)
{
    if (point_index >= ccs->cpt_points || is_done)return; // no more points to consider

    // Include the current point in the context iff it does commute with all points in the current context
    for (size_t i = 0; i < ctx_size; i++)
    {
        size_t point = ccs->geometries[ccs->cpt_geometries][i];
        //print("checking point %ld with index %ld\n", point, point_index);
        if (bit_matrix_get_bit(ccs->commutation_matrix, point, point_index))return; // if the point does not commute with the current context, we cannot include it
    }

    // Include the current point in the context
    ccs->geometries[ccs->cpt_geometries][ctx_size] = point_index; // store the point index (1-based)
    ctx_size++;


    // print ("context found: ");
    //for (size_t i = 0; i < ctx_size; i++)print("%ld ", ccs->geometries[ccs->cpt_geometries][i]);
    //print("\n");

    // check that the vector sum of the context is 0
    bit_vector context_sum = bit_set_create(ccs->cpt_points, NULL);
    for (size_t i = 0; i < ctx_size; i++)
    {
        size_t point = ccs->geometries[ccs->cpt_geometries][i];

        bit_set_apply_op(&context_sum, ccs->commutation_matrix.bit_sets[point], '^');
    }
    if (bit_set_is_empty(context_sum) && ctx_size > 0)
    {
        ccs->geometries[ccs->cpt_geometries][ctx_size] = 0; // terminate the context
        //copy into the next context in case another greater context is found
        for (size_t i = 0; i < ctx_size + 1; i++)
        {
            ccs->geometries[ccs->cpt_geometries + 1][i] = ccs->geometries[ccs->cpt_geometries][i];
        }
        ccs->cpt_geometries++;
        if (ccs->cpt_geometries+1 >= *capacity)print("capacity reached !\n");
        
        if (ctx_size > ccs->max_points_per_geometry)ccs->max_points_per_geometry = ctx_size;
    }
    bit_set_free(context_sum);

    for (size_t i = point_index + 1; i < ccs->cpt_points; i++)
    {
        // Recur with the next point
        iterate_contexts_rec(ccs, ctx_size, i, capacity);
    }
}

//same as old one but recursive
hypergram build_full_CCS_fast(bit_matrix bm)
{

    size_t capacity = bm.size * 10000; // initial capacity for geometries, can be adjusted

    hypergram ccs = {
        .geometries = (size_t **)init_matrix(capacity, bm.size + 1, sizeof(size_t)),
        .cpt_geometries = 0,
        .cpt_points = bm.size,
        .n_qubits = bm.size,//just in case
        .max_points_per_geometry = 0,
        .commutation_matrix = copy_bit_matrix(bm),
        .assignment = calloc(bm.size + 1, sizeof(bv)),

    };

    bit_matrix_print(bm);

    
    for (size_t i = 1; i <= bm.size; i++)
    {
        iterate_contexts_rec(&ccs, 0, i,&capacity);
    }

    // increment the every non-empty context by 1 to make it 1-based
    for (size_t i = 0; i < ccs.cpt_geometries; i++)
    {
        for (size_t j = 0; j < ccs.max_points_per_geometry && ccs.geometries[i][j] != 0; j++)
        {
            ccs.geometries[i][j]++;
        }
    }

    return ccs;
}


quantum_assignment gram_matrix_to_quantum_assignment(bit_matrix gram){
    hypergram ccs = build_full_CCS_fast(gram);

    for (size_t i = 0; i < ccs.cpt_geometries; i++){
        for (size_t j = 0; j < ccs.max_points_per_geometry; j++){
            if(ccs.geometries[i][j] == 0)break;
            //print("%ld ", ccs.geometries[i][j]);
        }
        //print("\n");
    }

    ccs.n_qubits = pauli_assignment_from_anticommutations(ccs.commutation_matrix, &(ccs.assignment[1]));
    
    
    
    quantum_assignment qa = hypergram_to_quantum_assignment(ccs);
    
    
    //remove any I from each context
    for (size_t i = 0; i < qa.cpt_geometries; i++){
        size_t shift = 0;
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            if (qa.geometries[qa.geometry_indices[i]][j] == I)shift++;
            if(j+shift >= qa.points_per_geometry){
                if(j < qa.points_per_geometry)qa.geometries[qa.geometry_indices[i]][j] = I;
                continue;
            }

            qa.geometries[qa.geometry_indices[i]][j] = qa.geometries[qa.geometry_indices[i]][j + shift];
        }

    }

    //remove empty contexts
    size_t shift = 0;
    for (size_t i = 0; i + shift < qa.cpt_geometries; i++){
        while (i + shift < qa.cpt_geometries && qa.geometries[i + shift][0] == I)
        {
            shift++;
        }
        if(i+shift < qa.cpt_geometries){
            qa.geometry_indices[i] = qa.geometry_indices[i + shift];
        }
        
        
        
        
        
        
    }
    qa.cpt_geometries -= shift;
    
    
    
    
    /*test all products*/
    for(size_t i = 0; i < qa.cpt_geometries; i++){
        for(size_t j = 0; j < qa.points_per_geometry; j++){
            bv bv1 = I;
            for(size_t k = 0; k < qa.points_per_geometry; k++){
                if(k == I)break;
                bv1 = bv1 Qplus qa.geometries[i][k];
            }
            if(bv1 != I){
                print("error : %d != I\n",bv1);
                exit(1);
            }
        }
    }
    //print(":");
    hypergram_free(ccs);
    return qa;
}


