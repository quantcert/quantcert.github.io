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

#ifndef _SUB_GEO
#define _SUB_GEO

#include "constants.c"
#include "contextuality_degree.c"
#include "quadrics.c"
#include "bitset.c"

/**
 * @brief Context-commutation structure
 * Represents a geometry from which the contextuality degree can be computed
 */
typedef struct hypergram {
    size_t** geometries;/*array of geometries 0 is the end of an array*/
    size_t cpt_geometries;
    size_t cpt_points;
    size_t max_points_per_geometry;
    bit_matrix commutation_matrix;
    bv* assignment;
    size_t n_qubits;
} hypergram;

/**
 * @return true if the gram matrix is valid
 * i.e. it is square and symmetric
 */
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

/**
 * @brief Returns the incidence matrix of the hypergraph of the contexts of a hypergram
 *
 * @param obss Set of observables, will be sorted as an output
 * @param size Number of observables
 * @param basis (output) the basis
 *
 * @return int Cardinal of the basis
 */
bit_matrix incidence_matrix(hypergram ccs)
{

    bit_matrix bm = bit_matrix_create(ccs.cpt_geometries, ccs.cpt_points+1);

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

/**
 * @brief Checks if a hypergram is valid
 * 
 * @return true If the hypergram is valid
 * i.e. the commutation matrix is valid and the geometries are valid
 */
bool is_hypergram_valid(hypergram ccs)
{
    print("number of points = %ld\nnumber of geometries = %ld\n", ccs.cpt_points, ccs.cpt_geometries);

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

/**
 * @brief Parses a gram matrix from a file (csv style)
 * 
 * ex:
 * 0,1,1
 * 1,0,1
 * 1,1,0
 * 
 * or
 * 
 * 011
 * 101
 * 110
 * 
 * @param f 
 * @return bit_matrix 
 */
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
            bit_matrix_set_bit(bm,cpt+1,index+1,c-'0');
            index++;
        }
        cpt++;
        if(cpt == bm.size)break;
    }
    free(line);
    return bm;
}

/**
 * @brief Parses a hypergraph from a file
 * 
 * ex:
 * 1,2,3
 * 4,5,6
 * 7,8,9
 * 
 * no zero allowed
 * 
 * @param f 
 * @return CCS 
 */
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

/**
 * @brief Frees the memory allocated for a hypergram
 *
 * @param ccs Hypergram to free
 */
void hypergram_free(hypergram ccs)
{
    free_matrix(ccs.geometries);
    bit_matrix_free(ccs.commutation_matrix);
    free(ccs.assignment);
}

/**
 * @brief Creates a hypergram from a two files
 * The first file contains the geometries
 * ex:
 * 1,2,3
 * 4,5,6
 * 7,8,9
 *
 * The second file contains the gram matrix
 * ex:
 * 011
 * 101
 * 110
 *
 * @param geometries_file
 * @param gram_file
 * @return hypergram or (hypergram){0} if the hypergram is not valid
 */
hypergram hypergram_create_from_file(FILE *geometries_file, FILE *gram_file){
    hypergram res = {
        .geometries = parse_geometries(geometries_file, &res.cpt_geometries, &res.max_points_per_geometry, &res.cpt_points),
        .cpt_points = 0,
        .n_qubits = 0,
        .commutation_matrix = parse_gram_matrix(gram_file)
    };

    if(!is_hypergram_valid(res)){
        print("error : the hypergram is not valid\n");
        hypergram_free(res);
        return (hypergram){0};
    }

    return res;
}

/**
 * @brief Prints a hypergram
 * 
 * @param ccs Hypergram to print
 */
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

/**
 * @brief Creates a list of all the observables present in a quantum assignment object
 * 
 * @param qa 
 * @param res 
 * @return size_t 
 */
size_t quantum_assignment_compute_assignment(quantum_assignment qa,bv** res)
{
    (*res) = calloc(BV_LIMIT_CUSTOM(qa.n_qubits),sizeof(bv));

    //we will use this bitset to keep track of the observables we have already seen
    bit_set bs = bit_set_create(BV_LIMIT_CUSTOM(qa.n_qubits),NULL);

    for (size_t i = 0; i < qa.cpt_geometries; i++){
        for (size_t j = 0; j < qa.points_per_geometry; j++){
            bit_set_set_bit(bs,qa.geometries[qa.geometry_indices[i]][j],true);
        }
    }
    //we will now fill the res array with the observables we have seen (set bits in the bitset)
    size_t index = 1;
    for (bv i = I; i < BV_LIMIT_CUSTOM(qa.n_qubits); i++){
        if (bit_set_get_bit(bs, i))(*res)[index++] = i;
    }

    bit_set_free(bs);

    return index;
}


/**
 * @brief Transforms a quantum assignment into a hypergram
 */
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



/**
 * @brief Finds a set bit in a gram matrix
 * 
 * @param bm Bit matrix
 * @param pivot_i (output) index of the row of the first set bit
 * @param pivot_j (output) index of the column of the first set bit
 * 
 * @return true if a set bit was found, false otherwise
 */
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


/**
 * @brief Finds the assignment of a CCS from its independent gram matrix
 * by following the algorithm described in the paper
 * 
 * Contrary to the paper, instead of using a mutable matrix alone, we use a static matrix 
 * and the assignment to keep track of the anticommutations
 * 
 * @param gram_matrix Independent gram matrix
 * @param assignment (output) Assignment
 * @return int Number of qubits
 */
int pauli_assignment_from_anticommutations(bit_matrix gram_matrix, bv *assignment)
{
    //clean the assignment
    for (size_t i = 0; i < gram_matrix.size; i++)assignment[i] = I;
    int n_qubits_max = MIN((gram_matrix.size / 2) + 1,15);

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
    }

    for (size_t i = 0; i < gram_matrix.size; i++){
        assignment[i] = extend_bv(assignment[i], n_qubits_max, ith_qubit);
        // print_BV_custom(assignment[i], ith_qubit + 1);
    }
    // print("\n");
    return ith_qubit;
}

/**
 * @brief Computes the assignment of a hypergram
 * 
 * @param ccs Hypergram
 */
void hypergram_compute_assignment(hypergram *ccs){
    ccs->assignment = calloc(ccs->cpt_points+1,sizeof(bv));
    ccs->n_qubits = pauli_assignment_from_anticommutations(ccs->commutation_matrix, ccs->assignment);
}

/**
 * @brief Transforms a hypergram into a quantum assignment
 */
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
    compute_negativity(&res);
    
    return res;
}



#endif // _SUB_GEO