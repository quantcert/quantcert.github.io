/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file quantum_assignment.c
 * @brief Data structure and functions to handle quantum assignments
 * 
 */


#ifndef QUANTUM_ASSIGNMENT
#define QUANTUM_ASSIGNMENT

#include "bv.c"
#include "complex.c"

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
typedef struct
{
    size_t *geometry_indices;
    bv **geometries;
    size_t cpt_geometries;
    size_t points_per_geometry;
    int n_qubits;

    bool *lines_negativity;
} quantum_assignment;

bool quantum_assignment_autofill_indices(quantum_assignment* qa){
    bool has_no_indices = qa->geometry_indices == NULL;
    if(has_no_indices){
        qa->geometry_indices = calloc(qa->cpt_geometries,sizeof(size_t));
        for (size_t i = 0; i < qa->cpt_geometries; i++)qa->geometry_indices[i] = i;
    }
    return has_no_indices;
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
bool is_negative_custom(bv geometry[], int size, int n_qubits, bool verbose, FILE *output)
{
    pauli_matrix mat = get_matrix(I);

    /*for each qubit*/
    for (int i = 0; i < n_qubits; i++)
    {
        /*for each observable*/

        pauli_matrix printed = get_matrix(I);

        for (int j = size - 1; j >= 0; j--)
        { /*multiplication is performed from right to left*/
            /*we compute the product of the nth qubits of each observable
            (it should be the identity modulo some phase)
            and then we compute the product of all these matrices for all qubits*/
            bv gate = get_gate(geometry[j], i, n_qubits);
            mat = matrix_mult(mat, get_matrix(gate));
            if (verbose || true)
                printed = matrix_mult(printed, get_matrix(gate));
        }

        W2_complex phase = get_id_matrix_phase(printed);
        if (verbose)
        {
            print_complex_to_file(phase, output);
            if (i != n_qubits - 1)
                fprintf(output, " * ");
        }
        // if(phase.imag != 0 || phase.real != 1)print("ERROR : QUBIT PRODUCT IS NOT +I : %d+%di; ",phase.real,phase.imag);//uncomment to check for total positivity
    }
    
    
    /*we check that the phase is correct*/
    W2_complex phase = get_id_matrix_phase(mat);

    if (verbose)
    {
        fprintf(output, " = ");
        print_complex_to_file(phase, output);
    }

    if (phase.imag != 0 || phase.real == 0)
    {
        print("geometry phase error");
    }
    return phase.real == -1;
}
bool is_negative(bv geometry[], int size, int n_qubits) { return is_negative_custom(geometry, size, n_qubits, false, stderr); }

/**
 * @brief Prints a quantum assignment
 * 
 * @param qa 
 */
void print_quantum_assignment(quantum_assignment* qa){
    quantum_assignment_autofill_indices(qa);
    
    print("geometries : \n");
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        for (size_t j = 0; j < qa->points_per_geometry; j++){
            print_BV_custom(qa->geometries[qa->geometry_indices[i]][j],qa->n_qubits);
            //print("%ld ",qa.geometries[qa.geometry_indices[i]][j]);
        }
        print("%c",is_negative_custom(qa->geometries[qa->geometry_indices[i]],qa->points_per_geometry,qa->n_qubits,false,NULL)?'-':'+');
        print("\n");
    }
    print("\n");

}

/**
 * @brief Prints a quantum assignment as a readable configuration to a file
 * 
 * @param qa 
 * @param output 
 */
void quantum_assignment_print_to_file(quantum_assignment* qa, FILE *output)
{
    quantum_assignment_autofill_indices(qa);
    
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        for (size_t j = 0; j < qa->points_per_geometry; j++){
            print_BV_to_file(qa->geometries[qa->geometry_indices[i]][j],qa->n_qubits,output);
            if(j != qa->points_per_geometry-1)fprintf(output,",");
        }
        fprintf(output,"\n");
    }
}

/**
 * @brief Computes the negativite contexts of a quantum assignment
 * 
 * @param qa 
 */
bool compute_negativity(quantum_assignment* qa){

    if(qa->lines_negativity != NULL)return false;

    quantum_assignment_autofill_indices(qa);

    qa->lines_negativity = calloc(qa->cpt_geometries,sizeof(bool));
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        qa->lines_negativity[i] = is_negative_custom(qa->geometries[qa->geometry_indices[i]],qa->points_per_geometry,qa->n_qubits,false,NULL);
    }

    return true;
}

/**
 * @brief Conts the number of negative lines in a quantum assignment
 * 
 * @param qa quantum assignment
 * @return int 
 */
int negative_lines_count(quantum_assignment* qa){
    quantum_assignment_autofill_indices(qa);
    compute_negativity(qa);
    
    int res = 0;
    for (size_t i = 0; i < qa->cpt_geometries; i++){
        if(qa->lines_negativity[i])res++;
    }

    return res;
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
 * @brief merges two quantum assignments
 * 
 * @param qa1 
 * @param qa2
 * @return quantum_assignment 
 */
quantum_assignment quantum_assignment_merge(quantum_assignment qa1,quantum_assignment qa2){
    quantum_assignment res = {0};

    res.points_per_geometry = MAX(qa1.points_per_geometry,qa2.points_per_geometry);
    res.cpt_geometries = qa1.cpt_geometries+qa2.cpt_geometries;
    res.geometries = (bv**)init_matrix(res.cpt_geometries,res.points_per_geometry,sizeof(bv));
    res.n_qubits = qa1.n_qubits;

    for (size_t i = 0; i < qa1.cpt_geometries; i++)for (size_t j = 0; j < qa1.points_per_geometry; j++)res.geometries[i][j] = qa1.geometries[qa1.geometry_indices[i]][j];
    for (size_t i = 0; i < qa2.cpt_geometries; i++)for (size_t j = 0; j < qa2.points_per_geometry; j++)res.geometries[qa1.cpt_geometries+i][j] = qa2.geometries[qa2.geometry_indices[i]][j];
    
    quantum_assignment_autofill_indices(&res);

    print_quantum_assignment(&res);


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

#endif //QUANTUM_ASSIGNMENT