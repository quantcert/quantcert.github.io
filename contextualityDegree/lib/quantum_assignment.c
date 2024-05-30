#ifndef QUANTUM_ASSIGNMENT
#define QUANTUM_ASSIGNMENT

#include "bv.c"
#include "complex.c"

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

#endif //QUANTUM_ASSIGNMENT