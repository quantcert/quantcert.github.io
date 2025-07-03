/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file cayley_hexagon.c
 * @brief
 *
 * Generates all Cayley hexagons, their skew embeddings and their complements,
 * and checks their contextuality degree (24 for the skew complements and 0 for all
 * the other ones)
 *
 * Uses equations from the paper [HBS22]. This code was used for the paper [SHKMGD23].
 *
 * References
 * [HBS22] Frédéric Holweck, Henri de Boutray, and Metod Saniga. Three‑qubit‑embedded split Cayley hexagon is contextuality sensitive. Scientific Reports 12 (2022), no. 1, 8915
 * [SHKMGD23] Metod Saniga, Frédéric Holweck, Colm Kelleher, Axel Muller, Alain Giorgetti and Henri de Boutray. Classically-embedded split Cayley hexagons rule three-qubit contextuality with three-element contexts. https://arxiv.org/abs/2312.07738. 2023.
 *

 *
 */

#include "cayley_hexagon.h"

#include "constants.h"
#include "complex_int.h"
#include "quadrics.h"
#include "hashset.h"
#include "bit_vector.h"
#include "config_checker.h"

bv add_seventh_coord(bv bv1){
    bool seventh_coord = bit_parity(get_X(bv1,N_QUBITS_HEX) Qtimes get_Z(bv1,N_QUBITS_HEX),N_QUBITS_HEX);
    return bv1 | (seventh_coord << (2*N_QUBITS_HEX));/*7th coord placed at the end*/ 
}

bool com(bv x,bv y,int i,int j){
    return (BGET(x,i) Qtimes BGET(y,j)) Qplus (BGET(x,j) Qtimes BGET(y,i));
}

/**
 * @brief checks if the 6 points of a hexagon are aligned
 * using the corresponding Plücker equation
 * 
 * equation 13 from [HBS22]
 * 
 * @param bv1 
 * @param bv2 
 * @return true 
 * @return false 
 */
bool points_aligned(bv bv1,bv bv2){
    bv bv71 = add_seventh_coord(bv1);
    bv bv72 = add_seventh_coord(bv2);
    return 
        com(bv71,bv72,5,1) == com(bv71,bv72,0,6) &&
        com(bv71,bv72,0,2) == com(bv71,bv72,6,1) &&
        com(bv71,bv72,1,3) == com(bv71,bv72,2,6) &&
        com(bv71,bv72,2,4) == com(bv71,bv72,6,3) &&
        com(bv71,bv72,3,5) == com(bv71,bv72,4,6) &&
        com(bv71,bv72,4,0) == com(bv71,bv72,6,5) &&
        (com(bv71,bv72,0,3) Qplus com(bv71,bv72,1,4) Qplus com(bv71,bv72,2,5)) == 0;
}

/**
 * @brief finds the index of the line passing through 3 points
 * 
 * @param a 
 * @param b 
 * @param c 
 * @param lines_indices 
 * @return size_t 
 */
size_t get_line_index(bv a,bv b,bv c,size_t** lines_indices){
    for (size_t i = 0; i < NB_LINES_PER_POINT_CUSTOM(N_QUBITS_HEX); i++){
        for (size_t j = 0; j < NB_LINES_PER_POINT_CUSTOM(N_QUBITS_HEX); j++){
            if(lines_indices[a][i] == lines_indices[b][j]){
                for (size_t k = 0; k < NB_LINES_PER_POINT_CUSTOM(N_QUBITS_HEX); k++){
                    if(lines_indices[b][j] == lines_indices[c][k]){
                        return lines_indices[c][k];
                    }
                }
                break;
            }
        }
    }
    print("\nerror line:");
    print_BV_custom(a,N_QUBITS_HEX);
    print_BV_custom(b,N_QUBITS_HEX);
    print_BV_custom(c,N_QUBITS_HEX);
    return 0;
}

/**
 * @brief computes the coordinate of a skew embedding of a split Cayley hexagon
 * from a classical one using this map
 * 
 * equation 15 from [HBS22]
 * 
 * @param lines_indices 
 * @return true 
 * @return false 
 */
bv epsilon(bv bv1){
    bv bv7 = add_seventh_coord(bv1);
    bool res_tab[7]= {0};
    res_tab[0] = BGET(bv7,0) Qplus BGET(bv7,5) Qplus (BGET(bv7,3) Qtimes BGET(bv7,5)) Qplus (BGET(bv7,6) Qtimes BGET(bv7,4));
    res_tab[1] = BGET(bv7,1) Qplus BGET(bv7,2) Qplus (BGET(bv7,2) Qtimes BGET(bv7,4)) Qplus (BGET(bv7,6) Qtimes BGET(bv7,3));
    res_tab[2] = BGET(bv7,2);
    res_tab[3] = BGET(bv7,3);
    res_tab[4] = BGET(bv7,4);
    res_tab[5] = BGET(bv7,5);
    //res_tab[6] = BGET(bv7,6);
    bv res = 0;
    for (size_t i = 0; i < 7; i++)if(res_tab[i])BSET(res,i);
    return res;
}

/**
 * @brief computes the lines indices of embeddings of a split Cayley hexagon and its complement
 * 
 * @param permut 
 * @param set 
 * @param complement 
 */
void indices_arrays_from_sets(hash_set_bitset permut,size_t* set,size_t* complement){
    size_t set_index = 0,complement_index = 0;

    for (size_t i = 1; i <= NB_LINES_CUSTOM(N_QUBITS_HEX); i++)
    {
        if (hash_set_bitset_get(permut, i))
        {
            set[set_index] = i;
            set_index++;
        }
        else
        {
            complement[complement_index] = i;
            complement_index++;
        }
    }
    if(set_index != NB_LINES_CAYLEY_HEXAGON || complement_index != NB_LINES_CUSTOM(N_QUBITS_HEX) - NB_LINES_CAYLEY_HEXAGON){
        print("error ars : %ld %ld\n",set_index,complement_index);
    }
}


hash_set _hexagon_set = {0};

bool next_classical_cayley_hexagon(quantum_assignment lines_qa,size_t **lines_indices, quantum_assignment *qa, quantum_assignment *complement_qa)
{
    static int cpt_copies = 0;

    static bv i = 0, j = 0;

    static bv hexagon_indices[NB_LINES_CAYLEY_HEXAGON];

    if (lines_qa.cpt_geometries != NB_LINES_CUSTOM(N_QUBITS_HEX))
    {
        print("incorrect lines for hexagons\n");
        return false;
    }

    /*if the function is called for the first time, we initialize the structures*/
    if (_hexagon_set.list == NULL){

        *qa = (quantum_assignment){
            .geometries = lines_qa.geometries,
            .cpt_geometries = NB_LINES_CAYLEY_HEXAGON,
            .points_per_geometry = NB_POINTS_PER_LINE,
            .n_qubits = N_QUBITS_HEX};
        
        qa->geometry_indices = calloc(NB_LINES_CAYLEY_HEXAGON,sizeof(size_t));

        *complement_qa = (quantum_assignment){
            .geometries = lines_qa.geometries,
            .cpt_geometries = NB_LINES_CUSTOM(N_QUBITS_HEX) - NB_LINES_CAYLEY_HEXAGON,
            .points_per_geometry = NB_POINTS_PER_LINE,
            .n_qubits = N_QUBITS_HEX};

        complement_qa->geometry_indices = calloc(NB_LINES_CUSTOM(N_QUBITS_HEX) - NB_LINES_CAYLEY_HEXAGON,sizeof(size_t));

        /*We first initialize the hexagon formed by the equation in [HBS22]*/
        size_t cpt = 0;
        for (size_t i = 1; i <= NB_LINES_CUSTOM(N_QUBITS_HEX); i++){
            bv a = lines_qa.geometries[i][0], b = lines_qa.geometries[i][1], c = lines_qa.geometries[i][2];
            if (points_aligned(a, b) && points_aligned(b, c) && points_aligned(a, c)){
                hexagon_indices[cpt] = i;

                cpt++;
            }
        }
        _hexagon_set = hash_set_init();
    }
    
    /*Using at least 2 transvections, it is possible to generate all the 120 classical embeddings of hexagons*/
    for (; i < BV_LIMIT_CUSTOM(N_QUBITS_HEX); i++){
        for (; j < BV_LIMIT_CUSTOM(N_QUBITS_HEX); j++){
            hash_set_bitset permut = {0};

            for (size_t k = 0; k < NB_LINES_CAYLEY_HEXAGON; k++)
            {
                bv a = lines_qa.geometries[hexagon_indices[k]][0], b = lines_qa.geometries[hexagon_indices[k]][1], c = lines_qa.geometries[hexagon_indices[k]][2];

                a = transvection(j, transvection(i, a, N_QUBITS_HEX), N_QUBITS_HEX);
                b = transvection(j, transvection(i, b, N_QUBITS_HEX), N_QUBITS_HEX);
                c = transvection(j, transvection(i, c, N_QUBITS_HEX), N_QUBITS_HEX);
                hash_set_bitset_add(&permut, get_line_index(a, b, c, lines_indices));
            }
            /*We check if the hexagon has already been generated*/
            if(hash_set_exists(&_hexagon_set, permut))continue;

            cpt_copies++;

            indices_arrays_from_sets(permut, qa->geometry_indices, complement_qa->geometry_indices);
            free(qa->lines_negativity);
            free(complement_qa->lines_negativity);
            qa->lines_negativity = NULL;
            complement_qa->lines_negativity = NULL;
            quantum_assignment_compute_negativity(qa);
            quantum_assignment_compute_negativity(complement_qa);

            print("%d,", cpt_copies);

            return true;
            
        }
        j = 0;
    }

    return false;
}

hash_set _skew_hexagon_set = {0};

bool next_skew_cayley_hexagon(quantum_assignment lines_qa, size_t **lines_indices, quantum_assignment *qa, quantum_assignment *complement_qa)
{
    static int skew_copies = 0;

    static bv i = 0, j = 0, k = 0, l = 0;

    static bv skex_hexagon_indices[NB_LINES_CAYLEY_HEXAGON];

    if(lines_qa.cpt_geometries != NB_LINES_CUSTOM(N_QUBITS_HEX)){
        print("incorrect lines for hexagons\n");
        return false;
    }

    /*if the function is called for the first time, we initialize the structures*/
    if (_skew_hexagon_set.list == NULL){

        *qa = (quantum_assignment){
            .geometries = lines_qa.geometries,
            .cpt_geometries = NB_LINES_CAYLEY_HEXAGON,
            .points_per_geometry = NB_POINTS_PER_LINE,
            .n_qubits = N_QUBITS_HEX};

        qa->geometry_indices = calloc(NB_LINES_CAYLEY_HEXAGON,sizeof(size_t));

        *complement_qa = (quantum_assignment){
            .geometries = lines_qa.geometries,
            .cpt_geometries = NB_LINES_CUSTOM(N_QUBITS_HEX) - NB_LINES_CAYLEY_HEXAGON,
            .points_per_geometry = NB_POINTS_PER_LINE,
            .n_qubits = N_QUBITS_HEX};

        complement_qa->geometry_indices = calloc(NB_LINES_CUSTOM(N_QUBITS_HEX) - NB_LINES_CAYLEY_HEXAGON,sizeof(size_t));

        /*We first initialize the skew embedding of the first hexagon formed by the equation in [HBS22]*/
        size_t cpt = 0;
        for (size_t i = 1; i <= NB_LINES_CUSTOM(N_QUBITS_HEX); i++){
            bv a = lines_qa.geometries[i][0], b = lines_qa.geometries[i][1], c = lines_qa.geometries[i][2];
            if (points_aligned(a, b) && points_aligned(b, c) && points_aligned(a, c)){

                bv sa = epsilon(a), sb = epsilon(b), sc = epsilon(c);
                skex_hexagon_indices[cpt] = get_line_index(sa, sb, sc, lines_indices);

                cpt++;
            }
        }

        print("cpt:%ld\n",cpt);

        _skew_hexagon_set = hash_set_init();
    }

    /*Using at least 4 transvections, it is possible to generate all the 7560 skew embeddings of hexagons*/
    for (/*i=0*/; i < BV_LIMIT_CUSTOM(N_QUBITS_HEX); i++){
        for (/*j=i*/; j < BV_LIMIT_CUSTOM(N_QUBITS_HEX); j++){
            for (/*k=j*/; k < BV_LIMIT_CUSTOM(N_QUBITS_HEX); k++){
                for (/*l=k*/; l < BV_LIMIT_CUSTOM(N_QUBITS_HEX); l++){

                    hash_set_bitset permut = {0};

                    for (size_t m = 0; m < NB_LINES_CAYLEY_HEXAGON; m++)
                    {
                        bv a = lines_qa.geometries[skex_hexagon_indices[m]][0], b = lines_qa.geometries[skex_hexagon_indices[m]][1], c = lines_qa.geometries[skex_hexagon_indices[m]][2];

                        a = transvection(l, transvection(k, transvection(j, transvection(i, a, N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX);
                        b = transvection(l, transvection(k, transvection(j, transvection(i, b, N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX);
                        c = transvection(l, transvection(k, transvection(j, transvection(i, c, N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX), N_QUBITS_HEX);
                        hash_set_bitset_add(&permut, get_line_index(a, b, c, lines_indices));
                    }
                    /*We check if the hexagon has already been generated*/
                    if (hash_set_exists(&_skew_hexagon_set, permut))continue;
                    
                    skew_copies++;
                    indices_arrays_from_sets(permut, qa->geometry_indices, complement_qa->geometry_indices);
                    free(qa->lines_negativity);
                    free(complement_qa->lines_negativity);
                    qa->lines_negativity = NULL;
                    complement_qa->lines_negativity = NULL;
                    quantum_assignment_compute_negativity(qa);
                    quantum_assignment_compute_negativity(complement_qa);
                    
                    return true;
                    
                    
                }
                l = k;
            }
            k = j;
        }
        j = i;
    }

    return false;
}

void free_cayley_hexagons()
{
    hash_set_free(&_hexagon_set);
    _hexagon_set.list = NULL;
    hash_set_free(&_skew_hexagon_set);
    _skew_hexagon_set.list = NULL;
}

