/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file cayley_hexagon.h
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
#ifndef CAYLEY_HEXAGON_C
#define CAYLEY_HEXAGON_C

#include "constants.h"

#include "quadrics.h"
#include "config_checker.h"

#define NB_LINES_CAYLEY_HEXAGON 63
#define N_QUBITS_HEX 3

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
bv epsilon(bv bv1);


/**
 * @brief computes the next classical embeddings of the split cayley hexagon of order 2
 * 
 * @param lines list of all the three qubit lines
 * @param lines_indices list of the indices of the lines
 * @param qa cayley hexagon generated
 * @param complement_qa complement of the cayley hexagon generated
 * @return true if there is a next hexagon
 * @return false when all hexagons have been generated
 */
bool next_classical_cayley_hexagon(quantum_assignment lines_qa, size_t **lines_indices, quantum_assignment *qa, quantum_assignment *complement_qa);

/**
 * @brief computes the next skew embeddings of the split cayley hexagon of order 2
 *
 * @param lines list of all the three qubit lines
 * @param lines_indices list of the indices of the lines
 * @param qa cayley hexagon generated
 * @param complement_qa complement of the cayley hexagon generated
 * @return true if there is a next hexagon
 * @return false when all hexagons have been generated
 */
bool next_skew_cayley_hexagon(quantum_assignment lines_qa, size_t **lines_indices, quantum_assignment *qa, quantum_assignment *complement_qa);

void free_cayley_hexagons();

#endif //CAYLEY_HEXAGON_C
