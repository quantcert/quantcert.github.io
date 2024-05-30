/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU Lesser                 */
/* General Public License version 2                                             */
/**********************************************************************************/

/**
 * @file quadrics.c
 * @author Axel Muller
 *
 * [dHG+21] Henri de Boutray, Frédéric Holweck, Alain Giorgetti, Pierre-Alain Masson and Metod Saniga. Contextuality
 * degree of quadrics in multi-qubit symplectic polar spaces. arXiv:2105.13798
 */

#ifndef QUADRICS
#define QUADRICS

#include "bv.c"
#include "contextuality_degree.c"




/**All these formula are described in Contextuality degree of quadrics
in multi-qubit symplectic polar spaces*/
#define NB_HYPERBOLICS(N) ((pow4(N)+pow2(N))/2)
#define NB_ELLIPTICS(N) ((pow4(N)-pow2(N))/2)
#define NB_QUADRICS(N) (NB_HYPERBOLICS(N)+NB_ELLIPTICS(N))
#define NB_PERPSETS(N) (pow4(N)-1)
    int NB_GENERATORS(int N){return (N<=0)?1:(pow2(N)+1)*NB_GENERATORS(N-1);}

/**
 * @brief Counts the number of totally isotropic subspaces
 * up to 8 qubits, before the numbers get too large 
 * (see QPL2023)
 * 
 * @param N number of qubits
 * @param K dimension of the subspace
 * @return size_t 
 */
size_t NB_SUBSPACES(int N,int K){
    size_t a = 1,b = 1;
    for (int i = 1; i <= K+1; i++){
        a *= (pow2(N-K-1+i)-1) * (pow2(N-i+1)+1);
        b *= pow2(i)-1;
    }
    return a/b;
}

#define NB_OBS_PER_HYPERBOLIC(N) (NB_HYPERBOLICS(N)-1)
#define NB_OBS_PER_ELLIPTIC(N) (NB_ELLIPTICS(N)-1)
#define NB_OBS_PER_QUADRIC(N) (NB_OBS_PER_HYPERBOLIC(N)+NB_OBS_PER_ELLIPTIC(N))
#define NB_OBS_PER_PERPSET(N) ((pow4(N)/2)-1)
#define NB_OBS_PER_GENERATOR(N) (pow2(N)-1)

#define NB_LINES_PER_HYPERBOLIC(N) ((NB_OBS_PER_HYPERBOLIC(N))*(NB_OBS_PER_HYPERBOLIC(N-1))/3)
#define NB_LINES_PER_ELLIPTIC(N) ((NB_OBS_PER_ELLIPTIC(N)*(NB_OBS_PER_ELLIPTIC(N-1)))/3)
#define NB_LINES_PER_QUADRIC(N) (NB_LINES_PER_HYPERBOLIC(N)+NB_LINES_PER_ELLIPTIC(N))
#define NB_LINES_PER_PERPSET(N) (pow4((N)-1)-1)

#define NO_LINE (0)

/**all 6 2-spreads as lists of indexes of the W2 doily lines*/
size_t index_two_spreads[NB_TWO_SPREADS_PER_DOILY][NB_LINES_TWO_SPREAD];



/**
 * @brief inserts a line index in the incidence map of observables
 * 
 * @param obs 
 * @param index 
 */
void insert_line_index(bv obs, size_t index,size_t** lines_indices,int n_qubits) {
    int i;

    for (i = 0; lines_indices[obs][i] != NO_LINE && i<NB_LINES_PER_POINT_CUSTOM(n_qubits); i++)
    

    
    if (obs >= BV_LIMIT_CUSTOM(n_qubits) || i >= NB_LINES_PER_POINT_CUSTOM(n_qubits))
        print("error out of bounds %d %d / %d %d [%d,%d]!\n", obs, i, BV_LIMIT_CUSTOM(n_qubits),
            NB_LINES_PER_POINT_CUSTOM(n_qubits),obs,i);

    lines_indices[obs][i] = index;

}

/**
 * @brief generates all the lines for a given number of qubits
 * 
 * @param lines resulting array of lines of observables(The index 0 is NOT a line !)
 * @param lines_indices indicates all the lines each observable belongs to (all lists of lines are sorted)
 * @param n_qubits number of qubits of the geometry
*/
quantum_assignment generate_total_lines(size_t*** lines_indices,int n_qubits) {

    bv ** lines = (bv **)init_matrix(NB_LINES_CUSTOM(n_qubits) + 1, NB_POINTS_PER_LINE, sizeof(bv));
    *lines_indices = (size_t**)init_matrix(BV_LIMIT_CUSTOM(n_qubits),NB_LINES_PER_POINT_CUSTOM(n_qubits),sizeof(size_t));//calloc(BV_LIMIT_CUSTOM(VARQ),sizeof(size_t[NB_LINES_PER_POINT_CUSTOM(VARQ)]));
    if(lines == NULL || *lines_indices == NULL){
        print("no memory !");return (quantum_assignment){0};
    }

    size_t index = 1;/*not 0 since this number represents the absence of lines*/

    /*for every pair of ascending observables*/
    for (bv i = 1; i < BV_LIMIT_CUSTOM(n_qubits); i++) {
        for (bv j = i+1; j < BV_LIMIT_CUSTOM(n_qubits); j++) {
            bv k = i Qplus j;
            if (i > j || j > k || innerProduct_custom(i, j, n_qubits) != 0)
                continue;
            /*if they commute, we build the third one (it automatically commutes with the other ones)*/
            bv list[] = { i, j, k };
            for (size_t j = 0; j < NB_POINTS_PER_LINE; j++) {
                lines[index][j] = list[j];
                insert_line_index(list[j], index,*lines_indices,n_qubits);
            }
            index++;            
        }
        if(PRINT_PROGRESSION && i%500 == 0)print(" %.2f%% ", 100*((float)index/(float)(NB_LINES_CUSTOM(n_qubits) + 1)));
    }
    /*if the number of lines generated is unexpected*/
    if (index != (size_t)NB_LINES_CUSTOM(n_qubits) + 1)
        print("INDEX INCOHERENT ! : %ld != %d\n", index,NB_LINES_CUSTOM(n_qubits) + 1);

    quantum_assignment qa = {
        .geometry_indices = calloc(NB_LINES_CUSTOM(n_qubits), sizeof(size_t)),
        .geometries = lines,
        .cpt_geometries = NB_LINES_CUSTOM(n_qubits),
        .points_per_geometry = NB_POINTS_PER_LINE,
        .n_qubits = n_qubits
    };

    for (size_t i = 0; i < (size_t)NB_LINES_CUSTOM(n_qubits); i++) qa.geometry_indices[i] = i+1;

    return qa;
}

/**
 * @brief returns the index of the leftmost set bit, or -1 if it doesn't exist
 * 
 * @param obs 
 * @return int 
 */
int bv_left_most(bv obs){
    if(obs == 0)return -1;
    for (size_t i = 0; i < 8*sizeof(bv); i++){
        obs >>= 1;
        if(obs == 0)return i;
    }
    return 8*sizeof(bv);//should not happen
}

/**
 * @brief Comparison function to sort a bv array in DESCENDING order
 */
int bv_anti_cmp(const void * first, const void * second ) {
    return * (const bv *) second - * (const bv *) first;
}

/**
 * @brief computes the basis of a set of observables
 * using Gaussian elimination
 * 
 * @param obss set of observables, will be sorted as an output
 * @param size number of observables
 * @param basis (output) the basis
 * 
 * @return int cardinal of the basis
 */
int find_basis(bv *obss,int size,bv** basis){

    bv copy[size];
    for (int i = 0; i < size; i++)copy[i] = obss[i];
    
    qsort(copy,size,sizeof(bv),bv_anti_cmp);

    int current_basis_index = 0;
    /*for each ith element of each vector*/
    while (copy[current_basis_index] != 0)
    {
        bv b = copy[current_basis_index];
        int index = bv_left_most(b);
        
        current_basis_index++;
        for (int j = current_basis_index; j < size; j++)/*we then unset this bit on all other vectors below by adding it*/
        {
            if(BGET(copy[j],index))copy[j] = copy[j] Qplus b;
        }
        int remaining_size = size-current_basis_index-2;
        if(remaining_size>0)qsort(&(copy[current_basis_index]),remaining_size,sizeof(bv),bv_anti_cmp);
    }

    if(basis != NULL){//we copy each non-null vector into the basis list
        *basis = calloc(current_basis_index,sizeof(bv));
        int cpt = 0;
        for (int i = 0; i < size; i++)if(copy[i] != I)(*basis)[cpt++] = copy[i];
    }
    return current_basis_index;
}

/**
 * @brief Get the basis combination of an observable given that 
 * the basis is sorted in descending order
 * 
 * @param obs 
 * @param basis 
 * @param size 
 * @return uint64_t 
 */
uint64_t get_basis_combination(bv obs,bv *basis,int size){
    uint64_t res = 0;
    for (int i = 0; i < size; i++)/*we find the possible basis by checking their leftmost part since they are all unique*/
    {
        if(bv_left_most(obs) == bv_left_most(basis[i])){/*we subtract the element of the matrix when its leftmost part is the same*/
            BSET(res,i);
            obs = obs Qplus basis[i];
            if(obs == 0)return res;
        }
    }
    print("UNEXPECTED BASIS OVERFLOW");
    return 0;
}


bv**    _subspaces_tab;
size_t  _subspaces_current_index;
size_t  _subspaces_limit;
size_t _subspaces_n_negative = 0;

/**
 * @brief Callback function used to get all generators from the commuting function
*/
bool generate_subspace(bv bv1[],int size){
    
    bv tab[NB_OBS_PER_GENERATOR(size)+1];//_subspaces_tab[_subspaces_current_index];

    
    for (bv i = 1; i <= (bv)NB_OBS_PER_GENERATOR(size); i++){
        bv sum = I;
        for (int j = 0; j < size; j++)
        {
            if(BGET(i,j))sum = sum Qplus bv1[j];
        }
        if(sum == I)return false;
        if(i > 1 && tab[i-1] > sum)return false;
        tab[i] = sum;
    }

    bool is_neg = is_negative(&(tab[1]),NB_OBS_PER_GENERATOR(size),N_QUBITS);
    
    #pragma omp critical
    {
        
        _subspaces_current_index++;
        if(is_neg)_subspaces_n_negative ++;
    }

    if(_subspaces_current_index % 2000000 == 1)print("%ld,",_subspaces_current_index);
    if(_subspaces_current_index > _subspaces_limit)print("error subspaces");


    if(DOUBLE_CHECK)return false;
    return _subspaces_current_index == _subspaces_limit;
}

int cptor = 0;
bool cpteur(bv bv1,int n_qubits){
    (void)bv1;
    (void)n_qubits;
    #pragma omp critical
    {cptor++;}
    return false;
}

/**
 * @brief internal function used to generate all mutually commuting elements
 * 
 * @param bv1 list of observable to generate
 * @param index current index to generate observables on
 * @param size size of the observable list
 * @param n_qubits number of qubits per observable
 * @param callback function to call when a list of observables is generated (return false for the program to end)
 * 
 * @return true iff the iteration has to stop
 */
bool commuting_rec(bv bv1[],uint32_t index,uint32_t size,int n_qubits,bool (*callback)(bv[],int)){
    if(is_done)return true;
    if(index == size){
        /*If every observable is generated we call the callback function with the right parameters*/
        bool res;
        res = (*callback)(bv1,size);
        //if(res)print("\n--stop--\n");
        return res;
    }
    /*starting observable*/
    bv start = index==0?I+1:bv1[index-1]+1;

    for (bv i = start; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        bool valid = true;
        for (uint32_t j = 0; j < index; j++)/*we filter out every observable that doesn't commute with all the other ones*/
            if(innerProduct_custom(i,bv1[j],n_qubits) != 0){valid = false;break;}
        if(!valid)continue;
        bv1[index] = i;
        /*if we stop the iterations*/
        if(commuting_rec(bv1,index+1,size,n_qubits,callback))return true;
    }
    return false;
}

/**
 * @brief Generates all possible sets of mutually commuting observables
 * 
 * @param size number of mutually commuting observables
 * @param n_qubits number of qubits per observable
 * @param callback (bool)(array,size) function called every time a set of observables is generated (as an array) and returns false when we want to stop iteration
 */
void commuting(uint32_t size,int n_qubits,bool (*callback)(bv[],int)){
    /*Splitting the first for loop between every thread*/
    int num_done = 0;
    bool stop_iters = false;
    #pragma omp parallel for schedule(dynamic,1) if(size > 2)
    for (size_t i = I+1; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        if(stop_iters)continue;
        bv bv1[size];
        bv1[0] = i;
        if(commuting_rec(bv1,1,size,n_qubits,callback))stop_iters = true;
        #pragma omp critical
        {
        num_done++;
        if (N_QUBITS > 3)
            print(" %.2f%% ", 100*((float)num_done/(float)BV_LIMIT));
        }
    }
}

/**
 * @brief Tests the conjecture stating that any subspace with k >= 4 has for each 
 * of its points an even number of negative lines passing through it
 * 
 * @param n_qubits 
 * @param k 
 * @param bv1 
 */
void subspace_neg_lines_count(int n_qubits,int k,bv bv1[]){
    
    for (bv i = 1; i < (bv)pow2(k+1); i++)//subspaces don't have a 0 index, because it is like the identity observable
    {
        int total_count = 0;
        int neg_count = 0;
        for (bv j = 1; j < (bv)pow2(k+1); j++)
        {
            if((i Qplus j) < j || i == j)continue;// to avoid copies / includes the i = j case
            bv line[] = {bv1[i],bv1[j],bv1[i Qplus j]};
            neg_count += is_negative(line,NB_POINTS_PER_LINE,n_qubits);
            total_count++;
        }
        print("(%d,%d)",total_count,neg_count);
        if(neg_count %2 !=0)print("EVEN POINTS PER LINE CONJECTURE WRONG!");
    }
}



void subspaces_rec(int n_qubits,int k,int l,bv bv1[]){
    if(l == k+1){
        bool is_neg = is_negative(&(bv1[1]),NB_OBS_PER_GENERATOR(k+1),n_qubits);
        if(CHECK_SUBSPACES_LINES_EVEN)subspace_neg_lines_count(n_qubits,k,bv1);
        #pragma omp critical
        {
            for (size_t i = 0; i < (size_t)NB_OBS_PER_GENERATOR(k+1); i++)_subspaces_tab[_subspaces_current_index][i] = bv1[i+1];
            _subspaces_current_index++;
            if(is_neg)_subspaces_n_negative++;
            int last_right_most = -1;
            for (int i = 0; i < k+1; i++){
                //print_BV_custom(bv1[pow2(i)],n_qubits);
                int right_most = 2*n_qubits;
                while(!BGET(bv1[pow2(i)],right_most))right_most--;
                //print("%d",right_most);
                if(right_most <= last_right_most)print("CONJECTURE WRONG");
                last_right_most = right_most;
            }
            
            
            
            

            
            
            
            
        
        }
        if(_subspaces_current_index % 200000 == 1)print("%ld,",_subspaces_current_index);
        return;
    }
    for (bv i = bv1[pow2(l)-1]+1; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        bv1[pow2(l)] = i;

        bool valid = true;
        for (int j = pow2(l)+1; j < pow2(l+1); j++)
        {
            bv1[j] = bv1[j&mask(l)] Qplus i;
            
            /*no need to check commutativity for the previous ones, they form a line with them*/
            if(bv1[j] < bv1[pow2(l)] || innerProduct_custom(i,bv1[j],n_qubits) != 0){
                valid = false;
                break;
            }
        }
        if(!valid)continue;
        
        subspaces_rec(n_qubits,k,l+1,bv1);
    }
}

quantum_assignment subspaces(int n_qubits,int k){

    int size = NB_OBS_PER_GENERATOR(k+1);

    //print("s:%ld %d %ld\n",NB_SUBSPACES(n_qubits,k),size,NB_SUBSPACES(n_qubits,k)*size);
    /*instanciate the global variables accessible in the callback function (+1 for the callback function)*/
    /*                                   nb lines is an upper bound because we don't know the real number*/
    _subspaces_tab = (bv**)init_matrix(NB_SUBSPACES(n_qubits,k)+1,size,sizeof(bv));
    _subspaces_current_index = 0;
    _subspaces_limit = NB_SUBSPACES(n_qubits,k);

    quantum_assignment qa = (quantum_assignment){0};
    qa.geometry_indices = calloc(NB_SUBSPACES(n_qubits,k), sizeof(size_t));
    qa.geometries = _subspaces_tab;
    qa.points_per_geometry = NB_OBS_PER_GENERATOR(k+1);
    qa.cpt_geometries = _subspaces_limit;
    qa.n_qubits = n_qubits;

    /*Splitting the first for loop between every thread*/
    int num_done = 0;
    bool stop_iters = false;
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i = I+1; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        if(stop_iters)continue;
        bv bv1[NB_OBS_PER_GENERATOR(k+1)];
        bv1[0] = I;bv1[1] = i;
        subspaces_rec(n_qubits,k,1,bv1);
        #pragma omp critical
        {num_done++;
        if(PRINT_PROGRESSION)print(" %.2f%% ", 100*((float)num_done/(float)BV_LIMIT));
        }
    }
    for (size_t i = 0; i < NB_SUBSPACES(n_qubits,k); i++) qa.geometry_indices[i] = i;

    return qa;
}

quantum_assignment affine_planes(quantum_assignment planes){
    quantum_assignment ap = {
        .geometries = (bv**)init_matrix(planes.cpt_geometries*7,NB_POINTS_PER_AFFINE,sizeof(bv)),
        .cpt_geometries = planes.cpt_geometries*7,
        .points_per_geometry = NB_POINTS_PER_AFFINE,
        .n_qubits = planes.n_qubits
    };

    quantum_assignment_autofill_indices(&ap);

    for (size_t i = 0; i < planes.cpt_geometries; i++)
    {
        bv* plane = planes.geometries[planes.geometry_indices[i]];
        
        int j = 0;

        
        for (size_t k = 0; k < 7; k++)
        {
            for (size_t l = k+1; l < 7; l++)
            {
                bv c = plane[k] Qplus plane[l];
                if(c > plane[k] && c > plane[l]){
                    bv *affine_plane = ap.geometries[i * 7 + j];

                    int ap_index = 0;
                    for (size_t m = 0; m < 7; m++){
                        if(m == k || m == l || plane[m] == c)continue;
                        affine_plane[ap_index++] = plane[m];
                    }

                    j++;
                }
            }
            
        }
        if (j != 7)print("INCORRECT NUMBER OF POINTS PER AFFINE PLANES! :%d\n",j);
    }
    return ap;
}

/**
 * @brief Computes the contextuality degree of generators 
 * 
 * @param n_qubits number of qubits
 *  
 * @return int contextuality degree of the geometry
 */
int subspaces_contextuality_degree(int n_qubits,int k){


    print("\nGenerating subspaces...(expected: at most %ld * %d)\n\n",NB_SUBSPACES(n_qubits,k),NB_OBS_PER_GENERATOR(k+1));

    quantum_assignment qa = subspaces(n_qubits,k);

    print("\nfirst : %ld\n",_subspaces_current_index);

    print("\n\nGeometries generated((%ld,%d,%d,n_neg=%ld))\n\nnow checking contextuality...\n\n",_subspaces_current_index,NB_OBS_PER_GENERATOR(k+1),n_qubits,_subspaces_n_negative);
    
    bool contextuality_only = (k != 1) && !(n_qubits == 2 && k == 1);
    
    int c_degree = geometry_contextuality_degree(qa,contextuality_only,false,true,NULL);

    free_matrix(_subspaces_tab);

    return c_degree;
}

/**
 * @brief function generating the lines following a given form from a source observable and a given form
 * (even though for perpsets an additionnal condition is that the source observable must belong to each line)
 * 
 * WARNING : the array geometries_indices is allocated on the heap and must be freed after use
 * 
 * @param obs source observable of the geometry
 * @param form function determining if a point belongs to a geometry
 * @param lines_indices array specifying for each point all the lines it belongs to
 * @param n_qubits number of qubits of the geometry
 * @param lines array of all the lines of n qubits
 * 
 * @param lines_res array of all the lines of the resulting geometry
*/
quantum_assignment zero_locus(bv obs, unsigned int (*form)(bv, bv, int),size_t** lines_indices,int n_qubits,bv** lines,size_t* lines_res,bool complement) {
    
    quantum_assignment qa = {
        .geometry_indices = lines_res,
        .geometries = lines,
        .cpt_geometries = 0,
        .points_per_geometry = NB_POINTS_PER_LINE,
        .n_qubits = n_qubits
    };

    int cpt = 0;/*number of points in the geometry*/

    /*number of points of each line satisfying the form*/
    uint_least8_t* line_tab = calloc(NB_LINES_CUSTOM(n_qubits)+1,sizeof(uint_least8_t));

    /*first we count how many times a line contains a point belonging to a geometry*/
    for (bv i = 1; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        if((*form)(obs,i,n_qubits) == 0){
            //nis[n_I(i)]++;
            cpt++;
            /*we increment each the line counter for each line the point is in*/
            for (int j = 0; j < NB_LINES_PER_POINT_CUSTOM(n_qubits); j++){
                line_tab[lines_indices[i][j]]++;
            }             
        }
    }

    /*Then for each line for which all points satisfy the form, we add it to the geometry*/
    for (int i = 1; i <= NB_LINES_CUSTOM(n_qubits); i++){
        int n = line_tab[i];
        if(n != 0 && n != 1 && n != 3)print("INCORRECT NUMBER OF OBS PER LINE! :%d\n",n);
        
        bool goal = (form != &innerProduct_custom);
        /*if all the 3 points of the line are in the geometry*/
        if((n == NB_POINTS_PER_LINE) != complement){
            /*perpsets needs the origin to belong to every line, unlike quadrics*/
            if(form == &innerProduct_custom){
                for (int j = 0; j < NB_POINTS_PER_LINE; j++)
                {
                    if(obs == lines[i][j])goal = true;
                }
                if(goal == complement)continue;
            }
            lines_res[qa.cpt_geometries] = i;
            qa.cpt_geometries++;
        }
    }
    /*we check that the number of observables and lines are the ones expected for quadrics and perpsets*/
    free(line_tab);

    return qa;
}

quantum_assignment perpset(bv obs,size_t** lines_indices,int n_qubits,bv** lines,bool complement) {
    size_t *perp_res = calloc(complement?(NB_LINES_CUSTOM(n_qubits)-NB_LINES_PER_PERPSET(n_qubits)):NB_LINES_PER_PERPSET(n_qubits),sizeof(size_t));
    return zero_locus(obs, &innerProduct_custom,lines_indices,n_qubits,lines,perp_res,complement);
}
quantum_assignment quadric(bv obs,size_t** lines_indices,int n_qubits,bv** lines,bool complement) {
    size_t *quad_res = calloc(complement?(NB_LINES_CUSTOM(n_qubits)-NB_LINES_PER_ELLIPTIC(n_qubits)):NB_LINES_PER_QUADRIC(n_qubits),sizeof(size_t));
    return zero_locus(obs, &quadraticForm_custom,lines_indices,n_qubits,lines,quad_res,complement);
}

size_t current_two_spread = 0;
/**
 * @brief Generates a 2-spread by removing the given spread from the doily
 * 
 * @param w2_doily_lines 
 * @param lines_indices 
 */
void generate_two_spread(size_t* lines_indices){
    size_t current_index = 0;
    for (size_t i = 0; i < NB_LINES_DOILY; i++)
    {
        bool valid = true;
        for (size_t j = 0; j < 5; j++)if(i == lines_indices[j]){valid = false;break;}
        if(!valid)continue;
        index_two_spreads[current_two_spread][current_index++] = i;
    }
    current_two_spread++;
}

/**
 * @brief Builds recursively a spread by checking every possible combination of 5 lines such that they 
 * don't contain a common point (and since 5*3 = 15 the 5 lines contain all the points of the doily)
 * 
 * @param w2_doily_lines 
 * @param incidence_list 
 * @param lines_indices 
 * @param depth 
 */
void doily_spreads_rec(bv w2_doily_lines[][NB_POINTS_PER_LINE],uint16_t incidence_list,size_t* lines_indices,int depth){
    if(depth <= 0){
        generate_two_spread(lines_indices);
        return;
    }
    size_t start = depth == 5?0:lines_indices[depth];
    for (size_t i = start; i < NB_LINES_DOILY; i++)
    {
        bool valid = true;
        uint16_t inc_copy = incidence_list;
        for (size_t j = 0; j < NB_POINTS_PER_LINE; j++)
        {
            if (BGET(inc_copy,w2_doily_lines[i][j])){valid = false;break;}
            BSET(inc_copy,w2_doily_lines[i][j]);
        }
        if(!valid)continue;
        lines_indices[depth-1] = i;
        doily_spreads_rec(w2_doily_lines,inc_copy,lines_indices,depth-1);
    }
    
}

/**
 * @brief generates the 6 spreads of the W2 doily in order to build the corresponding 2-spreads
 * 
 */
void doily_spreads(bv W2_doily_lines[][NB_POINTS_PER_LINE]){
    uint32_t incidence_list = 0;
    size_t lines_indices[5] = {0};
    doily_spreads_rec(W2_doily_lines,incidence_list,lines_indices,5);
}

/**
 * @brief generates a 2-spread from a doily and then computes its contextuality degree
 * 
 */
void two_spread_contextuality_degree(bv doily[],bv W2_doily_lines[][NB_POINTS_PER_LINE],int n_qubits){

    bv** cdegree_two_spread_lines = (bv**)init_matrix(NB_LINES_DOILY,NB_POINTS_PER_LINE,sizeof(bv));

    quantum_assignment qa;
    qa.cpt_geometries = NB_LINES_TWO_SPREAD;
    qa.points_per_geometry = NB_POINTS_PER_LINE;
    qa.n_qubits = n_qubits;

    for (size_t i = 0; i < NB_LINES_DOILY; i++){
        for (size_t j = 0; j < NB_POINTS_PER_LINE; j++)cdegree_two_spread_lines[i][j] = doily[W2_doily_lines[i][j]];
    }
    for (size_t i = 0; i < NB_TWO_SPREADS_PER_DOILY; i++){
        qa.geometry_indices = index_two_spreads[i];
        qa.geometries = cdegree_two_spread_lines;
        int deg = geometry_contextuality_degree(qa,false,false,false,NULL);
        if(deg != 1)print("UNEXPECTED DEGREE:%d, ",deg);
        print("-");
    }

    free_matrix(cdegree_two_spread_lines);
}

/**
 * @brief Contextualit degree of the two spread Metod described in a mail at 27/01/2023
 * 
 * WARNING : 
 */
void contextuality_degree_of_particular_two_spread_27_01_2023(bv W2_doily_lines[][NB_POINTS_PER_LINE]){
    static const int n_qubits = 5;
    bv particular_2_spread[] = {
        [XI]=str_to_bv_custom("ZXYYI",n_qubits),[IX]=str_to_bv_custom("YZXZX",n_qubits),[IZ]=str_to_bv_custom("XYZXY",n_qubits),
        [XX]=str_to_bv_custom("XYZXX",n_qubits),[YI]=str_to_bv_custom("ZXXIX",n_qubits),[ZI]=str_to_bv_custom("XYXYY",n_qubits),
        [XZ]=str_to_bv_custom("YZXZY",n_qubits),[YY]=str_to_bv_custom("XYXYX",n_qubits),[ZX]=str_to_bv_custom("ZXIXZ",n_qubits),
        [YZ]=str_to_bv_custom("YZYXZ",n_qubits),[ZY]=str_to_bv_custom("ZXXIY",n_qubits),[IY]=str_to_bv_custom("YZIYI",n_qubits),
        [ZZ]=str_to_bv_custom("IIYZI",n_qubits),[XY]=str_to_bv_custom("XYYII",n_qubits),[YX]=str_to_bv_custom("XYIZI",n_qubits)
    };
    two_spread_contextuality_degree(particular_2_spread,W2_doily_lines,5);
}

#endif //QUADRICS