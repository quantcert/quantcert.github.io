/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file quadrics.c
 * @brief C functions to generate various geometries, among which quadrics.
 */
#include "quadrics.h"

#include "bv.h"
#include "contextuality_degree.h"




int NB_GENERATORS(int N){return (N<=0)?1:(pow2(N)+1)*NB_GENERATORS(N-1);}

size_t NB_SUBSPACES(int N,int K){
    size_t a = 1,b = 1;
    for (int i = 1; i <= K+1; i++){
        a *= (pow2(N-K-1+i)-1) * (pow2(N-i+1)+1);
        b *= pow2(i)-1;
    }
    return a/b;
}

/**all 6 2-spreads as lists of indexes of the W2 doily lines*/
size_t index_two_spreads[NB_TWO_SPREADS_PER_DOILY][NB_LINES_TWO_SPREAD];


void insert_line_index(bv obs, size_t index,size_t** lines_indices,int n_qubits) {
    int i;

    for (i = 0; lines_indices[obs][i] != NO_LINE && i<NB_LINES_PER_POINT_CUSTOM(n_qubits); i++)
    

    
    if (obs >= BV_LIMIT_CUSTOM(n_qubits) || i >= NB_LINES_PER_POINT_CUSTOM(n_qubits))
        print("error out of bounds %d %d / %d %ld [%d,%d]!\n", obs, i, BV_LIMIT_CUSTOM(n_qubits),
            NB_LINES_PER_POINT_CUSTOM(n_qubits),obs,i);

    lines_indices[obs][i] = index;

}

quantum_assignment generate_total_lines(size_t*** lines_indices,int n_qubits) {

    print("generating lines...");

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
        print("INDEX INCOHERENT ! : %ld != %ld\n", index,NB_LINES_CUSTOM(n_qubits) + 1);

    quantum_assignment qa = {
        .geometry_indices = calloc(NB_LINES_CUSTOM(n_qubits), sizeof(size_t)),
        .geometries = lines,
        .cpt_geometries = NB_LINES_CUSTOM(n_qubits),
        .points_per_geometry = NB_POINTS_PER_LINE,
        .n_qubits = n_qubits
    };
    for (size_t i = 0; i < (size_t)NB_LINES_CUSTOM(n_qubits); i++) qa.geometry_indices[i] = i+1;
    quantum_assignment_compute_negativity(&qa);

    return qa;
}

int bv_left_most(bv obs){
    if(obs == 0)return -1;
    for (size_t i = 0; i < 8*sizeof(bv); i++){
        obs >>= 1;
        if(obs == 0)return i;
    }
    return 8*sizeof(bv);//should not happen
}

int bv_anti_cmp(const void * first, const void * second ) {
    return * (const bv *) second - * (const bv *) first;
}

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

bool commuting_rec(bv bv1[],uint32_t index,uint32_t size,int n_qubits,quantum_assignment* res){
    if(is_done)return true;
    if(index == size-1){
        /*If every observable is generated we call the callback function with the right parameters*/
        bv last = I;
        for (size_t i = 0; i < size-1; i++)
        {
            res->geometries[res->cpt_geometries][i] = bv1[i];
            last = last Qplus bv1[i];
        }
        //check if the last observable is not in the list
        for (size_t i = 0; i < size-1; i++)if(bv1[i] >= last)return false;
        res->geometries[res->cpt_geometries][size-1] = last;
        res->cpt_geometries++;
        return false;
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
        if(commuting_rec(bv1,index+1,size,n_qubits,res))return true;
    }
    return false;
}

quantum_assignment commuting(uint32_t size,int n_qubits){

    quantum_assignment res = {
        .geometries = (bv**)init_matrix(NB_LINES_CUSTOM(n_qubits)*1000,size,sizeof(bv)),
        .cpt_geometries = 0,
        .points_per_geometry = size,
        .n_qubits = n_qubits
    };
    
    /*Splitting the first for loop between every thread*/
    int num_done = 0;
    bool stop_iters = false;
    //#pragma omp parallel for schedule(dynamic,1) if(size > 2)
    for (size_t i = I+1; i < BV_LIMIT_CUSTOM(n_qubits); i++){
        if(stop_iters)continue;
        bv bv1[size];
        bv1[0] = i;
        if(commuting_rec(bv1,1,size,n_qubits,&res))stop_iters = true;
        #pragma omp critical
        {
        num_done++;
        if (N_QUBITS > 3)
            print(" %.2f%% ", 100*((float)num_done/(float)BV_LIMIT));
        }
    }
    quantum_assignment_autofill_indices(&res);
    quantum_assignment_compute_negativity(&res);
    return res;
}

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
        if(_subspaces_current_index % 200000 == 1 && PRINT_PROGRESSION)print("%ld,",_subspaces_current_index);
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

quantum_assignment subspaces(int n_qubits, int k)
{

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

void doily_spreads(bv W2_doily_lines[][NB_POINTS_PER_LINE]){
    uint32_t incidence_list = 0;
    size_t lines_indices[5] = {0};
    doily_spreads_rec(W2_doily_lines,incidence_list,lines_indices,5);
}

