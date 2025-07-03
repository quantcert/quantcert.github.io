/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file bv.h
 * @brief Contains functions about the bitvector used to represent observables
 *
 * Some functions are based on [DHGMS22] (https://doi.org/10.1088/1751-8121/aca36f)
 * and [HDS22] (https://doi.org/10.1038/s41598-022-13079-3)
 *
 */
#include "bv.h"

char *HEXAD[][DOILY_SIZE] = {
    {"IYY", "YYY", "XZX", "YII", "XYZ", "XIX", "XXZ", "YZI", "ZYZ", "ZXZ", "IXY", "YXY", "IZI", "ZZX", "ZIX"},
    {"IYI", "ZYX", "YYI", "ZIX", "ZXZ", "ZZZ", "YII", "YXY", "IXY", "XIX", "XZZ", "YZY", "XXZ", "XYX", "IZY"},
    {"YII", "YZI", "YYI", "IZI", "XXZ", "ZXZ", "IYI", "XIZ", "XYZ", "IXI", "ZZZ", "ZIZ", "XZZ", "YXI", "ZYZ"},
    {"IYI", "ZIZ", "ZYZ", "ZYX", "ZXZ", "ZZZ", "IXI", "ZXX", "ZIX", "IXY", "IYY", "IZY", "ZZX", "IYY", "IZI"},
    {"XZZ", "XXX", "IYY", "XIZ", "IXI", "XYZ", "XIX", "XZX", "IZI", "IYI", "IXY", "XXZ", "IIY", "XYX", "IZY"},
    {"YII", "YYY", "YZY", "IYY", "XIX", "ZIX", "IZY", "XXX", "XYZ", "IXI", "ZZZ", "ZXX", "XZZ", "YXI", "ZYZ"}};


unsigned int popcnt(bv bv1){
  return __builtin_popcount(bv1);
}


int factorial(int a){
  int res = 1;
  for (int i = 2; i <= a; i++)res*=i;
  return res;
}

word get_Z(bv bv1,int n_qubits){
  return bv1 >> n_qubits; /* another mask & mask(N_QUBITS) could be applied in case the bv is beyond limits*/
}

word get_X(bv bv1,int n_qubits){
  return bv1 & mask(n_qubits);
}

bv to_index_custom(word z,word x,int n_qubits){
  return (z << n_qubits) | x;
}

int bit_parity(word n,const int n_qubits){
  if(n_qubits > 16)n ^= n >> 16;/* <- up to W32*/
  if(n_qubits > 8)n ^= n >> 8;/* <- up to W16*/
  if(n_qubits > 4)n ^= n >> 4;/* <- up to W8*/
  if(n_qubits > 2)n ^= n >> 2;/* <- up to W4*/
  n ^= n >> 1;/* <- enough for W2*/
  return (n & 1);/* == 1*/;
}

word innerProductVector(bv i1,bv i2,int n_qubits){
    /*We perform the multiplications, plus one step of the bit parity check with the addition*/
    return (get_Z(i1,n_qubits) Qtimes get_X(i2,n_qubits)) Qplus
           (get_Z(i2,n_qubits) Qtimes get_X(i1,n_qubits));
}

unsigned int innerProduct_custom(bv i1,bv i2,int n_qubits){
  /*bit parity check with the addition*/
  return bit_parity(innerProductVector(i1,i2,n_qubits),n_qubits);
}
unsigned int innerProduct(bv i1,bv i2){return innerProduct_custom(i1,i2,N_QUBITS);}

unsigned int baseQuadraticFormVector(bv i, int n_qubits) {
  /*We perform the multiplications*/
  return get_Z(i, n_qubits) Qtimes get_X(i, n_qubits);
}

unsigned int quadraticForm_custom(bv base, bv bv1, int n_qubits) {
  /* The Xor from Qplus preserves the parity of the two vectors */
  word n = baseQuadraticFormVector(bv1, n_qubits)
      Qplus innerProductVector(base, bv1, n_qubits);

  return bit_parity(n, n_qubits);
}

unsigned int quadraticForm(bv base, bv bv1) {
  return quadraticForm_custom(base, bv1, N_QUBITS);
}

unsigned int innerProduct_W2(bv i1,bv i2){

  word n = (get_Z(i1,INDEX_DOILY_SIZE) Qtimes get_X(i2,INDEX_DOILY_SIZE)) Qplus
            (get_Z(i2,INDEX_DOILY_SIZE) Qtimes get_X(i1,INDEX_DOILY_SIZE));

  return bit_parity(n,INDEX_DOILY_SIZE);
}

bv transvection(bv p,bv q,int n_qubits){
  return q Qplus (innerProduct_custom(p,q,n_qubits) * p);
}

bool is_symmetric(bv bv1,int n_qubits){
  return !bit_parity(get_X(bv1,n_qubits) Qtimes get_Z(bv1,n_qubits),n_qubits);
}

void print_bits_custom(bv b,int n_bits){
  for (int i = 0; i < n_bits; i++)print("%d", (b >> (n_bits - i - 1)) & 1);
}

int n_I_custom(bv bv1,int n_qubits){
  /*for each bit gives 1 iff even and odd aren't I*/
  //bv tmp = bv1.even | bv1.odd;
  bv1 = (bv1 | (bv1 >> n_qubits)) & mask(n_qubits);
  int res = n_qubits - popcnt(bv1);
  if (DOUBLE_CHECK && (res < 0 || res >= n_qubits))
    print("NUMBER OF QUBITS ERROR");
  return res;
}
int n_I(bv bv1){return n_I_custom(bv1,N_QUBITS);}

unsigned int get_gate(bv i,unsigned int gate,int n_qubits){
  unsigned int res = 0;
  unsigned int shift = (n_qubits - gate - 1);/*how much we have to shift each word*/
  res |= (get_X(i,n_qubits) >> shift) & 0b01;/*we put the X part at the end*/
  res |= ((get_Z(i,n_qubits) >> shift) & 0b01) << 1;/*we put the Z part in front*/
  return res;
}

bv set_gate_custom(bv bv1,bv obs,int index,int n_qubits){
  int shift = n_qubits-index-1;//the bit shift is the opposite of the index
  bv1 &= ~(1<<shift) & ~(1<<(shift+n_qubits));//clearing even and odd values for the old value
  bv1 |= (obs>>1) << (shift+n_qubits);//set even on first bit, then placing it correctly
  bv1 |= (obs& 1) << shift;//same for odd
  return bv1;
}

bv set_gate(bv bv1,bv obs,int index){return set_gate_custom(bv1,obs,index,N_QUBITS);}

void print_BV_to_file(bv bv1,int n_qubits,FILE* output){
  for(size_t i = 0;i < (size_t)n_qubits;i++){
    switch (get_gate(bv1,i,n_qubits)){
    case I: fprintf(output,"I");break;
    case X: fprintf(output,"X");break;
    case Y: fprintf(output,"Y");break;
    case Z: fprintf(output,"Z");break;
    }
  }
}

void print_BV_custom(bv bv1,int n_qubits){

  print("[");
  for(size_t i = 0;i < (size_t)n_qubits;i++){
    switch (get_gate(bv1,i,n_qubits)){
    case I: print("I");break;
    case X: print("X");break;
    case Y: print("Y");break;
    case Z: print("Z");break;
    }
  }
  print("]");
}
void print_BV_W2(bv bv1){print_BV_custom(bv1,INDEX_DOILY_SIZE);}

unsigned int char_to_gate(char gate){
  switch (gate)
  {
    case 'I':return I;
    case 'X':return X;
    case 'Y':return Y;
    case 'Z':return Z;
    default:print("read error !");return I;
  }
}

bv str_to_bv_custom(char *str,int n_qubits){
  bv res = I;
  int i = 0;
  for (char c = str[i]; c >= 'A' && c <='Z'; c = str[++i])res = set_gate_custom(res,char_to_gate(c),i,n_qubits);
  return res;
}

bv str_to_bv(char *str){return str_to_bv_custom(str,N_QUBITS);}

bv extend_bv(bv bv1,size_t from_qubits,size_t to_qubits){
    word z = get_Z(bv1,from_qubits);
    word x = get_X(bv1,from_qubits);

    return to_index_custom(z,x,to_qubits);
}
