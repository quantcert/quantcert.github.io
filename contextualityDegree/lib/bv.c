/**********************************************************************************/
/* Copyright (C) 2024 Axel Muller and Alain Giorgetti                             */
/* Université de Franche-Comté, CNRS, institut FEMTO-ST, F-25000 Besançon, France */
/**********************************************************************************/
/* This software is distributed under the terms of the GNU General Public License */
/* version 2                                                                      */
/**********************************************************************************/
/**
 * @file bv.c
 * @brief Contains functions about the bitvector used to represent observables
 *
 * Some functions are based on [DHGMS22] (https://doi.org/10.1088/1751-8121/aca36f)
 * and [HDS22] (https://doi.org/10.1038/s41598-022-13079-3)
 *
 */

#ifndef MY_BV
#define MY_BV 1/*safe guard*/

#include <time.h>

#include "constants.c"

/*bit vector representation of an observable*/
typedef uint32_t bv;
/*even or odd part of an observable*/
typedef uint32_t word;

/* Binary representation of each 1-qubit observable */
#define I 0b00
#define X 0b01
#define Y 0b11
#define Z 0b10

/* In the bitvector 0bRSTU, RS is the Z part and TU is the X part. */
#define II 0b0000//0
#define IX 0b0001//1
#define XI 0b0010//2
#define XX 0b0011//3
#define IZ 0b0100//4
#define IY 0b0101//5
#define XZ 0b0110//6
#define XY 0b0111//7
#define ZI 0b1000//8
#define ZX 0b1001//9
#define YI 0b1010//10
#define YX 0b1011//11
#define ZZ 0b1100//12
#define ZY 0b1101//13
#define YZ 0b1110//14
#define YY 0b1111//15

#define BV_LIMIT_CUSTOM(N) ((bv)pow4(N))//represents the limit value (not included) for an observable bit vector
#define BV_LIMIT BV_LIMIT_CUSTOM(N_QUBITS)

//#define mult(a,b) ((a) & (b))//multiplies an observable by another
//#define add(a,b)  ((a) ^ (b))//applies an observable to the other

#define Qplus ^ //applies an observable to the other
#define Qtimes &

char *HEXAD[][DOILY_SIZE] = {
  {"IYY","YYY","XZX","YII","XYZ","XIX","XXZ","YZI","ZYZ","ZXZ","IXY","YXY","IZI","ZZX","ZIX"},
  {"IYI","ZYX","YYI","ZIX","ZXZ","ZZZ","YII","YXY","IXY","XIX","XZZ","YZY","XXZ","XYX","IZY"},
  {"YII","YZI","YYI","IZI","XXZ","ZXZ","IYI","XIZ","XYZ","IXI","ZZZ","ZIZ","XZZ","YXI","ZYZ"},
  {"IYI","ZIZ","ZYZ","ZYX","ZXZ","ZZZ","IXI","ZXX","ZIX","IXY","IYY","IZY","ZZX","IYY","IZI"},
  {"XZZ","XXX","IYY","XIZ","IXI","XYZ","XIX","XZX","IZI","IYI","IXY","XXZ","IIY","XYX","IZY"},
  {"YII","YYY","YZY","IYY","XIX","ZIX","IZY","XXX","XYZ","IXI","ZZZ","ZXX","XZZ","YXI","ZYZ"}
};

/**
 * @brief returns the number of bits set to one on a bit vector
 * 
 * @param bv1 
 * @return unsigned int 
 */
unsigned int popcnt(bv bv1){
  return __builtin_popcount(bv1);
}

/**
 * @brief returns the factorial of a given number
 * 
 * @param a 
 * @return int 
 */
int factorial(int a){
  int res = 1;
  for (int i = 2; i <= a; i++)res*=i;
  return res;
}

/**
 * @brief Get the even part of the observable
 * 
 * @param bv1 
 * @param n_qubits 
 * @return word 
 */
word get_Z(bv bv1,int n_qubits){
  return bv1 >> n_qubits; /* another mask & mask(N_QUBITS) could be applied in case the bv is beyond limits*/
}
/**
 * @brief Get the odd part of the observable
 * 
 * @param bv1 
 * @param n_qubits 
 * @return word 
 */
word get_X(bv bv1,int n_qubits){
  return bv1 & mask(n_qubits);
}

/**
 * @brief merge the x and z parts into an bv
 * 
 * @param bv1 
 * @param n_qubits number of qubits in the bit vector
 * @return bv
 */
bv to_index_custom(word z,word x,int n_qubits){
  return (z << n_qubits) | x;
}

/**
 * @brief returns 1 iff the number of bits in the word
 * is odd
 * (compilation optimization removes automatically unnecessary ifs)
 * bit parity trick found in http://graphics.stanford.edu/~seander/bithacks.html#ParityLookupTable
 * @param n 
 * @param n_qubits 
 * @return int 
 */
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

/**
 * @brief symplectic product of 2 observables
 * 
 * bitvector parity in O(log2(N))
*/
unsigned int innerProduct_custom(bv i1,bv i2,int n_qubits){
  /*bit parity check with the addition*/
  return bit_parity(innerProductVector(i1,i2,n_qubits),n_qubits);
}
unsigned int innerProduct(bv i1,bv i2){return innerProduct_custom(i1,i2,N_QUBITS);}

/**
 * @brief Base quadratic form of an observable
 *
 * Q_0(bv1) = x1x2 + x3x4 + ... + x_{N-1}x_N
 *
 * See [DHGMS22] Section 5.1
 *
 */
unsigned int baseQuadraticFormVector(bv i, int n_qubits) {
  /*We perform the multiplications*/
  return get_Z(i, n_qubits) Qtimes get_X(i, n_qubits);
}

/**
 * @brief Quadratic of bv1 given a base
 * 
 * Q_base(bv1) = Q_0(bv1) + <bv1,base>
 * 
 * See [DHGMS22] Section 5.1
 * 
 */
unsigned int quadraticForm_custom(bv base, bv bv1, int n_qubits) {
  /* The Xor from Qplus preserves the parity of the two vectors */
  word n = baseQuadraticFormVector(bv1, n_qubits)
      Qplus innerProductVector(base, bv1, n_qubits);

  return bit_parity(n, n_qubits);
}

unsigned int quadraticForm(bv base, bv bv1) {
  return quadraticForm_custom(base, bv1, N_QUBITS);
}

/**
 * @brief same as inner product but with 2 qubits (faster)
 * 
 * @param i1 
 * @param i2 
 * @return int 
 */
unsigned int innerProduct_W2(bv i1,bv i2){

  word n = (get_Z(i1,INDEX_DOILY_SIZE) Qtimes get_X(i2,INDEX_DOILY_SIZE)) Qplus
            (get_Z(i2,INDEX_DOILY_SIZE) Qtimes get_X(i1,INDEX_DOILY_SIZE));

  return bit_parity(n,INDEX_DOILY_SIZE);
}

/**
 * @brief Performs the transvection of the point q from the base p
 * 
 * T_p(q) = q + <p,q>p
 * 
 * See [HDS22] Section 1
 */
bv transvection(bv p,bv q,int n_qubits){
  return q Qplus (innerProduct_custom(p,q,n_qubits) * p);
}

/**
 * @brief Returns true if the observable is symmetric i.e. if its number of Y's is even
 * 
 * @param bv1 
 * @return int 
 */
bool is_symmetric(bv bv1,int n_qubits){
  return !bit_parity(get_X(bv1,n_qubits) Qtimes get_Z(bv1,n_qubits),n_qubits);
}

/**
 * @brief Prints every bit of a bitvector
 * 
 * @param b 
 */
void print_bits_custom(bv b,int n_bits){
  for (int i = 0; i < n_bits; i++)print("%d", (b >> (n_bits - i - 1)) & 1);
}

/**
 * @brief gives the number of Is in an observable of size N_QUBITS
 * 
 * @param bv1 
 * @return int 
 */
int n_I_custom(bv bv1,int n_qubits){
  /*for each bit gives 1 iff even and odd aren't I*/
  //bv tmp = bv1.even | bv1.odd;
  bv1 = (bv1 | (bv1 >> n_qubits)) & mask(n_qubits);
  int res = n_qubits - popcnt(bv1);
  if (DOUBLE_CHECK && (res < 0 || res >= n_qubits))
    print("!5!%d!!", res); // shouldn't occur
  return res;
}
int n_I(bv bv1){return n_I_custom(bv1,N_QUBITS);}

/**
 * @brief returns a 2bit long integer corresponding to the selected gate
 *                     
 * example : get_gate(IXZZ,1) = X
 *                     ^
*/
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

/**
 * @brief Sets the gate "obs" to the quantum gate "index" of the observable "bv1" 
 * 
 * @param bv1 
 * @param obs 
 * @param gate 
 * @return bv 
 */
bv set_gate(bv bv1,bv obs,int index){return set_gate_custom(bv1,obs,index,N_QUBITS);}

/**
 * @brief Prints a bitvector in the standard output (used to write files)
 * 
 * @param bv1 
 * @param N_QUBITS 
 */
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

/**
 * @brief prints a bit vector with its pauli gates in the console
 * 
 * @param bv1 bit vector to print
 */
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

/**
 * @brief gives for each character its corresponding qubit
 * 
 * @param gate 
 * @return unsigned int 
 */
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

/**
 * @brief Transforms a string into a qubit bit vector
 * 
 * @param str observable in string form
 * @param n_qubits number of qubits in the observable
 * @return bv 
 */
bv str_to_bv_custom(char *str,int n_qubits){
  bv res = I;
  int i = 0;
  for (char c = str[i]; c >= 'A' && c <='Z'; c = str[++i])res = set_gate_custom(res,char_to_gate(c),i,n_qubits);
  return res;
}

/**
 * @brief tTransforms a string into a N_QUBITS qubit bit vector
 * 
 * @param str observable in string form
 * @return bv 
 */
bv str_to_bv(char *str){return str_to_bv_custom(str,N_QUBITS);}

bv extend_bv(bv bv1,size_t from_qubits,size_t to_qubits){
    word z = get_Z(bv1,from_qubits);
    word x = get_X(bv1,from_qubits);

    return to_index_custom(z,x,to_qubits);
}

#endif //MY_BV
