#ifndef Complex_Math_hpp
#define Complex_Math_hpp

#include "config.hpp"

class CSRcomplex;

void genComplexOnes(int m, int n, double real_value, double imag_value, complex<double> * ones);
void genDiagonalDComplex(int n, double d_real, double d_imag, complex< double >* diag);

void make3DLaplace_complex(int nx, int ny, int nz, CSRcomplex& L, double imag);

void solveSystemComplex(CSRcomplex& A, complex< double >* X, complex< double >* B, int pardiso_mtype, int number_of_rhs);

void shiftIndicesComplex(CSRcomplex& A, int value);

void create2x2SymBlockMatrixComplex(CSRcomplex& A, CSRcomplex& B, CSRcomplex& T, CSRcomplex& C);

#endif
