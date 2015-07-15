#ifndef Real_Math_hpp
#define Real_Math_hpp

#include "config.hpp"

class CSRcomplex;

void genComplexOnes(int m, int n, double real_value, double imag_value, complex<double> * ones);
void genDiagonalDComplex(int n, double d_real, double d_imag, complex< double >* diag);

void make3DLaplace_complex(int nx, int ny, int nz, CSRcomplex& L, double imag);

#endif