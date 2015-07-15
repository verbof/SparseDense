#include "config.hpp"
#include "ComplexMath.hpp"
#include "timing.hpp"
#include "CSRcomplex.hpp"
#include "ParDiSO.hpp"
#include "shared_var.h"
#include <cassert>
#include <cstring>


#define mesh(i,j,k) mesh[nx*ny*((k)-1) + nx*((j)-1)+(i)-1]
#define nzvals(i) nzvals[(i)-1]
#define rowind(i) rowind[(i)-1]
#define colptr(i) colptr[(i)-1]
#define xsol(i)   xsol[(i)-1]
#define rhs(i)    rhs[(i)-1]
#define diag(i)   diag[(i)-1]
#define diag2(i)  diag2[(i)-1]


void genComplexOnes(int m, int n, double real_value, double imag_value, complex< double >* ones)
{

    for(int j = 0; j < m*n; j++)
    {
	ones[j] = complex<double>(real_value, imag_value);
    }

}

void genDiagonalDComplex(int n, double d_real, double d_imag, complex< double >* diag)
{
    for (int k = 0; k < n; k++)
    {
        diag[k*n + k] = complex<double>(d_real, d_imag);
    }

}


void make3DLaplace_complex(int nx, int ny, int nz, CSRcomplex& L, double imag)
{

   int i, j, k, nnodes,nedges, nnz;
   int count, node;
   int* mesh;
   int* rowind;
   int* colptr;
   complex< double >* nzvals;
   int nnz_full;
   double dval;


   nnodes = nx*ny*nz;
   cout << "Generating a " << nnodes << "-by-" << nnodes << " Laplacian matrix..." << endl;

   mesh = (int*) calloc(nnodes, sizeof(int));

   for (i = 0; i<nnodes; i++) mesh[i]=i+1;
   
   /* first pass to count the number of edges */
   /* Dirichlet BC */ 
   nedges = 0; 
  
   for (k = 1; k <= nz; k++) {
     for (j = 1; j <= ny; j++) {
       for (i = 1; i <= nx; i++) {
       if (k < nz) nedges++;
       if (k > 1)  nedges++;
       if (j < ny) nedges++;
       if (j > 1)  nedges++;
       if (i < nx) nedges++;
       if (i > 1)  nedges++;   
     }
   }
   }

   /* print the matrix dimension and number of nonzeros */
   nnz = nedges/2 + nnodes;
   //printf("%d  %d\n", nnodes, nnz);
   nnz_full = nnz;
   
   colptr  = (int*)malloc((nnodes+1)*sizeof(int));
   rowind  = (int*)malloc(nnz_full*sizeof(int));
   nzvals  = (complex< double >*)malloc(nnz_full*sizeof(complex< double >));

   colptr(1) = 1;
   count = 0;
   node  = 0;

  
   dval = 20.0;
   //printf(" Laplacian with Dirichlet boundary condition\n");

   for (k = 1; k <= nz; k++) {
     for (j = 1; j <= ny; j++) {
        for (i = 1; i <= nx; i++) {
	 /* diagonal */
         dval = nx*ny*nz + (i + j + k) / 100 ;
         //if (printa)
            //printf("%d %d  %8.2e\n", mesh(i,j,k), mesh(i,j,k), dval); 

         rowind[count] = mesh(i,j,k); 
         nzvals[count] = complex<double>(dval + node,imag); /* no pivots */
         count++;

         /* lower */
         if (i < nx) {
            //if (printa) 
	       //printf("%d %d -1.0\n", mesh(i+1,j,k), mesh(i,j,k));

            rowind[count] = mesh(i+1,j,k);
            dval = -1.0;
            dval = (i + j + k) / 100 ;
            nzvals[count] = complex<double>(dval,imag);
            count++;
         }


         /* right */
         if (j < ny) {
            //if (printa) 
	       //printf("%d %d -1.0\n", mesh(i,j+1,k), mesh(i,j,k)); 

            rowind[count] = mesh(i,j+1,k);
            dval = -1.0;
            dval = (i + j + k) / 100 ;
            nzvals[count] = complex<double>(dval,imag);
            count++;
         }       

         /* lower */
         if (k < nz) {
            //if (printa) 
	       //printf("%d %d -1.0\n", mesh(i,j,k+1), mesh(i,j,k));

            rowind[count] = mesh(i,j,k+1);
            dval = -1.0;
            nzvals[count] = complex<double>(dval,imag);
            count++; 
         }

         node++;
         colptr(node+1) = count+1; 
         
      } 
   }
   }
   if (count != nnz) {
       printf(" is this correct? count = %d, nnz = %d\n", count, nnz);  
       return;
   }

    //rowind = ja, nzvals = a  nnodes colptr(nnodes) colptr(j) = ia, 


   L.allocate(nnodes, nnodes, colptr(nnodes));
   L.make(nnodes, nnodes, colptr(nnodes), colptr, rowind, nzvals);

   /*
   if (printa) 
   { 
	printf("%d \n", nnodes );
	printf("%d \n", nnodes );
	printf("%d \n", colptr(nnodes) );
        for (j = 1; j <= nnodes+1; j++) 
		printf("%d \n", colptr(j)-1 );
        for (j = 1; j <= colptr(nnodes); j++) 
		printf("%d \n", rowind(j)-1 );
        for (j = 1; j <= colptr(nnodes); j++)
		printf("%e \n", nzvals(j) );
   }
   */

   free(mesh);
   free(colptr); 
   free(rowind); 
   free(nzvals); 

   cout << "- Done!" << endl;

}
