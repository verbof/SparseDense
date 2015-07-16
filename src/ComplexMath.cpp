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

void solveSystemComplex(CSRcomplex& A, complex< double >* X, complex< double >* B, int pardiso_mtype, int number_of_rhs)
{
    /*cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout << "@@@ S O L V I N G     A    L I N E A R    S Y S T E M  @@@" << endl;
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

    cout << "*** G e n e r a t i n g    # " << number_of_rhs << "   r h s *** " << endl;*/

    // initialize pardiso and forward to it minimum number of necessary parameters
    /*int pardiso_message_level = 0;

    ParDiSO pardiso(pardiso_mtype, pardiso_message_level);*/


    pardiso_mtype = pardiso_mtype;

    // Numbers of processors, value of OMP_NUM_THREADS
    int number_of_processors = 1;
    char* var = getenv("OMP_NUM_THREADS");

    if(var != NULL)
        sscanf( var, "%d", &number_of_processors );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }

    pardiso_var.iparm[3]  = number_of_processors;
    pardiso_var.iparm[8]  = 0;


    timing secs;
    double initializationTime = 0.0;
    double factorizationTime  = 0.0;
    double solutionTime       = 0.0;



    //cout << "S Y M B O L I C     V O O D O O" << endl;

    secs.tick(initializationTime);
    pardiso_var.init(A, number_of_rhs);
    secs.tack(initializationTime);



    //cout << "L U                 F A C T O R I Z A T I O N" << endl;

    secs.tick(factorizationTime);
    pardiso_var.factorize(A);
    secs.tack(factorizationTime);


    //cout << "L U                 B A C K - S U B S T I T U T I O N" << endl;

    secs.tick(solutionTime);
    pardiso_var.solve(A, X, B);
    secs.tack(solutionTime);


    //errorReport(number_of_rhs, A, B, X);
    // writeSolution(number_of_rhs, A.nrows, X);

    if(iam==0){
    cout << "-------------------------------" << endl;
    cout << "T I M I N G         R E P O R T" << endl;
    cout << "-------------------------------" << endl;
    cout.setf(ios::floatfield, ios::scientific);
    cout.precision(4);
    cout << "Initialization phase: " << initializationTime*0.001 << " sec" << endl;
    cout << "Factorization  phase: " << factorizationTime*0.001 << " sec" << endl;
    cout << "Solution       phase: " << solutionTime*0.001 << " sec" << endl;}
}

void shiftIndicesComplex(CSRcomplex& A, int value)
{
  int i;
  for (i = 0; i <= A.nrows; i++)
  {
    A.pRows[i] += value;
  }

  for (i = 0; i < A.nonzeros; i++)
  {
    A.pCols[i] += value;
  }    
}

void create2x2SymBlockMatrixComplex(CSRcomplex& A, CSRcomplex& B, CSRcomplex& T, // input
                             CSRcomplex& C)  // output
{
    assert(A.nrows==B.nrows);
    assert(A.ncols==B.nrows);
    assert(T.ncols==B.ncols);
    assert(T.nrows==B.ncols);

    int nrows    = A.nrows + T.nrows;
    int ncols    = A.ncols + T.ncols;
    int nonzeros = A.nonzeros + B.nonzeros + T.nonzeros;

    int* ic   = new int[nrows + 1];
    int* jc   = new int[nonzeros];
    complex< double >* c = new complex< double >[nonzeros];

    int nonzero_counter = 0;
    ic[0] = nonzero_counter;
    for (int i = 0; i < A.nrows; i++)
    {
        // push ith row of A
        for (int index = A.pRows[i]; index < A.pRows[i+1]; index++)
        {
            int& j              = A.pCols[index];
            if (j>=i)
            {
                complex< double >& a_ij        = A.pData[index];

                c[nonzero_counter] = a_ij;
                jc[nonzero_counter] = j;

                nonzero_counter++;
            }
        }

        // push ith row of B
        for (int index = B.pRows[i]; index < B.pRows[i+1]; index++)
        {
            int& j              = B.pCols[index];
            complex< double >& b_ij        = B.pData[index];

            c[nonzero_counter] = b_ij;
            jc[nonzero_counter] = A.ncols + j;

            nonzero_counter++;
        }

        ic[i+1] = nonzero_counter;
    }

    for (int i = 0; i < T.nrows; i++)
    {
        // push ith row of T
        for (int index = T.pRows[i]; index < T.pRows[i+1]; index++)
        {
            int& j              = T.pCols[index];
            complex< double >& t_ij        = T.pData[index];

            c[nonzero_counter] = t_ij;
            jc[nonzero_counter] = A.ncols + j;

            nonzero_counter++;
        }

        ic[A.nrows+i+1] = nonzero_counter;
    }
    if (nonzero_counter != nonzeros)
        cout << "Nonzeroes do not match, nonzero_counter= " << nonzero_counter << "; nonzeros= " << nonzeros <<endl;


    C.make(nrows, ncols, nonzero_counter, ic, jc, c);
    C.sortColumns();
    // C.writeToFile("C.csr");
}

