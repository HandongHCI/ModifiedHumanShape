/**
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Stefanie Wuhrer, Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#ifndef GAUSS_VECTOR
#define GAUSS_VECTOR

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>
#include "matrix.h"
//Mex 

#include <math.h>
#include "mex.h"
#include "blas.h"
//#include "stddef.h"
#include "lapack.h"
//#include  "blascompat32.h"

typedef double doublereal;

//includes for CLAPACK
// namespace clapack
// {
	// extern "C"
 //    {
		// #include "blaswrap.h"
		// #include "f2c.h"
	// 	extern int dgemm_(char *transa, char *transb, ptrdiff_t *m, ptrdiff_t *n, ptrdiff_t *k, doublereal *alpha, doublereal *a, ptrdiff_t *lda, doublereal *b, ptrdiff_t *ldb, doublereal *beta, doublereal *c, ptrdiff_t *ldc);
	// 	extern int dgemv_(char *trans, ptrdiff_t *m, ptrdiff_t *n, doublereal *alpha, doublereal *a, ptrdiff_t *lda, doublereal *x, ptrdiff_t *incx, doublereal *beta, doublereal *y, ptrdiff_t *incy);
	// 	extern int dsyevx_(char *jobz, char *range, char *uplo, ptrdiff_t *n, doublereal *a, ptrdiff_t *lda, doublereal *vl, doublereal *vu, ptrdiff_t *il, ptrdiff_t *iu, doublereal *abstol, ptrdiff_t *m, doublereal *w, doublereal *z__, ptrdiff_t *ldz, doublereal *work, ptrdiff_t *lwork, ptrdiff_t *iwork, ptrdiff_t *ifail, ptrdiff_t *info);
	// 	extern double dlamch_(char *cmach);
	// 	extern int dgetrf_(ptrdiff_t *m, ptrdiff_t *n, doublereal *a, ptrdiff_t *lda, ptrdiff_t *ipiv, ptrdiff_t *info);
	// 	extern int dgetri_(ptrdiff_t *n, doublereal *a, ptrdiff_t *lda, ptrdiff_t *ipiv, doublereal *work, ptrdiff_t *lwork, ptrdiff_t *info);
	// 	extern int dgesdd_(char *jobz, ptrdiff_t *m, ptrdiff_t *n, doublereal *a, ptrdiff_t *lda, doublereal *s, doublereal *u, ptrdiff_t *ldu, doublereal *vt, ptrdiff_t *ldvt, doublereal *work, ptrdiff_t *lwork, ptrdiff_t *iwork, ptrdiff_t *info);
	// 	extern int dgels_(char *trans, ptrdiff_t *m, ptrdiff_t *n, ptrdiff_t *nrhs, doublereal *a, ptrdiff_t *lda, doublereal *b, ptrdiff_t *ldb, doublereal *work, ptrdiff_t *lwork, ptrdiff_t *info);
	// }
// }


class GaussVector
{
public:
	//constructor:
	GaussVector();
	GaussVector(int dimension, double * mean, double * covMat); 
	~GaussVector();

	//generate a random vector of a Gaussian distribution:
    //Monica-Changed qualification error    
    //void GaussVector::generateRandomVector(double *& result);
    void generateRandomVector(double *& result);
	//if returnUnits = true, unitVec returns unit vectors. Otherwise, you must provide unit vectors!!!
	void generateRandomVector(double *& result, bool returnUnits, double *& unitVec);	
	//generate numberSamples vectors of a Gaussian distribution. If returnUnits = true, unitVec returns unit vectors.
	//Otherwise, you must provide unit vectors!!!
	void generateRandomVectors(int numberSamples, double *& results, bool returnUnits, double *& unitVecs);

	//Method to sample on a surface of equal likelihood (an ellipsoid)
	void generateRandomVector(double distFromMean, double *& result);

	//methods to return computed quantities:
	int getDimension(){return dimension;}
	double * getMean(){return mean;}
	double * getCovarianceMat(){return covarianceMat;}
	double * getEigenvalues(){return eigenvalues;}
	double * getEigenvectors(){return eigenvectors;}
	double * getSqrtEigenvalueMatrix(){return sqrtEigenvalueMatrix;}
	double * getMultMat(){return multMatrix;}
	double getDeterminantOfCov();

	//method to compute and print numberSamples random vectors to a file (filename) that is compatible with Maple 7 
	//for plotting:
	void printSampleFile(int numberSamples, char * filename);
	void printSampleFile(int numberSamples, char * filename, double *& samples);
	//print only selected dimensions:
	void printSampleFile(int numberSamples, char * filename, int numDimensions, int * dimensions);
	void printSampleFile(int numberSamples, char * filename, double *& samples, int numDimensions, int * dimensions);

private:
	//compute the eigenvalues and eigenvectors of the covariance matrix using the CLAPACK library
	bool computeEigen(double * covMatrix);
	//compute a one-dimensional unit normal distribution:
	double getUnitNormalDistribution();
	//member variables:
	int dimension;
	double * mean; 
	double * covarianceMat;
	double * eigenvalues;
	double * eigenvectors;
	double * sqrtEigenvalueMatrix;
	double * multMatrix;
};

#endif // GAUSS_VECTOR
