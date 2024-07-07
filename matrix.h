/*
 *
 * matrix.h
 * HAL_HAS
 *
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * All rights reserved. CSIRO is willing to grant you a license to HAL-HAS on the terms of the GNU General Public
 * License version 3 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.html), except
 * where otherwise indicated for third party material.
 * The following additional terms apply under clause 7 of that license:
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES
 * OR OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF
 * LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS;
 * OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL,
 * INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
 * THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS,
 * OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT
 * SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER
 * GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY
 * OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF FOLLOWING
 * REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
 * (a)               THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT
 *                   SERVICES AGAIN;
 * (b)               THE REPAIR OF THE SOFTWARE;
 * (c)               THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE
 *                   RELEVANT SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED
 * WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 * Third Party Components
 * The following third party components are distributed with the Software.  You agree to comply with the license
 * terms for these components as part of accessing the Software.  Other third party software may also be identified
 * in separate files distributed with the Software.
 * ___________________________________________________________________
 * 
 * R : A Computer Language for Statistical Data Analysis version 3.0.1 (http://cran.r-project.org/src/base/R-3/R-3.0.1.tar.gz)
 * Copyright (C) 2000-2004 The R Core Team
 * This software is licensed under GNU GPL
 * 
 * JACOBI_EIGENVALUE.C (http://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.c)
 * Copyright (C) 2003-2013 John Burkardt
 * This software is licensed under GNU LGPL (http://www.gnu.org/licenses/lgpl.html)
 * ___________________________________________________________________
 */


#ifndef __RAL_RAS__matrix__
#define __RAL_RAS__matrix__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi_eigenvalue.h"

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#endif

// #include "eigenvalue_4_by_4.h"

using namespace std;

// The following paramaters are designed for Jacobi algorithm
// to find the eigenvalues and eigenvectors of a real symmetrix square matrix
#define MAX_ROUND 1000
// #define VERY_SMALL_NUMBER 1.0e-20

// construct a matrix
// the values inside are not resetted
double* matrix(int row, int col);

// reset the matrix to zero
void resetMat(double* mat, int row, int col);

// set all entries inside the matrix to a specific value
void setMatrixToVal(double* mat, int row, int col, double value);

// reset a square matrix to an identity matrix
void ident(double* sqMat, int dim);

// to check whether the square matrix is an identity matrix
int isIdent(double* sqMat, int dim, double allowance);

// set the values of the particular row of the matrix same as inArray
void setMatrix(double* matrix, int col, int row, double inArray[]);

// to set the diagional of the square matrix same as inArray
void setDiag(double* sqMatrix, double* inArray, int dim);

// output a diagional matrix from an array
double* toDiagMat(double* inArray, int dim);

// to collect the diagional of a square matrix
void diag(double* outArray, double* inSqMat, int dim);
double* diag(double* inSqMat, int dim);

// duplicate
void duplicate(double* toMat, double* frMat, int row, int col);
double* duplicate(double* frMat, int row, int col);

// get the first non-zero item
// return the position of the item, or -1 if none
int firstNonZeroPos(int* matrix, int dim);

// print out a matrix
void printMatrix(double* matrix, int row, int col);

//=======================================//
// Matrix operations                     //
//=======================================//

// to compute all the eigenvectors and eigenvalues for a real symmetric square matrix
// void computeEigenVect(double* sqMatrix, int dim, double* eigenVals, double* eigenVects, int& numRound);

// to compute the transpose of a matrix with dimension row x col
void transpose(double* outMat, double* inMat, int row, int col);
double* transpose(double* inMat, int row, int col);

// to compute the multiplication of two matrice
// dimension of matrix 1: row x common_len
// dimension of matrix 2: common_len x col
void multiply(double *outMat, double* inMat1, double* inMat2, int row, int common_len, int col);
double* multiply(double * inMat1, double* inMat2, int row, int common_len, int col);
// outMat cannot be the same as inMat1 or inMat2
void multiplyQuick(double *outMat, const double* inMat1, const double* inMat2, int row, int common_len, int col);
void multiplyQuick44(double *r, const double* a, const double* b);
void multiplyQuick44_t(double *r, const double* a_t, const double* b);
void multiplyQuick41(double *r, const double* a, const double* b);
// to compute the Euclidean distance between two row matrix
double euclidDist(double* rowMat1, double* rowMat2, int size);

// compute the distance matrix for the Euclidean distance between any two rows of the matrix
// and return a distance matrix
void computeDistMatrix(double* distMat, double* inMat, int row, int column);

// perform hierarchical clustering
// input: distMat (i.e. distance matrix) and size
// output: merge and height
// dimension of merge: (size-1) * 2
// dimension of height: 1 * (size-1)
void hierCluster(int* merge, double* height, double* distMat, int size);

#endif /* defined(__RAL_RAS__matrix__) */
