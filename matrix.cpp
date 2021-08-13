/*
 *
 * matrix.cpp
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


#include "matrix.h"

// construct a matrix
// the values inside are not resetted
double* matrix(int row, int col) {
	return (double*) malloc (sizeof(double)*row*col);
}

// reset the matrix to zero
void resetMat(double* mat, int row, int col) {
	memset(mat,0,sizeof(double)*row*col);
}

// reset a square matrix to an identity matrix
void ident(double* sqMat, int dim) {
	memset(sqMat,0,sizeof(double)*dim*dim);
	for (int i=0; i<dim; i++)
		sqMat[i*dim+i] = 1.0;
}

// to check whether the square matrix is an identity matrix
// "allowance" is a small number, like 1.0e-3
int isIdent(double* sqMat, int dim, double allowance) {
	allowance = fabsl(allowance);
	for (int i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {
			if (i != j) {
				if (fabsl(sqMat[i*dim +j]) > allowance) {
					return 0;
				}
			} else {
				if (fabsl(sqMat[i*dim +j]-1.0) > allowance) {
					return 0;
				}
			}
		}
	}
	return 1;
}

// set the values of the particular row of the matrix same as inArray
void setMatrix(double* mat, int col, int row, double inArray[]) {
	double* fr_ptr = &(inArray[0]);
	double* to_ptr = &(mat[row*col]);
	memcpy(to_ptr,fr_ptr,col*sizeof(double));
}

// set all entries inside the matrix to a specific value
void setMatrixToVal(double* mat, int row, int col, double value) {
	for (int i=0; i<row*col; i++) {
		mat[i] = value;
	}
}

// to set the diagional of the square matrix same as inArray
void setDiag(double* sqMatrix, double* inArray, int dim) {
	for (int i=0; i<dim; i++)
		sqMatrix[i*dim+i]=inArray[i];
}

// output a diagional matrix from an array
double* toDiagMat(double* inArray, int dim) {
	double* result = matrix(dim,dim);
	resetMat(result, dim, dim);
	setDiag(result, inArray, dim);
	return result;
}

// to collect the diagional of a square matrix
void diag(double* outArray, double* inSqMat, int dim) {
	for (int i=0; i<dim; i++)
		outArray[i] = inSqMat[i*dim+i];
}
double* diag(double* inSqMat, int dim) {
	double* result = matrix(1, dim);
	diag(result,inSqMat,dim);
	return result;
}

// duplicate
void duplicate(double* toMat, double* frMat, int row, int col) {
	memcpy(toMat,frMat,row*col*sizeof(double));
}
double* duplicate(double* frMat, int row, int col) {
	double* toMat = matrix(row, col);
	memcpy(toMat,frMat,row*col*sizeof(double));
	return toMat;
}

// get the first non-zero item
// return the position of the item, or -1 if none
int firstNonZeroPos(int* matrix, int dim) {
	for (int i=0; i<dim; i++) {
		if (matrix[i] != 0)
			return i;
	}
	return -1;
}

// print out a matrix
void printMatrix(double* matrix, int row, int col){
	for (int i=0; i<row; i++) {
		for (int j=0; j<col; j++) {
			cout << matrix[i*col+j] << "\t";
		}
		cout << endl;
	}
}

/*
// compute eigen vector
void computeEigenVect(double* sqMatrix, int dim, double* eigenVals, double* eigenVects, int& numRound)
{
	int rot_num;
	// int success;
	// int i,j;

	jacobi_eigenvalue ( dim, sqMatrix, MAX_ROUND, eigenVects, eigenVals, &numRound, &rot_num );

	// transpose eigenVects
	double tmp;
	// int i, j;
	for (i=0; i<dim; i++) {
		for (j=i+1; j<dim; j++) {
			tmp = eigenVects[i*dim+j];
			eigenVects[i*dim+j] = eigenVects[j*dim+i];
			eigenVects[j*dim+i] = tmp;
		}
	}
}*/

// to compute the transpose of a matrix with dimension row x col
void transpose(double* outMat, double* inMat, int row, int col) {
	for (int i=0; i<row; i++) {
		for (int j=0; j<col; j++) {
			outMat[j*row+i] = inMat[i*col+j];
		}
	}
}
double* transpose(double* inMat, int row, int col) {
	double* result = matrix(col,row);
	transpose(result, inMat, row, col);
	return result;
}

// to compute the multiplication of two matrice
// dimension of matrix 1: row x common_len
// dimension of matrix 2: common_len x col
void multiply(double *outMat, double* inMat1, double* inMat2, int row, int common_len, int col){
	double* result = new double[row*col];
	memset(result, 0, sizeof(double)*row*col);
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			for (int k=0; k<common_len; k++)
				result[i*col+j] += inMat1[i*common_len+k]*inMat2[k*col+j];
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			outMat[i*col+j] = result[i*col+j];
	delete[] result;
}

// outMat cannot be the same as inMat1 or inMat2
void multiplyQuick(double *outMat, const double* inMat1, const double* inMat2, int row, int common_len, int col){
	memset(outMat, 0, sizeof(double)*row*col);
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			for (int k=0; k<common_len; k++)
				outMat[i*col+j] += inMat1[i*common_len+k]*inMat2[k*col+j];
}

void multiplyQuick44(double *r, const double* a, const double* b){
	for (int i=0; i<16; i+=4)
		for (int j=0; j<4; j++)
			r[i+j] = a[i]*b[j] + a[i+1]*b[j+4] + a[i+2]*b[j+8] + a[i+3]*b[j+12];
}

void multiplyQuick44_t(double *r, const double* a_t, const double* b){
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			r[i*4+j] = a_t[i]*b[j] + a_t[i+4]*b[j+4] + a_t[i+8]*b[j+8] + a_t[i+12]*b[j+12];
}

void multiplyQuick41(double *r, const double* a, const double* b){
	/*
	float rr[4];
	__m128 a_line, b_line, r_line;
	// b_line = _mm_set_ps1((float)b[0]);
	b_line = _mm_set_ps((float)b[0],(float)b[0],(float)b[0],(float)b[0]);
	a_line = _mm_set_ps((float)a[0],(float)a[4],(float)a[8],(float)a[12]);
	r_line = _mm_mul_ps(a_line, b_line);
	for (int i=1; i<4; i++) {
		// b_line = _mm_set_ps1((float)b[i]);
		b_line = _mm_set_ps((float)b[i],(float)b[i],(float)b[i],(float)b[i]);
		a_line = _mm_set_ps((float)a[i],(float)a[i+4],(float)a[i+8],(float)a[i+12]);
		r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line);
	}
	_mm_store_ps(rr, r_line);
	for (int i=0; i<4; i++)
		r[i] = rr[i];
	 */
	for (int i=0; i<16; i+=4)
		r[i/4] = a[i]*b[0] + a[i+1]*b[1] + a[i+2]*b[2] + a[i+3]*b[3];
}



double* multiply(double * inMat1, double* inMat2, int row, int common_len, int col) {
	double* result = matrix(row,col);
	multiplyQuick(result,inMat1,inMat2,row,common_len,col);
	return result;
}

// to compute the Euclidean distance between two row matrix
double euclidDist(double* rowMat1, double* rowMat2, int size) {
	double dist = 0;
	double diff;
	for (int i=0; i<size; i++) {
		diff = fabs(rowMat1[i] - rowMat2[i]);
		dist += diff*diff;
	}
	return sqrt(dist);
}

// compute the distance matrix for the Euclidean distance between any two rows of the matrix
// and return a distance matrix
void computeDistMatrix(double* distMat, double* inMat, int row, int column) {

	for (int i=0; i<row; i++)
		distMat[i*row+i]=0.0;

	for (int i=0; i<row-1; i++) {
		for (int j=i+1; j<row; j++) {
			distMat[i*row+j] = euclidDist(&(inMat[i*column]),&(inMat[j*column]),column);
			distMat[j*row+i] = distMat[i*row+j];
		}
	}
}


// create a new distance matrix when the i-th item and j-th item are combined
// assume size>=2
double* updateDistMat(double* distMat, int i, int j, int size) {
	int new_m=0;
	int new_n=0;
	int new_size=size-1;
	double *newDistMat = new double[new_size*new_size];
	int m,n;
	for (m=0; m<size; m++) {
		if (m==i || m==j)
			continue;
		new_n=0;
		for (n=0; n<size; n++) {
			if (n==i || n==j)
				continue;
			newDistMat[new_m*new_size+new_n] = distMat[m*size+n];
			new_n++;
		}
		// consider for the new item (grouping with i and j)
		double dist_i = distMat[m*size+i];
		double dist_j = distMat[m*size+j];
		if (dist_j > dist_i)
			newDistMat[new_m*new_size+(new_size-1)]=dist_j;
		else
			newDistMat[new_m*new_size+(new_size-1)]=dist_i;
		newDistMat[(new_size-1)*new_size+(new_m)]=newDistMat[new_m*new_size+(new_size-1)];
		newDistMat[(new_size-1)*new_size+(new_size-1)]=0.0;
		new_m++;
	}
	return newDistMat;
}

// create a new label when the i-th item and j-th item are combined
// assum size>=2
int* updateLabel(int* label, int i, int j, int size) {
	int new_size = size-1;
	int* newLabel = new int[new_size];
	int new_m=0;
	int largest_label=0;
	for (int m=0; m<size; m++) {
		if (largest_label < label[m])
			largest_label = label[m];
		if (m==i || m==j)
			continue;
		newLabel[new_m] = label[m];
		new_m++;
	}
	newLabel[new_size-1] = largest_label+1;
	return newLabel;
}

void sub_hierCluster(int* merge, double* height, double* distMat, int* label, int size) {

	if (size<=1)
		return;

	// get the smallest dist
	double small_dist = distMat[1];
	int small_i=0;
	int small_j=1;
	for (int i=0; i<size-1; i++) {
		for (int j=i+1; j<size; j++) {
			if (small_dist > distMat[i*size+j]) {
				small_dist = distMat[i*size+j];
				small_i = i;
				small_j = j;
			}
		}
	}
	merge[0] = label[small_i];
	merge[1] = label[small_j];
	height[0] = small_dist;

	double* newDistMat = updateDistMat(distMat, small_i, small_j, size);
	int *newLabel = updateLabel(label, small_i, small_j, size);

	sub_hierCluster(&(merge[2]), &(height[1]), newDistMat, newLabel, size-1);

	delete[] newDistMat;
	delete[] newLabel;
}

// perform hierarchical clustering
// input: distMat (i.e. distance matrix) and size
// output: merge and height
// dimension of merge: (size-1) * 2
// dimension of height: 1 * (size-1)
void hierCluster(int* merge, double* height, double* distMat, int size) {

	int* label = new int[size];
	for (int m=0; m<size; m++) {
		label[m] = -m-1;
	}

	sub_hierCluster(merge, height, distMat, label, size);

	delete[] label;
}
