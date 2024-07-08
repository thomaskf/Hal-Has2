/*
 *
 * parameters.h
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

#ifndef __RAL_RAS__parameters__
#define __RAL_RAS__parameters__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "alignment.h"
#include "matrix.h"
#include "tool_box.h"
#include "definitions.h"
#include "simd.h"


#define IC_Thres_groupAll1 1.005 // The IC threshold when the groupAllCandidates = 1
#define IC_Thres_groupAll2 1.0 // The IC threshold when the groupAllCandidates = 2 or (Original's and John's algorithm)
#define max_MAXIT 1000
#define MIN_PROB_VALUE 1.0e-10 // the minimum value for the probability matrix

class VariableSet {
public:
	double beta;
	double* alpha;
	double* probXGivenInv;
	int num_alpha;
	int num_chars;
	double*  rootNodeFreq;

	// the following variables are for optimization
	double beta_t;
	double* alpha_t;
	double* probXGivenInv_t;
	double* rootNodeFreq_t;
	double*  rootNodeFreq_alpha;
	double* invar_constant;
	
	bool isGTRUpilson;

	VariableSet(int num_chars); // constructor
	VariableSet(int num_alpha, int num_chars); // constructor
	~VariableSet(); //destructor

	void setGTRUpilson(); // set "isGTRUpilson" to TRUE

	void readVariableFile(char* fileName); // read the values of the variables from the file

	void resetAllVariables(Alignment& alignment); // reset the variables according to the constant sites

	// set up the parameters for optimisation method
	void set_for_opt();

	void randomInitAllVariables(); // randomly initialize the variables

	// copy the whole content (except rateIDs) from another instance of ParameterSet
	void copyFrom(VariableSet& vs);

	void showContent();
	void showContent(ofstream* fout);
	void showContent(string &outStr, int format); // output the content to the string

	// roundup all the values
	void roundContent();

	//=================================================================
	// for reading the new format of parameter and variable files
	//=================================================================
	// read the values of the variables from the site info file    
	void readSiteInfoFile(char* fileName, vector<string>* siteCatNames);
	
	void compute_rootNodeFreq_alpha();
	
	void compute_invar_constant();

};

class ParameterSet {

public:
	double*  w;  // variables in the S vector for all edges
	double*  pi; // stationary probabilities for all edges
	double*  pi_t; // this is for optimization only
	double*  t;  // time for all edges
	int num_w; // number of variables in S vector
	int num_edges; // number of edges
	int num_chars;

	int num_rate_matrices; // number of rate matrices
	vector<vector<int>* > rateIDs; // the ID of rate matrices for each edge

	// double* allEigenMat;  // the eigen matrix for all edges
	// i.e. sqrt(pi) * R * inverse of sqrt(pi), where R = S * Pi

	double*  allEigenVals; // the eigen values for all edges
	double*  allEigenVects_t; // the transpose of eigen vectors for all edges

	double*  allCondProb_t; // transpose of conditional probability matrices for all edges
	// double*  allCondProb; // conditional probability matrices for all edges

	ParameterSet(int num_chars); // constructor
	ParameterSet(int num_w, int num_edges, int num_chars); // constructor
	~ParameterSet(); // destructor

	// read the values of w and pi from the file
	void readParamFile(char* fileName, int num_w, int num_edges);

	// load the file of rate groups
	// return the number of rate groups
	void loadRateMat(char* rateGrpFile, int* topMatrix, vector<string>* leafList);

	// import the rate matrix
	// prerequisite: the value of "num_edges" has to be set
	void loadRateMat(vector<int>& rateMat);

    // import the rate matrix string
    // prerequisite: the value of "num_edges" has to be set
    void loadRateMat(string rateMatStr);

	// reset the parameters
	void reset();

	// set up the parameters for optimisation method
	void set_for_opt();

	void initialize(int num_w, int num_edges, int isSingleGrp);
	// initialize the parameters
	// isSingleGrp : 1 - all edges are initialized as the same group
	//               0 - all edges are initizlized as different groups

	void randomInit(int num_w, int num_edges, int initSeed, int isSingleGrp);
	// initialize the parameters randomly
	// but still the same set of values for each edge
	// isSingleGrp : 1 - all edges are initialized as the same group
	//               0 - all edges are initizlized as different groups

	// update the parameters such that
	// outFormat = 1 : Edge length is set to the rate of substitution; OR
	//             2 : S6 in the rate matrix is set to 1
	void updateContent(int outFormat);
	void showContent();
	void showContent(ofstream* fout);
	void showContent(string &outStr, int* topMatrix, vector<string>* leafList); // output the content to the string

	// compute the eigenvalues and transpose of eigenvectors for the edges
	// if "edge" is -1, then do for all edges
	// output: allEigenVals and allEigenVects_t
	void computeAllEigenValues_t(int edge);

	// compute the eigenvalues and transpose of eigenvectors for the edges of a particular rate group
	// output: allEigenVals and allEigenVects_t
	void computeAllEigenValuesForRateGrp_t(int rateGrp);

	// OBJECTIVE:             compute the conditional probabilities for all the edges.
	// OUTPUT:                allCondProb
	// SIMILAR FUNCTION IN R: getCondProb
	// PREREQUISITE:          computeAllEigenMatrix()
	void computeAllCondProb(int isReversible);

	// OBJECTIVE:             compute the conditional probabilities for a particular edge
	// OUTPUT:                allCondProb
	// SIMILAR FUNCTION IN R: getCondProb
	// PREREQUISITE:          computeAllEigenMatrix()
	void computeCondProbForEdge(int edge);
	
	// to normalize S and the edge length
	// so that the edge length equals to the number of substitution
	// only for number of w = 6
	void normalize();

	void printAllCondProb();

	void printAllEigens();

	void print_t();

	// output the rate matrix in single array form
	// 1-based
	void outputRateMat(vector<int>& outRateMatArray);

	// update the rate matrix
	void updateRateMat(vector<int>& rateMatArray, int numRateGrp);

	// roundup the parameters based on the format
	// which the edge length is set to the rate of substitution
	void roundContent();

	void printRateMatID();

	// copy the whole content (except rateIDs) from another instance of ParameterSet
	void copyFrom(ParameterSet& ps);


	//==================================
	// For new format of parameter files
	//==================================
	// read the values of w and pi from the parameter file
	void readParamFile(char* fileName, vector<string>& nodeList, vector<double>& edgeLens, 
			int catID, int num_edges, int edgeRepresent);



private:
	// Temporary memory allocation
	// for internal computation use only
	double*  eigenMat;  // the eigen matrix
	double*  s;
	double*  tmpMat1;
	double*  tmpMat2;
	double*  sqrtPi;

	// get the initialize ID of rate matrices for each edge
	// before: the values of we have been read from the file
	void initializeRateMatIDFrParam();

	// initialize the rate matrix for each edge
	// isSingleGrp : 1 - all edges are initialized as the same group
	//               0 - all edges are initizlized as different groups
	void initializeRateMatID(int isSingleGrp);

};


class AllParameterSet {
public:
	vector<ParameterSet*> ps;

	vector<double*> allCondProbSet_t;

	int numRateCat;

	// constructor
	AllParameterSet();
	AllParameterSet(int numRateCat, int num_chars);
	AllParameterSet(int numRateCat, int num_chars, int num_edges);
	
	// destructor
	~AllParameterSet();

	// read the values of w and pi from all the parameter files
	void readParamFile(char* fileName, int num_w, int num_edges);

	// load the file of rate groups
	// return the number of rate groups
	int loadRateMat(char* rateGrpFile, int* topMatrix, vector<string>* leafList);

	// update the rate matrix
	void updateRateMat(vector<int>& rateMatArray, int numRateGrp);    

	// update the parameters such that
	// outFormat = 1 : Edge length is set to the rate of substitution; OR
	//             2 : S6 in the rate matrix is set to 1
	void updateContent(int outFormat);

	void showContent();
	void showContent(ofstream* fout);
	void showContent(string &outStr, int* topMatrix, vector<string>* leafList); // output the content to the string

	// compute the eigen-values and eigen-vectors for all edges for all parameter sets
	void computeAllEigenValues_t();

	// OBJECTIVE:             compute the conditional probabilities for all the edges for all parameter sets
	// OUTPUT:                allCondProb
	// SIMILAR FUNCTION IN R: getCondProb
	// PREREQUISITE:          computeAllEigenMatrix()
	void computeAllCondProb(int isReversible);

	// to normalize S and the edge length
	// so that when the number of substitution equals to the edge length
	// void normalizeAll();

	void printAllCondProb();

	void printAllEigens();

	// reset the parameters
	void reset();

	// set up the parameters for optimisation method
	void set_for_opt();

	void initialize(int num_w, int num_edges); // initialize all the parameters

	void randomInit(int num_w, int num_edges); // initialize all the parameters randomly

	int size();

	// copy the whole content (except rateIDs) from another instance of AllParameterSet
	void copyFrom(AllParameterSet& all_ps);

	//==================================
	// For new format of parameter files
	//==================================

	// read the values of w and pi from the parameter file list
	void readParamFileList(char* fileName, vector<string>& categoryList, vector<string>& nodeList, 
			vector<double>& edgeLens, int num_edges, int edgeRepresent);


};

// ===================================================
// the other functions
// ===================================================

// obtain the s-matrix
void computeSMatrix(double* s_matrix, double* w, double* pi);

// compute the eigen-matrix
void computeEigenMatrix(double* eigenMat, double* w, double* pi, double* s);

// print all conditional probabilities for all rate categories
void printAllCatCondProb(vector<ParameterSet*>& ps_set);

// ===================================================
// the following functions are defined in core.cpp
// ===================================================

// get the arrangement of rate matrix from the tree format
vector<int>* treeFormatToRateMatrix(int* topMatrix, string rateMatrixTree, vector<string>* leafList);

// load the rateGrpFile and collect the list of rate groups
// Input: (1) rateGrpFile; (2) topology matrix; (3) leafList
// Output: array of the rate groups (i.e. array of integers) and number of rate group
vector<int>* collectRateGrp(char* rateGrpFile, int* topMatrix, vector<string>* leafList, int& numRateGrp);


#endif /* defined(__RAL_RAS__parameters__) */
