/*
 *
 * core.h
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


#ifndef __RAL_RAS__CORE__
#define __RAL_RAS__CORE__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tool_box.h"
#include "definitions.h"
#include "matrix.h"
#include "alignment.h"
#include "charSet.h"
#include "rateMatrixResult.h"
#include "parameters.h"
#include "simd.h"

// #define SHOW_LIKEHOOD_EACH_SITE_CAT

// Input: Topology file
// Output: 1. Topology matrix; 2. number of rows of the matrix, 3. List of leaf names,
//         4. corresponding rate matrix, 5. edge length
void genTopMatrix(char* topFile, int** topMatrix, int* rowNum, vector<string>** leafList, vector<int>* rateMatrix, vector<double>* edgeLens);

// for the tree file using for simulation
// Input: Trees file
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length (# of items = # of edges * # of site categories)
//         4. List of node names, 5. List of category
// return: Topology matrix
int* genTopMatrix2(char* treesFile, int& rowNum, vector<string>& leafList, vector<double>& edgeLens, vector<string>& nodeList, vector<string>& siteCatNames);

// Output: 1. the mapping between the leaf ID to the line number in the paramater file, and
//         2. the mapping between the internal node ID to the line number in the parameter file
// Input: the topology matrix
void getMapping(int* leafID2line, int* internalID2line, int* topologyMatrix, int numOfLines);

// Output: Print the newick tree format with edge length (i.e. the edge length = # of substituation in the edge)
// Input: the topology matrix, the leaf names, and the edge lengths
void outTreeWithEdgeLen(int* topMatrix, vector<string>* leafList, vector<double>* edgeLens);

// Output: Print the newick tree format with edge length (i.e. the edge length = # of substituation in the edge)
// Input: the topology matrix, the leaf names, and the edge lengths
void outTreeWithEdgeLen(int* topMatrix, vector<string>* leafList, double* edgeLens, ofstream& fout);

// Output: Print the newick tree format with rate group ID
// Input: the topology matrix, the leaf names, and the rate group IDs
void outTreeWithGrpID(int* topMatrix, vector<string>* leafList, vector<int>* rateGrps);

// compute the likelihood matrix (dimension: internal nodes # x num_chars) of all the internal nodes for the specific site
// i.e. internal nodes # = number of lines in the topology matrix
// output: likelihoodMat
void computeLikelihoodMatrix(double* likelihoodMat, double* allCondProb_t, int* topologyMatrix, char* specificSites,
		int numLineTopMat, int num_chars);

// compute the likelihood matrix (dimension: internal nodes # x num_chars) of all the internal nodes for the specific site
// i.e. internal nodes # = number of lines in the topology matrix
// output: likelihoodMat
// Note: leftLikelihoodMat (size: num_chars) and rightLikelihoodMat (size: num_chars) are dummy
void computeLikelihoodMatrix(double* likelihoodMat, double* allCondProb_t, int* topologyMatrix, char* specificSites,
		int numLineTopMat, char* changedNodes, double* leftLikelihoodMat, double* rightLikelihoodMat,int num_chars);

// compute the total log-likelihood for a phylogenetic tree
double getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat,
		Alignment& alignment, VariableSet& variables, int num_chars);

// compute the total log-likelihood for a phylogenetic tree
double getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* likelihoodMat, double* leftLikelihoodMat, double* rightLikelihoodMat); //, char* changedNodes);

// compute the total log-likelihood for a phylogenetic tree for different variables
void getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, vector<VariableSet*>* variables, int num_chars, double* log_Ls);

// compute the total log-likelihood for a phylogenetic tree
// for many conditional probabilities
void getLogLikelihood(double* allCondProbArr_t, int n, vector<int>* corr_edges, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* result_Log_Ls, int gr_type);

// compute the total log-likelihood for a phylogenetic tree
// (for more than one rate category)
double getLogLikelihood(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars);

// compute the total log-likelihood for a phylogenetic tree
// (for more than one rate category)
// NEWLY ADDED: siteLikelihoods - the likelihoods of each nucleotide of each site for each variable category
//                                dimension: num_chars * numUniqSites * numVarCategories
//
//              updatedCat      - if "updatedCat" >= 0, only update the category "updatedCat"
//                                else update all the categories
//
//              temporary use:  likelihoodMat, leftLikelihoodMat, rightLikelihoodMat

double getLogLikelihood(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* siteLikelihoods, int updatedCat, double* likelihoodMat, double* leftLikelihoodMat, double* rightLikelihoodMat);

void getLogLikelihoodDetails(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars);

// get the name of the IC
string getICName(int ICType);

// get the IC-type of the IC name
int getICType(string ICName);

// compute the IC value with no regularization
double getIC_no_reg(double logVal, int numRates, int nsites, int nTaxa, int ICType, int sameRootMat);

// compute the IC value
double getIC(double logVal, int numRates, int nsites, int nTaxa, int ICType, int numRateChanges, int sameRootMat);

// compute the IC value
double getIC(double logVal, int numRates, int nsites, int nTaxa, int& df, int ICType, 
		int numRateChanges, int sameRootMat, int use_reg);

// compute the IC value for different RAL-RAS model
// modelType  : 1 - 5
// numRateCat : number of rate categories
// numRateMat : number of rate matrices
// isUpsilonModel : whether it is an upsilon model
// output: 1. df - degree of freedom; 2. bic - IC value
double getIC(double logVal, int numRateMat, int nsites, int nTaxa, int numRateCat,
		int modelType, int& df, int isUpsilonModel, int ICType, 
		int numRateChanges, int sameRootMat, int use_reg);

// load the rateGrpFile and collect the list of rate groups
// Input: (1) rateGrpFile; (2) topology matrix; (3) leafList
// Output: array of the rate groups (i.e. array of integers) and number of rate group
vector<int>* collectRateGrp(char* rateGrpFile, int* topMatrix, vector<string>* leafList, int& numRateGrp);

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, vector<string>* leafList);

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, string rateMatrix, vector<string>* leafList);

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, vector<string>* leafList, vector<int>* untouchNodes, int currNode);

// get the arrangement of rate matrix from the tree format
vector<int>* treeFormatToRateMatrix(int* topMatrix, string rateMatrixTree, vector<string>* leafList);

// initialize the thread locks
void initThreadLocks(ThreadLocks* threadLocks);

// return the number of rate groups of a given rate matrix
int getNumRateGrp(vector<int>& rateMatrix);

// read the HAS result file
// output: 1. number of site categories; 2. HALHAS model number
void readHASResultfile(char* file, int& catNum, int& modelNum);

//===========================================
// for loading the files used for simulation
//===========================================

// Input: topology string (can also include the rate group and the edge lengths)
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length, 4. List of node names
// Return: Topology matrix
int* genTopMatrixFrStr2(string newickStr, int& rowNum, vector<string>& leafList,
		vector<double>& edgeLens, vector<string>& nodeList);


//============================================
// to check how many changes on the rate matrices (convergent or non-convergent)
// across the tree
//============================================
int numChanges(int* topMatrix, vector<int>* rateMatrix, int line, int parentRate);
// int numChanges(int* topMatrix, vector<int>* rateMatrix, int line, int parentRate, vector<string>* leafList, string& str);

#endif /* defined(__RAL_RAS__BUP__) */
