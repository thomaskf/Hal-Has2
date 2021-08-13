/*
 *
 * RAL_mpi.h
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

#ifndef __RAL_RAS__RAL_NEW__
#define __RAL_RAS__RAL_NEW__

#include <mpi.h>
#include <iostream>
#include <fstream>
#include "alignment.h"
#include "core.h"
#include "core_mpi.h"
#include "charSet.h"
#include "tool_box.h"
#include "matrix.h"
#include "optimization.h"
#include "rateMatrixResult.h"
#include "definitions.h"
#include "parameters.h"
#include "user_options.h"
#include "RAL_common.h"


// number of the best rate matrix arrangement to check
// whether the edges are assigned for the same rate matrix
#define CONSISTENCY_CHECK_NUM 10

using namespace std;


// perform RAL algorithm
//
// options              : the user options
// rateGrpFile          : the rate group to start with
// preLogFile           : The output file of the previous run.
//                        One wants to resume the process from the previous run.
// world_size           : total number of processes ( > 1)
// world_rank           : which processor is runnning the function
void performRAL(UserOptions* options, char* rateGrpFile, char* preLogFile, int num_chars,
				int world_size, int world_rank);



// ============================================================================================================================
// ITERATION FUNCTIONS for RAL
// ============================================================================================================================


// perform RAL iteration for the corresponding rate matrices
// groupAllCandidates  : For each of the best rate matrices in this iteration,
//                       generate the potential candidates for next iteration
//                       1 - All the potential candidate form the same group to be processed in next iteration.
//                       2 - The potential candidates form separate groups according to the rate matrix generated from.
//                           (i.e. Original design)
void RALIterationNew (TempRateMats* rateMatrixList, int* topMatrix, int numLineTopMat,
		Alignment* alignment, int isReversible, int numSpecies,
		OptimalRateMats* optRateMatrix, int* leafID2line, int* internalID2line,
		double thresHeight, ofstream* fout, TempRateMats* preInterResult,
		TempRateMats* bestHundredRateMatrices, UserOptions* options, int num_chars,
		int saveToInterResult, int tagID, ParameterSet& start_ps, VariableSet& start_vs);

// ============================================================================================================================
// FUNCTIONS for each CPU thread
// ============================================================================================================================

// the Wrapper to perform "optimThread" in each CPU thread
void * optimThreadWrapper ( void * arg );

// The optimization process run by each CPU thread
void optimThread(int* topMatrix, int numLineTopMat, Alignment* alignment, int isReversible,
		int numSpecies, UserOptions* userOptions, int num_chars, int tagID);

// ============================================================================================================================
// OTHER FUNCTIONS
// ============================================================================================================================

// rearrange all items inside the matrices
// such that the topK items (among all items) are with the least IC values
void findConsistentRates(TempRateMats* bestHundredRateMatrices, vector<int>* consistentRates, int topK);

// recalculate all the IC values in the matrix
void reCalculateAllICs(TempRateMats* interResults, int ICType, int num_sites, int num_species, double logLDeviation, int* topMatrix);

// Load previous intermediate results from the log file
TempRateMats* loadPreInterResult(char* logFile, int num_w, int num_edges, int num_alpha, int ICType, 
		int num_sites, int num_species, double logLDeviation, int num_chars, int* topMatrix);

// ==================================================================================================
// Maximize the likelihood value for the specific matrice arrangement
// ==================================================================================================

// inputFormat: 1 - FASTA; 2 - PHY format
double getLogL(char* alignFile, char* treeFile, int inputFormat, char* rateGrpFile, int mode, int maxIT, int listOutParamDetail, int constantSiteExist, int includeGap, int num_chars);

// ==================================================================================================
// Compute the likelihood values for the corresponding parameter and variable values
// ==================================================================================================

double getLogLFrFiles(char* alignFile, char* treeFile, char* siteInfoFile, char* paramFilelist, int edgeRepresent, int numRateMatrix, int ICType, double& IC, int& df, int includeGap, int num_chars);

// ==================================================================================================
// use the consensus method to come up an answer
// ==================================================================================================

void consensusMethod(ofstream* fout, TempRateMats* allResults, UserOptions* options, double IC_Thres,
		int numSpecies, int nEdge, int numSites, int* topMatrix, vector<string>* leafList);


#endif
