/*
 *
 * RAS_mpi.h
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

#ifndef __RAL_RAS__RAS_MPI__
#define __RAL_RAS__RAS_MPI__

#include <mpi.h>
#include <iostream>
#include <vector>
#include "alignment.h"
#include "core.h"
#include "charSet.h"
#include "tool_box.h"
#include "matrix.h"
#include "optimization.h"
#include "rateMatrixResult.h"
#include "user_options.h"
#include "core_mpi.h"

typedef struct Quad {
	int modelType;
	int numCat;
	int optMode;
	int optMaxIT;
} Quad;

struct classcomp {
	bool operator() (const Quad& q1, const Quad& q2) const {
		if (q1.modelType != q2.modelType)
			return (q1.modelType < q2.modelType);
		else if (q1.numCat != q2.numCat)
			return q1.numCat < q2.numCat;
		else if (q1.optMode != q2.optMode)
			return q1.optMode < q2.optMode;
		else
			return q1.optMaxIT < q2.optMaxIT;
	}
};


// The data structure for holding the previous results from the log file
class PreHASLogResult {
public:
	int num_chars;
	vector<int> modelType;
	vector<int> numCat;
	vector<int> optMode;
	vector<int> optMaxIT;
	vector<double> loglikelihood;
	vector<int> df;
	vector<double> ICs;
	vector<int> ICTypes; // icType: type of Information Criteria (1 - AIC; 2 - Adjusted BIC; 3 - BIC; 4 - CAIC)
	vector<AllParameterSet*> PSs;
	vector<VariableSet*> VSs;
	map<Quad,int,classcomp> index;

	PreHASLogResult(int num_chars); // Constructor
	void readLogFile(char* logFile, int numEdges); // read the log file
	int getRecord(int modelType, int numCat, int optMode, int optMaxIT); // get the record according to the modelType and numCat
};

//================================================================================================


// perform RAS algorithm (MPI version)
//
// alignFile     : the alignment file
// topFile       : the topology file
// rateGrpFile   : the rate group file
// minNumRateCat : the min number of rate categories allowed
// maxNumRateCat : the max number of rate categories allowed
// numCpuThreads : number of CPU threads is used
// randomInt     : 0 - the parameters are initialized as default values; 1: initialized randomly
// numIterations : number of iterations
// prefixOut     : the prefix of the output file
// preLogFile    : the previous HAS log file
void performRAS(char* alignFile, char* topFile, char* rateGrpFile, int minNumRateCat, int maxNumRateCat,
		int randomInit, int numIterations, string prefixOut, char* preLogFile, int ICType,
		UserOptions* userOptions, int num_chars, int isGTRUpilson, int world_size, int world_rank);

// get the likelihood value
//
// alignFile    : the alignment file
// topFile      : the topology file
// paramFile    : the prefix of the parameter files
//                the parameter files: "paramFile"1, "paramFile"2, ...
// variableFile : the variable file
// numRateCat   : the total number of parameter files
// modelType    : the type of RAL-RAS model 1 - 5
// void getLikelihoodMx(char* alignFile, char* topFile, char* paramFile, char* variableFile, int numRateCat, int modelType);


// compute the number of substitutions for each edge
//
// topFile      : the topology file
// paramFile    : the prefix of the parameter files
//                the parameter files: "paramFile"1, "paramFile"2, ...
// numRateCat   : the total number of parameter files
void computeSubsEachEdge(char* topFile, char* paramFile, int numRateCat);

// The procedure for RAS executed by each process
void RASIter(Alignment* alignment, int* topMatrix, int numLineTopMat, int numIterations,
		int frNumRateCat, int toNumRateCat, int frModelType, int toModelType,
		vector<int>* rateMatArray, int numRateGrp, int randomInit, int isGTRUpilson, int* jobAllocated,
		double* bestIC, AllParameterSet* bestPs, VariableSet* bestVs, string* bestModel, ofstream* foutChkpt,
		double* summaryHASLogL, double* summaryHASIC, int* summaryHASDf, int ICType,
		UserOptions* userOptions, int num_chars, int totProcs, int tagID);

// The optimization process run by each CPU thread
void optimThread(int* topMatrix, int numLineTopMat, Alignment* alignment, int isReversible,
		int numSpecies, UserOptions* userOptions, int num_chars, int tagID, int randomInit,
		int isGTRUpilson, int threadID);
		
#endif /* defined(__RAL_RAS__RAS_MPI__) */
