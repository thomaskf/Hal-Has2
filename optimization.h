/*
 *
 * optimization.h
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


#ifndef __RAL_RAS__optimization__
#define __RAL_RAS__optimization__

#include <iostream>
#include <iomanip>
#include "alignment.h"
#include "parameters.h"
#include "core.h"
#include "definitions.h"
#include "lbfgsb_new.h"
#include "tool_box.h"

#define MIN_values 1.0e-5
#define MAX_values 1.0e4
#define MIN_PI 1.0e-2
#define MAX_PI 1.0e2

#define NUMBER_OPT_MODES 8  // the number of different modes of the optimization method
#define RAS_OPT_MODE7    2  // For RAS, opt mode 7 will be the same as opt mode 2
#define RAS_OPT_MODE8    6  // For RAS, opt mode 8 will be the same as opt mode 6

class Optim {
public:
	VariableSet *vs;
	Alignment *algn;
	int* topMatrix;
	int isReversible;
	int paramType; // 1 : all edge lengths; 2 : s-matrix along an edge; 3: pi-vector along an endge
	int currRateMatID;
	int currRateCatID;

	int num_edges;
	int num_w;
	int num_rateCat;
	int num_rateMat;

	int num_chars;

	// for RAL (i.e. number of rate categories = 1)
	ParameterSet *ps;

	// for RAS (i.e. number of rate categories > 1)
	AllParameterSet *ps_set;
	int modelType; // type of RAL-RAS model, range of 1 - 5
	int isShare; //(for RAS) indicate whether it is shared among all the rate categories
	double* edgeScalars; // scalar multiplication for edges (for model 5)

	// constructor
	Optim(Alignment *algn, int *topMatrix, int isReversible, int num_chars);

	// constructor
	Optim(VariableSet *vs, ParameterSet *ps, Alignment *algn, int *topMatrix, int isReversible, int num_chars);

	// constructor
	Optim(VariableSet *vs, AllParameterSet *ps_set, Alignment *algn, int *topMatrix, int isReversible, int modelType, int num_chars);

	// destructor
	~Optim();

	// optimize all parameters
	// return the number of iterations processed
	int optimizeAllParam(int numIterations, double& loglikelihood, double& IC, ThreadLocks* threadLocks, int masterID, int threadID, int maxit, int mode, int ICType, int numRateChanges, int sameRootMat, int precise);

	// optimize all parameters for RAS
	// return the number of iterations processed
	// isGTRUpilson : true if it is GTR Upilson model
	// threadID : thread ID
	// Output: 1. loglikelihood; 2. IC; 3. df (degree of freedoms)
	int optimizeAllParamRAS(int numIterations, int numRateMat, double& loglikelihood,
		double& IC, int& df, int isGTRUpilson, int threadID, ThreadLocks* threadLocks,
		int maxit, int mode, int ICType, int sameRootMat, int numRateChanges, int precise);

	// set the variables and parameters
	void setVarParam(VariableSet *vs, ParameterSet *ps);

// for temparary use
	double* likelihoodMat;
	double* leftLikelihoodMat;
	double* rightLikelihoodMat;
	double* HASsiteLikelihoods; // the likelihoods of each site for each variable category (for HAS)
								// dimension: numUniqSites * numVarCategories
	
private:
	// optimizing the values
	void getOptValues(int type, int isShare, int maxit);
	// type -- 1 : all edge lengths; 2 : s-matrix along an edge; 3: pi-vector along an endge

	// optimize the root vector of each rate category
	// isShare : whether the root vector is shared
	void getOptRootVector(int isShare, int maxit);

	// optimize the alpha and invariant site
	void getOptAlphaInvariantSite(int maxit);

	// optimize the alpha and invariant site when beta = 0
	void getOptAlphaInvariantSiteZeroBeta(VariableSet *vs_zero_beta);
};


#endif /* defined(__RAL_RAS__optimization__) */
