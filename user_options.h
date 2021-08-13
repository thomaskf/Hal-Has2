/*
 *
 * user_options.h
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

#ifndef __RAL_RAS__user_options__
#define __RAL_RAS__user_options__

#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <cstdlib>
#include "tool_box.h"
#include "definitions.h"
#include "core.h"

using namespace std;

class GetOptions {
public:
	vector<char> flags;
	vector<string> values;

	int size();

	void read(int argc, char** argv);
};

class UserOptions {
public:

	// ======================================================
	// Variables
	// ======================================================

	// alignFile           : the alignment file
	//                       (compulsory)
	string alignFile;

	// topFile             : the topology file
	//                       (compulsory)
	string topFile;

	// rateMatrixFile      : the rate matrix file
	string rateMatrixFile;

	// inputFormat         : the format of the input multiple sequence alignment file
	//                       1 - FASTA format
	//                       2 - sequential PHYLIP format
	//                       (default = 1)
	int inputFormat;

	// numIterations       : the maximum number of iterations allowed
	//                       the iterations stop either when the loglikelihood value converages,
	//                       or the number exceeds the maximum number.
	//                       numIterations = -1, indicating no limit on the number of iterations.
	//                       (default = -1);
	int numIterations;

	// numSteps            : the maximum number of steps allowed for bottom-up Algorithm
	//                       numSteps = -1, indicating no limit on the number of iterations.
	//                       (default = -1);
	int numSteps;

	// numMaxCycles        : the maximum number of cycles for the new bottom-up algorithm
	//                       to indicate no limit on the maximum number of cycles, set numMaxCycles = -1
	//                       (default = -1);
	int numMaxCycles;

	// numCPUThreads       : the number of CPU threads is used (parallel processing)
	//                       (default = 4);
	int numCPUThreads;

	// outputPrefix        : prefix for output files
	//                       (default = <alignment file w/o .ext)
	string outputPrefix;

	// programType:
	//    1 - HAL-BU
	//    2 - HAL-TD
	//    3 - HAS
	//    4 - Upsilon
	//    5 - HAL-BU new version 1
	int programType;

	// HALMode:
	//    1 - BU
	//    2 - TD
	//    3 - BUTD
	//    4 - BU new version 1
	//    7 - John's
	int HALMode;

	// topK                : The number of best results will be proceeded in next iteration
	//                       (default = 3)
	int topK;

	// randInitParam       : 0 - the parameters are initialized by using default values;
	//                       1 - the parameters are randomly initialized.
	//                       (default = 0)
	int randInitParam;

	// groupAllCandidates  : For each of the best rate matrices in this iteration,
	//                       generate the potential candidates for next iteration
	//                       1 - All the potential candidate from the same-level iterations are grouped for the next-level iterations.
	//                       2 - All the potential candidate from the same-level iterations are NOT grouped for the next-level iterations.
	//                           (i.e. Original design)
	//                       (default = 1)
	int groupAllCandidates;

	// preLogFile[optional]: The output CHECKPOINT file of the previous run.
	//                       One wants to resume the process from the previous run.
	//                       (default = "")
	string preLogFile;
    
    // preLog[optional]: The output LOG file of the previous run.
    //                       One wants to resume the process from the previous run.
    //                       (default = "")
    string preLog;


	// minRateCat  : Maximum number of rate categories allowed
	//               (default = 2)
	int minRateCat;

	// maxRateCat  : Maximum number of rate categories allowed
	//               (default = 4)
	int maxRateCat;
	
	// stepRateCat : The step of rate categories
	//               (default = 1)
	// for example if minRateCat=3; maxRateCat=11; stepRateCat=2
	// then the program will check the number of rate categories: 3, 5, 7, 9, 11
	int stepRateCat;

	// optimization mode (default = 1)
	vector<int> opt_mode;

	// optimization maxIT (default = 1)
	vector<int> opt_maxIT;

	// range of deviations of resulting value of log-likelihood expected
	// default: 10.0
	double logLDeviation;

	// information criteria
	// 1 - AIC; 2 - Adjusted IC; 3 - BIC (default); 4 - CAIC
	int ic_list[TOTAL_NUM_IC]; // if AIC is selected, ic_list[0] = 1; if BIC is selected, then ic_list[2]=1; ....
	int info_criteria;

	// output the consensus
	int outputConsensus;

	// gap handling
	// 1 - ignore all columns with gaps; 2 - treat gaps as missing data (default)
	int gapHandling;
	
	// enable the checkpoint function
	int outChkPtFile;
	
	// the starting state of the optimization method
	// 0: from the original state
	// 1: from the tunned state of the parent configuration (obsoleted)
	// 2: optimization starts from the state using one model (default)
	int startOptState;
	
	// whether list out the values of the parameters for the best HAL model
	// 0: Not listing out the values of the parameters for the best HAL model
	// 1: listing out the values of the parameters for the best HAL model (default)
	int listParamBestHAL;
	
	// whether TD algorithm will continue after BU algorithm
	// 0: Not continue with TD algorithm after BU algorithm
	// 1: Continue with TD algorithm after BU algorithm (default)
	int continueTD;
	
	// whether it is unrooted tree
	int unrooted;

	vector<int>* startRateMatrix;

	vector<int>* startUntouchEdges;

	int maxNumRatMat;
	
	int precise; // More precise optimzation, but needs more time (default: 0)

	// =================================
	// Functions
	// =================================

	// constructor
	UserOptions(int programType);
	UserOptions(int programType, int HALMode);

	// set the default values
	void reset();
	
	// set the default values of mode and maxit
	void setDefaultModeMaxIT();

	// set the values of mode and maxit for precise option
	void setPreciseModeMaxIT();

	// To show the usage of the HAL program
	void outputHALBUUsage(char* progName, int isMPI);

	// To show the usage of the HAL TD program
	void outputHALTDUsage(char* progName, int isMPI);

	// To show the usage of the HAS program
	void outputHASUsage(char* progName, int isMPI);

	// To show the usage of the HAS program for Upsilon model
	void outputUpsilonUsage(char* progName);

	// read the arguments and collect the user options
	int readArguments(int argc, char** argv, string* errMsg);

	// To show the summary of the parameters
	void showSummary(int isMPI);

	// by default, prefixOut = <alignment file> w/o .ext
	void setDefaultOutputPrefix();
};


#endif /* defined(__RAL_RAS__user_options__) */
