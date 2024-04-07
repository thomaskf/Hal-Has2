/*
 *
 * user_options.cpp
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

#include "user_options.h"

// read the arguments and collect the options
void GetOptions::read(int argc, char** argv) {

	int i;
	char* currArg;
	int currLen;
	char flag;
	string value;
	for (i=1; i<argc; i++) {
		// skip the first argument
		currArg = argv[i];
		currLen = (int) strlen(currArg);
		if (currArg[0]=='-' && currLen>1) {
			// get the flag
			flag=currArg[1];
			// get the value
			if (currLen>2) {
				value = (&currArg[2]);
			} else if (argc > i+1 && argv[i+1][0]!='-') {
				i++;
				value = argv[i];
			} else {
				value = "";
			}
			flags.push_back(flag);
			values.push_back(value);
		}
	}

}

int GetOptions::size() {

	return (int) flags.size();
}

// constructor
UserOptions::UserOptions(int programType) {
	this->programType = programType;
	this->HALMode = 0;
	// output consensus
	this->outputConsensus = 0;
	// initialize the ic_list
	memset(ic_list, 0, 4*sizeof(int));
	this->outChkPtFile = 1;
	setDefaultModeMaxIT();
}

UserOptions::UserOptions(int programType, int HALMode) {
	this->programType = programType;
	this->HALMode = HALMode;
	// output consensus
	this->outputConsensus = 0;
	// initialize the ic_list
	memset(ic_list, 0, 4*sizeof(int));
	this->outChkPtFile = 1;
	setDefaultModeMaxIT();
}

// set the default values
void UserOptions::reset() {
	alignFile = "";
	topFile = "";
	rateMatrixFile = "";
	inputFormat = 2;
	numIterations = -1;
	numMaxCycles = -1;
	numCPUThreads = 4;
	topK = DEFAULT_TOP_K;
	if (HALMode == 1 || HALMode == 7 || HALMode == 8) {
		// original BU algorithm
		topK = DEFAULT_TOP_K_OLD;
	}
	randInitParam = 0;
	// All the potential candidate from the same-level iterations are grouped for the next-level iterations.
	groupAllCandidates = 1;
	if (HALMode == 1 || HALMode == 7) {
		groupAllCandidates = 2;
	}
	preLogFile = "";
    preLog = "";
	outputPrefix = "";
	minRateCat = 1;
	maxRateCat = 4;
	stepRateCat = 1;
	startRateMatrix=NULL;
	startUntouchEdges=NULL;
	maxNumRatMat = 0;
	logLDeviation = DEFAULT_LOG_L_THRES;
	info_criteria = DEFAULT_IC; // 1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC
	// output consensus
	this->outputConsensus = 0;
	// by default, gaps are treated as missing data
	gapHandling = 2;
	outChkPtFile = 1;
	if (HALMode == 4 || HALMode == 7 || HALMode == 8) {
		startOptState = DEFAULT_OPT_START;
	} else {
		// original BU algorithm / TD
		startOptState = DEFAULT_OPT_START_OLD;
	}
	listParamBestHAL = 1;
	continueTD = 1;
	unrooted = 0;
	numSteps = -1;
	precise = 0;
	setDefaultModeMaxIT();
}

// set the default values of mode and maxit
void UserOptions::setDefaultModeMaxIT() {
	vector<string> token;
	int i;
	// optimization mode
	string value = DEFAULT_OPT_MODE;
	tokenizer(value, ",", &token);
	opt_mode.clear();
	for (i=0; i<token.size(); i++)
		opt_mode.push_back(atoi(token[i].c_str()));
	// optimization maxIT
	value = DEFAULT_OPT_MAXIT;
	tokenizer(value, ",", &token);
	opt_maxIT.clear();
	for (i=0; i<token.size(); i++)
		opt_maxIT.push_back(atoi(token[i].c_str()));
}

// set the values of mode and maxit for precise option
void UserOptions::setPreciseModeMaxIT() {
	vector<string> token;
	int i;
	// optimization mode
	string value = PRECISE_MODE;
	tokenizer(value, ",", &token);
	opt_mode.clear();
	for (i=0; i<token.size(); i++)
		opt_mode.push_back(atoi(token[i].c_str()));
	// optimization maxIT
	value = PRECISE_MAXIT;
	tokenizer(value, ",", &token);
	opt_maxIT.clear();
	for (i=0; i<token.size(); i++)
		opt_maxIT.push_back(atoi(token[i].c_str()));
}


// To show the usage of the HAL BU program
void UserOptions::outputHALBUUsage(char* progName, int isMPI) {
	cout << "================================================================================" << endl;
    
    cout << "Welcome to";
    
    if (HALMode==1)
        cout << "                  Welcome to HAL";
    else
        cout << "                  Welcome to HAL-NEW";
    
    // show MPI version if so
    if (isMPI)
        cout << " MPI";
    
    // show regularization if so
    if (C_REG == 1)
        cout << " (regularization)";
    
    // show version
    cout << " Version " << VERSION << endl << endl;

	if (isMPI) {
		cout << "Syntax: mpirun -n <machines> " << progName << " <alignment file> <topology file> <other options>" << endl;
		cout << endl;
	} else {
		cout << "Syntax: " << progName << " <alignment file> <topology file> <other options>" << endl;
		cout << "        " << progName << " -h" << endl;
		cout << endl;
	}
	
	if (isMPI)
		cout << "  <machines>          : Number of machines ( > 1 )" << endl;
	cout << "  <alignment file>    : Multiple alignment file" << endl;
	cout << "  <topology file>     : Topology file in Newick format" << endl << endl;
	cout << "other options:" << endl;
	cout << "  -f <format of file> : The format of the multiple alignment file" << endl;
	cout << "                        1 - FASTA format" << endl;
	cout << "                        2 - Sequential PHYLIP format (default)" << endl;
	if (!isMPI)
		cout << "  -u <# of CPUs>      : Number of CPU threads used (default: 4)" << endl;
	cout << "  -o <output prefix>  : Prefix for output files" << endl;
	cout << "                         (default: <alignment file> w/o .ext)" << endl;
	cout << "  -r <checkpoint file>: Resume from the last execution;" << endl;
	cout << "                        <checkpoint file>, which was outputted from the last" << endl;
	cout << "                        execution of the program, stores all the immediate" << endl;
	cout << "                        results (with file extension: '.chkpt.txt')" << endl;
    cout << "  -z <previous log>   : Load the log file too (for testing only)" << endl;
	cout << "  -i <# of iterations>: Number of iterations performed for each-time" << endl;
	cout << "                        parameter tuning" << endl;
	cout << "                        (default: -1 [no limit, run until converge])" << endl;
	cout << "  -y <matrix file>    : Optimize the parameters according to the matrix file" << endl;
	cout << "                        (no search algorithm will be performed)" << endl;
	cout << "  -b <info criteria>  : Information Criteria (IC) to be used (default: " << DEFAULT_IC << ")" << endl;
	cout << "                        1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC" << endl;
	cout << "  -g <gap handling>   : 1 - ignore all columns with gaps or ambiguous chars;" << endl;
	cout << "                        2 - treat gaps as missing data (default)" << endl;
	cout << "  -precise            : More precise optimization, but needs more time." << endl;
	cout << "                        (default: disable)" << endl;

#ifndef RELEASE_VERSION

	cout << "  -v <grp candidates> : 1 - All the potential candidates from the same-level" << endl;
	cout << "                            iterations are grouped and top K candidates are" << endl;
	cout << "                            selected to proceed to the next level;" << endl;
	cout << "                        2 - Top K candidates from one iteration are selected" << endl;
	cout << "                            to proceed to the next level - original design." << endl;
	cout << "                            (this can slow down the program significantly)" << endl;
	if (HALMode == 1 || HALMode == 7) {
		cout << "                        (default: 2)" << endl;
	} else {
		cout << "                        (default: 1)" << endl;
	}
	if (HALMode==1) {
		// original BU
		cout << "  -k <topK value>     : 1 - 100, top K value (default: " << DEFAULT_TOP_K_OLD << ")" << endl;
		cout << "  -p <opt start mode> : 0 - optimization always starts from the original state" << endl;
		// cout << "                        1 - optimization starts from the previous iteration (obsoleted)" << endl;
		cout << "                        2 - optimization starts from the state using one model" << endl;
		cout << "                        (default: " << DEFAULT_OPT_START_OLD << ")" << endl;
	} else {
		// HALMode == 4 (new BU)
		cout << "  -c <max # of cycles>: Maximum number of cycles to be repeated for " << endl;
		cout << "                        the new HAL bottom-up algorithm" << endl;
		cout << "                        (default: -1 [no limit, run until converge])" << endl;
		cout << "  -e <optim mode>     : 1 - 9, Optimization mode (default: " << DEFAULT_OPT_MODE << ")" << endl;
		cout << "                        (for more than one, format: mode1,mode2,...)" << endl;
		cout << "  -j <optim maxit>    : 1 - 100, Optimization maxIT (default: " << DEFAULT_OPT_MAXIT << ")" << endl;
		cout << "                        (for more than one, format: maxit1,maxit2,...)" << endl;
		cout << "  -k <topK value>     : 1 - 100, top K value (default: " << DEFAULT_TOP_K << ")" << endl;
		// cout << "  -d <logL deviation>  : Deviation range of log-likelihood (default: " << DEFAULT_LOG_L_THRES << ")" << endl;
		cout << "  -p <opt start mode> : 0 - optimization always starts from the original state" << endl;
		// cout << "                        1 - optimization starts from the previous iteration (obsoleted)" << endl;
		cout << "                        2 - optimization starts from the state using one model" << endl;
		cout << "                        (default: " << DEFAULT_OPT_START << ")" << endl;
		// cout << "  -s                   : Use consensus approach (default: disable)" << endl;
	}
	
#endif
	
	/*
	cout << "  -l                   : 0 - Not listing out the parameters of the best HAL model" << endl;
	cout << "                         1 - Listing out the parameters of the best HAL model" << endl;
	cout << "                         (default: 1)" << endl;
	*/
	
	cout << "  -q                  : 0 - Not run top-down algorithm after bottom-up algorithm" << endl;
	cout << "                        1 - Running top-down algorithm after bottom-up algorithm" << endl;
	cout << "                        (default: 1)" << endl;

	/*
	cout << "  -x                   : 0 - rooted tree (default)" << endl;
	cout << "                         1 - unrooted tree" << endl;
	*/
	if (!isMPI)
		cout << "  -h                  : This help page" << endl << endl;

	if (HALMode==1) {
		// original BU
		cout << "Output files:" << endl;
		cout << "1. <output prefix>.BU.chkpt.txt - all the intermediate results (for" << endl;
		cout << "                                  resuming the program when necessary)" << endl;
		cout << "2. <output prefix>.BU.<IC>.result.txt - the resulting optimal rate matrices" << endl;
		cout << "                                        from the bottom-up algorithm" << endl;
		cout << "3. <output prefix>.BU.HAL.<IC>.result.txt - the resulting optimal rate matrices" << endl;
		cout << "                                            after both the bottom-up & the top-" << endl;
		cout << "                                            down algorithm" << endl<<endl;
	} else {
		// new BU
		cout << "Output files:" << endl;
		cout << "1. <output prefix>.NEW-BU.chkpt.txt - all the intermediate results (for" << endl;
		cout << "                                      resuming the program when necessary)" << endl;
		cout << "2. <output prefix>.NEW-BU.<IC>.result.txt - the resulting optimal rate matrices" << endl;
		cout << "                                            from the bottom-up algorithm" << endl;
		cout << "3. <output prefix>.NEW-BU.HAL.<IC>.result.txt - the resulting optimal rate" << endl;
		cout << "                                            matrices after both the bottom-up &" << endl;
		cout << "                                            the top-down algorithm" << endl<<endl;
	}

	cout << "Example:" << endl;
	cout << "  Alignment file: data.phy" << endl;
	cout << "  Topology file : tree.txt" << endl;
	if (isMPI)
		cout << "  # of machines: 4" << endl;

	cout << "  To execute the " << progName << " program using the default parameters:" << endl;
		
	if (isMPI)
		cout << "      $ mpirun -n 4 " << progName << " data.phy tree.txt" << endl << endl;
	else
		cout << "      $ " << progName << " data.phy tree.txt" << endl << endl;

	cout << "Contact: " << CONTACTPERSON << endl;
	cout << "================================================================================" << endl;
}


// To show the usage of the HAL TD program
void UserOptions::outputHALTDUsage(char* progName, int isMPI) {
	cout << "================================================================================" << endl;
	if (isMPI)
		cout << "                  Welcome to HAL-Top-Down MPI Version " << VERSION << endl << endl;
	else
		cout << "                  Welcome to HAL-Top-Down Version " << VERSION << endl << endl;

	if (isMPI) {
		cout << "Syntax: mpirun -n <machines> " << progName << " <alignment file> <topology file> <BU result file> <other options>" << endl;
		cout << endl;
	} else {
		cout << "Syntax: " << progName << " <alignment file> <topology file> <BU result file> <other options>" << endl;
		cout << "        " << progName << " -h" << endl << endl;
	}
	if (isMPI)
		cout << "  <machines>           : Number of machines ( > 1 )" << endl;
	cout << "  <alignment file>     : Multiple alignment file" << endl;
	cout << "  <topology file>      : Topology file in Newick format" << endl;
	cout << "  <BU result file>     : The resulting file from HAL-Bottom-Up program" << endl << endl;
	cout << "other options:" << endl;
	cout << "  -f <format of file>  : The format of the multiple alignment file" << endl;
	cout << "                         1 - FASTA format" << endl;
	cout << "                         2 - Sequential PHYLIP format (default)" << endl;
	if (!isMPI)
		cout << "  -u <# of CPUs>       : Number of CPU threads used (default: 4)" << endl;
	cout << "  -o <output prefix>   : Prefix for output files" << endl;
	cout << "                         (default: <alignment file> w/o .ext)" << endl;
	cout << "  -r <checkpoint file> : Resume from the last execution;" << endl;
	cout << "                         <checkpoint file>, which was outputted from the last" << endl;
	cout << "                         execution of the program, stores all the immediate" << endl;
	cout << "                         results (with file extension: '.chkpt.txt')" << endl;
	cout << "  -a                   : not outputting the checkpoint file" << endl;
	cout << "                         (this may increase the program speed, but the program" << endl;
	cout << "                          cannot be resumable)" << endl;
	cout << "  -b <info criteria>   : Information Criteria to be used (default: " << DEFAULT_IC << ")" << endl;
	cout << "                         1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC" << endl;
	cout << "  -i <# of iterations> : Number of iterations performed for each-time" << endl;
	cout << "                         parameter tuning" << endl;
	cout << "                         (default: -1 [no limit, run until converge])" << endl;
	if (!isMPI)
		cout << "  -h                   : This help page" << endl << endl;

	cout << "Output files:" << endl;
	cout << "  1. <output prefix>.TD.chkpt.txt  - all the intermediate results (for resuming" << endl;
	cout << "                                     the program when necessary)" << endl;
	cout << "  2. <output prefix>.TD.result.txt - the resulting optimal rate matrices" << endl << endl;

	cout << "Example:" << endl;
	cout << "  Example alignment file: data.phy" << endl;
	cout << "  Example topology file : tree.txt" << endl;
	cout << "  HAL_BU result file (which will appear after execution of HAL_BU):" << endl;
	cout << "      data.BU.result.txt" << endl;
	if (isMPI)
		cout << "  # of machines: 4" << endl;
	cout << "  To execute the HAL_TD program using the default parameters:" << endl;
	if (isMPI)
		cout << "      $ mpirun -n 4 " << progName << " data.phy tree.txt data.BU.result.txt" << endl << endl;
	else
		cout << "      $ " << progName << " data.phy tree.txt data.BU.result.txt" << endl << endl;

	cout << "Contact: " << CONTACTPERSON << endl;
	cout << "================================================================================" << endl;
}


// To show the usage of the HAS program
void UserOptions::outputHASUsage(char* progName, int isMPI) {
	cout << "================================================================================" << endl;
    cout << "                  Welcome to HAS";

    // show MPI version if so
    if (isMPI)
        cout << " MPI";
    
    // show regularization if so
    if (C_REG == 1)
        cout << " (regularization)";
    
    // show version
    cout << " Version " << VERSION << endl << endl;
    
	if (isMPI) {
		cout << "Syntax: mpirun -n <machines> " << progName << " <alignment file> <topology file> <other options>" << endl;
		cout << endl;
	} else {
		cout << "Syntax: " << progName << " <alignment file> <topology file> <other options>" << endl;
		cout << "        " << progName << " -h" << endl << endl;
	}
	if (isMPI)
		cout << "  <machines>           : Number of machines ( > 1 )" << endl;
	cout << "  <alignment file>     : Multiple alignment file" << endl;
	cout << "  <topology file>      : Topology file in Newick format" << endl;
	cout << endl;
	cout << "other options:" << endl;
	cout << "  -t <HAL result file> : The resulting file from HAL-NEW program" << endl;
	cout << "                         If no HAL result file is provided, then all edges" << endl;
	cout << "                         are assumed to be the same rate group" << endl;
	cout << "  -f <format of file>  : The format of the multiple alignment file" << endl;
	cout << "                         1 - FASTA format" << endl;
	cout << "                         2 - Sequential PHYLIP format (default)" << endl;
	if (!isMPI)
		cout << "  -u <# of CPUs>       : Number of CPU threads used (default: 4)" << endl;
	cout << "  -o <output prefix>   : Prefix for output files" << endl;
	cout << "                         (default: <alignment file> w/o .ext)" << endl;
	cout << "  -r <checkpoint file> : Resume from the last execution;" << endl;
	cout << "                         <checkpoint file>, which was outputted from the last" << endl;
	cout << "                         execution of the program, stores all the immediate" << endl;
	cout << "                         results (with file extension: '.chkpt.txt')" << endl;
	/*
	cout << "  -a                   : not outputting the checkpoint file" << endl;
	cout << "                         (this may increase the program speed, but the program" << endl;
	cout << "                          cannot be resumable)" << endl;
	*/
	cout << "  -i <# of iterations> : Number of iterations performed for each-time" << endl;
	cout << "                         parameter tuning" << endl;
	cout << "                         (default: -1 [no limit, run until converge])" << endl;
	cout << "  -m <min category #>  : Minimum number of site categories allowed" << endl;
	cout << "                         (default: 1)" << endl;
	cout << "  -n <max category #>  : Maximum number of site categories allowed" << endl;
	cout << "                         (default: 4)" << endl;
	cout << "  -w <step category #> : Step of change of site categories" << endl;
	cout << "                         (default: 1)" << endl;
	/*
	cout << "  -e <optim mode>      : 1 - 9, Optimization mode (default: " << DEFAULT_OPT_MODE << ")" << endl;
	cout << "                         (for more than one, format: mode1,mode2,...)" << endl;
	cout << "  -j <optim maxit>     : 1 - 100, Optimization maxIT (default: " << DEFAULT_OPT_MAXIT << ")" << endl;
	cout << "                         (for more than one, format: maxit1,maxit2,...)" << endl;
	*/
	cout << "  -precise             : More precise optimization, but needs more time." << endl;
	cout << "                         (default: disable)" << endl;

	cout << "  -b <info criteria>   : Information Criteria (IC) to be used (default: " << DEFAULT_IC << ")" << endl;
	cout << "                         1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC" << endl;
	if (!isMPI)
		cout << "  -h                   : This help page" << endl << endl;

	cout << "Output files:" << endl;
	cout << "  1. <output prefix>.HAS.chkpt.txt  - all the intermediate results (for resuming" << endl;
	cout << "                                      the program when necessary)" << endl;
	cout << "  2. <output prefix>.HAS.<IC>.result.txt - the resulting optimal number of site" << endl;
	cout << "                                      categories, the corresponding proportion" << endl;
	cout << "                                      and the details of parameters for each" << endl;
	cout << "                                      site category" << endl << endl;

	cout << "Example:" << endl;
	cout << "  Example alignment file: data.phy" << endl;
	cout << "  Example topology file : tree.txt" << endl;
	
	string hal_result_file = "data.NEW-BU.HAL.BIC.result.txt";
	string hal_program ="HAL-NEW";
	
#ifdef USE_OLD_PARAM
	hal_program = "HAL";
	hal_result_file = "data.BU.HAL.BIC.result.txt";
#endif
	
	cout << "  " << hal_program << " result file (which will appear after execution of " << hal_program << "):" << endl;
	cout << "      " << hal_result_file << endl;

	if (isMPI)
		cout << "  # of machines: 4" << endl;
	cout << "  To execute the " << progName << " program using the default parameters:" << endl;
	
	if (isMPI)
		cout << "      $ mpirun -n 4 " << progName << " data.phy tree.txt -t " << hal_result_file << endl << endl;
	else
		cout << "      $ " << progName << " data.phy tree.txt -t " << hal_result_file << endl << endl;

	cout << "Contact: " << CONTACTPERSON << endl;
	cout << "================================================================================" << endl;
}

// To show the usage of the Upsilon program
void UserOptions::outputUpsilonUsage(char* progName) {
	cout << "================================================================================" << endl;
	cout << "                  Welcome to HAS for Upsilon Model Version " << VERSION << endl << endl;

	cout << "Syntax: " << progName << " <alignment file> <topology file> <other options>" << endl;
	cout << "        " << progName << " -h" << endl << endl;
	cout << "  <alignment file>     : Multiple alignment file" << endl;
	cout << "  <topology file>      : Topology file in Newick format" << endl;
	cout << "other options:" << endl;
	cout << "  -f <format of file>  : The format of the multiple alignment file" << endl;
	cout << "                         1 - FASTA format" << endl;
	cout << "                         2 - Sequential PHYLIP format (default)" << endl;
	cout << "  -u <# of CPUs>       : Number of CPU threads used (default: 4)" << endl;
	cout << "  -o <output prefix>   : Prefix for output files" << endl;
	cout << "                         (default: <alignment file> w/o .ext)" << endl;
	cout << "  -r <checkpoint file> : Resume from the last execution;" << endl;
	cout << "                         <checkpoint file>, which was outputted from the last" << endl;
	cout << "                         execution of the program, stores all the immediate" << endl;
	cout << "                         results (with file extension: '.chkpt.txt')" << endl;
	cout << "  -a                   : not outputting the checkpoint file" << endl;
	cout << "                         (this may increase the program speed, but the program" << endl;
	cout << "                          cannot be resumable)" << endl;
	cout << "  -i <# of iterations> : Number of iterations performed for each-time" << endl;
	cout << "                         parameter tuning" << endl;
	cout << "                         (default: -1 [no limit, run until converge])" << endl;
	cout << "  -m <min category #>  : Minimum number of site categories allowed" << endl;
	cout << "                         (default: 2)" << endl;
	cout << "  -n <max category #>  : Maximum number of site categories allowed" << endl;
	cout << "                         (default: 4)" << endl;
	cout << "  -w <step category #> : Step of change of site categories" << endl;
	cout << "                         (default: 1)" << endl;
	cout << "  -b <info criteria>   : Information Criteria to be used (default: " << DEFAULT_IC << ")" << endl;
	cout << "                         1 - AIC; 2 - Adjusted IC; 3 - BIC; 4 - CAIC" << endl;
	cout << "  -h                   : This help page" << endl << endl;

	cout << "Output files:" << endl;
	cout << "  1. <output prefix>.Upsilon.chkpt.txt  - all the intermediate results (for resuming" << endl;
	cout << "                                      the program when necessary)" << endl;
	cout << "  2. <output prefix>.Upsilon.result.txt - the resulting optimal rate matrices" << endl << endl;
	cout << "Contact: " << CONTACTPERSON << endl;
	cout << "================================================================================" << endl;
}


// read the arguments and collect the user options
int UserOptions::readArguments(int argc, char** argv, string* errMsg) {

	// mode:
	//    1 - HAL-BU
	//    2 - HAL-TD
	//    3 - HAS
	//    4 - Upsilon

	// return:
	//    0 if the user options are valid
	//    1 if there exists error message
	//    2 if a help menu is called


	int i,j;
	bool duplicateOption = false;
	bool emptyOption = false;
	vector<string> token;
	int minArgNum;

	bool a_option_assigned = false;
	bool b_option_assigned = false;
	bool c_option_assigned = false;
	bool d_option_assigned = false;
	bool e_option_assigned = false;
	bool f_option_assigned = false;
	bool g_option_assigned = false;
	bool i_option_assigned = false;
	bool j_option_assigned = false;
	bool k_option_assigned = false;
	bool l_option_assigned = false;
	bool m_option_assigned = false;
	bool n_option_assigned = false;
	bool o_option_assigned = false;
	bool p_option_assigned = false;
	bool q_option_assigned = false;
	bool r_option_assigned = false;
	bool s_option_assigned = false;
	bool t_option_assigned = false;
	bool u_option_assigned = false;
	bool v_option_assigned = false;
	bool w_option_assigned = false;
	bool x_option_assigned = false;
	bool y_option_assigned = false;
    bool z_option_assigned = false;
	bool precise_option_assigned = false;
	int num_ic_selected = 0;

	*errMsg = "";

	switch(programType) {
	case 2:
		minArgNum = 4;
		break;
	default:
		minArgNum = 3;
		break;
	}


	if (argc < minArgNum) {
		// invoke the help menu
		return 2;
	}

	reset(); // set the options as default values

	alignFile = argv[1];
	topFile = argv[2];
	if (programType==2) {
		// for HAL-TD
		rateMatrixFile = argv[3];
	}
	// newly added, for those people who get used to the previous input format
	if (programType==3 && argc>=4 && argv[3][0]!='-') {
		// for HAS
		rateMatrixFile = argv[3];
	}


	// get the other options if there is any
	if (argc > minArgNum) {
		GetOptions options;
		options.read(argc, argv);

		for (i=0; i<options.size(); i++) {

			char flag = options.flags[i];
			string value = options.values[i];

			switch (flag) {

			case 'f':
			if (f_option_assigned)
				duplicateOption = true;
			else if (value.length() == 0)
				emptyOption = true;
			else {
				f_option_assigned = true;
				int valueInt = atoi(value.c_str());
				if (valueInt < 1 || valueInt > 2)
					*errMsg = "Error! The value for '-f' has to be 1 or 2";
				else
					inputFormat = valueInt;
			}
			break;

			case 'u':
				if (u_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					u_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt <= 0)
						*errMsg = "Error! The value for '-u' has to be greater than 0";
					else
						numCPUThreads = valueInt;
				}
				break;

			case 'k':
				if (k_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					k_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt <= 0)
						*errMsg = "Error! The value for '-k' has to be greater than 0";
					else
						topK = valueInt;
				}
				break;

			case 'o':
				if (o_option_assigned)
					duplicateOption = true;
				else {
					o_option_assigned = true;
					if (value.length()==0)
						emptyOption = true;
					else
						outputPrefix = value;
				}
				break;

			case 't':
				if (t_option_assigned)
					duplicateOption = true;
				else {
					t_option_assigned = true;
					if (value.length()==0)
						emptyOption = true;
					else
						rateMatrixFile = value;
				}
				break;

			case 'r':
				if (r_option_assigned)
					duplicateOption = true;
				else {
					r_option_assigned = true;
					if (value.length() == 0)
						emptyOption = true;
					else
						preLogFile = value;
				}
				break;

            case 'z':
                if (z_option_assigned)
                    duplicateOption = true;
                else {
                    z_option_assigned = true;
                    if (value.length() == 0)
                        emptyOption = true;
                    else
                        preLog = value;
                }
                break;

			case 'i':
				if (i_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					i_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt == 0)
						*errMsg = "Error! The value for '-i' has to be greater than 0 or -1 (i.e. unlimited)";
					else
						numIterations = valueInt;
				}
				break;

			case 'y':
				if (y_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					y_option_assigned = true;
					rateMatrixFile = value; // set the starting matrix
					numSteps = 1; // only do the first step for bottom-up
					continueTD = 0; // not proceeding to top-down
				}
				break;

			case 'c':
				if (c_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					c_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt == 0)
						*errMsg = "Error! The value for '-c' has to be greater than 0 or -1 (i.e. unlimited)";
					else
						numMaxCycles = valueInt;
				}
				break;

			case 's':
				if (s_option_assigned)
					duplicateOption = true;
				else {
					s_option_assigned = true;
					outputConsensus = 1;
				}
				break;

			case 'g':
				if (g_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					g_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 1 || valueInt > 2)
						*errMsg = "Error! The value for '-g' has to be 1 or 2";
					else
						gapHandling = valueInt;
				}
				break;

			case 'v':
				if (v_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					v_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 1 || valueInt > 2)
						*errMsg = "Error! The value for '-v' has to be 1 or 2";
					else
						groupAllCandidates = valueInt;
				}
				break;


			case 'p':
				if (value == "recise") {
					if (precise_option_assigned)
						duplicateOption = true;
					else {
						precise_option_assigned = true;
						precise = 1;
					}
				} else {
					if (p_option_assigned)
						duplicateOption = true;
					else if (value.length() == 0)
						emptyOption = true;
					else {
						p_option_assigned = true;
						int valueInt = atoi(value.c_str());
						if (valueInt != 0 && valueInt != 2)
							*errMsg = "Error! The value for '-p' has to be 0 or 2";
						else
							startOptState = valueInt;
					}
				}
				break;
					
			case 'm':
				if (programType!=3 && programType!=4)
					*errMsg = "Unknown option '-" + string(1,flag);
				else if (m_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					m_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 1)
						*errMsg = "Error! The value for '-m' has to be greater than 0";
					else
						minRateCat = valueInt;
				}
				break;

			case 'n':
				if (programType!=3 && programType!=4)
					*errMsg = "Unknown option '-" + string(1,flag);
				else if (n_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					n_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 1)
						*errMsg = "Error! The value for '-n' has to be greater than 0";
					else
						maxRateCat = valueInt;
				}
				break;

			case 'w':
				if (programType!=3 && programType!=4)
					*errMsg = "Unknown option '-" + string(1,flag);
				else if (w_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					w_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 1)
						*errMsg = "Error! The value for '-n' has to be greater than 0";
					else
						stepRateCat = valueInt;
				}
				break;

			case 'a':
			if (a_option_assigned)
				duplicateOption = true;
			else {
				a_option_assigned = true;
				outChkPtFile = 0;
			}
			break;

			case 'h':
				return 2;
				break;

			case 'e':
				opt_mode.clear();
				if (e_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					e_option_assigned = true;
					tokenizer(value,",",&token);
					for (j=0; j<(int)token.size(); j++) {
						opt_mode.push_back(atoi(token[j].c_str()));
					}
				}
				break;

			case 'j':
				opt_maxIT.clear();
				if (j_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					j_option_assigned = true;
					tokenizer(value,",",&token);
					for (j=0; j<(int)token.size(); j++) {
						opt_maxIT.push_back(atoi(token[j].c_str()));
					}
				}
				break;

			case 'd':
				if (d_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					d_option_assigned = true;
					logLDeviation = atof(value.c_str());
				}
				break;

			case 'b':
				if (b_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					b_option_assigned = true;
					tokenizer(value,",",&token);
					for (j=0; j<(int)token.size(); j++) {
						int valueInt = atoi(token[j].c_str());
						if (valueInt < 1 || valueInt > 4) {
							*errMsg = "Error! The value for '-b' has to be 1-4";
							break;
						} else {
							ic_list[valueInt-1] = 1;
							num_ic_selected++;
						}
					}
				}
				break;
					
			case 'l':
				if (l_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					l_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 0 || valueInt > 1)
						*errMsg = "Error! The value for '-l' has to be 0 or 1";
					else
						listParamBestHAL = valueInt;
				}
				break;

			case 'q':
				if (q_option_assigned)
					duplicateOption = true;
				else if (value.length() == 0)
					emptyOption = true;
				else {
					q_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 0 || valueInt > 1)
						*errMsg = "Error! The value for '-q' has to be 0 or 1";
					else
						continueTD = valueInt;
				}
				break;

			case 'x':
				if (x_option_assigned)
					duplicateOption = true;
				else {
					x_option_assigned = true;
					int valueInt = atoi(value.c_str());
					if (valueInt < 0 || valueInt > 1)
						*errMsg = "Error! The value for '-x' has to be 0 or 1";
					else
						unrooted = valueInt;
				}
				break;

			default:
				*errMsg = "Unknown option '-" + string(1,flag);
				break;

			} // case

			if (duplicateOption) {
				*errMsg = "Error! Duplicate option -" + string(1,flag);
			} else if (emptyOption) {
				*errMsg = "Error! Empty value for the option -" + string(1,flag);
			}

			if (*errMsg != "") {
				return 1;
			}

		}

		// check the values of minRateCat and maxRateCat
		if (maxRateCat < minRateCat)
			*errMsg = "Error! The value for '-m' cannot be greater than the value for '-n'";

		// check the values of opt_mode and opt_maxIT for optimization
		else if ((e_option_assigned && !j_option_assigned) || (!e_option_assigned && j_option_assigned))
			*errMsg = "Error! both options '-e' and '-j' should be used at the same time";

		else if (opt_mode.size() != opt_maxIT.size())
			*errMsg = "Error! The number of values for '-e' and for '-j' are not the same";

		if (*errMsg != "") {
			return 1;
		}

	}
	if (outputPrefix == "") {
		setDefaultOutputPrefix();
	}
	if (outputConsensus == 1 && (!d_option_assigned)) {
		logLDeviation = DEFAULT_LOG_L_THRES_FOR_CONSENSUS;
	}
	if (num_ic_selected == 0) {
		ic_list[DEFAULT_IC-1] = 1;
	}
	if (!e_option_assigned && !j_option_assigned && precise_option_assigned) {
		setPreciseModeMaxIT();
	}
	if (opt_mode.size()==0) {
		setDefaultModeMaxIT();
	}

	return 0;
}

// by default, prefixOut = <alignment file> w/o .ext
void UserOptions::setDefaultOutputPrefix() {

	outputPrefix = removeExtension(alignFile);

}

// To show the summary of the parameters
void UserOptions::showSummary(int isMPI) {

	int i;

	cout << "================================================================================" << endl;

	if (programType==1)
		cout << "                   Welcome to HAL ";
	else if (programType == 2)
		cout << "                   Welcome to HAL Top down ";
	else if (programType == 3)
		cout << "                        Welcome to HAS ";
	else // programType == 4
		cout << "              Welcome to HAS for Upsilon model ";
	
	if (isMPI)
		cout << "MPI ";

    // show regularization if so
    if (C_REG == 1)
        cout << "(regularization) ";

	cout << "Version " << VERSION << endl << endl;

	cout << "Input file ............................................. " << alignFile << endl;

	cout << "Format of input file ................................... ";
	switch (inputFormat) {
	case 1:
		cout << "FASTA" << endl;
		break;
	case 2:
		cout << "sequential PHYLIP" << endl;
		break;
	}

	cout << "Topology file .......................................... " << topFile << endl;

	if (rateMatrixFile != "") {
		if (programType==3) // HAS
			cout << "BU/TD result file ...................................... " << rateMatrixFile << endl;
		else if (programType==2) // HAL-TD
			cout << "BU result file ......................................... " << rateMatrixFile << endl;
		else if (programType==1) // HAL-BU
			cout << "User input matrix file ................................. " << rateMatrixFile << endl;
	}

	if (!isMPI) {
		cout << "Number of CPU threads .................................. " << numCPUThreads << endl;
	} else {
		cout << "Number of processors ................................... " << numCPUThreads << endl;
	}

	if (programType==1 || programType==2) // HAL-BU or HAL-TD
		cout << "Number of best results proceeded in next iteration ..... " << topK << endl;

	cout << "Prefix for output files ................................ " << outputPrefix << endl;

	if (preLogFile != "") {
		cout << "Resume from the checkpoint file ........................ " << preLogFile << endl;
	}

	cout << "Maximum # of iterations allowed for each-time tuning ... ";
	if (numIterations == -1) {
		cout << "unlimited" << endl;
	} else {
		cout << numIterations << endl;
	}

	cout << "C_REG value ............................................ " << C_REG << endl;

	if (programType==1) {
		if (HALMode==1)
			cout << "Algrithm for searching the optimal rate matrix.......... " << "original bottom-up" << endl;
		else if (HALMode==4)
			cout << "Algrithm for searching the optimal rate matrix.......... " << "new bottom-up approach" << endl;
		else if (HALMode==7)
			cout << "Algrithm for searching the optimal rate matrix.......... " << "John's bottom-up approach" << endl;
		else if (HALMode==8)
			cout << "Algrithm for searching the optimal rate matrix.......... " << "new2 bottom-up approach" << endl;
		else
			cout << "Algrithm for searching the optimal rate matrix.......... " << "bottom-up approach" << endl;
	} else if (programType==2)
		cout << "Algrithm for searching the optimal rate matrix.......... " << "top-down approach" << endl;

	if (programType==3 || programType==4) {
		// HAS or Upsilon model
		cout << "Minimum number of site categories allowed............... " << minRateCat << endl;
		cout << "Maximum number of site categories allowed............... " << maxRateCat << endl;
		cout << "Step of change of site categories ...................... " << stepRateCat << endl;
	}

	cout << "Group all candidates of same-level iterations........... ";
	if (groupAllCandidates==1) {
		cout << "YES" << endl;
	} else {
		cout << "NO" << endl;
	}

	cout << "Information criteria.................................... ";
	int firstIC = 1;
	for (i=0; i<4; i++) {
		if (ic_list[i]==1) {
			if (!firstIC)
				cout << ",";
			cout << getICName(i+1);
			firstIC = 0;
		}
	}
	cout << endl;

	// optimization mode (default = 1)
	cout << "optimization mode....................................... ";
	for (i=0; i<(int)opt_mode.size(); i++) {
		if (i>0)
			cout << ",";
		cout << opt_mode[i];
	}
	cout << endl;

	// optimization maxIT (default = 1)
	cout << "optimization maxIT...................................... ";
	for (i=0; i<(int)opt_maxIT.size(); i++) {
		if (i>0)
			cout << ",";
		cout << opt_maxIT[i];
	}
	cout << endl;

	// the starting state of the optimization method
	cout << "optimization start mode................................. ";
	if (startOptState==0)
		cout << "from initial state" << endl;
	else if (startOptState==1)
		cout << "from previous iteration" << endl;
	else if (startOptState==2)
		cout << "from the simplest model" << endl;
	else
		cout << "unknown" << endl;

/*
	if (programType == 1 || programType == 2) {
		
		if (logLDeviation > 0.1) {
			// range of deviations of resulting value of log-likelihood expected
			cout << "range of deviations of log-likelihood................... ";
			cout << logLDeviation << endl;
		}

		cout << "use concensus approach.................................. ";

		if (outputConsensus == 1)
			cout << "Yes";
		else
			cout << "No";
		cout << endl;
	}
	*/
	
	if (precise == 1) {
		cout << "Precise optimization process............................ Yes" << endl;
	}
	
	// running top-down algorithm after bottom-up algorithm
	if (programType == 1) {
		cout << "running top-down algorithm after bottom-up algorithm.... ";
		if (continueTD == 1)
			cout << "Yes";
		else
			cout << "No";
		cout << endl;
	}


	cout << "================================================================================" << endl;
}

