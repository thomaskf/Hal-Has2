/*
 *
 * RAL_new.cpp
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

#include "RAL_new.h"

// perform RAL algorithm
//
// options              : the user options
// rateGrpFile          : the rate group to start with
// preChkptFile         : The checkpoint file of the previous run.
// preLogFile           : The log file of the previous run.
//                        One wants to resume the process from the previous run.
void performRAL(UserOptions* options, char* rateGrpFile, char* preChkptFile, int num_chars, char* preLogFile) {

#ifdef GET_TIME_STAT
	pre_time = clock();
	pre_optim_time = pre_time;
#endif

	double thresHeight = 0.025;
	int isReversible=1; int num_w=6;
	// int numRateCat = 1; // number of rate category
	double IC;

	// load the topology
	int* topMatrix;
	int numLineTopMat, numSpecies, numEdges;
	vector<string>* leafList;
	vector<int> rateMatrix;
	vector<double> edgeLens;
	genTopMatrix((char*)options->topFile.c_str(), &topMatrix, &numLineTopMat, &leafList, &rateMatrix, &edgeLens);
	numSpecies=(int)leafList->size();
	if (numSpecies < 2) {
		cerr << "Error! The number of species is too small (i.e. less than 2) for this algorithm!" << endl;
		exit(1);
	}
	numEdges=numSpecies*2-2;
	map<string,int>* leafPosMap = genReverseMap(leafList);
	int* leafID2line = new int[numSpecies+1];
	int* internalID2line = new int[numLineTopMat+1];
	getMapping(leafID2line, internalID2line, topMatrix, numLineTopMat);

	// load the alignment file
	Alignment alignment;
	if (options->inputFormat==1)
		alignment.readFASTA((char*)options->alignFile.c_str(), leafPosMap);
	else
		alignment.readFile((char*)options->alignFile.c_str(), leafPosMap);
	if (options->gapHandling == 1) {
		// 1 - ignore all columns with gaps or ambiguous characters
		alignment.keepSitesWithChars(0); // remove the columns with invalid characters and gaps
	} else {
		alignment.keepSitesWithChars(1); // only remove the columns with invalid characters
	}
	if (alignment.numSites == 0) {
		cerr << "Error! The length of sequences are zero!" << endl;
		exit(1);
	}
	int* siteOrder = alignment.getSortedSitesID(); // get the order of the sites according to the topology
	// alignment.computeNumUniqueSites(siteOrder); // compute the number of unique columns
	alignment.keepUniqueSites(siteOrder); // keep only the unique columns
	// alignment.reorderUniqueSites(); // for efficient calculation of log-likelihood
#ifdef GET_TIME_STAT
	curr_time = clock();
	second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
	printf("Time used after reording the unique sites: %4.2f ms\n", second);
	pre_time = curr_time;
#endif
	alignment.buildChangedNodes(topMatrix, numLineTopMat);

	int num_chars_to_consider = 4;
	int isSingleGrp = 1;
	if (options->HALMode==2) {
		isSingleGrp = 0; // top-down approach
	}

#ifdef GTRIFO
	// GTR + I + FO model
	options->numSteps = 1;
	isSingleGrp = 1;
#endif

	// for starting state of the optimization
	int num_alpha = 1;
	ParameterSet start_ps(num_w, numEdges, num_chars);
	VariableSet start_vs(num_alpha, num_chars);
	start_ps.reset();
	start_vs.resetAllVariables(alignment);
	int initSeed = 1;
	if (options->randInitParam==1) {
		start_ps.randomInit(num_w, numEdges, initSeed, isSingleGrp);
	} else {
		start_ps.initialize(num_w, numEdges, isSingleGrp);
	}

	// initialize the starting rate matrix
	if (options->startRateMatrix!=NULL) {
		start_ps.loadRateMat(*(options->startRateMatrix));
	}
	
#ifndef GTRIFO	
	if (rateGrpFile != NULL) {
		start_ps.loadRateMat(rateGrpFile, topMatrix, leafList);
	}
#endif

	vector<int> startRateMatArray; // the starting rate matrix
	start_ps.outputRateMat(startRateMatArray);
	
	vector<int> untouchNodes(numEdges,0);

	// Load previous intermediate results if the checkpoint file exists
	TempRateMats* preInterResult = NULL;
	if (preChkptFile != NULL) {
		int num_alpha = 1; // for RAL
		options->info_criteria = firstNonZeroPos(options->ic_list, TOTAL_NUM_IC)+1;
		preInterResult = loadPreInterResult(preChkptFile, num_w, numEdges, num_alpha, options->info_criteria, 
				alignment.numSites, alignment.numSeqs, options->logLDeviation, num_chars, topMatrix);
	} else {
		preInterResult = new TempRateMats;
	}
    
    // Load previous log file (not checkpoint file) (added on Oct 12, 2020)
    if (preLogFile != NULL) {
        int num_alpha = 1; // for RAL
        options->info_criteria = firstNonZeroPos(options->ic_list, TOTAL_NUM_IC)+1;
        loadPreLogResult(preLogFile, num_w, numEdges, num_alpha, options->info_criteria,
                         alignment.numSites, alignment.numSeqs, options->logLDeviation, num_chars, topMatrix, preInterResult);
    }

	// term of the approach
	string approachTerm = "";
	switch(options->HALMode) {
	case 1:
		approachTerm = "BU";
		break;
	case 2:
		approachTerm = "TD";
		break;
	case 3:
		approachTerm = "BUTD";
		break;
	case 4:
		approachTerm = "NEW-BU";
		break;
	case 6:
		approachTerm = "BF";
		break;
	case 7:
		approachTerm = "JOHN";
		break;
	case 8:
		approachTerm = "NEW2";
		break;
	default:
		approachTerm = "OTHER";
		break;
	}

	// checkpoint file
	ofstream* chkptOut = NULL;

#ifndef GTRIFO	
	string outChkptFile;
	if (options->outChkPtFile && options->outputPrefix.length() > 0) {
		outChkptFile = options->outputPrefix + "." + approachTerm + ".chkpt.txt";
		chkptOut = new ofstream;
		if (preChkptFile != NULL) {
			if (outChkptFile != string(preChkptFile)) {
				filecopy((char*) outChkptFile.c_str(), preChkptFile);
			}
			chkptOut->open(outChkptFile.c_str(), std::ofstream::out | std::ofstream::app); // append to the log file
		} else {
			chkptOut->open(outChkptFile.c_str());
		}
	}
#endif

	// The optimal rate matrix for each number of rates allowed
	OptimalRateMats optRateMatrix;

	// The rate matrice waiting for process
	TempRateMats rateMatrixList;

	// The first one hundred best rate matrices
	TempRateMats bestHundredRateMatrices;

	// ============================================================
	// setup and initialize the locks for multiple CPU threads
	// ============================================================

	ThreadLocks* threadLocks = NULL;
	if (options->numCPUThreads > 1) {
		threadLocks = new ThreadLocks();
		initThreadLocks(threadLocks);
	}
	
	
	// =================
	// iteration start
	// =================
	// perform RAL iteration for each user defined information criteria
	int i,j;
	string ICName;
	int firstIC = 1;
	ofstream* resultOut = NULL;
	string resultFileName = "";
	for (i=0; i<TOTAL_NUM_IC; i++) {

		if (options->ic_list[i] == 0)
			continue;

		options->info_criteria = i+1;
		ICName = getICName(options->info_criteria);

		// initialize untouchNodes
		for (j=0; j<numEdges; j++) {
			untouchNodes[j] = 0;
		}
		if (options->startUntouchEdges!=NULL) {
			untouchNodes.assign(options->startUntouchEdges->begin(), options->startUntouchEdges->end());
		}

		// initialize the matrix list for processing
		IC = 1e60;
		rateMatrixList.clear();
		rateMatrixList.insert(startRateMatArray, untouchNodes, IC, -1, 1);

		// initialize the optial matrix for each number of rates allowed
		optRateMatrix.clear();

		// initialize the matrix which keeps the first one hundred best rate matrices
		bestHundredRateMatrices.clear();

		// recalculate the IC values for the array preInterResult
		if (firstIC==0) {
			reCalculateAllICs(preInterResult, options->info_criteria, alignment.numSites, numSpecies, options->logLDeviation, topMatrix);
		}

		// start RAL iterations
		int saveToInterResult = 1;
		
#ifdef GET_TIME_STAT
		curr_time = clock();
		second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
		printf("Time used before starting the iteration: %4.2f ms\n", second);
		pre_time = curr_time;
#endif
		
		RALIterationNew (&rateMatrixList, topMatrix, numLineTopMat,
				&alignment, isReversible, numSpecies,
				threadLocks, &optRateMatrix, leafID2line, internalID2line, thresHeight,
				chkptOut, preInterResult, &bestHundredRateMatrices, options, num_chars_to_consider,
				saveToInterResult, start_ps, start_vs, leafList);

		// result file
		if (options->outputPrefix.length() > 0) {
			resultOut = new ofstream;
#ifndef GTRIFO	
			resultFileName = options->outputPrefix + "." + approachTerm + "." + ICName + ".result.txt";
#else
			resultFileName = options->outputPrefix + ".GTRI." + ICName + ".result.txt";
#endif			
			resultOut->open(resultFileName.c_str());
		}

		// print the set of optimal rate matrices
		// optRateMatrix.print_optimal_rateMatrix(resultOut, topMatrix, leafList);
		// if (options->outputConsensus)
		//     bestHundredRateMatrices.print_consensus(resultOut, topMatrix, leafList, ICName, options->logLDeviation*2.0);
		// else
		// bestHundredRateMatrices.print_optimal_rateMatrices(resultOut, topMatrix, leafList, ICName);
		bestHundredRateMatrices.print_optimal_rateMatrices(resultOut, topMatrix, leafList, alignment.numSites, options->info_criteria);
		// bestHundredRateMatrices.print_content();
		// print out the content of the optimal rate matrices
		optRateMatrix.print_content(resultOut, topMatrix, leafList, ICName);
		// print out the content of the best hundred rate matrices
		string description = "The best 100 arrangements of the rate matrices:";
		bestHundredRateMatrices.print_content(resultOut, topMatrix, leafList, description, ICName);
		if (options->listParamBestHAL && bestHundredRateMatrices.size() > 0 && bestHundredRateMatrices.PSs[0]!=NULL) {
			// list out the parameters of the best HAL model
			(*resultOut) << "=====================================================================================" << endl;
			(*resultOut) << "The corresponding parameters of the best HAL model:" << endl;
            string currResult = "";
            // (bestHundredRateMatrices.PSs[0])->updateContent(2);
            (bestHundredRateMatrices.PSs[0])->showContent(currResult, topMatrix, leafList);
            (bestHundredRateMatrices.VSs[0])->showContent(currResult, 1);
            (*resultOut) << currResult << endl;
			// show the tree
			(bestHundredRateMatrices.PSs[0])->updateContent(1);
			(*resultOut) << "The tree (whose edge lengths represent average # of substitutions per site):" << endl;
			outTreeWithEdgeLen(topMatrix, leafList, (bestHundredRateMatrices.PSs[0])->t, *resultOut);
			(*resultOut) << endl;
		}

		if (resultOut != NULL) {
			resultOut->close();
			resultOut = NULL;
		}
		
		
		// start the top-down algorithm if necessary
#ifndef GTRIFO	
		if ((options->HALMode==1 || options->HALMode==4 || options->HALMode==7 || options->HALMode==8) && (options->continueTD==1)) {
			options->HALMode=2;

			// initialize untouchNodes
			for (j=0; j<numEdges; j++) {
				untouchNodes[j] = 0;
			}
			// initialize the matrix list for processing
			IC = 1e60;
			rateMatrixList.clear();
			// select the top one
			for (j=0; j<1 && j<bestHundredRateMatrices.size(); j++) {
				rateMatrixList.insert(*(bestHundredRateMatrices.matrices[j]), untouchNodes, IC, -1, 1);
			}
			// start
			cout << "Starting the top-down algorithm..." << endl;
			
			RALIterationNew (&rateMatrixList, topMatrix, numLineTopMat,
							 &alignment, isReversible, numSpecies,
							 threadLocks, &optRateMatrix, leafID2line, internalID2line, thresHeight,
							 chkptOut, preInterResult, &bestHundredRateMatrices, options, num_chars_to_consider,
							 saveToInterResult, start_ps, start_vs, leafList);
			
			// result file
			approachTerm += ".HAL";
			if (options->outputPrefix.length() > 0) {
				resultOut = new ofstream;
				resultFileName = options->outputPrefix + "." + approachTerm + "." + ICName + ".result.txt";
				resultOut->open(resultFileName.c_str());
			}
			
			// print the set of optimal rate matrices
			// optRateMatrix.print_optimal_rateMatrix(resultOut, topMatrix, leafList);
			// if (options->outputConsensus)
			//     bestHundredRateMatrices.print_consensus(resultOut, topMatrix, leafList, ICName, options->logLDeviation*2.0);
			// else
			// bestHundredRateMatrices.print_optimal_rateMatrices(resultOut, topMatrix, leafList, ICName);
			bestHundredRateMatrices.print_optimal_rateMatrices(resultOut, topMatrix, leafList, alignment.numSites, options->info_criteria);
			// bestHundredRateMatrices.print_content();

			// print out the content of the optimal rate matrices
			optRateMatrix.print_content(resultOut, topMatrix, leafList, ICName);
			
            // print out the content of the best hundred rate matrices
			string description = "The best 100 arrangements of the rate matrices:";
			bestHundredRateMatrices.print_content(resultOut, topMatrix, leafList, description, ICName);

			if (options->listParamBestHAL && bestHundredRateMatrices.size() > 0 && bestHundredRateMatrices.PSs[0] != NULL) {
				// list out the parameters of the best HAL model
				(*resultOut) << "=====================================================================================" << endl;
				(*resultOut) << "The corresponding parameters of the best HAL model:" << endl;
				string currResult = "";
				// (bestHundredRateMatrices.PSs[0])->updateContent(2);
				(bestHundredRateMatrices.PSs[0])->showContent(currResult, topMatrix, leafList);
				(bestHundredRateMatrices.VSs[0])->showContent(currResult, 1);
				(*resultOut) << currResult << endl;
				// show the tree
				(bestHundredRateMatrices.PSs[0])->updateContent(1);
				(*resultOut) << "The tree (whose edge lengths represent average # of substitutions per site):" << endl;
				outTreeWithEdgeLen(topMatrix, leafList, (bestHundredRateMatrices.PSs[0])->t, *resultOut);
				(*resultOut) << endl;
			}
			
			if (resultOut != NULL) {
				resultOut->close();
				resultOut = NULL;
			}

		}
#endif		
		
		firstIC = 0;

	}

	if (chkptOut != NULL)
		chkptOut->close();
		
	// clear all the memories
	if (threadLocks != NULL)
		delete threadLocks;

	// for consensus approach
	if (options->outputConsensus == 1) {
		resultOut = new ofstream;
		resultFileName = options->outputPrefix + "." + approachTerm + ".consensus.result.txt";
		resultOut->open(resultFileName.c_str());
		consensusMethod(resultOut, preInterResult, options, options->logLDeviation*2,
				numSpecies, numEdges, alignment.numSites, topMatrix, leafList);
		resultOut->close();
	}

	cout << endl;
	cout << "================================================================================" << endl;

	// list out the names of the resulting files
	for (i=0; i<TOTAL_NUM_IC; i++) {
		if (options->ic_list[i] == 0)
			continue;
		ICName = getICName(i+1);

#ifndef GTRIFO	
			resultFileName = options->outputPrefix + "." + approachTerm + "." + ICName + ".result.txt";
#else
			resultFileName = options->outputPrefix + ".GTRI." + ICName + ".result.txt";
#endif			

		if (resultFileName.length() > 0) {
			cout << "Result file: " << resultFileName << endl;
			// cout << "Result file according to "<< ICName << ": " << resultFileName << endl;
		}
	}
	if (options->outputConsensus == 1) {
		resultFileName = options->outputPrefix + "." + approachTerm + ".consensus.result.txt";
		if (resultFileName.length() > 0)
			cout << "Result file of consensus approach: " << resultFileName << endl;
	}

	cout << "HAL Version " << VERSION << " finishes" << endl;


}


// ============================================================================================================================
// ITERATION FUNCTIONS for RAL
// ============================================================================================================================

// perform RAL iteration for the corresponding rate matrices
// groupAllCandidates  : For each of the best rate matrices in this iteration,
//                       generate the potential candidates for next iteration
//                       1 - All the potential candidate from the same-level iterations are grouped for the next-level iterations.
//                       2 - All the potential candidate from the same-level iterations are NOT grouped for the next-level iterations.
//                           (i.e. Original design)
void RALIterationNew (TempRateMats* rateMatrixList, int* topMatrix, int numLineTopMat,
		Alignment* alignment, int isReversible, int numSpecies,
		ThreadLocks* threadLocks, OptimalRateMats* optRateMatrix, int* leafID2line, int* internalID2line,
		double thresHeight, ofstream* fout, TempRateMats* preInterResult,
		TempRateMats* bestHundredRateMatrices, UserOptions* options, int num_chars,
		int saveToInterResult, ParameterSet& start_ps, VariableSet& start_vs, vector<string>* leafList){

	int i,j,k;

	int numIterations = options->numIterations;
	int availCPUs = options->numCPUThreads;
	int topK = options->topK;
	int whichApproach = options->HALMode;
	int groupAllCandidates = options->groupAllCandidates;
	int maxNumCycles = options->numMaxCycles;

	vector<int> theRateMatArray;
	vector< int> currUntouchNodes;
	int threadId;
	int numEdges = numSpecies*2-2;
	vector<int> untouchzero(numEdges,0);

	// initialize a set of variables for threads
	vector<pthread_t> thread_set;
	vector<OptimThreadArguments> args_threads;
	thread_set.resize(availCPUs);
	args_threads.resize(availCPUs);

	// initialize the job allocation array
	int* jobAllocated = NULL;
	int jobAllocatedSize = 0;

	// indicate whether the job is needed to output to the check point file
	vector<int> jobNeedToOutChkptFile;

	// how many nodes can be considered
	int untouchedNodes;

	// grouping array
	vector<int> grouping_array;
	grouping_array.push_back(rateMatrixList->size());

	// job list for optimization
	TempRateMats jobList;

	// an array to check the consistency of the resulting arrangment of rate matrices
	vector<int>* consistentRates = NULL;
	if (whichApproach==4) {
		consistentRates = new vector<int> (numEdges, 0);
	}
	
	// the order of the edges to be examined
	int* edgeOrder = NULL;

	int masterThresID = 0;
	bool needToProcess;
	int numMatToCheck = 0;
	int startItem;
	int insertToSet;
	int numCycles = 1;
	int numStep = 0;
	double preBestIC = 1e77;
	double currBestIC;
	double IC_Thres;

	if (groupAllCandidates==1 && whichApproach!=1 && whichApproach!=7) {
		IC_Thres = IC_Thres_groupAll1;
	} else {
		IC_Thres = IC_Thres_groupAll2;
	}

	string ICName = getICName(options->info_criteria);

	int opt_method_num = (int) options->opt_mode.size();

	while (grouping_array.size() > 0 && (options->numSteps==-1 || numStep++ < options->numSteps)) {

		// the number of rate matrices are going to check
		numMatToCheck = grouping_array[grouping_array.size()-1];
		grouping_array.pop_back();
		cout << "[master " << ++masterThresID << "] Checking " << numMatToCheck << " rate matices...." << endl << flush;

		// perform optimization on the last 'numMatToCheck' rate matrix
		jobList.clear();
		jobNeedToOutChkptFile.clear();
		startItem = rateMatrixList->size()-numMatToCheck;
		insertToSet = 0;

		for (i=startItem; i<rateMatrixList->size(); i++) {

			// check whether it needs to process
			needToProcess = true;
			if (rateMatrixList->PSs[i]!=NULL) {
				needToProcess = false;
				jobNeedToOutChkptFile.push_back(0); // no need to output the result to checkpoint file
			} else {
				int interID = preInterResult->getID(*(rateMatrixList->matrices[i]));
				if (interID > -1 && (options->HALMode!=2 || preInterResult->PSs[interID] != NULL)) {
					needToProcess = false;
                    if (preInterResult->PSs[interID] != NULL) {
                        rateMatrixList->PSs[i] = new ParameterSet(num_chars);
                        rateMatrixList->PSs[i]->copyFrom(*(preInterResult->PSs[interID]));
                        rateMatrixList->VSs[i] = new VariableSet(num_chars);
                        rateMatrixList->VSs[i]->copyFrom(*(preInterResult->VSs[interID]));
                    } else {
                        rateMatrixList->PSs[i] = NULL;
                        rateMatrixList->VSs[i] = NULL;
                    }
					rateMatrixList->loglikelihood[i] = preInterResult->loglikelihood[interID];
					rateMatrixList->ICs[i] = preInterResult->ICs[interID];
					rateMatrixList->numIters[i] = preInterResult->numIters[interID];
					rateMatrixList->lowPossibleICs[i] = preInterResult->lowPossibleICs[interID];
					jobNeedToOutChkptFile.push_back(0);
				} else {
					jobNeedToOutChkptFile.push_back(1); // need to output the result to checkpoint file
				}
			}

			if (needToProcess) {
				// duplicate the jobs for different optimization methods
				for (j=0; j<opt_method_num; j++) {
					if (options->startOptState==2)
						jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->parentICs[i], -1, insertToSet,start_ps, start_vs);
					else if (options->startOptState==0 || rateMatrixList->parentPSs[i]==NULL || rateMatrixList->parentVSs[i]==NULL)
						jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->parentICs[i], -1, insertToSet);
					else
						jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->parentICs[i], -1, insertToSet, *(rateMatrixList->parentPSs[i]), *(rateMatrixList->parentVSs[i]));
				}
			} else {
				// duplicate the empty jobs for different optimization methods
				for (j=0; j<opt_method_num; j++)
					jobList.insert();
			}
		}

		// initialize the "jobAllocated"
		if (jobAllocatedSize < jobList.size()) {
			// increase the size of the "jobAllocated"
			if (jobAllocated != NULL)
				delete[] jobAllocated;
			jobAllocatedSize = jobList.size();
			jobAllocated = new int[jobAllocatedSize];
		}

		// indicate which job is needed to process
		for (i=0; i<jobList.size(); i++) {
			if (jobList.matrices[i]==NULL)
				jobAllocated[i]=1; // processed
			else
				jobAllocated[i]=0; // waiting to process
		}
		
		if (availCPUs > 1) {
			for (i=1; i<availCPUs; i++) {
				thread_set[i] = 0;
				args_threads[i].num_chars = num_chars;
				args_threads[i].rateMatChkLst = &jobList;
				args_threads[i].topMatrix = topMatrix;
				args_threads[i].numLineTopMat = numLineTopMat;
				args_threads[i].alignment = alignment;
				args_threads[i].isReversible = isReversible;
				args_threads[i].numSpecies = numSpecies;
				args_threads[i].numIterations = numIterations;
				args_threads[i].masterThresID = masterThresID;
				args_threads[i].threadID = i;
				args_threads[i].availCPUs = availCPUs;
				args_threads[i].threadLocks = threadLocks;
				args_threads[i].jobAllocate = jobAllocated;
				args_threads[i].preInterResult = preInterResult;
				args_threads[i].fout = fout;
				args_threads[i].userOptions = options;
			}

			// optimizing all parameters and get the optimal likelihood and IC value
			// The optimization process run by each CPU thread (in parallel)
			for (threadId=1; threadId<availCPUs; threadId++) {
				if ( pthread_create ( & ( thread_set[threadId] ), NULL, optimThreadWrapper, ( void * ) & ( args_threads[threadId] ) ) ) {
					cerr << "[performRAL: ThreadId " << masterThresID << "." << threadId << "] Cannot create ""optimThreadWrapper""\n" << endl << flush;
					pthread_mutex_unlock(&(threadLocks->lock_createThread));
					exit(1);
				}
			}
		}
		
#ifdef GET_TIME_STAT
		curr_time = clock();
		second = (curr_time - pre_time)* 1000.0 / CLOCKS_PER_SEC;
		printf("Before calling optimThread, time used: %4.2f ms\n", second);
		pre_time = curr_time;
#endif
#ifdef GET_TIME_STAT_DETAIL
		pre_optim_time = curr_time;
		compute_eigen = 0.0;
		compute_prob = 0.0;
		compute_L = 0.0;
#endif

		// for this processor
		optimThread(&jobList, topMatrix, numLineTopMat, alignment, isReversible, numSpecies, numIterations, 0, masterThresID, availCPUs, threadLocks, jobAllocated, preInterResult, options, num_chars);

#ifdef GET_TIME_STAT
		curr_time = clock();
		second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
		printf("After calling optimThread, total time used: %4.2f ms\n", second);
		pre_time = curr_time;
#endif
#ifdef GET_TIME_STAT_DETAIL
		printf("(eigen: %4.2f ms;", compute_eigen);
		printf(" cond prob: %4.2f ms;", compute_prob);
		printf(" log-L: %4.2f ms)\n", compute_L);
#endif

		// wait for all the processes to finish
		if (availCPUs > 1) {
			for (threadId=1; threadId<availCPUs; threadId++) {
				if ( thread_set[threadId] != 0 ) {
					if ( pthread_join ( thread_set[threadId], NULL ) ) {
						pthread_mutex_lock(&(threadLocks->lock_output));
						cerr << "[performRAL: ThreadId " << masterThresID << "." << threadId << "] Crash!\n" << endl << flush;
						pthread_mutex_unlock(&(threadLocks->lock_output));
						exit(1);
					}
					// pthread_detach(thread_set[threadId]);
					thread_set[threadId] = 0;
				}
			}
		}

		// consolidate all the answers
		for (i=0; i<jobList.size(); i+=opt_method_num) {
			k = i/opt_method_num;
			if (jobList.matrices[i]!=NULL) {
				int opt_i = i;
				double opt_IC = jobList.ICs[i];
				for (j=1; j<opt_method_num; j++) {
					if (jobList.ICs[i+j] < opt_IC) {
						opt_IC = jobList.ICs[i+j];
						opt_i = i+j;
					}
				}

				// all numbers are rounded up, so that
				// it is consistent with the intermediate results loaded from chkpt file

                if (jobList.PSs[opt_i] != NULL) {
                    rateMatrixList->PSs[startItem + k] = new ParameterSet(num_chars);;
                    rateMatrixList->PSs[startItem + k]->copyFrom(*(jobList.PSs[opt_i]));
                    rateMatrixList->PSs[startItem + k]->roundContent();

                    rateMatrixList->VSs[startItem + k] = new VariableSet(num_chars);;
                    rateMatrixList->VSs[startItem + k]->copyFrom(*(jobList.VSs[opt_i]));
                    rateMatrixList->VSs[startItem + k]->roundContent();
                } else {
                    rateMatrixList->PSs[startItem + k] = NULL;
                    rateMatrixList->VSs[startItem + k] = NULL;
                }

				rateMatrixList->loglikelihood[startItem + k] = roundNumber(jobList.loglikelihood[opt_i], ANS_DECI);
				rateMatrixList->ICs[startItem + k] = roundNumber(opt_IC, ANS_DECI);
				rateMatrixList->numIters[startItem + k] = jobList.numIters[opt_i];
				rateMatrixList->lowPossibleICs[startItem + k] = roundNumber(jobList.lowPossibleICs[opt_i], ANS_DECI);
			}
			
			// for startOptState = 2 (i.e. optimization starts from the state using one model)
			// if this is the state where only one rate group is involved,
			// then set the starting parameters and starting variables
			if (options->startOptState==2 && rateMatrixList->numRates[startItem + k]==1 && rateMatrixList->PSs[startItem + k]!=NULL) {
				start_ps.copyFrom(*rateMatrixList->PSs[startItem + k]);
				start_vs.copyFrom(*rateMatrixList->VSs[startItem + k]);
			}
			
			// if this is the state where only one rate group is involved
			// then compute the order of the edges
			if (rateMatrixList->numRates[startItem + k]==1 && whichApproach==4 && rateMatrixList->PSs[startItem + k] != NULL) {
				edgeOrder = getEdgeOrder(rateMatrixList->PSs[startItem + k], options->unrooted);
			}
		}
		// output the intermediate results
		for (i=startItem; i<rateMatrixList->size(); i++) {
			if (jobNeedToOutChkptFile[i-startItem]==1) {

				// save the intermediate result to the array if necessary
				if (saveToInterResult) {
                    if (rateMatrixList->PSs[i] == NULL)
                        preInterResult->insertwMap(*(rateMatrixList->matrices[i]), rateMatrixList->loglikelihood[i], rateMatrixList->ICs[i],
                                -1, untouchzero, rateMatrixList->numIters[i], rateMatrixList->lowPossibleICs[i]);
                    else
                        preInterResult->insertwMap(*(rateMatrixList->matrices[i]), rateMatrixList->loglikelihood[i], rateMatrixList->ICs[i],
							    -1, *(rateMatrixList->PSs[i]), *(rateMatrixList->VSs[i]), untouchzero,
							    rateMatrixList->numIters[i], rateMatrixList->lowPossibleICs[i]);
				}

				// need to output the result
				// update the parameters such that the edge length equals to the average # of substitutions 
				// rateMatrixList->PSs[i]->updateContent(1);
				if (fout != NULL) {
					// output to the checkout file
					(*fout) << "===============================================================================================" << endl;
					(*fout) << "[# of iterations:" << rateMatrixList->numIters[i] << "] loglikelihood:" << longDoublToStr(rateMatrixList->loglikelihood[i],ANS_DECI);
					(*fout) << " " << ICName << ":" << longDoublToStr(rateMatrixList->ICs[i],ANS_DECI);
					// (*fout) << " lowest_possible_IC:" << longDoublToStr(rateMatrixList->lowPossibleICs[i],ANS_DECI)
					(*fout) << endl;
                    if (rateMatrixList->PSs[i] != NULL) {
                        rateMatrixList->PSs[i]->showContent(fout);
                        rateMatrixList->VSs[i]->showContent(fout);
                    }
					for (int k=0; k<(int)rateMatrixList->matrices[i]->size(); k++) {
						if (k>0)
							(*fout) << ",";
						(*fout) << rateMatrixList->matrices[i]->at(k);
					}
					(*fout) << " [numRateGrps: " << rateMatrixList->numRates[i] << "] " << endl << flush;
				}
			}
		}
		// update the optimal rate matrix
		for (i=startItem; i<rateMatrixList->size(); i++) {
			optRateMatrix->insertOpt(rateMatrixList->numRates[i], *(rateMatrixList->matrices[i]), rateMatrixList->loglikelihood[i],
					rateMatrixList->ICs[i]);
		}

		// rearrange the last 'numMatToCheck' items inside the matrices
		// such that the top K items (among the last 'numMatToCheck' items) are with the least IC values
		rateMatrixList->rearrangeLastNItems(topK, numMatToCheck);

		// print the last 'numMatToCheck' items inside rateMatrixList
		rateMatrixList->print_last_N(numMatToCheck, ICName);

		// save the top K items (among the last 'numMatToCheck' items) into the joblist
		jobList.clear();
		for (i=startItem; i<rateMatrixList->size(); i++) {
			if (groupAllCandidates==1) {
                if (i < startItem+topK) {
                    if (rateMatrixList->PSs[i] != NULL)
                        jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], *(rateMatrixList->PSs[i]), *(rateMatrixList->VSs[i]), rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
                    else
                        jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
                }
			} else {
                if (i < startItem+topK && rateMatrixList->ICs[i]<=rateMatrixList->parentICs[i]) {
                    if (rateMatrixList->PSs[i] == NULL)
                        jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], *(rateMatrixList->PSs[i]), *(rateMatrixList->VSs[i]), rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
                    else
                        jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
                }
			}
		}
        
		// move the last 'numMatToCheck' items inside rateMatrixList to bestHundredRateMatrices
		rateMatrixList->move_last_k_items(bestHundredRateMatrices, numMatToCheck, 1);

		// for new bottom-up algorithm
		// check how many nodes can be considered
		// (same for all, thus only need to check the first one)
		untouchedNodes = 0;
		if (whichApproach>=4 && jobList.size()> 0) {
			for (j=0; j<numEdges; j++) {
				if (jobList.untouches[0]->at(j)==0)
					untouchedNodes++;
			}
			// for groupAllCandidates = 1, it allows for more than one cycles
			if (untouchedNodes==0 && groupAllCandidates==1) {

				// show the best IC score for this cycle
				bestHundredRateMatrices->rearrangeFirstNItems(1, bestHundredRateMatrices->size());
				currBestIC = bestHundredRateMatrices->ICs[0];
				cout << "Best result for cycle " << numCycles << " is: ";
				for (k=0; k<(int)bestHundredRateMatrices->matrices[0]->size(); k++) {
					if (k>0)
						cout << ",";
					cout << bestHundredRateMatrices->matrices[0]->at(k);
				}
				cout << " [" << ICName << ": " << longDoublToStr(currBestIC,ANS_DECI) << "]" << endl;

				/*
				if (whichApproach==4) {
					// rearrange all items inside the matrices
					// such that the CONSISTENCY_CHECK_NUM items (among all items) are with the least IC values
					findConsistentRates(bestHundredRateMatrices, consistentRates, CONSISTENCY_CHECK_NUM);

					cout << "Consistent rate matrice assignment: ";
					for (k=0; k<(int)consistentRates->size(); k++) {
						if (k>0)
							cout << ",";
						cout << consistentRates->at(k);
					}
					cout << endl;

				}
				 */

				if ((maxNumCycles == -1 || numCycles < maxNumCycles) && currBestIC < preBestIC) {
					preBestIC = currBestIC;
					numCycles++;


					/*
					// increase the value of topK
					topK = origTopK + 2 * (numEdges - untouchedNodes);

					// replace all the job lists to the topK items
					bestHundredRateMatrices->rearrangeFirstNItems(topK, bestHundredRateMatrices->size());
					jobList.clear();
					for (i=0; i<topK && i<bestHundredRateMatrices->size(); i++) {
						jobList.insert(*(bestHundredRateMatrices->matrices[i]), *(bestHundredRateMatrices->untouches[i]),
						bestHundredRateMatrices->loglikelihood[i], bestHundredRateMatrices->parentICs[i], *(bestHundredRateMatrices->PSs[i]), 
					 *(bestHundredRateMatrices->VSs[i]), bestHundredRateMatrices->ICs[i], -1, insertToSet);
					}
					 */

					// reset all the untouchedNodes for those edges without the consistent rate matrix
					untouchedNodes = jobList.reset_untouches(consistentRates);
				}
			}
		}

		int new_candidates = 0;
		if (options->numSteps==-1 || numStep++ < options->numSteps) {
			for (i=0; i<jobList.size(); i++) {
				if (jobList.ICs[i] <= jobList.parentICs[i] * IC_Thres) {
					// the IC is small enough
					int pre_new_candidates = new_candidates;
					if ((whichApproach==1 || whichApproach==3) && jobList.numRates[i] < numEdges) {
						// generate all the potential rate matrices into the rateMatrixList (using bottom-up approach)
						new_candidates += generateRateMat(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], options->unrooted);
					}
					if ((whichApproach==2 || whichApproach==3) && jobList.numRates[i] > 1) {
						// generate all the potential rate matrices into the rateMatrixList (using top-down approach)
						new_candidates += generateRateMatTD(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, internalID2line, jobList.untouches[i], jobList.numRates[i], thresHeight, options->unrooted);
					}
					if (whichApproach==4 && untouchedNodes>0) {
						// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 1 approach)
						new_candidates += generateRateMatNewBU1(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], jobList.lowPossibleICs[i], edgeOrder, options->unrooted);
					}
					if (whichApproach==6 && untouchedNodes>0) {
						// generate all the potential rate matrices into the rateMatrixList (using bute-force approach)
						new_candidates += generateRateMatButeForce(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, options->maxNumRatMat, jobList.lowPossibleICs[i], options->unrooted);
					}
					if ((whichApproach==7) && jobList.numRates[i] < numEdges) {
						// generate all the potential rate matrices into the rateMatrixList (using John's bottom-up approach)
						new_candidates += generateRateMatJohn(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], options->unrooted);
					}
					if (whichApproach==8 && untouchedNodes>0) {
						// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 2 approach)
						new_candidates += generateRateMatNewBU2(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], jobList.lowPossibleICs[i], options->unrooted);
					}
					if (new_candidates > pre_new_candidates) {
						// print out the matrice
						cout << "[";
						for (int n=0; n<(int)jobList.matrices[i]->size(); n++) {
							if (n > 0)
								cout << ",";
							cout << jobList.matrices[i]->at(n);
						}
						cout << "] produces " << new_candidates-pre_new_candidates << " new candidates" << endl;
					}
					if (groupAllCandidates==2) {
						// The potential candidates form separate groups according to the rate matrix generated from. (i.e. Original design)
						if (new_candidates > 0) {
							grouping_array.push_back(new_candidates);
							new_candidates = 0;
						}
					}
				}
			}
		}

		if (groupAllCandidates!=2) {
			// All the potential candidate form the same group to be processed in next iteration.
			if (new_candidates > 0) {
				grouping_array.push_back(new_candidates);
				new_candidates = 0;
			}
		}
	}

	// rearrange all items inside the matrices
	// such that the top 100 items (among all items) are with the least IC values
	bestHundredRateMatrices->rearrangeFirstNItems(100, bestHundredRateMatrices->size());
	// remove the remaining item
	int numRemainItem = bestHundredRateMatrices->size() - 100;
	if (numRemainItem > 0)
		bestHundredRateMatrices->remove_last_k_items(numRemainItem);

	// clear the job allocation array
	delete[] jobAllocated;
	if (edgeOrder != NULL)
		delete[] edgeOrder;
}


// the Wrapper to perform "optimThread" in each CPU thread
void * optimThreadWrapper ( void * arg ) {

	OptimThreadArguments* currArg = (OptimThreadArguments*) arg;

	optimThread(currArg->rateMatChkLst, currArg->topMatrix, currArg->numLineTopMat, currArg->alignment, currArg->isReversible,
			currArg->numSpecies, currArg->
			numIterations, currArg->threadID, currArg->masterThresID, currArg->availCPUs,
			currArg->threadLocks, currArg->jobAllocate, currArg->preInterResult, currArg->userOptions, currArg->num_chars);

	pthread_exit ( 0 );
	return 0;

}

// The optimization process run by each CPU thread
void optimThread(TempRateMats* rateMatChkLst, int* topMatrix, int numLineTopMat, Alignment* alignment, int isReversible,
		int numSpecies, int numIterations, int threadID, int masterID, int numCPUThreads,
		ThreadLocks* threadLocks, int* jobAllocated, TempRateMats* preInterResult, UserOptions* userOptions,
		int num_chars) {


	int num_w = 6;
	int num_edges = 2 * numSpecies - 2;
	int num_alpha = 1;
	Optim* curr_opt = new Optim(alignment, topMatrix, isReversible, num_chars);
	double curr_loglikelihood, curr_IC;
	int numIter;
	int currMode;
	VariableSet vs(num_alpha, num_chars);

	ParameterSet ps(num_w, num_edges, num_chars);

	int totalJobs = (int) rateMatChkLst->size();

	// int optMethods[] = {3,5,8}; // optimal optimization methods are 3, 5 and 8 with MaxIt = max_MAXIT
	// int currMaxIt = max_MAXIT;
	int currMaxIt, num_rateMat; //, optMethod;
	int numOptMethods = (int)userOptions->opt_mode.size();

	string ICName = getICName(userOptions->info_criteria);

	for (int rateMatID=0; rateMatID<totalJobs; rateMatID++) {

		bool proceed = false;
		// check which job has not been processed
		if (threadLocks != NULL)
			pthread_mutex_lock(&(threadLocks->lock_jobAllocate));
		if (jobAllocated[rateMatID]==0) {
			jobAllocated[rateMatID]=1;
			proceed = true;
		}
		if (threadLocks != NULL)
			pthread_mutex_unlock(&(threadLocks->lock_jobAllocate));

		if (proceed) {

			// optMethod = userOptions->optMethods[rateMatID % (userOptions->optMethods.size())]-1;
			// currMode = optMethod/2+1;
			// currMaxIt = (optMethod%2==0)?1:max_MAXIT;
			currMode = userOptions->opt_mode[rateMatID%numOptMethods];
			currMaxIt = userOptions->opt_maxIT[rateMatID%numOptMethods];

			// print out the potential rate matrix
			vector<int>* curr_matrix = rateMatChkLst->matrices[rateMatID];
			num_rateMat = rateMatChkLst->numRates[rateMatID];
			
			// compute the number of change of rate matrices
			int numRateChanges = numChanges(topMatrix, curr_matrix, numSpecies-2, 0);
			
			// check whether the rate matrices are the same for two root edges
			int sameRootMat = (curr_matrix->at(num_edges-1)==curr_matrix->at(num_edges-2));

			if (rateMatChkLst->parentPSs[rateMatID]!=NULL && rateMatChkLst->parentVSs[rateMatID]!=NULL) {
				vs.copyFrom(*(rateMatChkLst->parentVSs[rateMatID]));
				ps.copyFrom(*(rateMatChkLst->parentPSs[rateMatID]));
			} else {
				vs.resetAllVariables(*alignment);
				ps.reset();
			}
			// update the rate matrix in the parameter
			ps.updateRateMat(*curr_matrix, num_rateMat);
			curr_opt->setVarParam(&vs, &ps);
			
			// optimizing all parameters for this rate matrices
			numIter = curr_opt->optimizeAllParam(numIterations, curr_loglikelihood, curr_IC, threadLocks, masterID, threadID, currMaxIt,
					currMode, userOptions->info_criteria, numRateChanges, sameRootMat, userOptions->precise);


			// print out the result
			/*
			pthread_mutex_lock(&(threadLocks->lock_output));
			cout << "[master " << masterID << " thread " << threadID << "] ";
			for (int k=0; k<(int)curr_matrix->size(); k++) {
				if (k>0)
					cout << ",";
				cout << curr_matrix->at(k);
			}
			cout << " numRateGrps:" << num_rateMat << " maxIT:" << currMaxIt << " mode:" << currMode;
			cout << " Log-L: " << curr_loglikelihood << " " << ICName << ":" << curr_IC << " numIter:" << numIter << endl << flush;
			pthread_mutex_unlock(&(threadLocks->lock_output));*/

			rateMatChkLst->VSs[rateMatID] = new VariableSet(num_chars);
			rateMatChkLst->PSs[rateMatID] = new ParameterSet(num_chars);
			rateMatChkLst->VSs[rateMatID]->copyFrom(vs);
			rateMatChkLst->PSs[rateMatID]->copyFrom(ps);

			rateMatChkLst->loglikelihood[rateMatID] = curr_loglikelihood;
			rateMatChkLst->ICs[rateMatID] = curr_IC;
			rateMatChkLst->numIters[rateMatID] = numIter;

			// compute the lowest possible IC value due to the deviation of
			// the resulting log likelihood value from the optimization engine
			rateMatChkLst->lowPossibleICs[rateMatID] = getIC(curr_loglikelihood + userOptions->logLDeviation, 
					num_rateMat, alignment->numSites, alignment->numSeqs,
					userOptions->info_criteria, numRateChanges, sameRootMat);


		}
	}
	delete(curr_opt);
}


// ============================================================================================================================
// OTHER FUNCTIONS
// ============================================================================================================================


void findConsistentRates(TempRateMats* bestHundredRateMatrices, vector<int>* consistentRates, int topK) {
	// rearrange all items inside the matrices
	// such that the topK items (among all items) are with the least IC values

	if (topK == 0)
		return;

	bestHundredRateMatrices->rearrangeFirstNItems(topK, bestHundredRateMatrices->size());

	int i,j;
	// for each edge, check whether those of the topK items are the same
	int isConsistent;
	int rateGrp;
	for (i=0; i<(int)consistentRates->size(); i++) {
		isConsistent=1;
		rateGrp = bestHundredRateMatrices->matrices[0]->at(i);
		for (j=1; j<topK && isConsistent; j++) {
			if (rateGrp != bestHundredRateMatrices->matrices[j]->at(i))
				isConsistent = 0;
		}
		consistentRates->at(i) = isConsistent;
	}
}

/*
void printCurrentMatrix(vector<int>& theRateMatArray, double logL, double IC, vector<int>& untouchNodes, int numRateGrp, string& ICName) {
	// print out the current rate matrix
	cout << "current rate matrix:" << endl;
	for (int k=0; k<(int)theRateMatArray.size(); k++) {
		if (k>0)
			cout << ",";
		cout << theRateMatArray[k];
	}
	cout << " [log-L: " << logL << "] ";
	cout << " [" << ICName << ": " << longDoublToStr(IC,ANS_DECI) << "] ";
	cout << " [untouchNode: ";
	for (int k=0; k<(int)untouchNodes.size(); k++) {
		if (k>0)
			cout << ",";
		cout << untouchNodes[k];
	}
	cout << "]";
	cout << " [numRateGrp:" << numRateGrp << "]" << endl;
}
 */

bool getDoubleValues(vector<string>& token, string desc, double* values, int num) {
	// get the values for desc
	// return true if the values are obtained completely

	int i=0;
	bool found = false;
	while ((i<(int)token.size()) && (!found)) {
		if (token[i]==desc)
			found = true;
		i++;
	}
	if (found) {
		int n=0;
		while (((int)token.size()) && (n<num)) {
			values[n] = atof(token[i].c_str());
			i++;n++;
		}
		if (n==num)
			return true;
	}
	return false;
}

bool getDoubleValues(vector<string>& token, double* values, int start, int num) {
	// get the values from start up to start+num
	// return true if the values are obtained completely

	int i=0;
	while (i<num && i+start<(int)token.size()) {
		values[i] = atof(token[i+start].c_str());
		i++;
	}
	if (i==num)
		return true;

	return false;
}

bool getIntValues(vector<string>& token, string desc, int* values, int num) {
	// get the values for desc
	// return true if the values are obtained completely

	int i=0;
	bool found = false;
	while ((i<(int)token.size()) && (!found)) {
		if (token[i]==desc)
			found = true;
		i++;
	}
	if (found) {
		int n=0;
		while ((i<(int)token.size()) && (n<num)) {
			values[n] = atoi(token[i].c_str());
			i++;n++;
		}
		if (n==num)
			return true;
	}
	return false;
}

bool getIntValues(vector<string>& token, vector<int>& values, int start, int num) {
	// get the values from start up to start+num
	// return true if the values are obtained completely

	int i=0;
	while (i<num && i+start<(int)token.size()) {
		values[i] = atoi(token[i+start].c_str());
		i++;
	}
	if (i==num)
		return true;

	return false;
}

// recalculate all the IC values in the matrix
void reCalculateAllICs(TempRateMats* interResults, int ICType, int num_sites, int num_species, double logLDeviation, int* topMatrix) {
	int i;
	double logL;
	double numR;
	int numRateChanges;
	int sameRootMat;
	int numEdges = 2*num_species-2;
	vector<int>* currMatrix;
	for (i=0; i<interResults->size(); i++) {
		logL = interResults->loglikelihood[i];
		numR = interResults->numRates[i];
		numRateChanges = numChanges(topMatrix, interResults->matrices[i], num_species-2, 0);
		currMatrix = interResults->matrices[i];
		sameRootMat = (currMatrix->at(numEdges-1)==currMatrix->at(numEdges-2));
		interResults->ICs[i] = getIC(logL, numR, num_sites, num_species, ICType, numRateChanges, sameRootMat);
		interResults->lowPossibleICs[i] = interResults->ICs[i] - 2.0 * logLDeviation;
	}
}

// Load previous intermediate results from the checkpoint file
TempRateMats* loadPreInterResult(char* chkptFile, int num_w, int num_edges, int num_alpha, int ICType, 
		int num_sites, int num_species, double logLDeviation, int num_chars, int* topMatrix) {
	ifstream fin;
	fin.open(chkptFile);
	string aline;
	vector<string> token;
	size_t pos;

	vector<int> matrix(num_edges);
	double loglikelihood = 0.0;
	double IC = 0.0;
	double lowest_IC = LARGE_NUMBER;
	int numRates = 0;
	ParameterSet PS(num_w,num_edges,num_chars);
	VariableSet VS(num_alpha, num_chars);
	vector<int> untouch(num_edges,0);
	int numIter = 0;
	bool need_recalculate_IC;
	int numRateChanges;

	TempRateMats* interResults = new TempRateMats;
	int numResults = 0;

	// get the name of the IC
	string ICName = getICName(ICType);

	cout << "loading the checkpoint file: " << chkptFile << endl << flush;

	while (getline(fin,aline)) {
		// search for the word "iterations"
		if (aline.find("iterations") != string::npos) {
			// appear

			trim(aline); // first trim the spaces/invisible characters at the end of the line
			tokenizer(aline, " :[],\t", &token);

			// get number of iterations
			if (!getIntValues(token, "iterations", &numIter, 1)) {
				continue;
			}

			// get loglikelihood
			if (!getDoubleValues(token, "loglikelihood", &loglikelihood, 1)) {
				continue;
			}

			// get IC
			need_recalculate_IC = true;
            if (numIter == 0)
                need_recalculate_IC = false;
			/*
			if (!getDoubleValues(token, ICName, &IC, 1)) {
				// The IC name does not match, then the corresponding IC and the lowest possible IC needs to recalculate
				need_recalculate_IC = true;
				// continue;
			} else {
				need_recalculate_IC = false;
			}*/

            if (numIter > 0) {
                // get the details of PS
                int edge = 0;
                while ((edge<num_edges) && (getline(fin,aline))) {
                    trim(aline);
                    tokenizer(aline, " :[],\t", &token);
                    double* curr_w = &(PS.w[edge*num_w]);
                    double* curr_pi = &(PS.pi[edge*4]);
                    double* curr_t = &(PS.t[edge]);
                    if ((!getDoubleValues(token, curr_w, 0, num_w)) ||
                            (!getDoubleValues(token, curr_pi, num_w, 4)) ||
                            (!getDoubleValues(token, curr_t, num_w+4, 1))) {
                        break;
                    }
                    edge++;
                }
                if (edge!=num_edges) {
                    continue;
                }


                // get the details of VS
                // beta
                if (getline(fin,aline)) {
                    trim(aline);
                    tokenizer(aline, " :[],\t", &token);
                    if (!getDoubleValues(token, &(VS.beta), 1, 1)) {
                        continue;
                    }
                } else {
                    continue;
                }

                // alpha
                if (getline(fin,aline)) {
                    trim(aline);
                    tokenizer(aline, " :[],\t", &token);
                    if (!getDoubleValues(token, VS.alpha, 1, num_alpha)) {
                        continue;
                    }
                    VS.num_alpha=num_alpha;
                } else {
                    continue;
                }

                // probXGivenInv
                if (getline(fin,aline)) {
                    trim(aline);
                    tokenizer(aline, " :[],\t", &token);
                    if (!getDoubleValues(token, VS.probXGivenInv, 1, 4)) {
                        continue;
                    }
                } else {
                    continue;
                }

                // rootNodeFreq
                if (getline(fin,aline)) {
                    trim(aline);
                    tokenizer(aline, " :[],\t", &token);
                    if (!getDoubleValues(token, VS.rootNodeFreq, 1, 4)) {
                        continue;
                    }
                } else {
                    continue;
                }
            }
            
			// rate matrix
			if (getline(fin,aline)) {
				trim(aline);
				// remove the "[...]" in the beginning of the line if there is any
				if (aline.length() > 0 && aline[0]=='[') {
					pos = aline.find("]");
					if (pos==string::npos) {
						continue;
					} else {
						aline = aline.substr(pos+1);
					}
				}
				tokenizer(aline, " :[],\t", &token);
				if (!getIntValues(token, matrix, 0, num_edges)) {
					continue;
				}

				// get number of rate groups
				if (!getIntValues(token, "numRateGrps", &numRates, 1)) {
					continue;
				}

				// recalculate the IC value if necessary
				if (need_recalculate_IC) {
					numRateChanges = numChanges(topMatrix, &matrix, num_species-2, 0);
					int sameRootMat = (matrix[num_edges-1]==matrix[num_edges-2]);
					IC = getIC(loglikelihood, numRates, num_sites, num_species, ICType, numRateChanges, sameRootMat);
				}

				// calculate the lowest IC
				lowest_IC = IC - 2.0 * logLDeviation;

			} else {
				continue;
			}

			// load to the interResults array
            if (numIter > 0)
                interResults->insertwMap(matrix, loglikelihood, IC, -1, PS, VS, untouch, numIter, lowest_IC);
            else
                interResults->insertwMap(matrix, loglikelihood, IC, untouch, lowest_IC);
            
			numResults++;
		}
	}
	return interResults;
}


// Load previous log file (not checkpoint file)
void loadPreLogResult(char* logFile, int num_w, int num_edges, int num_alpha, int ICType,
        int num_sites, int num_species, double logLDeviation, int num_chars, int* topMatrix, TempRateMats* interResults) {
    ifstream fin;
    fin.open(logFile);
    string aline;
    vector<string> token;
    size_t s_pos, t_pos;
    string bic_str;
    string rate_arrange_str;
    string loglike_str;

    vector<int> matrix(num_edges);
    vector<int> untouch(num_edges,0);
    double loglikelihood = 0.0;
    double IC = 0.0;
    double lowest_IC = LARGE_NUMBER;
    int numRates;
    int i;

    // get the name of the IC
    string ICName = getICName(ICType);

    cout << "loading the log file: " << logFile << endl << flush;

    while (getline(fin,aline)) {
        if (aline.length() > 5 && aline[aline.length()-1]==']') {
            // IC
            s_pos = aline.find("IC");
            if (s_pos==string::npos)
                continue;
            t_pos = aline.find_first_of("]", s_pos+4);
            if (t_pos==string::npos)
                continue;
            bic_str = aline.substr(s_pos+4, t_pos-s_pos-4);
            IC = atof(bic_str.c_str());

            // get the rate matrix
            s_pos = aline.find_first_of(" ");
            if (s_pos==string::npos)
                continue;
            rate_arrange_str = aline.substr(0,s_pos);
            tokenizer(rate_arrange_str.c_str(), " :[],\t", &token);
            if (!getIntValues(token, matrix, 0, num_edges))
                continue;
            s_pos = aline.find("loglikelihood:");
            if (s_pos==string::npos)
                continue;

            // log-likelihood
            t_pos = aline.find_first_of("]", s_pos+15);
            if (t_pos==string::npos)
                continue;
            loglike_str = aline.substr(s_pos+15, t_pos-s_pos-15);
            loglikelihood = atof(loglike_str.c_str());

            // number of rate matrix
            // minimum 2
            numRates = getNumRateGrp(matrix);
            if (numRates == 1)
                continue;
            
            // calculate the lowest IC
            lowest_IC = IC - 2.0 * logLDeviation;

            // load to the interResults array
            interResults->insertwMap(matrix, loglikelihood, IC, untouch, lowest_IC);
        }
    }
    // interResults->print_content();
    // return interResults;
}


// ==================================================================================================
// Maximize the likelihood value for the specific matrice arrangement
// ==================================================================================================

double getLogL(char* alignFile, char* treeFile, int inputFormat, char* rateGrpFile, int mode, int maxIT, int listOutParamDetail, int constantSiteExist, int includeGap, int num_chars) {
	// inputFormat: 1 - FASTA; 2 - PHY format

	// load the topology
	int* topMatrix;
	int numLineTopMat, numSpecies, numEdges;
	vector<string>* leafList;
	genTopMatrix(treeFile, &topMatrix, &numLineTopMat, &leafList, NULL, NULL);
	numSpecies=(int)leafList->size();
	if (numSpecies < 2) {
		cerr << "Error! The number of species is too small (i.e. less than 2) for this algorithm!" << endl;
		exit(1);
	}
	numEdges=numSpecies*2-2;
	map<string,int>* leafPosMap = genReverseMap(leafList);
	int* leafID2line = new int[numSpecies+1];
	int* internalID2line = new int[numLineTopMat+1];
	getMapping(leafID2line, internalID2line, topMatrix, numLineTopMat);

	// load the alignment file
	Alignment alignment;
	if (constantSiteExist==0)
		alignment.setGTRUpilson(); // no site is regarded as constant
	if (inputFormat==1)
		alignment.readFASTA(alignFile, leafPosMap);
	else
		alignment.readFile(alignFile, leafPosMap);
	alignment.keepSitesWithChars(includeGap); // remove the columns with invalid characters
	int* siteOrder = alignment.getSortedSitesID(); // get the order of the sites according to the topology
	alignment.computeNumUniqueSites(siteOrder); // compute the number of unique columns
	alignment.keepUniqueSites(siteOrder); // keep only the unique columns


	// initialize the parameters
	int isSingleGrp = 1;
	int num_w = 6;
	ParameterSet ps(num_chars);
	ps.initialize(num_w, numEdges, isSingleGrp);

	// initialize the variables
	int numRateCat = 1;
	VariableSet variables(numRateCat, num_chars);
	variables.resetAllVariables(alignment);

	// load the rate matrix
	ps.loadRateMat(rateGrpFile, topMatrix, leafList);

	// refine the rate matrix
	vector<int> theRateMatArray; // the current rate matrix
	ps.outputRateMat(theRateMatArray);
	refineRateMatrix(theRateMatArray);

	// load the refined rate matrix into the parameter
	ps.loadRateMat(theRateMatArray);
	
	// get the number of changes of rate matrices
	int numRateChanges = numChanges(topMatrix, &theRateMatArray, numSpecies-2, 0);

	// show the current rate matrix
	cout << "Rate matrix arrangement: ";
	int i;
	for (i=0; i<(int)theRateMatArray.size(); i++) {
		if (i>0)
			cout << ",";
		cout << theRateMatArray[i];
	}
	cout << endl;

	// optimization
	int isReversible = 1;
	int num_chars_to_consider = 4;
	Optim* curr_opt = new Optim(&alignment, topMatrix, isReversible, num_chars_to_consider);
	curr_opt->setVarParam(&variables, &ps);

	// optimizing all parameters for this rate matrices
	int numIterations = -1; // no limit on the number of iterations
	double loglikelihood;
	double bicValue;
	int infoCriteria = 3; // BIC
	curr_opt->optimizeAllParam(numIterations, loglikelihood, bicValue, NULL, 0, 0, maxIT, mode, infoCriteria, numRateChanges, 0, 1);

	if (listOutParamDetail==1) {
		ps.showContent();
		variables.showContent();
	}

	// release the memory
	delete curr_opt;

	return loglikelihood;
}


// ==================================================================================================
// Compute the likelihood values for the corresponding parameter and variable values
// ==================================================================================================

double getLogLFrFiles(char* alignFile, char* treeFile, char* siteInfoFile, char* paramFilelist, int edgeRepresent, int numRateMatrix, int ICType, double& IC, int& df, int includeGap, int num_chars) {

	// edgeRepresent: representation of the edge length
	//                1 - Average number of substitutions per site
	//                2 - Time
	// ICType: type of Information Criteria
	//                1 - AIC; 2 - Adjusted BIC; 3 - BIC; 4 - CAIC

	// load the topology
	int* topMatrix;
	int numLineTopMat, numSpecies, numEdges, numCategories;
	vector<string> leafList;
	vector<double> edgeLens;
	vector<string> nodeList;
	vector<string> siteCatList;
	topMatrix = genTopMatrix2(treeFile, numLineTopMat, leafList, edgeLens, nodeList, siteCatList);

	// print out the topMatrix
	for (int i=0;i<numLineTopMat; i++) {
		cout << topMatrix[i*2] << "," << topMatrix[i*2+1] << endl << flush;
	}

	numSpecies=(int)leafList.size();
	numEdges=numSpecies*2-2;
	numCategories=(int)siteCatList.size();

	// load the variables
	VariableSet vs(numCategories, num_chars);
	vs.readSiteInfoFile(siteInfoFile, &siteCatList);

	// print out the content of variables
	vs.showContent();

	// load the parameter list
	AllParameterSet allPS(numCategories, num_chars);

	allPS.readParamFileList(paramFilelist, siteCatList, nodeList, edgeLens, numEdges, edgeRepresent);

	// print out the content of variables
	allPS.showContent();

	// load the alignment file
	Alignment alignment;
	map<string,int>* leafPosMap = genReverseMap(&leafList);
	alignment.readFile(alignFile, leafPosMap);
	alignment.keepSitesWithChars(includeGap); // remove the columns with invalid characters
	int* siteOrder = alignment.getSortedSitesID(); // get the order of the sites according to the topology
	alignment.computeNumUniqueSites(siteOrder); // compute the number of unique columns
	alignment.keepUniqueSites(siteOrder); // keep only the unique columns


	int isReversible = 1;
	ParameterSet* ps = allPS.ps[0];
	int edge;
	for (edge=0; edge<numEdges; edge++)
		ps->computeAllEigenValues_t(edge);
	ps->computeAllCondProb(isReversible);

	// print out the conditional probabilities
	ps->printAllCondProb();

	int num_chars_to_consider = 4;    
	double logL = getLogLikelihood(ps->allCondProb_t, topMatrix, numLineTopMat, alignment, vs, num_chars_to_consider);
	IC = getIC(logL, numRateMatrix, alignment.numSites, numSpecies, df, ICType, 0);

	return logL;
}

// ==================================================================================================
// use the consensus method to come up an answer
// ==================================================================================================

void consensusMethod(ofstream* fout, TempRateMats* allResults, UserOptions* options, double IC_Thres,
		int numSpecies, int nEdge, int numSites, int* topMatrix, vector<string>* leafList) {

	vector<int> numRateGrps;
	vector<double> logLs;
	vector<string> rateArranges;
	double optICs[TOTAL_NUM_IC];
	int optICRsts[TOTAL_NUM_IC];
	double* allICValues = new double[allResults->size() * TOTAL_NUM_IC];
	int numRateChanges;
	int sameRootRate;
	vector<int>* currMatrix;
	int i,j,k;

	// initialization
	for (i=0; i<TOTAL_NUM_IC; i++) {
		optICs[i] = optICRsts[i] = 0;
	}

	// how many different IC types have been selected
	int ICTypeNum = 0;
	for (i=0; i<TOTAL_NUM_IC; i++)
		ICTypeNum += options->ic_list[i];

	// first for each selected IC type, get the lowest value
	for (i=0; i<TOTAL_NUM_IC; i++) {
		optICs[i] = -1;
	}
	for (i=0; i<allResults->size(); i++) {
		numRateChanges = numChanges(topMatrix, allResults->matrices[i], numSpecies-2, 0);
		for (j=0; j<TOTAL_NUM_IC; j++) {
			if (options->ic_list[j]==1) {
				// this IC type is selected
				currMatrix=allResults->matrices[i];
				sameRootRate=(currMatrix->at(nEdge-1)==currMatrix->at(nEdge-2));
				allICValues[i*TOTAL_NUM_IC+j] = getIC(allResults->loglikelihood[i], allResults->numRates[i], 
						numSites, numSpecies, j+1, numRateChanges, sameRootRate);
				if (optICs[j]==-1 || allICValues[i*TOTAL_NUM_IC+j] < optICs[j]) {
					optICs[j] = allICValues[i*TOTAL_NUM_IC+j];
					optICRsts[j] = i;
				}
			}
		}
	}

	// compute the weight for each result
	double* weight = new double[allResults->size()];
	memset(weight, 0, allResults->size()*sizeof(double));
	double curr_weight;
	for (i=0; i<allResults->size(); i++) {
		for (j=0; j<TOTAL_NUM_IC; j++) {
			if (options->ic_list[j]==1) {
				// this IC type is selected
				curr_weight = (optICs[j] + IC_Thres - allICValues[i*4+j]) / (IC_Thres * ICTypeNum);
				if (curr_weight > 0.0)
					weight[i] += curr_weight;
			}
		}
	}

	// get the range of number of rate groups to be considered
	int numRateGrpMin = -1;
	int numRateGrpMax = -1;
	int currNumRateGrp;
	if (ICTypeNum > 1) {
		for (i=0; i<TOTAL_NUM_IC; i++) {
			if (options->ic_list[i]==1) {
				// this IC type is selected
				currNumRateGrp = allResults->numRates[optICRsts[i]];
				if (numRateGrpMin == -1 || currNumRateGrp < numRateGrpMin) {
					numRateGrpMin = currNumRateGrp;
				}
				if (numRateGrpMax == -1 || currNumRateGrp > numRateGrpMax) {
					numRateGrpMax = currNumRateGrp;
				}
			}
		}
	}

	// compute the scores for all relationships
	// but only consider those results with weight > 0
	double* match_scores = new double[nEdge * nEdge];
	double* unmatch_scores = new double[nEdge * nEdge];
	memset(match_scores, 0, nEdge*nEdge*sizeof(double));
	memset(unmatch_scores, 0, nEdge*nEdge*sizeof(double));
	for (i=0; i<allResults->size(); i++) {
		if (weight[i] > 0.0 && allResults->numRates[i] >= numRateGrpMin
				&& (numRateGrpMax==-1 || allResults->numRates[i] <= numRateGrpMax)) {
			for (j=0; j<nEdge-1; j++) {
				for (k=j+1; k<nEdge; k++) {
					if (allResults->matrices[i]->at(j)==allResults->matrices[i]->at(k)) {
						// the same group
						match_scores[j*nEdge + k] += weight[i];
					} else {
						// not the same group
						unmatch_scores[j*nEdge + k] += weight[i];
					}
				}
			}
		}
	}

	// check each result and see which has the highest
	double* consensusScores = new double[allResults->size()];
	memset(consensusScores, 0, allResults->size()*sizeof(double));
	double bestConsensusScore = 0.0;
	int optResult = -1;
	for (i=0; i<allResults->size(); i++) {
		if (weight[i] > 0.0 && allResults->numRates[i] >= numRateGrpMin
				&& (numRateGrpMax==-1 || allResults->numRates[i] <= numRateGrpMax)) {
			for (j=0; j<nEdge-1; j++) {
				for (k=j+1; k<nEdge; k++) {
					if (allResults->matrices[i]->at(j)==allResults->matrices[i]->at(k)) {
						// the same group
						consensusScores[i] += match_scores[j*nEdge + k];
					} else {
						// not the same group
						consensusScores[i] += unmatch_scores[j*nEdge + k];
					}
				}
			}
			if (optResult == -1 || bestConsensusScore < consensusScores[i]) {
				optResult = i;
				bestConsensusScore = consensusScores[i];
			}
		}
	}

	if (optResult != -1) {
		// list out the result
		(*fout) << "[Best rate matrix arrangement by using consensus approach]" << endl;
		(*fout) << rateMatrixToTreeFormat(topMatrix, allResults->matrices[optResult], leafList);
		for (i=0; i<TOTAL_NUM_IC; i++) {
			if (options->ic_list[i]) {
				(*fout) << " " << getICName(i+1) << ": " << longDoublToStr(allICValues[optResult*4+i],ANS_DECI) << ";";
			}
		}
		(*fout) << endl;
	}
}
