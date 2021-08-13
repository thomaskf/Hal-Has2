/*
 *
 * RAL_mpi.cpp
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
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITeH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
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

#include "RAL_mpi.h"
// #define DEBUG_MD

// perform RAL algorithm
//
// options              : the user options
// rateGrpFile          : the rate group to start with
// preLogFile           : The output file of the previous run.
//                        One wants to resume the process from the previous run.
// world_size           : total number of processes ( > 1)
// world_rank           : which processor is runnning the function
void performRAL(UserOptions* options, char* rateGrpFile, char* preLogFile, int num_chars,
				int world_size, int world_rank) {

	double thresHeight = 0.025;
	int isReversible=1; int num_w=6;
	double IC;

	// ======================
	// for ALL PROCESSORS
	// ======================

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
		// 1 - ignore all columns with gaps
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
	//alignment.reorderUniqueSites(); // for efficient calculation of log-likelihood
	alignment.buildChangedNodes(topMatrix, numLineTopMat);

	int num_chars_to_consider = 4;

	int tagID = DEFAULT_RAL_TAGID;

	// =================================
	// FOR OTHER PROCESSORS EXCEPT 0
	// =================================

	if (world_rank > 0) {

		optimThread(topMatrix, numLineTopMat, &alignment, isReversible,
				numSpecies, options, num_chars_to_consider, tagID);
		return;
	}

	int isSingleGrp = 1;
	if (options->HALMode==2) {
		isSingleGrp = 0; // top-down approach
	}

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
	if (rateGrpFile != NULL) {
		start_ps.loadRateMat(rateGrpFile, topMatrix, leafList);
	}
	vector<int> startRateMatArray; // the starting rate matrix
	start_ps.outputRateMat(startRateMatArray);

	vector<int> untouchNodes(numEdges,0);

	// Load previous intermediate results from the log file if necessary
	TempRateMats* preInterResult = NULL;
	if (preLogFile != NULL) {
		int num_alpha = 1; // for RAL
		options->info_criteria = firstNonZeroPos(options->ic_list, TOTAL_NUM_IC)+1;
		preInterResult = loadPreInterResult(preLogFile, num_w, numEdges, num_alpha, options->info_criteria, 
				alignment.numSites, alignment.numSeqs, options->logLDeviation, num_chars, topMatrix);
	} else {
		preInterResult = new TempRateMats;
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
	string outChkptFile;
	if (options->outChkPtFile && options->outputPrefix.length() > 0) {
		outChkptFile = options->outputPrefix + "." + approachTerm + ".chkpt.txt";
		chkptOut = new ofstream;
		if (preLogFile != NULL) {
			if (outChkptFile != string(preLogFile)) {
				filecopy((char*) outChkptFile.c_str(), preLogFile);
			}
			chkptOut->open(outChkptFile.c_str(), std::ofstream::out | std::ofstream::app); // append to the log file
		} else {
			chkptOut->open(outChkptFile.c_str());
		}
	}

	// The optimal rate matrix for each number of rates allowed
	OptimalRateMats optRateMatrix;

	// The rate matrix waiting for process
	TempRateMats rateMatrixList;

	// The first one hundred best rate matrices
	TempRateMats bestHundredRateMatrices;

	// perform RAL iteration for each user defined information criteria
	int i,j,k;
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
		for (j=0;j<(int)options->opt_mode.size();j++)
			rateMatrixList.insert(startRateMatArray, untouchNodes, IC, -1, 1, options->opt_mode[j], options->opt_maxIT[j]);

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

		RALIterationNew (&rateMatrixList, topMatrix, numLineTopMat,
				&alignment, isReversible, numSpecies,
				&optRateMatrix, leafID2line, internalID2line, thresHeight,
				chkptOut, preInterResult, &bestHundredRateMatrices, options, num_chars_to_consider, saveToInterResult,
				tagID, start_ps, start_vs);

		// sort the set of optimal rate matrices
		bestHundredRateMatrices.rearrangeFirstNItems(100, bestHundredRateMatrices.size());
		

		// result file
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
		bestHundredRateMatrices.print_content(resultOut, topMatrix, leafList, description, ICName, 100);

		if (options->listParamBestHAL && bestHundredRateMatrices.size() > 0) {
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
		if ((options->HALMode==1 || options->HALMode==4 || options->HALMode==7 || options->HALMode==8) && (options->continueTD==1)) {
			options->HALMode=2;
			
			// initialize untouchNodes
			for (j=0; j<numEdges; j++) {
				untouchNodes[j] = 0;
			}
			// initialize the matrix list for processing
			IC = 1e60;
			rateMatrixList.clear();
			for (j=0; j<options->topK && j<(int)bestHundredRateMatrices.size(); j++) {
				for (k=0;k<(int)options->opt_mode.size();k++) {
					rateMatrixList.insert(*(bestHundredRateMatrices.matrices[j]), untouchNodes, IC, -1, 1,
										  options->opt_mode[k], options->opt_maxIT[k]);
				}
			}
			// start
			cout << "Starting the top-down algorithm..." << endl;
			
			RALIterationNew (&rateMatrixList, topMatrix, numLineTopMat,
							 &alignment, isReversible, numSpecies,
							 &optRateMatrix, leafID2line, internalID2line, thresHeight,
							 chkptOut, preInterResult, &bestHundredRateMatrices, options, num_chars_to_consider, saveToInterResult,
							 tagID, start_ps, start_vs);
			
			// sort the set of optimal rate matrices
			bestHundredRateMatrices.rearrangeFirstNItems(100, bestHundredRateMatrices.size());

			// result file
			approachTerm += ".HAL";
			if (options->outputPrefix.length() > 0) {
				resultOut = new ofstream;
				resultFileName = options->outputPrefix + "." + approachTerm + "." + ICName + ".result.txt";
				resultOut->open(resultFileName.c_str());
			}
			
			// optRateMatrix.print_optimal_rateMatrix(resultOut, topMatrix, leafList);
			// if (options->outputConsensus)
			//     bestHundredRateMatrices.print_consensus(resultOut, topMatrix, leafList, ICName, options->logLDeviation*2.0);
			// else
			bestHundredRateMatrices.print_optimal_rateMatrices(resultOut, topMatrix, leafList, alignment.numSites, options->info_criteria);
			// bestHundredRateMatrices.print_content();
			
			// print out the content of the optimal rate matrices
			optRateMatrix.print_content(resultOut, topMatrix, leafList, ICName);
			
			// print out the content of the best hundred rate matrices
			string description = "The best 100 arrangements of the rate matrices:";
			bestHundredRateMatrices.print_content(resultOut, topMatrix, leafList, description, ICName, 100);
			
			if (options->listParamBestHAL && bestHundredRateMatrices.size() > 0) {
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
		

		firstIC = 0;

	}
	
	// terminate all the processes
	terminateAllProcs(numEdges, num_w, num_chars, tagID);

	if (chkptOut != NULL)
		chkptOut->close();

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
	cout << "====================================================================" << endl;

	// list out the names of the resulting files
	for (i=0; i<TOTAL_NUM_IC; i++) {
		if (options->ic_list[i] == 0)
			continue;
		ICName = getICName(i+1);
		resultFileName = options->outputPrefix + "." + approachTerm + "." + ICName + ".result.txt";
		if (resultFileName.length() > 0)
			cout << "Result file according to "<< ICName << ": " << resultFileName << endl;
	}
	if (options->outputConsensus == 1) {
		resultFileName = options->outputPrefix + "." + approachTerm + ".consensus.result.txt";
		if (resultFileName.length() > 0)
			cout << "Result file of consensus approach: " << resultFileName << endl;
	}

	cout << "HAL-" << approachTerm << " Version " << VERSION << " finishes" << endl;


}


// ============================================================================================================================
// ITERATION FUNCTIONS for RAL
// To GENERATE NEW MODELS for MPI
// ============================================================================================================================
// check the number of items in resultList
// if the number of items > currLevelNumResults - topK
// generates some new models and can start some next iterations
// ============================================================================================================================

#define LAST_STAGE 0


int generateNewModelsMPI(TempRateMats* jobList, TempRateMats* resultList, int* leafID2line, int* internalID2line,
						 UserOptions* options, vector<int>& firstJobID, vector<int>& levelID, vector<int>& grouping_array,
						 int currLevelNumResults, int* topMatrix, int topK, double thresHeight, int whichApproach,
						 int& nextAvailableLevel, int numSpecies, int numEdges, int untouchedNodes, int opt_method_num, string& outMsg,
						 ParameterSet& start_ps, VariableSet& start_vs, int* edgeOrder) {

	int new_candidates = 0;
	if (resultList->size() > currLevelNumResults - topK) {
		int j,k,l,m;

		double IC_Thres;
		int groupAllCandidates = options->groupAllCandidates;
		if (groupAllCandidates==1 && whichApproach!=1 && whichApproach!=7) {
			IC_Thres = IC_Thres_groupAll1;
		} else {
			IC_Thres = IC_Thres_groupAll2;
		}
		
		
		int topN = topK + resultList->size() - currLevelNumResults;
		if (topN > resultList->size())
			topN = resultList->size();
		
		// rearrange the last 'numMatToCheck' items inside the matrices
		// such that the topN items (among the last 'resultList->size()' items) are with the least IC values
		resultList->createTopN(topK);
		// resultList->rearrangeLastNItems(topN, resultList->size());
		
		// generate the new models
		for (m=0; m<topN; m++) {
		// for (j=0; j<topN; j++) {
			j = resultList->sortedNItems[m];
			if (resultList->selected[j]==0 && resultList->ICs[j] <= resultList->parentICs[j] * IC_Thres) {
				// the IC is small enough
				resultList->selected[j] = 1;
				int curr_new_candidates = 0;
				if ((whichApproach==1 || whichApproach==3) && resultList->numRates[j]<numEdges) {
					// generate all the potential rate matrices into the rateMatrixList (using bottom-up approach)
					curr_new_candidates += generateRateMat(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->ICs[j], topMatrix, leafID2line, internalID2line, resultList->untouches[j], numSpecies, resultList->numRates[j], &(options->opt_mode), &(options->opt_maxIT));
				}
				if ((whichApproach==2 || whichApproach==3) && resultList->numRates[j]>1) {
					// generate all the potential rate matrices into the rateMatrixList (using top-down approach)
					curr_new_candidates += generateRateMatTD(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->ICs[j], topMatrix, internalID2line, resultList->untouches[j], resultList->numRates[j], thresHeight, &(options->opt_mode), &(options->opt_maxIT));
				}
				if (whichApproach==4 && untouchedNodes > LAST_STAGE) {
					// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 1 approach)
					curr_new_candidates += generateRateMatNewBU1(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->loglikelihood[j], resultList->parentICs[j], resultList->ICs[j], topMatrix, leafID2line, internalID2line, resultList->untouches[j], numSpecies, resultList->numRates[j], resultList->lowPossibleICs[j],
																 &(options->opt_mode), &(options->opt_maxIT), edgeOrder);
				}
				if (whichApproach==6 && nextAvailableLevel<numEdges) {
					// generate all the potential rate matrices into the rateMatrixList (using bute-force approach)
					curr_new_candidates += generateRateMatButeForce(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->loglikelihood[j], resultList->parentICs[j], resultList->ICs[j], topMatrix, leafID2line, internalID2line, resultList->untouches[j], numSpecies, options->maxNumRatMat, resultList->lowPossibleICs[j], &(options->opt_mode), &(options->opt_maxIT));
				}
				if ((whichApproach==7) && resultList->numRates[j]<numEdges) {
					// generate all the potential rate matrices into the rateMatrixList (using John's bottom-up approach)
					curr_new_candidates += generateRateMatJohn(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->ICs[j], topMatrix, leafID2line, internalID2line, resultList->untouches[j], numSpecies, resultList->numRates[j], &(options->opt_mode), &(options->opt_maxIT));
				}
				if (whichApproach==8 && nextAvailableLevel<numEdges-1) {
					// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 2 approach)
					curr_new_candidates += generateRateMatNewBU2(jobList, resultList->matrices[j], resultList->PSs[j], resultList->VSs[j], resultList->loglikelihood[j], resultList->parentICs[j], resultList->ICs[j], topMatrix, leafID2line, internalID2line, resultList->untouches[j], numSpecies, resultList->numRates[j], resultList->lowPossibleICs[j], &(options->opt_mode), &(options->opt_maxIT));
				}
				if (curr_new_candidates > 0) {
					// print out the new matrice to the messge box
					outMsg.append("[");
					for (int n=0; n<(int)resultList->matrices[j]->size(); n++) {
						if (n > 0)
							outMsg.append(",");
						outMsg.append(intToStr(resultList->matrices[j]->at(n)));
					}
					outMsg.append("] produces " + intToStr(curr_new_candidates) + " new candidates\n");
					
					// set up the array firstJobID
					int firstNewJob = jobList->size() - curr_new_candidates;
					if (curr_new_candidates % opt_method_num == 0) {
						for (k=0; k<curr_new_candidates/opt_method_num; k++) {
							for (l=0; l<opt_method_num; l++) {
								firstJobID.push_back(firstNewJob + k*opt_method_num);
							}
						}
					}
					
					// set up the array levelID
					for (k=0; k<curr_new_candidates; k++) {
						levelID.push_back(nextAvailableLevel);
						// cout << "levelID.push_back(" << nextLevelID << ")" << endl;
					}
					
					// set up the parent PS and VS
					if (options->startOptState==1 || options->startOptState==2) {
						for (k=0; k<curr_new_candidates; k++) {
							if (options->startOptState==2) // 2: optimization starts from the state using one model (default)
								jobList->assignParent(start_ps, start_vs, k+firstNewJob);
							else if (options->startOptState==1) // 1: from the tunned state of the parent configuration (obsoleted)
								jobList->assignParent(*(resultList->PSs[j]),
													  *(resultList->VSs[j]), k+firstNewJob);
						}
					}
				}
				if (groupAllCandidates==2 && curr_new_candidates>0) {
					// The potential candidates form separate groups according to the rate matrix generated from. (i.e. Original design)
					grouping_array.push_back(curr_new_candidates);
					nextAvailableLevel++;
				}
				new_candidates += curr_new_candidates;
			}
		}
	}
	return new_candidates;
}



// ============================================================================================================================
// ITERATION FUNCTIONS for RAL
// Higher flexibility on parallel computation
// ============================================================================================================================

void RALIterationNew (TempRateMats* jobList, int* topMatrix, int numLineTopMat,
					  Alignment* alignment, int isReversible, int numSpecies,
					  OptimalRateMats* optRateMatrix, int* leafID2line, int* internalID2line,
					  double thresHeight, ofstream* fout, TempRateMats* preInterResult,
					  TempRateMats* bestHundredRateMatrices, UserOptions* options, int num_chars,
					  int saveToInterResult, int tagID, ParameterSet& start_ps, VariableSet& start_vs) {
	
	
	int j,k,l;
	int num_w = 6;
	int numEdges = numSpecies*2-2;
	int num_alpha = 1;
	int topK = options->topK;
	int whichApproach = options->HALMode;
	int groupAllCandidates = options->groupAllCandidates;
	int maxNumRounds = options->numMaxCycles;
	double preBestIC = LARGE_NUMBER;
	double currBestIC = LARGE_NUMBER;
	vector<int> untouchzero(numEdges,0);
	
	// store the results of the current level
	TempRateMats* resultList = new TempRateMats();
	
	// the appending results of the level > the current level
	TempRateMats* appendResultList = new TempRateMats();
	// the corresponding level ID
	vector<int> appendResultLevelID;
	
	// indicate whether the job is needed to output to the checkpoint file
	vector<int> jobNeedToOutChkptFile;
	
	// number of processors available
	int world_size = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int totalP = world_size - 1;
	int availP = totalP;
	
	// number of optimization methods
	int opt_method_num = (int) options->opt_mode.size();

	// indicate which processors are available
	vector<int> availProcessors(totalP, 1);
	
	// inidcate whether the parameters and variables have been initialized for each processor
	vector<int> psvsInitialized(totalP, 0);
	
	// indicate the levelID for each job
	vector<int> levelID(jobList->size(), 0);
	
	// indicate the jobID of the first job of the same matrix arrangement
	// i.e. if there are more than one optmization methods, then there will be more than one job for the same matrix arrangement
	// vector<int> firstJobID(jobList->size(), 0);
	vector<int> firstJobID;
	if (jobList->size() % opt_method_num == 0) {
		for (k=0; k<jobList->size()/opt_method_num; k++) {
			for (l=0; l<opt_method_num; l++) {
				firstJobID.push_back(k*opt_method_num);
			}
		}
	} else {
		cerr << "[internal error] jobList->size() % opt_method_num != 0, please contact the software developer." << endl;
		exit(1);
	}
	
	// grouping array, which indicates the next group of rate matrices to be processed
	vector<int> grouping_array;
	grouping_array.push_back(jobList->size());
	
	// indicate the starting item inside resultList for the current level
	// int startItem = 0;
	
	// number of new models generated for next level
	int new_candidates = 0;
	
	// indicate whether it reaches the end of the round
	int reachRdEnd = 0;
	int nextReachRdEnd = 0;
	
	// indicate whether it can proceed to model generation
	int proceedModelGen = 0;
	
	// indicate whether the procedure should stop
	int procedureEnd = 0;
	
	// indicate whether the start_ps and start_vs has been set
	int start_ps_set = 0;
	
	// next available level (for groupAllCandidates = 1)
	int nextAvailableLevel = 1;
	
	// the order of the edges to be examined
	int* edgeOrder = NULL;

	string ICName = getICName(options->info_criteria);

	int currRoundID = 1; // current round ID
	int currLevelID = 0; // current level ID
	int currLevelNumResults = grouping_array[currLevelID] / opt_method_num;
	int iterationID = 1;
	
	// int outStandFinishJobs = 0;
	int currJob = 0;
	int readyJobExist = 0;
	int readyJobID = 0;
	string outMsg = "";
	
	if (jobList->size() > 0)
		cout << "[master " << iterationID << "] Checking " << jobList->size() << " rate matices...." << endl;
	
	while ((!procedureEnd && currJob < jobList->size() && availP > 0) || (!procedureEnd && readyJobExist) || (availP < totalP)) {
		// either some processors are still running or some jobs are not yet assigned to some processors
		
		// ==========================================
		// Send the jobs to available processors
		// ==========================================

		// there are jobs waiting for processing and there is available processor
		if (!procedureEnd && currJob < jobList->size() && availP > 0) {
			
			while (!readyJobExist && currJob < jobList->size() && availP > 0) {
				if (jobList->PSs[currJob]!=NULL) {
					readyJobExist = 1;
					readyJobID = currJob;
					jobNeedToOutChkptFile.push_back(0); // no need to output the result to checkpoint file
				} else {
					int interID = preInterResult->getID(*(jobList->matrices[currJob]));
					if (interID >= 0) {
						readyJobExist = 1;
						readyJobID = currJob;
						jobList->PSs[currJob] = new ParameterSet(num_chars);
						jobList->PSs[currJob]->copyFrom(*(preInterResult->PSs[interID]));
						jobList->VSs[currJob] = new VariableSet(num_chars);
						jobList->VSs[currJob]->copyFrom(*(preInterResult->VSs[interID]));
						jobList->loglikelihood[currJob] = preInterResult->loglikelihood[interID];
						jobList->ICs[currJob] = preInterResult->ICs[interID];
						jobList->numIters[currJob] = preInterResult->numIters[interID];
						jobList->lowPossibleICs[currJob] = preInterResult->lowPossibleICs[interID];
						jobNeedToOutChkptFile.push_back(0);
					} else {
						jobNeedToOutChkptFile.push_back(1);
					}
				}
				if (!readyJobExist) {
					// send the job to a processor
					for (j=0;j<totalP;j++) {
						if (availProcessors[j] == 1) {
							availProcessors[j] = 0;
							availP--;
							// check for the case when options->startOptState==2
							if (options->startOptState == 2) {
								if (start_ps_set==1 && psvsInitialized[j]==1) {
									sendInfo(jobList->matrices[currJob], currJob, jobList->modes[currJob], jobList->maxITs[currJob], 1, numEdges, num_w, num_chars, j+1, tagID, NULL, NULL);
								} else {
									if (start_ps_set==1)
										psvsInitialized[j]=1;
									sendInfo(jobList->matrices[currJob], currJob, jobList->modes[currJob], jobList->maxITs[currJob], 1, numEdges, num_w, num_chars, j+1, tagID, jobList->parentPSs[currJob], jobList->parentVSs[currJob]);
								}
							} else {
								sendInfo(jobList->matrices[currJob], currJob, jobList->modes[currJob], jobList->maxITs[currJob], 1, numEdges, num_w, num_chars, j+1, tagID, jobList->parentPSs[currJob], jobList->parentVSs[currJob]);
							}
							break;
						}
					}
				}
				currJob++;
			}
		}
		
		// ==========================================
		// process the finished job list
		// ==========================================
		
		if (!procedureEnd && readyJobExist) {
			// process the finished job list
			int canProcessForTheJob = 1;
			int frID = firstJobID[readyJobID];
			int toID = frID + opt_method_num - 1;
			int opt_j = -1;
			double opt_IC = LARGE_NUMBER;
			for (j=frID; j<=toID; j++) {
				if (jobList->PSs[j]==NULL) {
					canProcessForTheJob = 0;
					break;
				} else {
					if (jobList->ICs[j] < opt_IC) {
						opt_j = j;
						opt_IC = jobList->ICs[j];
					}
				}
			}
			
			if (canProcessForTheJob && opt_j != -1) {
				// output the intermediate results
				if (jobNeedToOutChkptFile[opt_j]==1) {
					
					// all numbers are rounded up, so that
					// it is consistent with the intermediate results loaded from chkpt file
					jobList->PSs[opt_j]->roundContent();
					jobList->VSs[opt_j]->roundContent();
					jobList->loglikelihood[opt_j] = roundNumber(jobList->loglikelihood[opt_j], ANS_DECI);
					jobList->ICs[opt_j] = roundNumber(jobList->ICs[opt_j], ANS_DECI);
					jobList->lowPossibleICs[opt_j] = roundNumber(jobList->lowPossibleICs[opt_j], ANS_DECI);
					
					// save the intermediate result to the array if necessary
					if (saveToInterResult) {
						preInterResult->insertwMap(*(jobList->matrices[opt_j]), jobList->loglikelihood[opt_j], jobList->ICs[opt_j],
												   -1, *(jobList->PSs[opt_j]), *(jobList->VSs[opt_j]), untouchzero,
												   jobList->numIters[opt_j], jobList->lowPossibleICs[opt_j]);
					}
					
					// need to output the result
					// update the parameters such that the edge length equals to the average # of substitutions
					// rateMatrixList->PSs[i]->updateContent(1);
					if (fout != NULL) {
						// output to the checkout file
						(*fout) << "===============================================================================================" << endl;
						(*fout) << "[# of iterations:" << jobList->numIters[opt_j] << "] loglikelihood:" << longDoublToStr(jobList->loglikelihood[opt_j],ANS_DECI);
						(*fout) << " " << ICName << ":" << longDoublToStr(jobList->ICs[opt_j],ANS_DECI);
						(*fout) << endl;
						jobList->PSs[opt_j]->showContent(fout);
						jobList->VSs[opt_j]->showContent(fout);
						for (int k=0; k<(int)jobList->matrices[opt_j]->size(); k++) {
							if (k>0)
								(*fout) << ",";
							(*fout) << jobList->matrices[opt_j]->at(k);
						}
						(*fout) << " [numRateGrps: " << jobList->numRates[opt_j] << "] " << endl << flush;
					}
				}
				
				// update the optimal rate matrix
				optRateMatrix->insertOpt(jobList->numRates[opt_j], *(jobList->matrices[opt_j]), jobList->loglikelihood[opt_j],
										 jobList->ICs[opt_j]);
				
				// for startOptState = 2 (i.e. optimization starts from the state using one model)
				// if this is the state where only one rate group is involved,
				// then set the starting parameters and starting variables
				if (options->startOptState==2 && jobList->numRates[opt_j]==1 && start_ps_set==0) {
					start_ps.copyFrom(*jobList->PSs[opt_j]);
					start_vs.copyFrom(*jobList->VSs[opt_j]);
					start_ps_set = 1;
				}
				
				if (jobList->numRates[opt_j]==1 && whichApproach==4) {
					// if this is the state where only one rate group is involved
					// then compute the order of the edges
					edgeOrder = getEdgeOrder(jobList->PSs[opt_j]);
				}
				
				
				// if this is the end of the round, then reset the untouchedNodes when necessary
				// (for "new" and "new2" algorithm only)
				int untouchedNodes=0;
				if (whichApproach==4 || whichApproach==8) {
					for (j=0; j<numEdges; j++) {
						if (jobList->untouches[opt_j]->at(j)==0)
							untouchedNodes++;
					}
					if ((untouchedNodes==LAST_STAGE) && (groupAllCandidates==1)) {
						// reach to the end of the round for the "new" and the "new2" algorithms,
						if (maxNumRounds == -1 || currRoundID < maxNumRounds) {
							// reset all the untouchedNodes for those edges
							for (j=0; j<numEdges; j++)
								jobList->untouches[opt_j]->at(j)=0;
							untouchedNodes = numEdges;
						}
						if (levelID[opt_j] == currLevelID)
							reachRdEnd = 1;
						else
							nextReachRdEnd = 1;
					}
				}
				
				if (levelID[opt_j] == currLevelID) {
					resultList->insert(jobList, opt_j, 1);
					proceedModelGen = 1;
				} else {
					appendResultList->insert(jobList, opt_j, 1);
					appendResultLevelID.push_back(levelID[opt_j]);
					proceedModelGen = 0;
				}
				
				while (proceedModelGen) {
					
					proceedModelGen = 0;
					
					// generates some new models
					new_candidates += generateNewModelsMPI(jobList, resultList, leafID2line, internalID2line,
														   options, firstJobID, levelID, grouping_array,
														   currLevelNumResults, topMatrix, topK, thresHeight, whichApproach,
														   nextAvailableLevel, numSpecies, numEdges, untouchedNodes, opt_method_num, outMsg,
														   start_ps, start_vs, edgeOrder);
					
					// if the number of items = currLevelNumResults
					// this round finishes and new round starts
					if (resultList->size() == currLevelNumResults) {
						
						// print the last 'numMatToCheck' items inside rateMatrixList
						resultList->print_last_N(resultList->size(), ICName);
						
						// move all the items inside resultList to bestHundredRateMatrices
						resultList->move_last_k_items(bestHundredRateMatrices, resultList->size(),1);
						
						// for groupAllCandidates = 1
						if (groupAllCandidates==1) {
							grouping_array.push_back(new_candidates);
							nextAvailableLevel++;
							new_candidates=0;
							
							if (reachRdEnd) {
								// show the best IC score for this round
								bestHundredRateMatrices->rearrangeFirstNItems(1, bestHundredRateMatrices->size());
								currBestIC = bestHundredRateMatrices->ICs[0];
								cout << "Best result for cycle " << currRoundID << " is: ";
								for (k=0; k<(int)bestHundredRateMatrices->matrices[0]->size(); k++) {
									if (k>0)
										cout << ",";
									cout << bestHundredRateMatrices->matrices[0]->at(k);
								}
								cout << " [" << ICName << ": " << longDoublToStr(currBestIC,ANS_DECI) << "]" << endl;
								
								if ((maxNumRounds == -1 || currRoundID < maxNumRounds) && currBestIC < preBestIC) {
									preBestIC = currBestIC;
									currRoundID++;
								} else {
									// the procedure should finish
									procedureEnd = 1;
								}
							}
						}
						
						// update the current level ID
						currLevelID++;
						if (currLevelID >= (int)grouping_array.size()) {
							procedureEnd = 1;
						} else {
							currLevelNumResults = grouping_array[currLevelID] / opt_method_num;
							
							// get the message of number of the candidates going to be checked
							if (currLevelNumResults > 0)
								outMsg.append("[master " + intToStr(++iterationID) + "] Checking " + intToStr(currLevelNumResults) + " rate matices....\n");
							
							// move the result of this level from appendResultList to resultList
							if (groupAllCandidates==1) {
								TempRateMats* tmp = resultList;
								resultList = appendResultList;
								appendResultList = tmp;
								// appendResultList->move_last_k_items(resultList, appendResultList->size(),0);
								appendResultLevelID.clear();
							} else {
								for (k=0; k<appendResultList->size(); k++) {
									if (appendResultLevelID[k] == currLevelID) {
										resultList->insert(appendResultList, k, 1);
										appendResultLevelID[k] = -1;
									}
								}
								// remove the item k inside appendResultlist if appendResultLevelID[k] == -1
								appendResultList->remove_items(&appendResultLevelID, -1);
								// remove the item inside appendResultLevelID if the item == -1
								int new_size = 0;
								for (k=0; k<(int)appendResultLevelID.size(); k++) {
									if (appendResultLevelID[k]!=-1) {
										// keep this item
										if (k > new_size) {
											appendResultLevelID[new_size] = appendResultLevelID[k];
										}
										new_size++;
									}
								}
								if (new_size < (int)appendResultLevelID.size()) {
									if (new_size > 0)
										appendResultLevelID.resize(new_size);
									else
										appendResultLevelID.clear();
								}
							}
							if (resultList->size() > 0)
								proceedModelGen = 1;
							
							
							reachRdEnd = nextReachRdEnd;
							nextReachRdEnd = 0;
						}
						
						// show the message
						if (!procedureEnd) {
							cout << outMsg;
							outMsg = "";
						}
					}
				}
			}
			readyJobExist = 0;
			
		} else {
		
			// =============================================
			// wait for the job result from other processors
			// =============================================
			
			if (availP < totalP) {
				
				// wait for the job result from other processors
				ParameterSet* ps = new ParameterSet(num_w, numEdges, num_chars);
				VariableSet* vs = new VariableSet(num_alpha, num_chars);
				double logL; double ic; int numIter; int sender;
				receiveResult(readyJobID, logL, ic, *ps, *vs, numIter, tagID, sender);
				
				// processor sender-1 becomes available now
				availProcessors[sender-1] = 1;
				availP++;
				
				// save to the job list
				jobList->loglikelihood[readyJobID] = logL;
				jobList->ICs[readyJobID] = ic;
				jobList->PSs[readyJobID] = ps;
				jobList->VSs[readyJobID] = vs;
				jobList->numIters[readyJobID] = numIter;
				jobList->lowPossibleICs[readyJobID] = ic - options->logLDeviation * 2.0;
				readyJobExist = 1;
				
			}
		}
		
	}
	
	// clear memory
	resultList->clear();
	appendResultList->clear();
}

/*
// ============================================================================================================================
// ITERATION FUNCTIONS for RAL
// ============================================================================================================================

// perform RAL iteration for the corresponding rate matrices (for processor ZERO only)
// groupAllCandidates  : For each of the best rate matrices in this iteration,
//                       generate the potential candidates for next iteration
//                       1 - All the potential candidate from the same-level iterations are grouped for the next-level iterations.
//                       2 - All the potential candidate from the same-level iterations are NOT grouped for the next-level iterations.
//                           (i.e. Original design)
void RALIterationNew (TempRateMats* rateMatrixList, int* topMatrix, int numLineTopMat,
		Alignment* alignment, int isReversible, int numSpecies,
		OptimalRateMats* optRateMatrix, int* leafID2line, int* internalID2line,
		double thresHeight, ofstream* fout, TempRateMats* preInterResult,
		TempRateMats* bestHundredRateMatrices, UserOptions* options, int num_chars,
		int saveToInterResult, int tagID, ParameterSet& start_ps, VariableSet& start_vs) {

	int i,j,k;

	int num_w = 6;
	int num_alpha = 1;
	// int numIterations = options->numIterations;
	int topK = options->topK;
	int whichApproach = options->HALMode;
	int groupAllCandidates = options->groupAllCandidates;
	int maxNumCycles = options->numMaxCycles;

	vector<int> theRateMatArray;
	vector< int> currUntouchNodes;
	// int threadId;
	int numEdges = numSpecies*2-2;
	vector<int> untouchzero(numEdges,0);

	// initialize the job allocation array
	int* jobAllocated = NULL;
	int jobAllocatedSize = 0;

	// indicate whether the job is needed to output to the check point file
	vector<int> jobNeedToOutChkptFile;

	// how many nodes can be considered
	int untouchedNodes;

	// grouping array, which indicates the next group of rate matrices to be processed
	vector<int> grouping_array;
	grouping_array.push_back(rateMatrixList->size());

	// job list for optimization
	TempRateMats jobList;

	// an array to check the consistency of the resulting arrangment of rate matrices
	vector<int>* consistentRates = NULL;
	if (whichApproach==4) {
		consistentRates = new vector<int> (numEdges, 0);
	}

	int masterThresID = 0;
	bool needToProcess;
	int numMatToCheck = 0;
	int startItem;
	int insertToSet;
	int numCycles = 1;
	double preBestIC = LARGE_NUMBER;
	double currBestIC;
	double IC_Thres;

	if (groupAllCandidates==1) {
		IC_Thres = IC_Thres_groupAll1;
	} else {
		IC_Thres = IC_Thres_groupAll2;
	}

	string ICName = getICName(options->info_criteria);

	int opt_method_num = (int) options->opt_mode.size();

	while (grouping_array.size() > 0) {

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
				if (interID > -1) {
					needToProcess = false;
					rateMatrixList->PSs[i] = new ParameterSet(num_chars);
					rateMatrixList->PSs[i]->copyFrom(*(preInterResult->PSs[interID]));
					rateMatrixList->VSs[i] = new VariableSet(num_chars);
					rateMatrixList->VSs[i]->copyFrom(*(preInterResult->VSs[interID]));
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
		
		// sending the jobs to other processors
		int nextJob = 0;
		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		int totalProcess = world_size-1;
		int numBusyProcess = 0;
		while (nextJob < jobList.size() && numBusyProcess < totalProcess) {
			if (jobAllocated[nextJob] == 0) {
				jobAllocated[nextJob] = 1;
				numBusyProcess++;
				// sending the job to the processor #<numBusyProcess>
				int currMode = options->opt_mode[nextJob%opt_method_num];
				int currMaxIT = options->opt_maxIT[nextJob%opt_method_num];
				sendInfo(jobList.matrices[nextJob], nextJob, currMode, currMaxIT, 1,
					numEdges, num_w, num_chars, numBusyProcess, tagID, jobList.parentPSs[nextJob], jobList.parentVSs[nextJob]);
			}
			nextJob++;
		}
		
		while (numBusyProcess > 0) {
			
			// receiving the job results from other processors
			ParameterSet* ps = new ParameterSet(num_w, numEdges, num_chars);
			VariableSet* vs = new VariableSet(num_alpha, num_chars);
			int jobID; double logL; double ic; int numIter; int sender;
			receiveResult(jobID, logL, ic, *ps, *vs, numIter, tagID, sender);
			
			// save to the job list
			jobList.loglikelihood[jobID] = logL;
			jobList.ICs[jobID] = ic;
			jobList.PSs[jobID] = ps;
			jobList.VSs[jobID] = vs;
			jobList.numIters[jobID] = numIter;
			jobList.lowPossibleICs[jobID] = ic - options->logLDeviation * 2.0;
			
			while (nextJob<jobList.size() && jobAllocated[nextJob]==1) {
				nextJob++;
			}
			
			if (nextJob<jobList.size()) {
				jobAllocated[nextJob] = 1;
				// sending the job to the processor #<numBusyProcess>
				int currMode = options->opt_mode[nextJob%opt_method_num];
				int currMaxIT = options->opt_maxIT[nextJob%opt_method_num];
				sendInfo(jobList.matrices[nextJob], nextJob, currMode, currMaxIT, 1,
					numEdges, num_w, num_chars, sender, tagID, jobList.parentPSs[nextJob], jobList.parentVSs[nextJob]);
				nextJob++;
			} else {
				numBusyProcess--;
			}
		}
		
		// consolidate all the answers
		for (i=0; i<jobList.size(); i+=opt_method_num) {
			if (jobList.matrices[i]!=NULL) {
				int opt_i = i;
				double opt_IC = jobList.ICs[i];
				for (j=1; j<opt_method_num; j++) {
					if (jobList.ICs[i+j] < opt_IC) {
						opt_IC = jobList.ICs[i+j];
						opt_i = i+j;
					}
				}
				k = i/opt_method_num;

				// all numbers are rounded up, so that
				// it is consistent with the intermediate results loaded from chkpt file

				rateMatrixList->PSs[startItem + k] = new ParameterSet(num_chars);;
				rateMatrixList->PSs[startItem + k]->copyFrom(*(jobList.PSs[opt_i]));
				rateMatrixList->PSs[startItem + k]->roundContent();

				rateMatrixList->VSs[startItem + k] = new VariableSet(num_chars);;
				rateMatrixList->VSs[startItem + k]->copyFrom(*(jobList.VSs[opt_i]));
				rateMatrixList->VSs[startItem + k]->roundContent();

				rateMatrixList->loglikelihood[startItem + k] = roundNumber(jobList.loglikelihood[opt_i], ANS_DECI);
				rateMatrixList->ICs[startItem + k] = roundNumber(opt_IC, ANS_DECI);
				rateMatrixList->numIters[startItem + k] = jobList.numIters[opt_i];
				rateMatrixList->lowPossibleICs[startItem + k] = roundNumber(jobList.lowPossibleICs[opt_i], ANS_DECI);

				// for startOptState = 2 (i.e. optimization starts from the state using one model)
				// if this is the state where only one rate group is involved,
				// then set the starting parameters and starting variables
				if (options->startOptState==2 && rateMatrixList->numRates[startItem + k]==1) {
					start_ps.copyFrom(*rateMatrixList->PSs[startItem + k]);
					start_vs.copyFrom(*rateMatrixList->VSs[startItem + k]);
				}

			}
		}
		// output the intermediate results
		for (i=startItem; i<rateMatrixList->size(); i++) {
			if (jobNeedToOutChkptFile[i-startItem]==1) {

				// save the intermediate result to the array if necessary
				if (saveToInterResult) {
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
					rateMatrixList->PSs[i]->showContent(fout);
					rateMatrixList->VSs[i]->showContent(fout);
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
				if (i < startItem+topK)
					jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], *(rateMatrixList->PSs[i]), *(rateMatrixList->VSs[i]), rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
			} else {
				if (i < startItem+topK && rateMatrixList->ICs[i]<=rateMatrixList->parentICs[i])
					jobList.insert(*(rateMatrixList->matrices[i]), *(rateMatrixList->untouches[i]), rateMatrixList->loglikelihood[i], rateMatrixList->parentICs[i], *(rateMatrixList->PSs[i]), *(rateMatrixList->VSs[i]), rateMatrixList->ICs[i], -1, insertToSet, rateMatrixList->lowPossibleICs[i]);
			}
		}

		// move the last 'numMatToCheck' items inside rateMatrixList to bestHundredRateMatrices
		rateMatrixList->move_last_k_items(bestHundredRateMatrices, numMatToCheck);

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

				if ((maxNumCycles == -1 || numCycles < maxNumCycles) && currBestIC < preBestIC) {
					preBestIC = currBestIC;
					numCycles++;

					// reset all the untouchedNodes for those edges without the consistent rate matrix
					untouchedNodes = jobList.reset_untouches(consistentRates);
				}
			}
		}


		int new_candidates = 0;
		for (i=0; i<jobList.size(); i++) {

			if (jobList.ICs[i] <= jobList.parentICs[i] * IC_Thres) {
				// the IC is small enough
				int pre_new_candidates = new_candidates;
				if ((whichApproach==1 || whichApproach==3) && jobList.numRates[i] < numEdges) {
					// generate all the potential rate matrices into the rateMatrixList (using bottom-up approach)
					new_candidates += generateRateMat(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i]);
				}
				if ((whichApproach==2 || whichApproach==3) && jobList.numRates[i] > 1) {
					// generate all the potential rate matrices into the rateMatrixList (using top-down approach)
					new_candidates += generateRateMatTD(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, internalID2line, jobList.untouches[i], jobList.numRates[i], thresHeight);
				}
				if (whichApproach==4 && untouchedNodes>0) {
					// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 1 approach)
					new_candidates += generateRateMatNewBU1(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], jobList.lowPossibleICs[i]);
				}
				if (whichApproach==6 && untouchedNodes>0) {
					// generate all the potential rate matrices into the rateMatrixList (using bute-force approach)
					new_candidates += generateRateMatButeForce(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, options->maxNumRatMat, jobList.lowPossibleICs[i]);
				}
				if ((whichApproach==7) && jobList.numRates[i] < numEdges) {
					// generate all the potential rate matrices into the rateMatrixList (using John's bottom-up approach)
					new_candidates += generateRateMatJohn(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i]);
				}
				if (whichApproach==8 && untouchedNodes>0) {
					// generate all the potential rate matrices into the rateMatrixList (using new bottom-up version 2 approach)
					new_candidates += generateRateMatNewBU2(rateMatrixList, jobList.matrices[i], jobList.PSs[i], jobList.VSs[i], jobList.loglikelihood[i], jobList.parentICs[i], jobList.ICs[i], topMatrix, leafID2line, internalID2line, jobList.untouches[i], numSpecies, jobList.numRates[i], jobList.lowPossibleICs[i]);
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
}
*/

// The optimization process run by each CPU thread
void optimThread(int* topMatrix, int numLineTopMat, Alignment* alignment, int isReversible,
		int numSpecies, UserOptions* userOptions, int num_chars, int tagID) {


	int numIter;
	int currMode;
	int currMaxIT;
	int num_rateMat;
	int num_w = 6;
	int num_edges = 2 * numSpecies - 2;
	int num_alpha = 1;
	int jobID;
	double curr_loglikelihood;
	double curr_IC;
	int isPsVsUpdated;
	// int i;
	int numIterations = userOptions->numIterations;
	Optim* curr_opt = new Optim(alignment, topMatrix, isReversible, num_chars);
	// string ICName = getICName(userOptions->info_criteria);
	VariableSet vs(num_alpha, num_chars);
	ParameterSet ps(num_w, num_edges, num_chars);
	VariableSet start_vs(num_alpha, num_chars);
	ParameterSet start_ps(num_w, num_edges, num_chars);
	start_vs.resetAllVariables(*alignment);
	start_ps.reset();

	int proceed = 1;
	vector<int> rateMatrix(num_edges);

	while (true) {
		
		// receive information from processor 0
		receiveInfo(rateMatrix, jobID, currMode, currMaxIT, num_rateMat, proceed, 
					num_edges, num_w, num_chars, tagID, &ps, &vs, isPsVsUpdated);
		
		if (!isPsVsUpdated) {
			// set vs and ps the same as the starting ps and the starting vs
			vs.copyFrom(start_vs);
			ps.copyFrom(start_ps);
		} else {
			// update the starting ps and the starting vs
			start_vs.copyFrom(vs);
			start_ps.copyFrom(ps);
		}
					
		if (proceed == 0)
			break;
		
		// update the rate matrix in the parameter
		ps.updateRateMat(rateMatrix, num_rateMat);
		curr_opt->setVarParam(&vs, &ps);

		// get the number of changes of rate matrices
		int numRateChanges = numChanges(topMatrix, &rateMatrix, numSpecies-2, 0);

		// optimizing all parameters for this rate matrices
		numIter = curr_opt->optimizeAllParam(numIterations, curr_loglikelihood, curr_IC, NULL, 0, 0, currMaxIT,
				currMode, userOptions->info_criteria, numRateChanges);

		// send the result to processor 0
		sendResult(jobID, curr_loglikelihood, curr_IC, *curr_opt->ps, *curr_opt->vs, numIter, tagID);
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
	for (i=0; i<interResults->size(); i++) {
		logL = interResults->loglikelihood[i];
		numR = interResults->numRates[i];
		numRateChanges = numChanges(topMatrix, interResults->matrices[i], num_species-2, 0);
		interResults->ICs[i] = getIC(logL, numR, num_sites, num_species, ICType, numRateChanges);
		interResults->lowPossibleICs[i] = interResults->ICs[i] - 2.0 * logLDeviation;
	}
}

// Load previous intermediate results from the log file
TempRateMats* loadPreInterResult(char* logFile, int num_w, int num_edges, int num_alpha, int ICType, 
		int num_sites, int num_species, double logLDeviation, int num_chars, int* topMatrix) {
	ifstream fin;
	fin.open(logFile);
	string aline;
	vector<string> token;
	size_t pos;

	vector<int> matrix(num_edges);
	double loglikelihood = 0.0; double IC = 0.0;
	double lowest_IC;
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

	cout << "loading the log file: " << logFile << endl << flush;

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
			if (!getDoubleValues(token, ICName, &IC, 1)) {
				// The IC name does not match, then the corresponding IC and the lowest possible IC needs to recalculate
				need_recalculate_IC = true;
				// continue;
			} else {
				need_recalculate_IC = false;
			}

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
					IC = getIC(loglikelihood, numRates, num_sites, num_species, ICType, numRateChanges);
				}

				// calculate the lowest IC
				lowest_IC = IC - 2.0 * logLDeviation;

			} else {
				continue;
			}

			// load to the interResults array
			interResults->insertwMap(matrix, loglikelihood, IC, -1, PS, VS, untouch, numIter, lowest_IC); 
			numResults++;
		}
	}
	return interResults;
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
	curr_opt->optimizeAllParam(numIterations, loglikelihood, bicValue, NULL, 0, 0, maxIT, mode, infoCriteria, numRateChanges);

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
	IC = getIC(logL, numRateMatrix, alignment.numSites, numSpecies, df, ICType);

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
	int i,j,k;
	int numRateChanges;

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
				allICValues[i*TOTAL_NUM_IC+j] = getIC(allResults->loglikelihood[i], allResults->numRates[i], 
						numSites, numSpecies, j+1, numRateChanges);
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
