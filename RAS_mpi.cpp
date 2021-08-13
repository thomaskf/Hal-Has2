/*
 *
 * RAS_mpi.cpp
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

#include "RAS_mpi.h"

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
		UserOptions* userOptions, int num_chars, int isGTRUpilson, int world_size, int world_rank) {

	// =================================
	// FOR ALL PROCESSORS
	// =================================

	// load the topology
	int* topMatrix;
	int numLineTopMat;
	int numEdges;
	vector<string>* leafList;
	vector<int> rateMatrix;
	vector<double> edgeLens;
	genTopMatrix(topFile, &topMatrix, &numLineTopMat, &leafList, &rateMatrix, &edgeLens);
	map<string,int>* leafPosMap = genReverseMap(leafList);
	numEdges = numLineTopMat * 2;
	int isReversible = 1;

	// load the alignment file
	Alignment alignment;
	if (userOptions->inputFormat==1)
		alignment.readFASTA(alignFile, leafPosMap);
	else
		alignment.readFile(alignFile, leafPosMap);
	if (userOptions->gapHandling == 1) {
		// 1 - ignore all columns with gaps
		alignment.keepSitesWithChars(0); // remove the columns with invalid characters and gaps
	} else {
		alignment.keepSitesWithChars(1); // only remove the columns with invalid characters
	}
	int* siteOrder = alignment.getSortedSitesID(); // get the order of the sites according to the topology
	alignment.computeNumUniqueSites(siteOrder); // compute the number of unique columns
	alignment.keepUniqueSites(siteOrder); // keep only the unique columns
	// alignment.reorderUniqueSites(); // for efficient calculation of log-likelihood
	alignment.buildChangedNodes(topMatrix, numLineTopMat);
	
	int tagID = DEFAULT_RAS_TAGID;

	// =================================
	// FOR OTHER PROCESSORS EXCEPT 0
	// =================================

	if (world_rank > 0) {

		// The optimization process run by each CPU thread
		optimThread(topMatrix, numLineTopMat, &alignment, isReversible,
					alignment.numSeqs, userOptions, num_chars, tagID, randomInit, isGTRUpilson, world_rank);

		return;
	}

	// =================================
	// FOR ZERO PROCESSOR ONLY
	// =================================

	// Load previous intermediate results from the log file if necessary
	PreHASLogResult preHASLogResult(num_chars);
	if (preLogFile != NULL) {
		preHASLogResult.readLogFile(preLogFile, numEdges);
	}

	int frModelType;
	if (isGTRUpilson)
		frModelType = 5;
	else
		frModelType = 1;
	int toModelType = 5;
	int frNumRateCat = minNumRateCat;
	int toNumRateCat = maxNumRateCat;
	int numUniqModelType = toModelType - frModelType + 1;
	int numUniqRateCat = toNumRateCat - frNumRateCat + 1;
	int opt_method_num = userOptions->opt_mode.size();
	int modelType, numRateCat, optMethod, optMode, optMaxIT;
	int i, s;

	// get the rate matrix arrangement
	int numRateGrp = 0;
	vector<int>* rateMatArray;
	if (rateGrpFile != NULL && strlen(rateGrpFile)>0) {
		rateMatArray = collectRateGrp(rateGrpFile, topMatrix, leafList, numRateGrp);
	} else {
		rateMatArray = new vector<int>(numEdges,1);
		numRateGrp=1;
	}

	// initialize the jobAllocated array
	int* jobAllocated  = new int[numUniqModelType * numUniqRateCat * opt_method_num];
	memset(jobAllocated, 0, sizeof(int)*numUniqModelType*numUniqRateCat*opt_method_num);

	// the HAS summary tables
	double* summaryHASLogL = new double[numUniqModelType * numUniqRateCat * opt_method_num];
	double* summaryHASIC = new double[numUniqModelType * numUniqRateCat * opt_method_num];
	int* summaryHASDf = new int[numUniqModelType * numUniqRateCat * opt_method_num];

	// record the best IC value so far
	// and the corresponding parameters and variables
	double bestIC = 9.99e30;
	AllParameterSet bestPs;
	VariableSet bestVs(num_chars);
	string bestModel = "";
	int rec_no;

	// check whether some results exist in the previous log file
	for (numRateCat=frNumRateCat; numRateCat<=toNumRateCat; numRateCat+=userOptions->stepRateCat) {
		for (modelType=frModelType; modelType<=toModelType; modelType++) {
			s = (numRateCat-frNumRateCat)*numUniqModelType*opt_method_num + (modelType-frModelType)*opt_method_num;
			// s=((modelType-frModelType)*numUniqRateCat+(numRateCat-frNumRateCat))*opt_method_num;
			for (optMethod=0; optMethod<opt_method_num; optMethod++) {
				optMaxIT=userOptions->opt_maxIT[optMethod];
				optMode=userOptions->opt_mode[optMethod];
				rec_no = preHASLogResult.getRecord(modelType, numRateCat, optMode, optMaxIT);
				if (rec_no >= 0) {
					summaryHASLogL[s+optMethod] = preHASLogResult.loglikelihood[rec_no];
					summaryHASDf[s+optMethod] = preHASLogResult.df[rec_no];
					if (preHASLogResult.ICTypes[rec_no]!=ICType) {
						summaryHASIC[s+optMethod] = getIC(preHASLogResult.loglikelihood[rec_no], numRateGrp, alignment.numSites,
								alignment.numSeqs, numRateCat, modelType, preHASLogResult.df[rec_no],
								isGTRUpilson, ICType);
					} else {
						summaryHASIC[s+optMethod] = preHASLogResult.ICs[rec_no];
					}
					if (preHASLogResult.ICs[rec_no] < bestIC) {
						bestPs.copyFrom(*(preHASLogResult.PSs[rec_no]));
						bestVs.copyFrom(*(preHASLogResult.VSs[rec_no]));
						bestIC = preHASLogResult.ICs[rec_no];

						bestModel = "";
						if (isGTRUpilson==0) {
							bestModel.append("HAL-HAS model type : " + intToStr(modelType) + "\n");
						} else {
							bestModel.append("Upsilon model\n");
						}
						bestModel.append("Number of site categories (variable) : " + intToStr(numRateCat) + "\n");
						bestModel.append("Log-likelihood : " + longDoublToStr(preHASLogResult.loglikelihood[rec_no],ANS_DECI) + "\n");
						bestModel.append("Degree of freedom : " + intToStr(preHASLogResult.df[rec_no]) + "\n");
						bestModel.append(getICName(ICType) + " : " + longDoublToStr(summaryHASIC[s+optMethod],ANS_DECI) + "\n");

					}
					jobAllocated[s+optMethod] = 1;
				}
			}
		}
	}

	// the output file for checkpoint
	string outChkptFile;
	if (isGTRUpilson)
		outChkptFile = prefixOut + ".Upsilon.chkpt.txt";
	else
		outChkptFile = prefixOut + ".HAS.chkpt.txt";
	ofstream foutChkpt;
	if (preLogFile != NULL) {
		if (outChkptFile != string(preLogFile)) {
			filecopy((char*) outChkptFile.c_str(), preLogFile);
		}
		foutChkpt.open(outChkptFile.c_str(), std::ofstream::out | std::ofstream::app); // append to the log file
	} else {
		foutChkpt.open(outChkptFile.c_str());
	}
	
	RASIter(&alignment, topMatrix, numLineTopMat, numIterations,
			frNumRateCat, toNumRateCat, frModelType, toModelType,
			rateMatArray, numRateGrp, randomInit, isGTRUpilson,
			jobAllocated, &bestIC, &bestPs, &bestVs, &bestModel, &foutChkpt,
			summaryHASLogL, summaryHASIC, summaryHASDf, ICType,
			userOptions, num_chars, world_size, tagID);

	foutChkpt.close();

	// terminate all the processes
	terminateAllProcsHAS(numEdges, tagID);

	string ICTerm = getICName(ICType);

	// the output file for result
	string outResultFile;
	if (isGTRUpilson)
		outResultFile = prefixOut + ".Upsilon."+ICTerm+".result.txt";
	else
		outResultFile = prefixOut + ".HAS."+ICTerm+".result.txt";
	
	ofstream foutResult;
	foutResult.open(outResultFile.c_str());
	// list out the summary
	foutResult << "=====================================================================================" << endl;
	foutResult << "Summary of the result" << endl;
	if (isGTRUpilson) {
		foutResult << "Upsilon model" << endl;
	} else {
		foutResult << "HAL-HAS model\t";
	}
	foutResult << "# of site categories\tLog-L\t" << ICTerm << "\tdegree of freedom" << endl;
	for (numRateCat=frNumRateCat; numRateCat<=toNumRateCat; numRateCat+=userOptions->stepRateCat) {
		for (modelType=frModelType; modelType<=toModelType; modelType++) {
			// s=((modelType-frModelType)*numUniqRateCat+(numRateCat-frNumRateCat))*opt_method_num;
			s = (numRateCat-frNumRateCat)*numUniqModelType*opt_method_num + (modelType-frModelType)*opt_method_num;
			double currBestIC = summaryHASIC[s];
			int currBestOptMethod = 0;
			for (i=1; i<opt_method_num; i++) {
				if (summaryHASIC[s+i]<currBestIC) {
					currBestIC = summaryHASIC[s+i];
					currBestOptMethod = i;
				}
			}
			if (!isGTRUpilson)
				foutResult << modelType << "\t";
			foutResult << numRateCat << "\t";
			foutResult << longDoublToStr(summaryHASLogL[s+currBestOptMethod],ANS_DECI) << "\t";
			foutResult << longDoublToStr(summaryHASIC[s+currBestOptMethod],ANS_DECI) << "\t";
			foutResult << summaryHASDf[s+currBestOptMethod] << endl;
		}
	}
	foutResult << endl;
	// list out the parameters of the best HAL-HAS model
	foutResult << "=====================================================================================" << endl;
	if (isGTRUpilson)
		foutResult << "The resulting parameters:" << endl;
	else
		foutResult << "The best HAL-HAS model and the corresponding parameters:" << endl;
	foutResult << endl;
	foutResult << bestModel << endl;

	// show the tree
	bestPs.updateContent(1);
	for (i=0; i<bestPs.numRateCat; i++) {
		foutResult << "[Tree of site category " << i+1 << "]" << endl;
		outTreeWithEdgeLen(topMatrix, leafList, bestPs.ps[i]->t, foutResult);
		foutResult << endl;
	}

	// show the content of all parameter sets and variable sets
	// bestPs.updateContent(2);
	string currResult = "";
	bestPs.showContent(currResult, topMatrix, leafList);
	bestVs.showContent(currResult, 1);
	foutResult << currResult << endl;
	foutResult.close();

	// output finish message
	if (outResultFile.length() > 0)
		cout << "Result file: " << outResultFile << endl;

	if (isGTRUpilson)
		cout << "HAS for Upsilon Version " << VERSION << " finishes" << endl;
	else
		cout << "HAS Version " << VERSION << " finishes" << endl;
}


// The optimization process run by each CPU thread
void optimThread(int* topMatrix, int numLineTopMat, Alignment* alignment, int isReversible,
		int numSpecies, UserOptions* userOptions, int num_chars, int tagID, int randomInit,
		int isGTRUpilson, int threadID) {
			
	int jobID;
	int currMode;
	int currMaxIT;
	int num_rateMat;
	int proceed;
	int numRateCat;
	int modelType;
	int num_w = 6;
	int num_edges = 2 * numSpecies - 2;
	double loglikelihood;
	double IC;
	int df;
	vector<int> rateMatrix(num_edges);
	int numIterations = userOptions->numIterations; // maximum number of iterations allowed
	int numIter; // number of iterations really used for optimisation
	Optim* opt;

	while (true) {

		// receive information from processor 0
		receiveInfoHAS(rateMatrix, jobID, currMode, currMaxIT, num_rateMat, proceed,
			numRateCat, modelType, num_edges, tagID);

		if (proceed != 1) {
			break;
		}

		// parameters
		AllParameterSet *all_ps = new AllParameterSet(numRateCat, num_chars);

		if (randomInit==1) {
			initialRandSeed();
			all_ps->randomInit(num_w, num_edges);
		} else {
			all_ps->initialize(num_w, num_edges);
		}

		// set the rate matrixcomputeEigenMatrix
		all_ps->updateRateMat(rateMatrix, num_rateMat);

		// Variables
		VariableSet *variables = new VariableSet(numRateCat, num_chars);
		if (isGTRUpilson) {
			variables->setGTRUpilson();
		}

		if (randomInit==1) {
			variables->randomInitAllVariables();
		} else {
			variables->resetAllVariables(*alignment);
		}

		// compute the corresponding loglikelihood and IC
		all_ps->computeAllEigenValues_t();
		all_ps->computeAllCondProb(isReversible);

		//================
		// optimization
		//================
		opt = new Optim(variables, all_ps, alignment, topMatrix, isReversible, modelType, num_chars);
		numIter = opt->optimizeAllParamRAS(numIterations, num_rateMat, loglikelihood, IC, df, isGTRUpilson, 
											threadID, NULL, currMaxIT, currMode, userOptions->info_criteria);

		// send the result to processor 0
		sendResultHAS(jobID, loglikelihood, IC, df, *all_ps, *variables, numIter, tagID);
		
		// clear memory
		delete all_ps;
		delete variables;
		delete opt;

	}
}

			
void RASIter(Alignment* alignment, int* topMatrix, int numLineTopMat, int numIterations,
		int frNumRateCat, int toNumRateCat, int frModelType, int toModelType,
		vector<int>* rateMatArray, int numRateGrp, int randomInit, int isGTRUpilson, int* jobAllocated,
		double* bestIC, AllParameterSet* bestPs, VariableSet* bestVs, string* bestModel, ofstream* foutChkpt,
		double* summaryHASLogL, double* summaryHASIC, int* summaryHASDf, int ICType,
		UserOptions* userOptions, int num_chars, int totProcs, int tagID) {

	double loglikelihood,IC;
	int df;

	// int isReversible=1; int num_w = 6;
	int numEdges = 2 * numLineTopMat;

	int numRateCat, modelType;

	AllParameterSet* all_ps = new AllParameterSet(toNumRateCat, num_chars, numEdges);
	VariableSet* variables = new VariableSet(toNumRateCat, num_chars);
	int numIter;
	int currMaxIT;
	int currMode;
	int optMethod;

	int numModels = toModelType - frModelType + 1;
	int numUniqRate = toNumRateCat - frNumRateCat + 1;
	int opt_method_num = userOptions->opt_mode.size();
	string ICTerm = getICName(ICType);
	int availProcs = totProcs - 1;
	int proceed = 1;
	int jobID;
	int nextJob=0;
	int totJobs = numModels*numUniqRate*opt_method_num;
	int sender;
	int x;

	//==========================================
	// send jobs to the processors
	//==========================================
	for (numRateCat=frNumRateCat; numRateCat<=toNumRateCat && availProcs>0; numRateCat+=userOptions->stepRateCat) {
		for (modelType=frModelType; modelType<=toModelType && availProcs>0; modelType++) {
			for (optMethod=0; optMethod<opt_method_num && availProcs>0; optMethod++) {
				
				// jobID = (modelType-frModelType)*numUniqRate*opt_method_num + (numRateCat-frNumRateCat)*opt_method_num + optMethod;
				jobID = (numRateCat-frNumRateCat)*numModels*opt_method_num + (modelType-frModelType)*opt_method_num + optMethod;

				// cout << "model:" << modelType << "; # of site categories:" << numRateCat << "; jobID: " << jobID << "; jobAllocated[" << jobID << "]" << jobAllocated[jobID] << endl;
				
				nextJob = jobID+1;
				// check which job has not been processed
				if (jobAllocated[jobID]==0) {
					jobAllocated[jobID]=1;
					proceed = true;
				} else {
					continue;
				}

				currMaxIT=userOptions->opt_maxIT[optMethod];
				currMode=userOptions->opt_mode[optMethod];
					
				// send the job to the processor #availProcs
				sendInfoHAS(rateMatArray, jobID, currMode, currMaxIT, numRateGrp, proceed, 
						 numRateCat, modelType, numEdges, availProcs--, tagID);
			}
		}
	}
	
	while (availProcs < totProcs-1) {
		
		// some processors are still busy
		
		// receiving the job results from other processors
		receiveResultHAS(jobID, loglikelihood, IC, df, *all_ps, *variables,
					  numIter, tagID, sender);

		optMethod = (jobID % opt_method_num);
		x = (jobID / opt_method_num);
		modelType= (x % numModels) + frModelType;
		numRateCat = (x / numModels) + frNumRateCat;
		currMaxIT=userOptions->opt_maxIT[optMethod];
		currMode=userOptions->opt_mode[optMethod];

		// output system message
		if (isGTRUpilson==0)
			cout << "[threadID:" << sender << "] model:" << modelType << "; # of site categories:" << numRateCat << "; maxIT:" << currMaxIT << "; mode:" << currMode;
		else
			cout << "[threadID:" << sender << "] Upsilon model; # of site categories:" << numRateCat << "; maxIT:" << currMaxIT << "; mode:" << currMode;
		if (numIter>0) {
			cout << " # of iterations:" << numIter << ";";
		}
		cout << " loglikelihood:" << longDoublToStr(loglikelihood,ANS_DECI) << "; " << ICTerm << ":" << longDoublToStr(IC,ANS_DECI) << "; df:" << df  << endl;

		// output to the checkpoint file
		*foutChkpt << "=====================================================================================" << endl;
		*foutChkpt << "HAL-HAS model type: " << modelType << endl;
		*foutChkpt << "Number of site categories: " << numRateCat << endl;
		*foutChkpt << "Mode of optimization method: " << currMode << endl;
		*foutChkpt << "MaxIT of optimization method: " << currMaxIT << endl;
		*foutChkpt << "Log-likelihood: " << longDoublToStr(loglikelihood,ANS_DECI) << endl;
		*foutChkpt << "Degree of freedom: " << df << endl;
		*foutChkpt << ICTerm << " : " << longDoublToStr(IC,ANS_DECI) << endl;
		// show the content of all parameter sets and variable sets
		all_ps->showContent(foutChkpt);
		variables->showContent(foutChkpt);
		*foutChkpt << flush;

		// store in the summary tables
		summaryHASLogL[jobID] = loglikelihood;
		summaryHASIC[jobID] = IC;
		summaryHASDf[jobID] = df;

		// check whether it is the optimal result
		if (IC < *bestIC) {
			*bestIC = IC;
			bestPs->copyFrom(*all_ps);
			bestVs->copyFrom(*variables);
			*bestModel = "";
			if (isGTRUpilson==0) {
				bestModel->append("HAL-HAS model type : " + intToStr(modelType) + "\n");
			} else {
				bestModel->append("Upsilon model\n");
			}
			bestModel->append("Number of site categories (variable) : " + intToStr(numRateCat) + "\n");
			bestModel->append("Log-likelihood : " + longDoublToStr(loglikelihood,ANS_DECI) + "\n");
			bestModel->append("Degree of freedom : " + intToStr(df) + "\n");
			bestModel->append(ICTerm + " : " + longDoublToStr(IC,ANS_DECI) + "\n");
		}
		
		// get the next job which has not been processed
		while (nextJob < totJobs && jobAllocated[nextJob]==1)
			nextJob++;
		
		if (nextJob < totJobs) {
			jobAllocated[nextJob] = 1;
			optMethod = (nextJob % opt_method_num);
			x = (nextJob / opt_method_num);
			modelType= (x % numModels) + frModelType;
			numRateCat = (x / numModels) + frNumRateCat;
			// numRateCat= (x % numUniqRate) + frNumRateCat;
			// modelType = (x / numUniqRate) + frModelType;
			currMaxIT=userOptions->opt_maxIT[optMethod];
			currMode=userOptions->opt_mode[optMethod];
			// send the job to the processor #sender
			sendInfoHAS(rateMatArray, nextJob, currMode, currMaxIT, numRateGrp, proceed, 
					 numRateCat, modelType, numEdges, sender, tagID);
			nextJob++;
		} else {
			// all jobs are allocated
			availProcs++;
		}
		
	}
}

// compute the number of substitutions for each edge
//
// topFile      : the topology file
// paramFile    : the prefix of the parameter files
//                the parameter files: "paramFile"1, "paramFile"2, ...
// numRateCat   : the total number of parameter files
void computeSubsEachEdge(char* topFile, char* paramFile, int numRateCat) {

	int num_w = 6;
	int num_chars = 4;

	// load the topology
	int* topMatrix;
	int numLineTopMat, numSpecies, numEdges;
	vector<string>* leafList;
	genTopMatrix(topFile, &topMatrix, &numLineTopMat, &leafList, NULL, NULL);
	numSpecies=(int)leafList->size();
	if (numSpecies < 2) {
		cerr << "Error! The number of species is too small (i.e. less than 2) for this algorithm!" << endl;
		exit(1);
	}
	numEdges=numSpecies*2-2;

	/*
	// print out the topMatrix
	cout << "topMatrix:" << endl;
	for (int i=0; i<numLineTopMat; i++) {
		for (int j=0; j<2; j++) {
			cout << topMatrix[i*2 + j] << ",";
		}
		cout << endl;
	}

	// print out the leaftList
	cout << "leafList:" << endl;
	for (int i=0; i<leafList->size(); i++) {
		cout << leafList->at(i) << ",";
	}
	cout << endl;
	 */

	// load all the parameter files
	AllParameterSet all_ps(numRateCat, num_chars);
	all_ps.readParamFile(paramFile, num_w, numEdges);

	// compute the S and R matrix
	double* s = matrix(4, 4);
	double* pim = matrix(4,4); // pi matrix
	double* r = matrix(4, 4);

	// compute the substituation for each edge and print out the corresponding newick tree
	vector<double> dist;
	dist.resize(numEdges);

	for (int rateCat=0; rateCat<numRateCat; rateCat++) {
		ParameterSet* ps = all_ps.ps[rateCat];
		cout << "Category " << rateCat+1 << endl;
		for (int edge=0; edge<numEdges; edge++) {
			resetMat(s, 4, 4);
			resetMat(pim, 4, 4);
			resetMat(r, 4, 4);
			dist[edge] = 0.0;
			double* curr_w = &((ps->w)[edge*num_w]);
			double* curr_pi = &((ps->pi)[edge*4]);
			computeSMatrix(s, curr_w, curr_pi);
			setDiag(pim, curr_pi, 4);
			// multiplyQuick(r, s, pim, 4, 4, 4);
			multiplyQuick44(r, s, pim);
			for (int k=0; k<4; k++) {
				dist[edge] += curr_pi[k] * r[4*k+k];
			}
			dist[edge] = (-dist[edge]) * (ps->t)[edge];
		}

		/*
		// show the distances
		for (int edge=0; edge<numEdges; edge++) {
			cout << dist[edge] << ",";
		}
		cout << endl;
		 */

		// print out the corresponding newick tree
		outTreeWithEdgeLen(topMatrix, leafList, &dist);
	}
}


PreHASLogResult::PreHASLogResult(int num_chars) {
	// Constructor
	this->num_chars = num_chars;
}

// read the HAS log file
void PreHASLogResult::readLogFile(char* logFile, int numEdges) {
	ifstream fin;
	fin.open(logFile);
	string aline;
	vector<string> token;
	int model, cat, deg, mode, maxIT, icType;
	double logL, ic;
	// icType: type of Information Criteria (1 - AIC; 2 - Adjusted BIC; 3 - BIC; 4 - CAIC)
	int num_w = 6;
	int step = 0;
	int paramSet = 0;
	AllParameterSet* all_ps = NULL;
	VariableSet* vs = NULL;
	int i;
	while(getline(fin,aline)) {
		if (aline.length() > 0 && aline[0]!='=') {
			// not begin with "="
			tokenizer(aline, ":,", &token);
			// HAL-HAS model type
			trim(token[0]);
			if (token[0] == "HAL-HAS model type") {
				model = atoi(token[1].c_str());
				step=1;
				// reset some variables
				paramSet = 0;
			} else if (token[0] == "Number of site categories" && step==1) {
				cat = atoi(token[1].c_str());
				step++;
				if (all_ps != NULL) {
					delete (all_ps);
					all_ps = NULL;
				}
				if (vs != NULL) {
					delete (vs);
					vs = NULL;
				}
				if (cat > 0) {
					all_ps = new AllParameterSet(cat, num_chars);
					vs = new VariableSet(cat, num_chars);
				}
			} else if (token[0] == "Mode of optimization method" && step==2) {
				mode = atoi(token[1].c_str());
				step++;
			} else if (token[0] == "MaxIT of optimization method" && step==3) {
				maxIT = atoi(token[1].c_str());
				step++;
			} else if (token[0] == "Log-likelihood" && step==4) {
				logL = atof(token[1].c_str());
				step++;
			} else if (token[0] == "Degree of freedom" && step==5) {
				deg = atoi(token[1].c_str());
				step++;
			} else if (step==6) {
				// icType (1 - AIC; 2 - Adjusted BIC; 3 - BIC; 4 - CAIC)
				int t = getICType(token[0]);
				if (t > 0) {
					ic = atof(token[1].c_str());
					icType = t;
				}
				step++;
			} else if (token[0] == "beta" && step==7 && paramSet==cat && vs!=NULL && token.size()==2) {
				vs->beta = atof(token[1].c_str());
				step++;
			} else if (token[0] == "alpha" && step==8 && vs!=NULL && (int)token.size()==cat+1) {
				for (i=0; i<cat; i++)
					vs->alpha[i] = atof(token[i+1].c_str());
				step++;
			} else if (token[0] == "probXGivenInv" && step==9 && vs!=NULL && (int)token.size()==num_chars+1) {
				for (i=0; i<num_chars; i++)
					vs->probXGivenInv[i] = atof(token[i+1].c_str());
				step++;
			} else if (token[0] == "rootNodeFreq" && step==10 && vs!=NULL && (int)token.size()==num_chars*cat+1) {
				for (i=0; i<num_chars*cat; i++)
					vs->rootNodeFreq[i] = atof(token[i+1].c_str());
				// save the record
				Quad quad;
				quad.modelType=model; quad.numCat=cat;
				quad.optMode=mode;quad.optMaxIT=maxIT;
				index.insert(pair<Quad,int>(quad,modelType.size()));
				modelType.push_back(model);
				numCat.push_back(cat);
				optMode.push_back(mode);
				optMaxIT.push_back(maxIT);
				loglikelihood.push_back(logL);
				df.push_back(deg);
				ICs.push_back(ic);
				ICTypes.push_back(icType);
				PSs.push_back(all_ps);
				VSs.push_back(vs);
				all_ps = NULL;
				vs = NULL;
				step++;
			} else if (aline.length()>10 && aline.substr(1,9)=="parameter" && step==7 && all_ps!=NULL && paramSet<cat) {
				// load one parameter set
				all_ps->ps[paramSet] = new ParameterSet(num_w, numEdges, num_chars);
				int i;
				for (i=0; i<numEdges; i++) {
					if(getline(fin,aline)) {
						tokenizer(aline," \t", &token);
						if (token.size()<11)
							break;
						for (int k=0; k<num_w; k++) {
							all_ps->ps[paramSet]->w[i*num_w+k] = atof(token[k].c_str());
						}
						for (int k=0; k<4; k++) {
							all_ps->ps[paramSet]->pi[i*4+k] = atof(token[num_w+k].c_str());
						}
						all_ps->ps[paramSet]->t[i] = atof(token[num_w+4].c_str());
					} else {
						break;
					}
				}
				if (i!=numEdges)
					continue;
				paramSet++;
			}
		}
	}
}

int PreHASLogResult::getRecord(int modelType, int numCat, int optMode, int optMaxIT) {
	// get the record according to the modelType, numCat, optMode and optMaxIT

	map<Quad,int>::iterator itr;
	Quad quad;
	quad.modelType=modelType; quad.numCat=numCat;
	quad.optMode=optMode; quad.optMaxIT=optMaxIT;
	itr = index.find(quad);
	if (itr != index.end()) {
		return itr->second;
	}
	return -1;
} 
