/*
 *
 * rateMatrixResult.cpp
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

#include "rateMatrixResult.h"

// constructor
TempRateMats::TempRateMats() {
	N = 0;
	checkedItems = 0;
	topN_IC = LARGE_NUMBER;
}

// keep the top N items and keep them in a sorted order
// int N;
// vector<int> sortedNItems;

// get the top N items and arrange them in a sorted order in accending order of IC values
void TempRateMats::createTopN(int the_N) {

	if (checkedItems > 0 && (N==the_N || N==checkedItems)) {
		maintainTopN(the_N);
	} else {
		if (the_N > (int)ICs.size())
			the_N = (int)ICs.size();
		
		N = the_N;
		
		vector<int> selected (ICs.size(), 0);
		
		if ((int)sortedNItems.size() < the_N)
			sortedNItems.resize(the_N);
		if ((int)sortedICs.size() < the_N)
			sortedICs.resize(the_N);
		
		for (int i=0; i<the_N; i++) {
			double opt_IC = LARGE_NUMBER;
			int opt_pos = -1;
			for (int j=0; j<(int)ICs.size(); j++) {
				if (!selected[j] && opt_IC > ICs[j]) {
					opt_IC = ICs[j];
					opt_pos = j;
				}
			}
			if (opt_pos == -1) {
				cerr << "Error! opt_pos shold not be -1" << endl << flush;
				exit(1);
			}
			selected[opt_pos]=1;
			sortedNItems[i] = opt_pos;
			sortedICs[i] = opt_IC;
		}
		checkedItems = (int) ICs.size();
		topN_IC = sortedICs[the_N-1];
	}
}

// since new items have been added, this function is to maintain the top N items
void TempRateMats::maintainTopN(int the_N) {
	if (the_N > (int)ICs.size())
		the_N = (int)ICs.size();
	
	N = the_N;

	for (int i=checkedItems; i<(int)ICs.size(); i++) {
		int inserted = 0;
		int j;
		if (ICs[i] < topN_IC) {
			for (j=0; j<(int)sortedICs.size(); j++) {
				if (ICs[i] < sortedICs[j]) {
					sortedNItems.insert(sortedNItems.begin()+j,i);
					sortedICs.insert(sortedICs.begin()+j,ICs[i]);
					inserted = 1;
					break;
				}
			}
		}
		if (!inserted && (int)sortedICs.size() < the_N) {
			sortedNItems.push_back(i);
			sortedICs.push_back(ICs[i]);
		} else if ((int)sortedICs.size() > the_N) {
			sortedNItems.resize(the_N);
			sortedICs.resize(the_N);
		}
		topN_IC = sortedICs[sortedICs.size()-1];
		
	}
	
	if (the_N != (int)sortedNItems.size()) {
		cerr << "[inside maintainTopN] Error! N != sortedNItems.size()" << endl << flush;
		exit(1);
	}
	checkedItems = (int) ICs.size();
}


// insert an item
void TempRateMats::insert() {
	matrices.push_back(NULL);
	loglikelihood.push_back(0.0);
	ICs.push_back(0.0);
	parentICs.push_back(0.0);
	numRates.push_back(0);
	PSs.push_back(NULL);
	VSs.push_back(NULL);
	untouches.push_back(NULL);
	numIters.push_back(-1);
	lowPossibleICs.push_back(0.0);
	parentPSs.push_back(NULL);
	parentVSs.push_back(NULL);
	modes.push_back(1);
	maxITs.push_back(1);
	selected.push_back(0);
}


// insert an item
void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(0.0);
		ICs.push_back(0.0);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		PSs.push_back(NULL);
		VSs.push_back(NULL);
		parentPSs.push_back(NULL);
		parentVSs.push_back(NULL);
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		if (node!=-1)
			new_untouch->at(node) = 1;
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(0.0);
		modes.push_back(1);
		maxITs.push_back(1);
		selected.push_back(0);
	}
}

// insert an item
void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet,
						  int theMode, int theMaxIT) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(0.0);
		ICs.push_back(0.0);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		PSs.push_back(NULL);
		VSs.push_back(NULL);
		parentPSs.push_back(NULL);
		parentVSs.push_back(NULL);
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		if (node!=-1)
			new_untouch->at(node) = 1;
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(0.0);
		modes.push_back(theMode);
		maxITs.push_back(theMaxIT);
		selected.push_back(0);
	}
}

// insert an item
void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet, ParameterSet& parent_ps, VariableSet& parent_vs) {
	
	insert(matrix, untouch, parentIC, node, insertToSet, parent_ps, parent_vs, 1, 1);
	
}

// insert an item
void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet, ParameterSet& parent_ps, VariableSet& parent_vs, int mode, int maxIT) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(0.0);
		ICs.push_back(0.0);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		PSs.push_back(NULL);
		VSs.push_back(NULL);

		ParameterSet* new_ps = new ParameterSet(parent_ps.num_chars);
		new_ps->copyFrom(parent_ps);
		parentPSs.push_back(new_ps);

		VariableSet* new_vs = new VariableSet(parent_vs.num_chars);
		new_vs->copyFrom(parent_vs);
		parentVSs.push_back(new_vs);

		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		if (node!=-1)
			new_untouch->at(node) = 1;
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(0.0);
		modes.push_back(mode);
		maxITs.push_back(maxIT);
		selected.push_back(0);
	}
}


void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double parentIC, ParameterSet& ps, VariableSet& vs, int insertToSet) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(0.0);
		ICs.push_back(0.0);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		ParameterSet* new_ps = new ParameterSet(ps.num_chars);
		new_ps->copyFrom(ps);
		PSs.push_back(new_ps);
		VariableSet* new_vs = new VariableSet(vs.num_chars);
		new_vs->copyFrom(vs);
		VSs.push_back(new_vs);
		parentPSs.push_back(NULL);
		parentVSs.push_back(NULL);
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(0.0);
		modes.push_back(1);
		maxITs.push_back(1);
		selected.push_back(0);
	}
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC) {
	
	return insert(matrix, untouch, logL, parentIC, ps, vs, IC, node, insertToSet, lowPossibleIC, 1, 1);
	
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, int mode, int maxIT) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(logL);
		ICs.push_back(IC);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		ParameterSet* new_ps = new ParameterSet(ps.num_chars);
		new_ps->copyFrom(ps);
		PSs.push_back(new_ps);
		VariableSet* new_vs = new VariableSet(vs.num_chars);
		new_vs->copyFrom(vs);
		VSs.push_back(new_vs);
		parentPSs.push_back(NULL);
		parentVSs.push_back(NULL);
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		if (node!=-1) {
			new_untouch->at(node) = 1;
		}
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(lowPossibleIC);
		modes.push_back(mode);
		maxITs.push_back(maxIT);
		selected.push_back(0);
	}
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, double IC, int node, int insertToSet, double lowPossibleIC) {
    
    insert(matrix, untouch, logL, parentIC, IC, node, insertToSet, lowPossibleIC, 1, 1);
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, double IC, int node, int insertToSet, double lowPossibleIC, int mode, int maxIT) {
    if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
        vector<int>* new_matrix = new vector<int>;
        new_matrix->assign(matrix.begin(),matrix.end());
        matrices.push_back(new_matrix);
        loglikelihood.push_back(logL);
        ICs.push_back(IC);
        parentICs.push_back(parentIC);
        numRates.push_back(getNumRateGrp(matrix));
        
        PSs.push_back(NULL);
        VSs.push_back(NULL);
        parentPSs.push_back(NULL);
        parentVSs.push_back(NULL);
        
        vector<int>* new_untouch = new vector<int>;
        new_untouch->assign(untouch.begin(), untouch.end());
        if (node!=-1) {
            new_untouch->at(node) = 1;
        }
        untouches.push_back(new_untouch);
        numIters.push_back(-1);
        if (insertToSet==1)
            matrices_set.insert(matrix);
        lowPossibleICs.push_back(lowPossibleIC);
        modes.push_back(mode);
        maxITs.push_back(maxIT);
        selected.push_back(0);
    }
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, ParameterSet& parent_ps, VariableSet& parent_vs) {

	insert(matrix, untouch, logL, parentIC, ps, vs, IC, node, insertToSet, lowPossibleIC, parent_ps, parent_vs, 1, 1);
	
}

void TempRateMats::insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, ParameterSet& parent_ps, VariableSet& parent_vs, int mode, int maxIT) {
	if (insertToSet==0 || matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(logL);
		ICs.push_back(IC);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		
		ParameterSet* new_ps = new ParameterSet(ps.num_chars);
		new_ps->copyFrom(ps);
		PSs.push_back(new_ps);
		
		VariableSet* new_vs = new VariableSet(vs.num_chars);
		new_vs->copyFrom(vs);
		VSs.push_back(new_vs);
		
		ParameterSet* new_ps2 = new ParameterSet(parent_ps.num_chars);
		new_ps2->copyFrom(parent_ps);
		parentPSs.push_back(new_ps2);
		
		VariableSet* new_vs2 = new VariableSet(parent_vs.num_chars);
		new_vs2->copyFrom(parent_vs);
		parentVSs.push_back(new_vs2);
		
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		if (node!=-1) {
			new_untouch->at(node) = 1;
		}
		untouches.push_back(new_untouch);
		numIters.push_back(-1);
		if (insertToSet==1)
			matrices_set.insert(matrix);
		lowPossibleICs.push_back(lowPossibleIC);
		modes.push_back(mode);
		maxITs.push_back(maxIT);
		selected.push_back(0);
	}
}
// map<vector<int>,int,compForIntArray> matrices_map;

// insert an item
void TempRateMats::insertwMap(vector<int>& matrix, double logL, double IC, double parentIC,
		vector<int>& untouch, int numIter, double lowPossibleIC) {

	if (matrices_set.find(matrix)==matrices_set.end()) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign(matrix.begin(),matrix.end());
		int currID = (int) matrices.size();
		matrices.push_back(new_matrix);
		loglikelihood.push_back(logL);
		ICs.push_back(IC);
		parentICs.push_back(parentIC);
		numRates.push_back(getNumRateGrp(matrix));
		PSs.push_back(NULL);
		VSs.push_back(NULL);
		parentPSs.push_back(NULL);
		parentVSs.push_back(NULL);
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign(untouch.begin(), untouch.end());
		untouches.push_back(new_untouch);
		numIters.push_back(numIter);
		matrices_set.insert(matrix);
		matrices_map.insert(pair<vector<int>,int>(matrix,currID));
		lowPossibleICs.push_back(lowPossibleIC);
		modes.push_back(1);
		maxITs.push_back(1);
		selected.push_back(0);
	}
}

// insert an item
void TempRateMats::insertwMap(vector<int>& matrix, double logL, double IC, double parentIC, ParameterSet& ps, VariableSet& vs,
        vector<int>& untouch, int numIter, double lowPossibleIC) {

    if (matrices_set.find(matrix)==matrices_set.end()) {
        vector<int>* new_matrix = new vector<int>;
        new_matrix->assign(matrix.begin(),matrix.end());
        int currID = (int) matrices.size();
        matrices.push_back(new_matrix);
        loglikelihood.push_back(logL);
        ICs.push_back(IC);
        parentICs.push_back(parentIC);
        numRates.push_back(getNumRateGrp(matrix));
        ParameterSet* new_ps = new ParameterSet(ps.num_chars);
        new_ps->copyFrom(ps);
        PSs.push_back(new_ps);
        VariableSet* new_vs = new VariableSet(vs.num_chars);
        new_vs->copyFrom(vs);
        VSs.push_back(new_vs);
        parentPSs.push_back(NULL);
        parentVSs.push_back(NULL);
        vector<int>* new_untouch = new vector<int>;
        new_untouch->assign(untouch.begin(), untouch.end());
        untouches.push_back(new_untouch);
        numIters.push_back(numIter);
        matrices_set.insert(matrix);
        matrices_map.insert(pair<vector<int>,int>(matrix,currID));
        lowPossibleICs.push_back(lowPossibleIC);
        modes.push_back(1);
        maxITs.push_back(1);
        selected.push_back(0);
    }
}

// insert an item
void TempRateMats::insertwMap(vector<int>& matrix, double logL, double IC, vector<int>& untouch, double lowPossibleIC) {

    if (matrices_set.find(matrix)==matrices_set.end()) {
        vector<int>* new_matrix = new vector<int>;
        new_matrix->assign(matrix.begin(),matrix.end());
        int currID = (int) matrices.size();
        matrices.push_back(new_matrix);
        loglikelihood.push_back(logL);
        ICs.push_back(IC);
        parentICs.push_back(0);
        numRates.push_back(getNumRateGrp(matrix));
        PSs.push_back(NULL);
        VSs.push_back(NULL);
        parentPSs.push_back(NULL);
        parentVSs.push_back(NULL);
        vector<int>* new_untouch = new vector<int>;
        new_untouch->assign(untouch.begin(), untouch.end());
        untouches.push_back(new_untouch);
        numIters.push_back(0);
        matrices_set.insert(matrix);
        matrices_map.insert(pair<vector<int>,int>(matrix,currID));
        lowPossibleICs.push_back(lowPossibleIC);
        modes.push_back(1);
        maxITs.push_back(1);
        selected.push_back(0);
    }
}

// find ID by matrix
int TempRateMats::getID(vector<int>& matrix) {
	map<vector<int>,int>::iterator itr;
	itr = matrices_map.find(matrix);
	if (itr != matrices_map.end()) {
		// found
		return itr->second;
	} else {
		// not found
		return -1;
	}
}


// rearrange the first N items inside the matrices
// such that the top K items (among the first N items) are with the least IC values
void TempRateMats::rearrangeFirstNItems(int topK, int N) {

	if (N > (int)ICs.size())
		N = (int)ICs.size();

	if (topK > N)
		topK = N;

	for (int i=0; i<topK; i++) {
		double opt_IC = ICs[i];
		double opt_pos = i;
		for (int j=i+1; j<N; j++) {
			if (opt_IC > ICs[j]) {
				opt_IC = ICs[j];
				opt_pos = j;
			}
		}
		if (opt_pos > i) {
			// swap between position i and position op_pos
			vector<int>* ptr = matrices[i];
			matrices[i] = matrices[opt_pos];
			matrices[opt_pos] = ptr;

			ptr = untouches[i];
			untouches[i] = untouches[opt_pos];
			untouches[opt_pos] = ptr;

			double temp = loglikelihood[i];
			loglikelihood[i] = loglikelihood[opt_pos];
			loglikelihood[opt_pos] = temp;

			temp = ICs[i];
			ICs[i] = ICs[opt_pos];
			ICs[opt_pos] = temp;

			temp = parentICs[i];
			parentICs[i] = parentICs[opt_pos];
			parentICs[opt_pos] = temp;

			int tmp = numRates[i];
			numRates[i] = numRates[opt_pos];
			numRates[opt_pos] = tmp;

			ParameterSet* tps = PSs[i];
			PSs[i] = PSs[opt_pos];
			PSs[opt_pos] = tps;

			tps = parentPSs[i];
			parentPSs[i] = parentPSs[opt_pos];
			parentPSs[opt_pos] = tps;

			VariableSet* tpv = VSs[i];
			VSs[i] = VSs[opt_pos];
			VSs[opt_pos] = tpv;

			tpv = parentVSs[i];
			parentVSs[i] = parentVSs[opt_pos];
			parentVSs[opt_pos] = tpv;

			tmp = numIters[i];
			numIters[i] = numIters[opt_pos];
			numIters[opt_pos] = tmp;

			temp = lowPossibleICs[i];
			lowPossibleICs[i] = lowPossibleICs[opt_pos];
			lowPossibleICs[opt_pos] = temp;
			
			tmp = modes[i];
			modes[i] = modes[opt_pos];
			modes[opt_pos] = tmp;

			tmp = maxITs[i];
			maxITs[i] = maxITs[opt_pos];
			maxITs[opt_pos] = tmp;
			
			tmp = selected[i];
			selected[i] = selected[opt_pos];
			selected[opt_pos] = tmp;
		}
	}
}

// rearrange the last N items inside the matrices
// such that the top K items (among the last N items) are with the least IC values
void TempRateMats::rearrangeLastNItems(int topK, int N) {

	if (N > (int)ICs.size())
		N = (int)ICs.size();

	if (topK > N)
		topK = N;

	int i,j;
	int startItem = (int)ICs.size()- N;
	for (i=startItem; i<topK+startItem; i++) {
		double opt_IC = ICs[i];
		double opt_pos = i;
		for (j=i+1; j<(int)ICs.size(); j++) {
			if (opt_IC > ICs[j]) {
				opt_IC = ICs[j];
				opt_pos = j;
			}
		}
		if (opt_pos > i) {
			// swap between position i and position op_pos
			vector<int>* ptr = matrices[i];
			matrices[i] = matrices[opt_pos];
			matrices[opt_pos] = ptr;

			ptr = untouches[i];
			untouches[i] = untouches[opt_pos];
			untouches[opt_pos] = ptr;

			double temp = loglikelihood[i];
			loglikelihood[i] = loglikelihood[opt_pos];
			loglikelihood[opt_pos] = temp;

			temp = ICs[i];
			ICs[i] = ICs[opt_pos]; 
			ICs[opt_pos] = temp;

			temp = parentICs[i];
			parentICs[i] = parentICs[opt_pos];
			parentICs[opt_pos] = temp;

			int tmp = numRates[i];
			numRates[i] = numRates[opt_pos];
			numRates[opt_pos] = tmp;

			ParameterSet* tps = PSs[i];
			PSs[i] = PSs[opt_pos];
			PSs[opt_pos] = tps;

			tps = parentPSs[i];
			parentPSs[i] = parentPSs[opt_pos];
			parentPSs[opt_pos] = tps;

			VariableSet* tpv = VSs[i];
			VSs[i] = VSs[opt_pos];
			VSs[opt_pos] = tpv;

			tpv = parentVSs[i];
			parentVSs[i] = parentVSs[opt_pos];
			parentVSs[opt_pos] = tpv;

			tmp = numIters[i];
			numIters[i] = numIters[opt_pos];
			numIters[opt_pos] = tmp;

			temp = lowPossibleICs[i];
			lowPossibleICs[i] = lowPossibleICs[opt_pos]; 
			lowPossibleICs[opt_pos] = temp;

			tmp = modes[i];
			modes[i] = modes[opt_pos];
			modes[opt_pos] = tmp;
			
			tmp = maxITs[i];
			maxITs[i] = maxITs[opt_pos];
			maxITs[opt_pos] = tmp;

			tmp = selected[i];
			selected[i] = selected[opt_pos];
			selected[opt_pos] = tmp;
		}
	}
}

// clear all items
void TempRateMats::clear() {
	for (int i=0; i<(int)matrices.size(); i++)
		delete(matrices[i]); // checked
	matrices.clear();
	loglikelihood.clear();
	ICs.clear();
	parentICs.clear();
	numRates.clear();
	matrices_set.clear();
	for (int i=0; i<(int)PSs.size(); i++)
		if (PSs[i] != NULL)
			delete PSs[i];
	PSs.clear();
	for (int i=0; i<(int)VSs.size(); i++)
		if (VSs[i] != NULL)
			delete VSs[i];
	VSs.clear();
	untouches.clear();
	numIters.clear();
	matrices_map.clear();
	lowPossibleICs.clear();
	
	for (int i=0; i<(int)parentPSs.size(); i++)
		if (parentPSs[i] != NULL)
			delete parentPSs[i];
	parentPSs.clear();
	
	for (int i=0; i<(int)parentVSs.size(); i++)
		if (parentVSs[i] != NULL)
			delete parentVSs[i];
	parentVSs.clear();
	
	modes.clear();
	maxITs.clear();
	selected.clear();
	
	N = 0;
	checkedItems = 0;
	sortedNItems.clear();
	sortedICs.clear();
}

// check whether the rate matrix is inside
bool TempRateMats::isInside(vector<int>& matrix) {
	return (matrices_set.find(matrix)!=matrices_set.end());
}

// print content
void TempRateMats::print_content() {
	for (int i=0; i<(int)matrices.size(); i++) {
		for (int k=0; k<(int)matrices[i]->size(); k++) {
			if (k>0)
				cout << ",";
			cout << matrices[i]->at(k);
		}
		// cout << " [numRateGrps: " << numRates[i] << "] ";
		cout << " [loglikelihood: " << longDoublToStr(loglikelihood[i],ANS_DECI) << "] ";
		if (C_REG==0)
			cout << " [IC: " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		else
			cout << " [regularized IC: " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		// cout << " [best IC: " << longDoublToStr(lowPossibleICs[i],ANS_DECI) << "]";
		cout << endl;
	}
}

// print content
void TempRateMats::print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string description, string ICName) {
	if (fout==NULL) {
		cout << "------------------------------------------------------" << endl;
		cout << description << endl;
		if (C_REG==0)
			cout << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			cout << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<(int)matrices.size(); i++) {
			if (ICs[i] >= 9e30)
				continue;
			cout << numRates[i] << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
			cout << topStr << endl;
		}
	} else {
		(*fout) << "------------------------------------------------------" << endl;
		(*fout) << description << endl;
		if (C_REG==0)
			(*fout) << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			(*fout) << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<(int)matrices.size(); i++) {
			if (ICs[i] >= 9e30)
				continue;
			(*fout) << numRates[i] << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
			(*fout) << topStr << endl;
		}
	}
}

// print content for the first N items
void TempRateMats::print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string description, string ICName, int firstN) {
	if (firstN > (int) matrices.size())
		firstN = (int) matrices.size();
	if (fout==NULL) {
		cout << "------------------------------------------------------" << endl;
		cout << description << endl;
		if (C_REG==0)
			cout << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			cout << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<firstN; i++) {
			if (ICs[i] >= 9e30)
				continue;
			cout << numRates[i] << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
			cout << topStr << endl;
		}
	} else {
		(*fout) << "------------------------------------------------------" << endl;
		(*fout) << description << endl;
		if (C_REG==0)
			(*fout) << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			(*fout) << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<firstN; i++) {
			if (ICs[i] >= 9e30)
				continue;
			(*fout) << numRates[i] << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
			(*fout) << topStr << endl;
		}
	}
}


// print the first N items
void TempRateMats::print_first_N(int N, string& ICName) {
	if (N > (int)matrices.size())
		N = (int)matrices.size();
	for (int i=0; i<N; i++) {
		for (int k=0; k<(int)matrices[i]->size(); k++) {
			if (k>0)
				cout << ",";
			cout << matrices[i]->at(k);
		}
		// cout << " [numRateGrps: " << numRates[i] << "] ";
		cout << " [loglikelihood: " << longDoublToStr(loglikelihood[i],ANS_DECI) << "] ";
		if (C_REG==0)
			cout << " [" << ICName << ": " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		else
			cout << " [regularized " << ICName << ": " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		cout << endl;
	}
}

// print the last N items
void TempRateMats::print_last_N(int N, string& ICName) {
	if (N > (int)matrices.size())
		N = (int)matrices.size();
	int startItem = (int)matrices.size()-N;
	for (int i=startItem; i<(int)matrices.size(); i++) {
		for (int k=0; k<(int)matrices[i]->size(); k++) {
			if (k>0)
				cout << ",";
			cout << matrices[i]->at(k);
		}
		// cout << " [numRateGrps: " << numRates[i] << "] ";
		cout << " [loglikelihood: " << longDoublToStr(loglikelihood[i],ANS_DECI) << "] ";
		if (C_REG==0)
			cout << " [" << ICName << ": " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		else
			cout << " [regularized " << ICName << ": " << longDoublToStr(ICs[i],ANS_DECI) << "]";
		// cout << "  [Lowest possible IC: " << longDoublToStr(lowPossibleICs[i],ANS_DECI) << "]";
		cout << endl;
	}
}

// assign from other matrix
void TempRateMats::assign(TempRateMats& otherMats, int fromRows, int numRows) {
	clear();
	for (int i=fromRows; i<fromRows+numRows; i++) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign((otherMats.matrices[i])->begin(), (otherMats.matrices[i])->end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(otherMats.loglikelihood[i]);
		ICs.push_back(otherMats.ICs[i]);
		parentICs.push_back(otherMats.parentICs[i]);
		numRates.push_back(otherMats.numRates[i]);
		if (otherMats.PSs[i] == NULL) {
			PSs.push_back(NULL);
		} else {
			ParameterSet* newPS = new ParameterSet(otherMats.PSs[i]->num_chars);
			newPS->copyFrom(*(otherMats.PSs[i]));
			PSs.push_back(newPS);
		}
		if (otherMats.VSs[i] == NULL) {
			VSs.push_back(NULL);
		} else {
			VariableSet* newVS = new VariableSet(otherMats.VSs[i]->num_chars);
			newVS->copyFrom(*(otherMats.VSs[i]));
			VSs.push_back(newVS);
		}
		
		if (otherMats.parentPSs[i] == NULL) {
			parentPSs.push_back(NULL);
		} else {
			ParameterSet* newPS = new ParameterSet(otherMats.parentPSs[i]->num_chars);
			newPS->copyFrom(*(otherMats.parentPSs[i]));
			parentPSs.push_back(newPS);
		}
		
		if (otherMats.parentVSs[i] == NULL) {
			parentVSs.push_back(NULL);
		} else {
			VariableSet* newVS = new VariableSet(otherMats.parentVSs[i]->num_chars);
			newVS->copyFrom(*(otherMats.parentVSs[i]));
			parentVSs.push_back(newVS);
		}
		
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign((otherMats.untouches[i])->begin(), (otherMats.untouches[i])->end());
		untouches.push_back(new_untouch);
		numIters.push_back(otherMats.numIters[i]);
		lowPossibleICs.push_back(otherMats.lowPossibleICs[i]);
		modes.push_back(otherMats.modes[i]);
		maxITs.push_back(otherMats.maxITs[i]);
		selected.push_back(otherMats.selected[i]);
	}
}

// insert items from other matrix
void TempRateMats::insert(TempRateMats* otherMats, int fromRows, int numRows) {
	for (int i=fromRows; i<fromRows+numRows; i++) {
		vector<int>* new_matrix = new vector<int>;
		new_matrix->assign((otherMats->matrices[i])->begin(), (otherMats->matrices[i])->end());
		matrices.push_back(new_matrix);
		loglikelihood.push_back(otherMats->loglikelihood[i]);
		ICs.push_back(otherMats->ICs[i]);
		parentICs.push_back(otherMats->parentICs[i]);
		numRates.push_back(otherMats->numRates[i]);
		if (otherMats->PSs[i] == NULL) {
			PSs.push_back(NULL);
		} else {
			ParameterSet* newPS = new ParameterSet(otherMats->PSs[i]->num_chars);
			newPS->copyFrom(*(otherMats->PSs[i]));
			PSs.push_back(newPS);
		}
		if (otherMats->VSs[i] == NULL) {
			VSs.push_back(NULL);
		} else {
			VariableSet* newVS = new VariableSet(otherMats->VSs[i]->num_chars);
			newVS->copyFrom(*(otherMats->VSs[i]));
			VSs.push_back(newVS);
		}
		
		if (otherMats->parentPSs[i] == NULL) {
			parentPSs.push_back(NULL);
		} else {
			ParameterSet* newPS = new ParameterSet(otherMats->parentPSs[i]->num_chars);
			newPS->copyFrom(*(otherMats->parentPSs[i]));
			parentPSs.push_back(newPS);
		}
		
		if (otherMats->parentVSs[i] == NULL) {
			parentVSs.push_back(NULL);
		} else {
			VariableSet* newVS = new VariableSet(otherMats->parentVSs[i]->num_chars);
			newVS->copyFrom(*(otherMats->parentVSs[i]));
			parentVSs.push_back(newVS);
		}
		
		vector<int>* new_untouch = new vector<int>;
		new_untouch->assign((otherMats->untouches[i])->begin(), (otherMats->untouches[i])->end());
		untouches.push_back(new_untouch);
		numIters.push_back(otherMats->numIters[i]);
		lowPossibleICs.push_back(otherMats->lowPossibleICs[i]);
		modes.push_back(otherMats->modes[i]);
		maxITs.push_back(otherMats->maxITs[i]);
		selected.push_back(otherMats->selected[i]);
	}
}

// assign parent PS and parent VS
void TempRateMats::assignParent(ParameterSet& ps, VariableSet& vs, int row) {
	if (parentPSs[row]==NULL)
		parentPSs[row] = new ParameterSet(ps.num_chars);
	if (parentVSs[row]==NULL)
		parentVSs[row] = new VariableSet(vs.num_chars);
	parentPSs[row]->copyFrom(ps);
	parentVSs[row]->copyFrom(vs);
}


// number of items
int TempRateMats::size() {
	return (int)ICs.size();
}

// remove the first k item
void TempRateMats::remove_first_k_items(int k) {

	int i;
	for (i=0; i<k; i++) {
		delete matrices[i];
		if (PSs[i]!=NULL)
			delete PSs[i];
		if (VSs[i]!=NULL)
			delete VSs[i];
		delete untouches[i];
		if (parentPSs[i]!=NULL)
			delete parentPSs[i];
		if (parentVSs[i]!=NULL)
			delete parentVSs[i];
	}
	matrices.erase(matrices.begin(),matrices.begin()+k);
	loglikelihood.erase(loglikelihood.begin(),loglikelihood.begin()+k);
	ICs.erase(ICs.begin(),ICs.begin()+k);
	parentICs.erase(parentICs.begin(),parentICs.begin()+k);
	numRates.erase(numRates.begin(),numRates.begin()+k);
	PSs.erase(PSs.begin(),PSs.begin()+k);
	VSs.erase(VSs.begin(),VSs.begin()+k);
	untouches.erase(untouches.begin(),untouches.begin()+k);
	numIters.erase(numIters.begin(),numIters.begin()+k);
	lowPossibleICs.erase(lowPossibleICs.begin(),lowPossibleICs.begin()+k);
	parentPSs.erase(parentPSs.begin(),parentPSs.begin()+k);
	parentVSs.erase(parentVSs.begin(),parentVSs.begin()+k);
	modes.erase(modes.begin(),modes.begin()+k);
	maxITs.erase(maxITs.begin(),maxITs.begin()+k);
	selected.erase(selected.begin(),selected.begin()+k);
	
	N = 0;
	checkedItems = 0;
	sortedNItems.clear();
	sortedICs.clear();

}

// remove the last k item
void TempRateMats::remove_last_k_items(int k) {

	int i;
	if (k > (int)matrices.size())
		k = (int)matrices.size();

	int startItem = (int)matrices.size() - k;
	for (i=(int)matrices.size()-1; i>=startItem; i--) {
		delete matrices[i];
		if (PSs[i]!=NULL)
			delete PSs[i];
		if (VSs[i]!=NULL)
			delete VSs[i];
		delete untouches[i];
		if (parentPSs[i]!=NULL)
			delete parentPSs[i];
		if (parentVSs[i]!=NULL)
			delete parentVSs[i];
		matrices.pop_back();
		loglikelihood.pop_back();
		ICs.pop_back();
		parentICs.pop_back();
		numRates.pop_back();
		PSs.pop_back();
		VSs.pop_back();
		untouches.pop_back();
		numIters.pop_back();
		lowPossibleICs.pop_back();
		parentPSs.pop_back();
		parentVSs.pop_back();
		modes.pop_back();
		maxITs.pop_back();
		selected.pop_back();
	}
	
	N = 0;
	checkedItems = 0;
	sortedNItems.clear();
	sortedICs.clear();

}

// move the last k items to another array
void TempRateMats::move_last_k_items(TempRateMats* anotherArray, int k, int insertToSet) {

	int i;
	if (k > (int)matrices.size())
		k = (int)matrices.size();

	int startItem = (int)matrices.size() - k;
	for (i=(int)matrices.size()-1; i>=startItem; i--) {
		if (!insertToSet || anotherArray->matrices_set.find(*(matrices[i]))==anotherArray->matrices_set.end()) {
			if (insertToSet)
				anotherArray->matrices_set.insert(*(matrices[i]));
			anotherArray->matrices.push_back(matrices[i]);
			anotherArray->loglikelihood.push_back(loglikelihood[i]);
			anotherArray->ICs.push_back(ICs[i]);
			anotherArray->parentICs.push_back(parentICs[i]);
			anotherArray->numRates.push_back(numRates[i]);
			anotherArray->PSs.push_back(PSs[i]);
			anotherArray->VSs.push_back(VSs[i]);
			anotherArray->untouches.push_back(untouches[i]);
			anotherArray->numIters.push_back(numIters[i]);
			anotherArray->lowPossibleICs.push_back(lowPossibleICs[i]);
			anotherArray->parentPSs.push_back(parentPSs[i]);
			anotherArray->parentVSs.push_back(parentVSs[i]);
			anotherArray->modes.push_back(modes[i]);
			anotherArray->maxITs.push_back(maxITs[i]);
			anotherArray->selected.push_back(selected[i]);
		} else {
			delete matrices[i];
			if (PSs[i]!=NULL)
				delete PSs[i];
			if (VSs[i]!=NULL)
				delete VSs[i];
			delete untouches[i];
			if (parentPSs[i]!=NULL)
				delete parentPSs[i];
			if (parentVSs[i]!=NULL)
				delete parentVSs[i];
		}
		matrices.pop_back();
		loglikelihood.pop_back();
		ICs.pop_back();
		parentICs.pop_back();
		numRates.pop_back();
		PSs.pop_back();
		VSs.pop_back();
		untouches.pop_back();
		numIters.pop_back();
		lowPossibleICs.pop_back();
		parentPSs.pop_back();
		parentVSs.pop_back();
		modes.pop_back();
		maxITs.pop_back();
		selected.pop_back();
	}
	
	N = 0;
	checkedItems = 0;
	sortedNItems.clear();
	sortedICs.clear();

}

// reset all the untouches to zeros
int TempRateMats::reset_untouches() {
	int i,j;
	for (i=0; i<(int)untouches.size(); i++) {
		for (j=0; j<(int)untouches[i]->size(); j++) {
			untouches[i]->at(j)=0;
		}
	}
	if (untouches.size() > 0)
		return (int)untouches[0]->size();
	else
		return 0;
}

// reset all the untouchedNodes for those edges without the consistent rate matrix
int TempRateMats::reset_untouches(vector<int>* consistentRates) {
	int i,j;
	int num = 0;
	for (i=0; i<(int)untouches.size(); i++) {
		for (j=0; j<(int)untouches[i]->size(); j++) {
			if (consistentRates!=NULL && j < (int)consistentRates->size() && consistentRates->at(j)==1) {
				untouches[i]->at(j)=1;
			} else {
				if (i==0)
					num++;
				untouches[i]->at(j)=0;
			}
		}
	}
	return num;
}


// print the set of optimal rate matrices
void TempRateMats::print_optimal_rateMatrices(ofstream* fout, int* topMatrix, vector<string>* leafList, int nsites, int ICType) {
	double bestIC = ICs[0];
	int i;
	int sameRootMat;
	int numEdges = 2*leafList->size()-2;
	string ICName = getICName(ICType);
	if (fout==NULL) {
		cout << "[Set of best rate matrix arrangements]" << endl;
		for (i=0; i<(int)lowPossibleICs.size(); i++) {
			if (lowPossibleICs[i] <= bestIC) {
				cout << rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
				sameRootMat = (matrices[i]->at(numEdges-1)==matrices[i]->at(numEdges-2));
				if (C_REG==0) {
					cout << " " << ICName << " " << longDoublToStr(ICs[i],ANS_DECI) << endl;
				} else {
					cout << " " << ICName << "(regularized): " << longDoublToStr(ICs[i],ANS_DECI);
					double logL = loglikelihood[i];
					// compute the IC value with no regularization
					double IC_no_reg = getIC_no_reg(logL, numRates[i], nsites, (int) leafList->size(), ICType, sameRootMat);
					cout << " " << ICName << "(no regularization): " << longDoublToStr(IC_no_reg,ANS_DECI) << endl;
				}
			}
		}
	} else {
		(*fout) << "[Set of best rate matrix arrangements]" << endl;
		for (i=0; i<(int)lowPossibleICs.size(); i++) {
			if (lowPossibleICs[i] <= bestIC) {
				(*fout) << rateMatrixToTreeFormat(topMatrix, matrices[i], leafList);
				sameRootMat = (matrices[i]->at(numEdges-1)==matrices[i]->at(numEdges-2));
				if (C_REG==0) {
					(*fout) << " " << ICName << " " << longDoublToStr(ICs[i],ANS_DECI) << endl;
				} else {
					(*fout) << " " << ICName << "(regularized): " << longDoublToStr(ICs[i],ANS_DECI) << endl;
					double logL = loglikelihood[i];
					// compute the IC value with no regularization
					double IC_no_reg = getIC_no_reg(logL, numRates[i], nsites, (int) leafList->size(), ICType, sameRootMat);
					(*fout) << " " << ICName << "(no regularization): " << longDoublToStr(IC_no_reg,ANS_DECI) << endl;
				}
			}
		}
	}
}

// print out the consensus of the optimal rate matrices
void TempRateMats::print_consensus(ofstream* fout, int* topMatrix, vector<string>* leafList, string& ICName, double IC_Thres) {
	int numEdges = (int) matrices[0]->size();
	vector<int> consensus;
	vector<double> confidence_percent;
	int i,j;
	double* dist = new double[numEdges];
	double bestIC = ICs[0];
	double highestIC = bestIC + IC_Thres;
	double curr_weight;
	int curr_rate_grp;
	double curr_total;
	int highest_rate_grp;
	double highest_weight;
	for (i=0; i<numEdges; i++) {
		// reset the "dist"
		memset(dist,0,sizeof(double)*numEdges);
		curr_total = 0.0;
		for (j=0; j<(int)lowPossibleICs.size(); j++) {
			if (lowPossibleICs[j] <= bestIC) {
				curr_weight = (highestIC - ICs[j]) / IC_Thres;
				curr_rate_grp = matrices[j]->at(i);
				dist[curr_rate_grp-1] += curr_weight;
				curr_total += curr_weight;
			}
		}
		// check which rate grp has the highest total weight
		highest_rate_grp = matrices[0]->at(i);
		highest_weight = dist[highest_rate_grp-1];
		for (j=0; j<numEdges; j++) {
			if (dist[j] > highest_weight) {
				highest_weight = dist[j];
				highest_rate_grp = j+1;
			}
		}
		consensus.push_back(highest_rate_grp);
		confidence_percent.push_back(highest_weight / curr_total);
	}
	(*fout) << "[Consensus of best rate matrix arrangements]" << endl;
	(*fout) << rateMatrixToTreeFormat(topMatrix, &consensus, leafList) << endl;
}


// remove the item i if aList->at(i) == remove_value
void TempRateMats::remove_items(vector<int>* aList, int remove_value) {
	int new_size = 0;
	int i;

	for (i=0; i<size(); i++) {
		if (aList->at(i)!=remove_value) {
			// item i should be placed into the position "new_size"
			if (i > new_size) {
				matrices[new_size] = matrices[i];
				loglikelihood[new_size] = loglikelihood[i];
				ICs[new_size] = ICs[i];
				parentICs[new_size] = parentICs[i];
				numRates[new_size] = numRates[i];
				PSs[new_size] = PSs[i];
				VSs[new_size] = VSs[i];
				untouches[new_size] = untouches[i];
				numIters[new_size] = numIters[i];
				lowPossibleICs[new_size] = lowPossibleICs[i];
				parentPSs[new_size] = parentPSs[i];
				parentVSs[new_size] = parentVSs[i];
				modes[new_size] = modes[i];
				maxITs[new_size] = maxITs[i];
				selected[new_size] = selected[i];
			}
			new_size++;
		} else {
			// item i should be removed
			if (matrices[i]!=NULL)
				delete matrices[i];
			if (PSs[i]!=NULL)
				delete PSs[i];
			if (VSs[i]!=NULL)
				delete VSs[i];
			if (untouches[i]!=NULL)
				delete untouches[i];
			if (parentPSs[i]!=NULL)
				delete parentPSs[i];
			if (parentVSs[i]!=NULL)
				delete parentVSs[i];
		}
	}
	if (new_size < size()) {
		if (new_size > 0) {
			matrices.resize(new_size);
			loglikelihood.resize(new_size);
			ICs.resize(new_size);
			parentICs.resize(new_size);
			numRates.resize(new_size);
			PSs.resize(new_size);
			VSs.resize(new_size);
			untouches.resize(new_size);
			numIters.resize(new_size);
			lowPossibleICs.resize(new_size);
			parentPSs.resize(new_size);
			parentVSs.resize(new_size);
			modes.resize(new_size);
			maxITs.resize(new_size);
			selected.resize(new_size);
		} else {
			matrices.clear();
			loglikelihood.clear();
			ICs.clear();
			parentICs.clear();
			numRates.clear();
			PSs.clear();
			VSs.clear();
			untouches.clear();
			numIters.clear();
			lowPossibleICs.clear();
			parentPSs.clear();
			parentVSs.clear();
			modes.clear();
			maxITs.clear();
			selected.clear();
		}
		N = 0;
		checkedItems = 0;
		sortedNItems.clear();
		sortedICs.clear();
	}
}


//=======================================
// OptimalRateMats
//=======================================

// insert
void OptimalRateMats::insertOpt(int numRate, vector<int>& matrix, double loglike, double IC) {
	if ((int)matrices.size() < numRate) {
		matrices.resize(numRate);
		loglikelihood.resize(numRate,-9.99e30);
		ICs.resize(numRate,9.99e30);
		isValid.resize(numRate,false);
	}
	if ((!isValid[numRate-1]) || (IC < ICs[numRate-1])) {
		matrices[numRate-1].assign(matrix.begin(),matrix.end());
		loglikelihood[numRate-1] = loglike;
		ICs[numRate-1] = IC;
		isValid[numRate-1] = true;
	}
}

// obtain the optimal values for the specific number of rate matrix allowed
// return false if it is not valid
bool OptimalRateMats::getOpt(int numRate, vector<int>& matrix, double& loglike, double& IC) {
	if (numRate > (int)matrices.size()) {
		return false;
	}
	if (isValid[numRate-1]) {
		matrix.assign(matrices[numRate-1].begin(),matrices[numRate-1].end());
		loglike=loglikelihood[numRate-1];
		IC=ICs[numRate-1];
	}
	return isValid[numRate-1];
}

// print
void OptimalRateMats::print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string& ICName) {
	if (fout==NULL) {
		cout << "------------------------------------------------------" << endl;
		cout << "Optimal ICs for different # of rate matrices:" << endl;
		if (C_REG == 0)
			cout << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			cout << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<(int)matrices.size(); i++) {
			if (ICs[i] >= 9e30)
				continue;
			cout << i+1 << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, &(matrices[i]), leafList);
			cout << topStr << endl;
		}
	} else {
		(*fout) << "------------------------------------------------------" << endl;
		(*fout) << "Optimal ICs for different # of rate matrices:" << endl;
		if (C_REG == 0)
			(*fout) << "# of rate matrices \t loglikelihood \t " << ICName << " \t rate matrix arrangement" << endl;
		else
			(*fout) << "# of rate matrices \t loglikelihood \t regularized " << ICName << " \t rate matrix arrangement" << endl;
		for (int i=0; i<(int)matrices.size(); i++) {
			if (ICs[i] >= 9e30)
				continue;
			(*fout) << i+1 << "\t" << longDoublToStr(loglikelihood[i],ANS_DECI) << "\t" << longDoublToStr(ICs[i],ANS_DECI) << "\t";
			// the rate matrix arrangement
			string topStr = rateMatrixToTreeFormat(topMatrix, &(matrices[i]), leafList);
			(*fout) << topStr << endl;
		}
	}
}

// print the optimal rate matrix
void OptimalRateMats::print_optimal_rateMatrix(ofstream* fout, int* topMatrix, vector<string>* leafList) {
	double bestIC = ICs[0];
	int best_i = 0;
	for (int i=1; i<(int)matrices.size(); i++) {
		if (bestIC > ICs[i]) {
			bestIC = ICs[i];
			best_i = i;
		}
	}
	if (fout==NULL) {
		cout << "Best rate matrix arrangement:" << endl;
		cout << rateMatrixToTreeFormat(topMatrix, &(matrices[best_i]), leafList) << endl;
		if (C_REG==0)
			cout << "IC: " << longDoublToStr(bestIC,ANS_DECI) << endl;
		else
			cout << "regularized IC: " << longDoublToStr(bestIC,ANS_DECI) << endl;
	} else {
		(*fout) << "[Best rate matrix arrangement]" << endl;
		(*fout) << rateMatrixToTreeFormat(topMatrix, &(matrices[best_i]), leafList) << endl;
		if (C_REG==0)
			(*fout) << "IC: " << longDoublToStr(bestIC,ANS_DECI) << endl;
		else
			(*fout) << "regularized IC: " << longDoublToStr(bestIC,ANS_DECI) << endl;
	}
}

// clear the optimal rate matrix
void OptimalRateMats::clear() {
	int i;
	for (i=0; i<(int)matrices.size(); i++)
		matrices[i].clear();
	matrices.clear();
	loglikelihood.clear();
	ICs.clear();
	isValid.clear();
}
