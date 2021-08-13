/*
 *
 * rateMatrixResult.h
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

#ifndef __RAL_RAS__rateMatrixResult__
#define __RAL_RAS__rateMatrixResult__

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "parameters.h"
#include "definitions.h"
#include "core.h"

using namespace std;

// generate the set of rate matrices
struct compForIntArray {
	bool operator() (const vector<int>& lhs, const vector<int>& rhs) const {
		for (int i=0;i<(int)lhs.size();i++) {
			if (lhs[i] != rhs[i])
				return (lhs[i] < rhs[i]);
		}
		// both exactly the same
		return false;
	}
};

// Array of rate matrices used inside the iterations
class TempRateMats {
public:
	vector<vector<int>* > matrices;
	vector<double> loglikelihood;
	vector<double> ICs;
	vector<double> parentICs;
	vector<int> numRates;
	vector<ParameterSet*> PSs;
	vector<VariableSet*> VSs;
	vector<vector<int>* > untouches;
	vector<int> numIters;
	set<vector<int>,compForIntArray> matrices_set;
	map<vector<int>,int,compForIntArray> matrices_map;
	vector<double> lowPossibleICs; // the lowest possible IC value due to the deviation from the optimization engine
	vector<ParameterSet*> parentPSs;
	vector<VariableSet*> parentVSs;
	vector<int> modes;
	vector<int> maxITs;
	vector<int> selected;
	
	// keep the top N items and keep them in a sorted order
	int N;
	int checkedItems;
	vector<int> sortedNItems;
	vector<double> sortedICs;
	void createTopN(int N); // get the top N items and arrange them in a sorted order in accending order of IC values
	void maintainTopN(int N); // since new items have been added, this function is to maintain the top N items
	double topN_IC; // the IC score of the last item of the top N items

	// constructor
	TempRateMats();
	
	// clear all items
	void clear();

	// insert an empty item
	void insert();

	// insert an item
	void insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet);
	void insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet, int mode, int maxIT);
	void insert(vector<int>& matrix, vector<int>& untouch, double parentIC, ParameterSet& ps, VariableSet& vs, int insertToSet);
	void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC);
    void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, double IC, int node, int insertToSet, double lowPossibleIC);
    
	void insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet, ParameterSet& parent_ps, VariableSet& parent_vs);

	void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, ParameterSet& parent_ps, VariableSet& parent_vs);
	

	void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, int mode, int maxIT);
	
	void insert(vector<int>& matrix, vector<int>& untouch, double parentIC, int node, int insertToSet, ParameterSet& parent_ps, VariableSet& parent_vs, int mode, int maxIT);
	
    void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, double IC, int node, int insertToSet, double lowPossibleIC, int mode, int maxIT);

    void insert(vector<int>& matrix, vector<int>& untouch, double logL, double parentIC, ParameterSet& ps, VariableSet& vs, double IC, int node, int insertToSet, double lowPossibleIC, ParameterSet& parent_ps, VariableSet& parent_vs, int mode, int maxIT);

	// insert an item
	void insertwMap(vector<int>& matrix, double logL, double IC, double parentIC, ParameterSet& ps, VariableSet& vs, 
			vector<int>& untouch, int numIter, double lowPossibleIC);

    // insert an item
    void insertwMap(vector<int>& matrix, double logL, double IC, double parentIC,
            vector<int>& untouch, int numIter, double lowPossibleIC);

    // insert an item
    void insertwMap(vector<int>& matrix, double logL, double IC, vector<int>& untouch, double lowPossibleIC);
    
	// find ID by matrix
	int getID(vector<int>& matrix);

	// rearrange the first N items inside the matrices
	// such that the top K items (among the first N items) are with the least IC values
	void rearrangeFirstNItems(int topK, int N);

	// rearrange the last N items inside the matrices
	// such that the top K items (among the last N items) are with the least IC values
	void rearrangeLastNItems(int topK, int N);

	// check whether the rate matrix is inside
	bool isInside(vector<int>& matrix);

	// assign from other matrix
	void assign(TempRateMats& otherMats, int fromRows, int numRows);

	// insert items from other matrix
	void insert(TempRateMats* otherMats, int fromRows, int numRows);
	
	// assign parent PS and parent VS
	void assignParent(ParameterSet& ps, VariableSet& vs, int row);
	
	// number of items
	int size();

	// print content
	void print_content();
	void print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string description, string ICName);
	// print content for the first N items
	void print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string description, string ICName, int firstN);
	
	// print the first N items
	void print_first_N(int N, string& ICName);

	// print the last N items
	void print_last_N(int N, string& ICName);

	// remove the first k item
	void remove_first_k_items(int k);

	// remove the last k item
	void remove_last_k_items(int k);

	// move the last k items to another array
	void move_last_k_items(TempRateMats* anotherArray, int k, int insertToSet);

	// reset all the untouches to zeros
	int reset_untouches();

	// reset all the untouchedNodes for those edges without the consistent rate matrix
	int reset_untouches(vector<int>* consistentRates);

	// print the set of optimal rate matrices
	// void print_optimal_rateMatrices(ofstream* fout, int* topMatrix, vector<string>* leafList, string& ICName);
	// print the set of optimal rate matrices
	void print_optimal_rateMatrices(ofstream* fout, int* topMatrix, vector<string>* leafList, int nsites, int ICType);

	// print out the consensus of the optimal rate matrices
	void print_consensus(ofstream* fout, int* topMatrix, vector<string>* leafList, string& ICName, double IC_Thres);
	
	// remove the item i if aList->at(i) == remove_value
	void remove_items(vector<int>* aList, int remove_value);
};


// The optimal rate matrice for each number of rates allowed
class OptimalRateMats {
public:
	vector<vector<int> > matrices;
	vector<double> loglikelihood;
	vector<double> ICs;
	vector<bool> isValid;

	// insert
	void insertOpt(int numRate, vector<int>& matrix, double loglike, double IC);

	// print
	void print_content(ofstream* fout, int* topMatrix, vector<string>* leafList, string& ICName);

	// print the optimal rate matrix
	void print_optimal_rateMatrix(ofstream* fout, int* topMatrix, vector<string>* leafList);

	// obtain the optimal values for the specific number of rate matrix allowed
	// return false if it is not valid
	bool getOpt(int numRate, vector<int>& matrix, double& loglike, double& IC);

	// clear
	void clear();

};


#endif /* defined(__RAL_RAS__rateMatrixResult__) */


