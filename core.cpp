/*
 *
 * core.cpp
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


#include "core.h"


int howManyCharInside(string str, char c) {
	int n=0;
	for (int i=0; i<(int)str.length(); i++) {
		if (str[i]==c)
			n++;
	}
	return n;
}



// Input: topology string (can also include the rate group and the edge lengths)
// Output: 1. Topology matrix; 2. number of rows of the matrix, 3. List of leaf names,
//         4. corresponding rate matrix, 5. edge length
void genTopMatrixFrStr(string newickStr, int** topMatrix, int* rowNum, vector<string>** leafList, vector<int>* rateMatrix, vector<double>* edgeLens) {

	trim(newickStr);
	if (newickStr.length()==0) {
		cerr << "Error! The topology string is empty" << endl;
		exit(1);
	}

	// get the number of rows in the matrix
	(*rowNum) = howManyCharInside(newickStr, ',');

	if ((*rowNum) == 0) {
		*topMatrix = NULL;
		*leafList = NULL;
		return;
	}

	// initialize the topMatrix and leafList
	(*topMatrix) = new int[(*rowNum) * 2];
	(*leafList) = new vector<string>;

	// create another array to store the nodes waiting to be processed
	vector<int> waitNodeList;
	vector<int> waitNodeGrp;
	vector<double> waitNodeDist;

	// the number of new nodes formed
	int newNodeNum = 0;

	// name of current leaf reading
	string currLeaf = "";
	string currGrp = "";
	string currDist = "";
	int grpMode = 0;
	int distMode = 0;

	for (int i=0; i<(int)newickStr.length(); i++) {
		char c = newickStr.at(i);
		if (c == ';')
			break;
		if (c == '(' || c == ' ') {
			continue; // do nothing
		} else if (c == ',') {
			if (currGrp.length() > 0 && !(isdigit(currGrp[0]))) {
				// currGrp is not an integer
				currLeaf.append("_");
				currLeaf.append(currGrp);
				currGrp = "";
			}
			// add the current leaf to leafList if necessary
			if (currLeaf != "") {
				(*leafList)->push_back(currLeaf);
				waitNodeList.push_back(-1 * (int)((*leafList)->size()));
				currLeaf = "";
			}
			if (currGrp != "") {
				waitNodeGrp.push_back((int) atoi(currGrp.c_str()));
				currGrp = "";
			}
			if (currDist != "") {
				waitNodeDist.push_back((double) atof(currDist.c_str()));
				currDist = "";
			}
			grpMode = 0;
			distMode = 0;
		} else if (c == ')') {
			if (currGrp.length() > 0 && !(isdigit(currGrp[0]))) {
				// currGrp is not an integer
				currLeaf.append("_");
				currLeaf.append(currGrp);
				currGrp = "";
			}
			if (currLeaf != "") {
				// add the current leaf to leafList
				(*leafList)->push_back(currLeaf);
				waitNodeList.push_back(-1 * (int)((*leafList)->size()));
				currLeaf = "";
			}
			if (currGrp != "") {
				waitNodeGrp.push_back((int) atoi(currGrp.c_str()));
				currGrp = "";
			}
			if (currDist != "") {
				waitNodeDist.push_back((double) atof(currDist.c_str()));
				currDist = "";
			}
			grpMode = 1;
			distMode = 0;
			// get the last two nodes which are waiting to be processed
			int last1 = waitNodeList.back();
			waitNodeList.pop_back();
			int last2 = waitNodeList.back();
			waitNodeList.pop_back();
			// add a new row into the topMatrix
			(*topMatrix)[newNodeNum*2] = last2;
			(*topMatrix)[newNodeNum*2+1] = last1;
			if (rateMatrix != NULL && waitNodeGrp.size() >= 2) {
				int tmp = waitNodeGrp.back();
				waitNodeGrp.pop_back();
				rateMatrix->push_back(waitNodeGrp.back());
				waitNodeGrp.pop_back();
				rateMatrix->push_back(tmp);
			}
			if (edgeLens != NULL && waitNodeDist.size() >= 2) {
				double temp = waitNodeDist.back();
				waitNodeDist.pop_back();
				edgeLens->push_back(waitNodeDist.back());
				waitNodeDist.pop_back();
				edgeLens->push_back(temp);
			}
			waitNodeList.push_back(++newNodeNum);
		} else if (c == '_' && rateMatrix!=NULL) {
			if (currGrp!="") {
				// the one in "currGrp" is not the rate group, but a part of the name
				currLeaf.append("_");
				currLeaf.append(currGrp);
				currGrp = "";
			}
			grpMode = 1;
		} else if (c == ':' && edgeLens!=NULL) {
			if (currDist!="") {
				// the one in "currDist" is not the distance, but a part of the name
				currLeaf.append(":");
				currLeaf.append(currDist);
				currDist = "";
			}
			grpMode = 0;
			distMode = 1;
		} else if (grpMode==0 && distMode==0) {
			currLeaf.append(1,c); // append the character to currLeaf
		} else if (grpMode==1) {
			currGrp.append(1,c);
		} else if (distMode==1) {
			currDist.append(1,c);
		}
	}
}

// Input: topology string (NO RATE GROUP)
// Output: 1. Topology matrix; 2. number of rows of the matrix, 3. List of leaf names,
//         4. corresponding rate matrix, 5. edge length
void genTopMatrixFrStr0(string newickStr, int** topMatrix, int* rowNum, vector<string>** leafList, vector<int>* rateMatrix, vector<double>* edgeLens) {

	trim(newickStr);
	if (newickStr.length()==0) {
		cerr << "Error! The topology string is empty" << endl;
		exit(1);
	}

	// get the number of rows in the matrix
	(*rowNum) = howManyCharInside(newickStr, ',');

	if ((*rowNum) == 0) {
		*topMatrix = NULL;
		*leafList = NULL;
		return;
	}

	// initialize the topMatrix and leafList
	(*topMatrix) = new int[(*rowNum) * 2];
	(*leafList) = new vector<string>;

	// create another array to store the nodes waiting to be processed
	vector<int> waitNodeList;
	vector<double> waitNodeDist;

	// the number of new nodes formed
	int newNodeNum = 0;

	// name of current leaf reading
	string currLeaf = "";
	string currDist = "";
	int distMode = 0;
	
	// skip because of '['
	int skip = 0;

	for (int i=0; i<(int)newickStr.length(); i++) {
		char c = newickStr.at(i);
		if (c=='[') {
			skip = 1;
			continue;
		}
		if (c==']') {
			skip = 0;
			continue;
		}
		if (skip == 1)
			continue;
		if (c == ';')
			break;
		if (c == '(' || c == ' ') {
			continue; // do nothing
		} else if (c == ',') {
			// add the current leaf to leafList if necessary
			if (currLeaf != "") {
				(*leafList)->push_back(currLeaf);
				waitNodeList.push_back(-1 * (int)((*leafList)->size()));
				currLeaf = "";
			}
			if (currDist != "") {
				waitNodeDist.push_back((double) atof(currDist.c_str()));
				currDist = "";
			}
			distMode = 0;
		} else if (c == ')') {
			if (currLeaf != "") {
				// add the current leaf to leafList
				(*leafList)->push_back(currLeaf);
				waitNodeList.push_back(-1 * (int)((*leafList)->size()));
				currLeaf = "";
			}
			if (currDist != "") {
				waitNodeDist.push_back((double) atof(currDist.c_str()));
				currDist = "";
			}
			distMode = 0;
			// get the last two nodes which are waiting to be processed
			int last1 = waitNodeList.back();
			waitNodeList.pop_back();
			int last2 = waitNodeList.back();
			waitNodeList.pop_back();
			// add a new row into the topMatrix
			(*topMatrix)[newNodeNum*2] = last2;
			(*topMatrix)[newNodeNum*2+1] = last1;
			if (edgeLens != NULL && waitNodeDist.size() >= 2) {
				double temp = waitNodeDist.back();
				waitNodeDist.pop_back();
				edgeLens->push_back(waitNodeDist.back());
				waitNodeDist.pop_back();
				edgeLens->push_back(temp);
			}
			waitNodeList.push_back(++newNodeNum);
		} else if (c == ':' && edgeLens!=NULL) {
			if (currDist!="") {
				// the one in "currDist" is not the distance, but a part of the name
				currLeaf.append(":");
				currLeaf.append(currDist);
				currDist = "";
			}
			distMode = 1;
		} else if (distMode==0) {
			currLeaf.append(1,c); // append the character to currLeaf
		} else if (distMode==1) {
			currDist.append(1,c);
		}
	}
}

// Input: Topology file
// Output: 1. Topology matrix; 2. number of rows of the matrix, 3. List of leaf names,
//         4. corresponding rate matrix, 5. edge length
void genTopMatrix(char* topFile, int** topMatrix, int* rowNum, vector<string>** leafList, vector<int>* rateMatrix, vector<double>* edgeLens) {

	// read the topology file
	ifstream fin;
	fin.open(topFile); // open the file
	if (!fin.is_open()) {
		cerr << "Error opening the topology file :" << topFile << endl;
		exit(1);
	}
	string aline;
	// get the topology from the first line of the file
	if (!getline(fin, aline)) {
		cerr << "Error! The topology file " << topFile << " is empty" << endl;
		exit(1);
	}
	fin.close(); // close the file
	trim(aline);
	if (aline.length()==0) {
		cerr << "Error! The first line of the topology file " << topFile << " is empty" << endl;
		exit(1);
	}
	genTopMatrixFrStr0(aline, topMatrix, rowNum, leafList, rateMatrix, edgeLens);
}

// get the name of the IC
string getICName(int ICType) {
	// information criteria name
	// 1 - AIC; 2 - AICc; 3 - BIC
	switch (ICType) {
	case 1:
		return "AIC";
		break;
	case 2:
		return "AICc";
		break;
	case 3:
		return "BIC";
		break;
	default:
		return "IC";
		break;
	}
}

// get the IC-type of the IC name
int getICType(string ICName) {
	// information criteria name
	// 1 - AIC; 2 - AICc; 3 - BIC
	trim(ICName);
	if (ICName=="AIC")
		return 1;
	else if (ICName == "AICc")
		return 2;
	else if (ICName == "BIC")
		return 3;
	else
		return 0;
}

// compute the IC value under HAL model with regularization
double getIC_no_reg(double logVal, int numRates, int nsites, int nTaxa, int ICType, int sameRootMat) {
	int df, numRateChanges;
	int use_reg = 0;
	return getIC(logVal, numRates, nsites, nTaxa, df, ICType, numRateChanges, sameRootMat, use_reg);
}

// compute the IC value under HAL model with regularization
double getIC(double logVal, int numRates, int nsites, int nTaxa, int ICType, int numRateChanges, int sameRootMat) {
	int df;
	int use_reg = 1;
	return getIC(logVal, numRates, nsites, nTaxa, df, ICType, numRateChanges, sameRootMat, use_reg);
}

// f(n,k) =
//    2k,            for AIC
//    (2nk)/(n-k-1), for AICc
//    ln(n)k,        for BIC
double f(int n, double k, int ICType) {
	switch(ICType) {
	case 1:
		return (2.0 * k); // AIC
		break;
	case 2:
		return (2.0 * n * k) / (n - k - 1.0); // AICc
		break;
	default:
		return (log(n) * k); // BIC
		break;
	}
}

// core function for computing the IC value under HAL model
// using regularization : use_reg = 1
// no regularization : use_reg = 0
double getIC(double logVal, int numRates, int nsites, int nTaxa, int& df, int ICType, 
		int numRateChanges, int sameRootMat, int use_reg) {

	// for HAL model, number of parameters:
	// edge lengths: (for rooted tree) 2 * nTaxa - 2
	//               (for unrooted tree) 2 * nTaxa - 2 - 1
	// root frequencies: 3
	// proportion of variable site: 1
	// char distribution in invariable sites: 3
	// for each rate category:
	//     marginal distribution: 3
	//     variables in s-matrix: 5
	df = (2 * nTaxa - 2) + 3 + 1 + 3 + numRates * (3 + 5);
	if (sameRootMat)
		df--;

#ifdef GTRIFO
	df -= 6; // deduce: (root frequencies) and (char distribution in invariable sites)
#endif

	// information criteria formula
	// IC = -2 ln(L) + f(n, df)
	
	// regularization for IC
	// IC = -2 ln(L) + f(n, df + mc)

	// where f(n,k) =
	//    2k,            for AIC
	//    (2nk)/(n-k-1), for AICc
	//    ln(n)k,        for BIC

	double ic;
	if (use_reg) { // regularization
		int numConvergentChanges = numRateChanges - numRates + 1;
		if (numConvergentChanges == 0 || C_REG == 0)
			ic = -2.0 * logVal + f(nsites, df, ICType);
		else
			ic = -2.0 * logVal + f(nsites, df + numConvergentChanges * C_REG, ICType);
	} else {
		ic = -2.0 * logVal + f(nsites, df, ICType);
	}

	return ic;
}

// compute the IC value for different RAL-RAS model
// modelType  : 1 - 5
// numRateCat : number of rate categories
// numRateMat : number of rate matrices
// isUpsilonModel : whether it is an upsilon model
// output: 1. df - degree of freedom; 2. ic - IC value
double getIC(double logVal, int numRateMat, int nsites, int nTaxa, int numRateCat,
		int modelType, int& df, int isUpsilonModel, int ICType, 
		int numRateChanges, int sameRootMat, int use_reg) {

	int numEdges;
	if (isUpsilonModel==0) {
		numEdges = 2 * nTaxa - 2; // for rooted tree
		if (sameRootMat)
			numEdges--; // should be treated as unrooted tree
		switch(modelType) {
		case 1:
			// (1) edge lengths: numRateCat*numEdges, (2) root freq: numRateCat*3;
			// (3) S: numRateMat*numRateCat*5, (4) Pi: numRateMat*numRateCat*3 ,(5) alpha: numRateCat
			// (6) invariable sites distribution: 3
			df = (numEdges + 4) * numRateCat + 8 * numRateMat * numRateCat + 3;
			break;
		case 2:
			// share S
			// (1) edge lengths: numRateCat*numEdges, (2) root freq: numRateCat*3;
			// (3) S: numRateMat*5, (4) Pi: numRateMat*numRateCat*3 ,(5) alpha: numRateCat
			// (6) invariable sites distribution: 3
			df = (numEdges + 4) * numRateCat + 3 * numRateMat * numRateCat + 5 * numRateMat + 3;
			break;
		case 3:
			// share S and Pi
			// (1) edge lengths: numRateCat*numEdges, (2) root freq: numRateCat*3;
			// (3) S: numRateMat*5, (4) Pi: numRateMat*3 ,(5) alpha: numRateCat
			// (6) invariable sites distribution: 3
			df = (numEdges + 4) * numRateCat + 8 * numRateMat + 3;
			break;
		case 4:
			// share S, Pi and root freq
			// (1) edge lengths: numRateCat*numEdges, (2) root freq: 3;
			// (3) S: numRateMat*5, (4) Pi: numRateMat*3 ,(5) alpha: numRateCat
			// (6) invariable sites distribution: 3
			df = (numEdges + 1) * numRateCat + 8 * numRateMat + 6;
			break;
		case 5:
			// share S, Pi, root freq and scalar edge length
			// (1) edge lengths: numEdges, (2) root freq: 3;
			// (3) S: numRateMat*5, (4) Pi: numRateMat*3 ,(5) alpha: numRateCat
			// (6) invariable sites distribution: 3, (7) scalar: numRateCat-1
			df = 2 * numRateCat + 8 * numRateMat + numEdges + 5;
#ifdef GTRIFO
			df -= 6; // deduce: (2) and (6)
#endif
			break;
		default:
			cerr << "[getIC] Error! The type of model (which is " << modelType << ") is not between 1 and 5" << endl;
			exit(1);
		}
	} else {
		numEdges = 2 * nTaxa - 3; // for unrooted tree
		// share S, Pi, scalar edge length
		// root freq is same as Pi
		// no constant site
		// numRateMat = 1
		// (1) edge lengths: numEdges, (2) S: 5, (3) Pi: 3 , (4) alpha: numRateCat-1, (5) scalar: numRateCat-1
		df = numEdges + 6 + 2 * numRateCat;
	}
	
	double ic;
	if (use_reg) { // regularization
		int numConvergentChanges = numRateChanges - numRateMat + 1;
		if (numConvergentChanges == 0 || C_REG == 0)
			ic = -2.0 * logVal + f(nsites, df, ICType);
		else
			ic = -2.0 * logVal + f(nsites, df + numConvergentChanges * C_REG, ICType);
	} else {
		ic = -2.0 * logVal + f(nsites, df, ICType);
	}

	return ic;
}

// Input: the topology matrix
// Output: 1. the mapping between the leaf ID to the line number in the paramater file, and
//         2. the mapping between the internal node ID to the line number in the parameter file
void getMapping(int* leafID2line, int* internalID2line, int* topologyMatrix, int numOfLines) {
	int lineNum = 0;
	for (int i=0; i<numOfLines; i++) {
		for (int j=0; j<2; j++) {
			int nodeID = topologyMatrix[i*2+j];
			if (nodeID < 0) {
				// leaf
				leafID2line[-nodeID] = lineNum;
			} else {
				// internal node
				internalID2line[nodeID] = lineNum;
			}
			lineNum++;
		}
	}
}

void subOutTreeWithGrpID(int* topMatrix, vector<string>* leafList, vector<int>* rateGrps, int line) {
	int node1 = topMatrix[line*2];
	int node2 = topMatrix[line*2+1];
	int grp1 = rateGrps->at(line*2);
	int grp2 = rateGrps->at(line*2+1);
	cout << "(";
	if (node1 < 0) {
		// a leaf
		cout << leafList->at(-node1-1);
	} else {
		// an internal node
		subOutTreeWithGrpID(topMatrix,leafList,rateGrps,node1-1);
	}
	cout << "_" << grp1 << ",";
	if (node2 < 0) {
		// a leaf
		cout << leafList->at(-node2-1);
	} else {
		// an internal node
		subOutTreeWithGrpID(topMatrix,leafList,rateGrps,node2-1);
	}
	cout << "_" << grp2 << ")";

}

// Output: Print the newick tree format with rate group ID
// Input: the topology matrix, the leaf names, and the rate group IDs
void outTreeWithGrpID(int* topMatrix, vector<string>* leafList, vector<int>* rateGrps) {

	int numSpecies=(int)leafList->size();
	int numLineTopMatrix=numSpecies-1;

	cout << "(";
	subOutTreeWithGrpID(topMatrix, leafList, rateGrps, numLineTopMatrix-1);
	cout << ");" << endl;
}

void subOutTreeWithEdgeLen(int* topMatrix, vector<string>* leafList, double* edgeLens, int line, string& str) {
	int node1 = topMatrix[line*2];
	int node2 = topMatrix[line*2+1];
	double dist1 = edgeLens[line*2];
	double dist2 = edgeLens[line*2+1];
	str.append("(");
	if (node1 < 0) {
		// a leaf
		str.append(leafList->at(-node1-1));
	} else {
		// an internal node
		subOutTreeWithEdgeLen(topMatrix,leafList,edgeLens,node1-1, str);
		str.append(intToStr(node1));
	}
	str.append(":" + longDoublToStr(dist1,PARAM_DECI) + ",");
	if (node2 < 0) {
		// a leaf
		str.append(leafList->at(-node2-1));
	} else {
		// an internal node
		subOutTreeWithEdgeLen(topMatrix,leafList,edgeLens,node2-1, str);
		str.append(intToStr(node2));
	}
	str.append(":" + longDoublToStr(dist2,PARAM_DECI));
	str.append(")");
}

// Output: Print the newick tree format with edge length (i.e. the edge length = # of substituation in the edge)
// Input: the topology matrix, the leaf names, and the edge lengths
void outTreeWithEdgeLen(int* topMatrix, vector<string>* leafList, vector<double>* edgeLens) {

	int numSpecies=(int)leafList->size();
	int numLineTopMatrix=numSpecies-1;
	string str = "";
	double* edgeLensArray = new double[edgeLens->size()];
	int i;
	for (i=0; i<(int)edgeLens->size(); i++)
		edgeLensArray[i] = edgeLens->at(i);

	subOutTreeWithEdgeLen(topMatrix, leafList, edgeLensArray, numLineTopMatrix-1, str);
	cout << str << ";" << endl;
	delete[] edgeLensArray;
}

// Output: Print the newick tree format with edge length (i.e. the edge length = # of substituation in the edge)
// Input: the topology matrix, the leaf names, and the edge lengths
void outTreeWithEdgeLen(int* topMatrix, vector<string>* leafList, double* edgeLens, ofstream& fout) {

	int numSpecies=(int)leafList->size();
	int numLineTopMatrix=numSpecies-1;
	string str = "";

	subOutTreeWithEdgeLen(topMatrix, leafList, edgeLens, numLineTopMatrix-1, str);
	fout << str << ";" << endl;
}


// compute the likelihood matrix (dimension: internal nodes # x num_chars) of all the internal nodes for the specific site
// i.e. internal nodes # = number of lines in the topology matrix
// output: likelihoodMat
void computeLikelihoodMatrix(double* likelihoodMat, double* allCondProb_t, int* topologyMatrix, char* specificSites, int numLineTopMat, int num_chars) {
	
	double  leftLikelihoodMat[4];
	double  rightLikelihoodMat[4];
	
	int nodeID = 0;
	double *baseDist, *prob_t, *result;
	int i, j, k;
	for (i=0; i<numLineTopMat; i++) {
		int leftNode = topologyMatrix[i*2];
		int rightNode = topologyMatrix[i*2+1];
		if (leftNode < 0) {
			// a leaf node
			baseDist = &(base_arrays[specificSites[-leftNode-1]*4]);
		} else {
			// an iternal node
			baseDist = &(likelihoodMat[(leftNode-1)*num_chars]);
		}
		prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
		multiply41_t(leftLikelihoodMat,prob_t,baseDist);
		nodeID++;
		if (rightNode < 0) {
			// a leaf node
			baseDist = &(base_arrays[specificSites[-rightNode-1]*4]);
		} else {
			// an iternal node
			baseDist = &(likelihoodMat[(rightNode-1)*num_chars]);
		}
		prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
		multiply41_t(rightLikelihoodMat,prob_t,baseDist);
		nodeID++;
		// compute the likelihood for this new internal node
		result = &(likelihoodMat[i*num_chars]);
		mul_decimals(result, leftLikelihoodMat, rightLikelihoodMat, num_chars);
	}
}

// compute the likelihood matrix (dimension: internal nodes # x num_chars) of all the internal nodes for the specific site
// i.e. internal nodes # = number of lines in the topology matrix
// output: likelihoodMat
// Note: leftLikelihoodMat (size: num_chars) and rightLikelihoodMat (size: num_chars) are dummy
void computeLikelihoodMatrix(double* likelihoodMat, double* allCondProb_t, int* topologyMatrix, char* specificSites, int numLineTopMat, char* changedNodes,
		double* leftLikelihoodMatArr, double* rightLikelihoodMatArr, int num_chars) {


	int nodeID = 0;
	int numSpecies = numLineTopMat+1;
	double *baseDist, *prob_t, *result;
	double *leftLikelihoodMat, *rightLikelihoodMat;
	int leftNode, rightNode;
	int l,r;
	int i,j,k;
	for (i=0; i<numLineTopMat; i++) {
		leftNode = topologyMatrix[i*2];
		rightNode = topologyMatrix[i*2+1];
		// l = (leftNode < 0)?(-leftNode-1):(leftNode-1+numSpecies);
		l = leftNode+numSpecies;
		// r = (rightNode < 0)?(-rightNode-1):(rightNode-1+numSpecies);
		r = rightNode+numSpecies;
		leftLikelihoodMat = &leftLikelihoodMatArr[i*num_chars];
		rightLikelihoodMat = &rightLikelihoodMatArr[i*num_chars];

		if (changedNodes==NULL || changedNodes[l]==1) {
			if (leftNode < 0) {
				// a leaf node
				baseDist = &(base_arrays[specificSites[-leftNode-1]*4]);
			} else {
				// an iternal node
				baseDist = &(likelihoodMat[(leftNode-1)*num_chars]);
			}
			prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
			multiply41_t(leftLikelihoodMat,prob_t,baseDist);
		}
		nodeID++;
		if (changedNodes==NULL || changedNodes[r]==1) {
			if (rightNode < 0) {
				// a leaf node
				baseDist = &(base_arrays[specificSites[-rightNode-1]*4]);
			} else {
				// an iternal node
				baseDist = &(likelihoodMat[(rightNode-1)*num_chars]);
			}
			prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
			multiply41_t(rightLikelihoodMat,prob_t,baseDist);
		}
		nodeID++;

		if (changedNodes==NULL || changedNodes[l]==1 || changedNodes[r]==1) {
			// compute the likelihood for this new internal node
			result = &(likelihoodMat[i*num_chars]);
			mul_decimals(result, leftLikelihoodMat, rightLikelihoodMat, num_chars);

		}
	}
}


// compute the likelihood matrix (dimension: internal nodes # x num_chars) of all the internal nodes for the specific site
// i.e. internal nodes # = number of lines in the topology matrix
// output: likelihoodMat
// Note: leftLikelihoodMat (size: num_chars) and rightLikelihoodMat (size: num_chars) are dummy
void computeLikelihoodMatrix(double* likelihoodMat, double* preLikelihoodMat,
							 double* allCondProb_t, int* topologyMatrix, char* specificSites,
							 int numLineTopMat, char* changedNodes, int* changedEdges,
							 double* leftLikelihoodMatArr, double* rightLikelihoodMatArr,
							 double* preLeftLikelihoodMatArr, double* preRightLikelihoodMatArr,
							 int num_chars) {
	
	
	int nodeID = 0;
	int numSpecies = numLineTopMat+1;
	double *baseDist, *prob_t, *result, * preResult;
	double *leftLikelihoodMat, *rightLikelihoodMat;
	int leftNode, rightNode;
	int l,r;
	int i,j,k;
	
	for (i=0; i<numLineTopMat; i++) {
		leftNode = topologyMatrix[i*2];
		rightNode = topologyMatrix[i*2+1];
		// l = (leftNode < 0)?(-leftNode-1):(leftNode-1+numSpecies);
		l = leftNode+numSpecies;
		// r = (rightNode < 0)?(-rightNode-1):(rightNode-1+numSpecies);
		r = rightNode+numSpecies;
		leftLikelihoodMat = &leftLikelihoodMatArr[i*num_chars];
		rightLikelihoodMat = &rightLikelihoodMatArr[i*num_chars];
		
		if (changedNodes[l]) {
			
			if (changedEdges[i*2]==0) {
				// same as the previous one
				copy_decimals(leftLikelihoodMat, &preLeftLikelihoodMatArr[i*num_chars], num_chars);
			} else {
				if (leftNode < 0) {
					// a leaf node
					baseDist = &(base_arrays[specificSites[-leftNode-1]*4]);
				} else {
					// an iternal node
					baseDist = &(likelihoodMat[(leftNode-1)*num_chars]);
				}
				prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
				multiply41_t(leftLikelihoodMat,prob_t,baseDist);
			}
		}
		nodeID++;
		if (changedNodes[r]) {
			
			if (changedEdges[i*2+1]==0) {
				// same as the previous one
				copy_decimals(rightLikelihoodMat, &preRightLikelihoodMatArr[i*num_chars], num_chars);
			} else {
				if (rightNode < 0) {
					// a leaf node
					baseDist = &(base_arrays[specificSites[-rightNode-1]*4]);
				} else {
					// an iternal node
					baseDist = &(likelihoodMat[(rightNode-1)*num_chars]);
				}
				prob_t =&(allCondProb_t[nodeID*num_chars*num_chars]);
				multiply41_t(rightLikelihoodMat,prob_t,baseDist);
			}
		}
		nodeID++;

		if (changedNodes[l] | changedNodes[r]) {
			
			// compute the likelihood for this new internal node
			result = &(likelihoodMat[i*num_chars]);

			if (changedEdges[i*2]==0 && changedEdges[i*2+1]==0) {
				// same as the previous result
				preResult = &(preLikelihoodMat[i*num_chars]);
				copy_decimals(result, preResult, num_chars);
			} else {
				mul_decimals(result, leftLikelihoodMat, rightLikelihoodMat, num_chars);
			}
			
		}
	}
}



// compute the total log-likelihood for a phylogenetic tree for different variables
void getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, vector<VariableSet*>* variables, int num_chars, double* log_Ls) {
	
	int numInternalNodes = alignment.numSeqs-1;
	double*  likelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	double currLogLikelihood;
	int numSeqs = alignment.numSeqs;
	int numUniqSites = alignment.numUniqueSites;
	// char* preSites;
	char* currSites;
	int i,j,k,numOfSites;
	double var_inv;
	
	if (variables==NULL || variables->size() == 0)
		return;
	
	int n = (int) variables->size();
	memset(log_Ls, 0, n*sizeof(double));
	
	// for dummy use
	double*  leftLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	double*  rightLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	char* currChangedNodes;
	// int prePosMinDiff;
	
	for (i=0; i<numUniqSites; i++) {
		currSites = &(alignment.uniqueSites[i*numSeqs]);
		currChangedNodes = &(alignment.changedNodes[i*(numInternalNodes+numSeqs+1)]);
		computeLikelihoodMatrix(likelihoodMat, allCondProb_t, topologyMatrix, currSites, numLineTopMat, currChangedNodes, leftLikelihoodMat, rightLikelihoodMat, num_chars);
		
		numOfSites = alignment.numEachUniqueSites[i];
		
		for (k=0; k<n; k++) {
			currLogLikelihood = 0.0;
			for (j=0; j<num_chars; j++) {
				currLogLikelihood += likelihoodMat[(numLineTopMat-1)*num_chars+j]*variables->at(k)->rootNodeFreq[j];
			}
			if (alignment.isConstantSites[i]==0) {
				currLogLikelihood = log(currLogLikelihood * (variables->at(k)->alpha[0]))*numOfSites;
			} else {
				var_inv = 0.0;
				for (j=0; j<num_chars; j++) {
					if ((alignment.isConstantSites[i] >> j) & 1)
						var_inv += variables->at(k)->probXGivenInv[j];
				}
				currLogLikelihood = log(currLogLikelihood * (variables->at(k)->alpha[0]) + variables->at(k)->beta * var_inv)*numOfSites;
			}
			log_Ls[k] += currLogLikelihood;
		}
	}
	
	free_decimals(likelihoodMat);
	free_decimals(leftLikelihoodMat);
	free_decimals(rightLikelihoodMat);
	
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_L += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
	
}



// compute the total log-likelihood for a phylogenetic tree
double getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars) {

	int numInternalNodes = alignment.numSeqs-1;
	double*  likelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	double logLikelihood = 0.0;
	double currLogLikelihood;
	int numSeqs = alignment.numSeqs;
	int numUniqSites = alignment.numUniqueSites;
	// char* preSites;
	char* currSites;
	int i,j;
	double var_inv;

	// for dummy use
	double*  leftLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	double*  rightLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	char* currChangedNodes;

	for (i=0; i<numUniqSites; i++) {
		currSites = &(alignment.uniqueSites[i*numSeqs]);
		/*
		memset(changedNodes,0,sizeof(int)*(numInternalNodes+numSeqs));
		if (i>0) {
			preSites = &(alignment.uniqueSites[(i-1)*numSeqs]);
			for (j=0; j<numSeqs; j++)
				if (preSites[j]!=currSites[j])
					changedNodes[j]=1;
		} else {
			for (j=0; j<numSeqs; j++)
				changedNodes[j]=1;
		}*/
		currChangedNodes = &(alignment.changedNodes[i*(numInternalNodes+numSeqs+1)]);
		computeLikelihoodMatrix(likelihoodMat, allCondProb_t, topologyMatrix, currSites, numLineTopMat, currChangedNodes, leftLikelihoodMat, rightLikelihoodMat, num_chars);
		// computeLikelihoodMatrix(likelihoodMat, allCondProb, topologyMatrix, currSites, numLineTopMat, NULL, leftLikelihoodMat, rightLikelihoodMat, num_chars);
		currLogLikelihood = 0.0;
		for (j=0; j<num_chars; j++) {
			currLogLikelihood += likelihoodMat[(numLineTopMat-1)*num_chars+j]*variables.rootNodeFreq[j];
		}
#ifdef SHOW_LIKEHOOD_EACH_SITE_CAT
		cerr << "likelihood of this site (before multiply by alpha and numOfSites): " << currLogLikelihood << endl;
#endif
		int numOfSites = alignment.numEachUniqueSites[i];
		if (alignment.isConstantSites[i]==0) {
			currLogLikelihood = log(currLogLikelihood * (variables.alpha[0]))*numOfSites;
		} else {
			var_inv = 0.0;
			for (j=0; j<num_chars; j++) {
				if ((alignment.isConstantSites[i] >> j) & 1)
					var_inv += variables.probXGivenInv[j];
			}
			currLogLikelihood = log(currLogLikelihood * variables.alpha[0] + variables.beta * var_inv)*numOfSites;
		}
		logLikelihood += currLogLikelihood;
#ifdef SHOW_LIKEHOOD_EACH_SITE_CAT
		cerr << "logL of this site: " << currLogLikelihood << " numOfSites: " << numOfSites << " variables.alpha[0]: " << variables.alpha[0] << " alignment.isConstantSites[i]: " << (int) alignment.isConstantSites[i] << endl;
#endif
	}

	free_decimals(likelihoodMat);
	free_decimals(leftLikelihoodMat);
	free_decimals(rightLikelihoodMat);

#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_L += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif

	return logLikelihood;
}

// compute the total log-likelihood for a phylogenetic tree
double getLogLikelihood(double* allCondProb_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* likelihoodMat, double* leftLikelihoodMat, double* rightLikelihoodMat) {
	
	int numInternalNodes = alignment.numSeqs-1;
	double logLikelihood = 0.0;
	double currLogLikelihood;
	int numSeqs = alignment.numSeqs;
	int numUniqSites = alignment.numUniqueSites;
	// char* preSites;
	char* currSites;
	char* currChangedNodes;
	int i,j;
	double var_inv;

	/*
	if (likelihoodMat==NULL)
		likelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	if (leftLikelihoodMat==NULL)
		leftLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	if (rightLikelihoodMat==NULL)
		rightLikelihoodMat = malloc_decimals(numInternalNodes*num_chars);
	*/
	for (i=0; i<numUniqSites; i++) {
		currSites = &(alignment.uniqueSites[i*numSeqs]);
		currChangedNodes = &(alignment.changedNodes[i*(numInternalNodes+numSeqs+1)]);
		computeLikelihoodMatrix(likelihoodMat, allCondProb_t, topologyMatrix, currSites, numLineTopMat, currChangedNodes, leftLikelihoodMat, rightLikelihoodMat, num_chars);
		currLogLikelihood = 0.0;
		for (j=0; j<num_chars; j++) {
			currLogLikelihood += likelihoodMat[(numLineTopMat-1)*num_chars+j]*variables.rootNodeFreq[j];
		}
#ifdef SHOW_LIKEHOOD_EACH_SITE_CAT
		cerr << "likelihood of this site (before multiply by alpha and numOfSites): " << currLogLikelihood << endl;
#endif
		int numOfSites = alignment.numEachUniqueSites[i];
		if (alignment.isConstantSites[i]==0) {
			currLogLikelihood = log(currLogLikelihood * (variables.alpha[0]))*numOfSites;
		} else {
			var_inv = 0.0;
			for (j=0; j<num_chars; j++) {
				if ((alignment.isConstantSites[i] >> j) & 1)
					var_inv += variables.probXGivenInv[j];
			}
			currLogLikelihood = log(currLogLikelihood * (variables.alpha[0]) + variables.beta * var_inv)*numOfSites;
		}
		logLikelihood += currLogLikelihood;
#ifdef SHOW_LIKEHOOD_EACH_SITE_CAT
		cerr << "logL of this site: " << currLogLikelihood << " numOfSites: " << numOfSites << " variables.alpha[0]: " << variables.alpha[0] << " alignment.isConstantSites[i]: " << (int) alignment.isConstantSites[i] << endl;
#endif
	}
	
	// need to further address, potential memory leak here!!!
	// free(likelihoodMat);
	// free(leftLikelihoodMat);
	// free(rightLikelihoodMat);
	
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_L += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
	
	return logLikelihood;
}


// compute the total log-likelihood for a phylogenetic tree
// for many conditional probabilities
void getLogLikelihood(double* allCondProbArr_t, int n, vector<int>* corr_edges, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* result_Log_Ls, int gr_type) {
	
	int numInternalNodes = alignment.numSeqs-1;
	double currLogLikelihood;
	int numSeqs = alignment.numSeqs;
	int numEdges = 2 * numSeqs - 2; // for a rooted tree
	int numUniqSites = alignment.numUniqueSites;
	// char* preSites;
	char* currSites;
	int i,j,k;
	double var_inv;
	
	double*  likelihoodMatArr = malloc_decimals(n * numInternalNodes * num_chars);
	double*  leftLikelihoodMatArr = malloc_decimals(n * numInternalNodes * num_chars);
	double*  rightLikelihoodMatArr = malloc_decimals(n * numInternalNodes * num_chars);
	char* currChangedNodes;
	double *allCondProb_t, *likelihoodMat, *leftLikelihoodMat, * rightLikelihoodMat;
	double *preLikelihoodMat = NULL;
	double *preLeftLikelihoodMat = NULL;
	double *preRightLikelihoodMat = NULL;
	int *changedEdges = NULL;
	int *changedEdgesAllOne = (int*) malloc (numEdges*sizeof(int));
	int grandchild1, grandchild2;
	
	memset(result_Log_Ls, 0, sizeof(double)* n);
	
	if (corr_edges != NULL) {
		// initialize the changedEdges array
		changedEdges = (int*) malloc (numEdges*sizeof(int));
		memset(changedEdges, 0, sizeof(int)*numEdges);
		for (i=0; i<(int)corr_edges->size(); i++) {
			changedEdges[corr_edges->at(i)]=1;
		}
		for (i=0; i<numLineTopMat*2; i++) {
			if (topologyMatrix[i] > 0) {
				// internal node
				grandchild1 = (topologyMatrix[i] - 1) * 2;
				grandchild2 = (topologyMatrix[i] - 1) * 2 + 1;
				changedEdges[i] |= changedEdges[grandchild1] | changedEdges[grandchild2];
			}
		}
	} else {
		// the case when optimizing the edge length
		changedEdges = (int*) malloc (n*numEdges*sizeof(int));
		memset(changedEdges, 0, sizeof(int)*n*numEdges);
		int* curr;
		for (i=1; i<n; i++) {
			curr = &changedEdges[i*numEdges];
			if (gr_type == 1) {
				if (i%2==0) {
					// even number
					curr[i/2-1]=1;
					curr[i/2]=1;
				} else {
					// odd number
					curr[(i-1)/2]=1;
				}
			} else {
				// gr_type = 2
				curr[i-1] = 1;
				curr[i] = 1;
			}
			for (j=0; j<numLineTopMat*2; j++) {
				if (topologyMatrix[j] > 0) {
					// internal node
					grandchild1 = (topologyMatrix[j] - 1) * 2;
					grandchild2 = (topologyMatrix[j] - 1) * 2 + 1;
					curr[j] |= curr[grandchild1] | curr[grandchild2];
				}
			}
		}
	}
	for (i=0; i<numEdges; i++)
		changedEdgesAllOne[i] = 1;

	// int prePosMinDiff;
	for (i=0; i<numUniqSites; i++) {
		currSites = &(alignment.uniqueSites[i*numSeqs]);
		currChangedNodes = &(alignment.changedNodes[i*(numInternalNodes+numSeqs+1)]);
		
		for (k=0; k<n; k++) {
			
			allCondProb_t = &(allCondProbArr_t[k*numEdges*num_chars*num_chars]);
			likelihoodMat = &(likelihoodMatArr[k*numInternalNodes*num_chars]);
			leftLikelihoodMat = &(leftLikelihoodMatArr[k*numInternalNodes*num_chars]);
			rightLikelihoodMat = &(rightLikelihoodMatArr[k*numInternalNodes*num_chars]);
			
			if (k==0) {
				computeLikelihoodMatrix(likelihoodMat, preLikelihoodMat, allCondProb_t, topologyMatrix, currSites,
										numLineTopMat, currChangedNodes, changedEdgesAllOne, leftLikelihoodMat, rightLikelihoodMat,
										preLeftLikelihoodMat, preRightLikelihoodMat, num_chars);
			} else {
				
				if (corr_edges == NULL) {
					// set up the changedEdges for the case when optimizing the edge length
					computeLikelihoodMatrix(likelihoodMat, preLikelihoodMat, allCondProb_t, topologyMatrix, currSites,
											numLineTopMat, currChangedNodes, &changedEdges[k*numEdges], leftLikelihoodMat, rightLikelihoodMat,
											preLeftLikelihoodMat, preRightLikelihoodMat, num_chars);
				} else {
					computeLikelihoodMatrix(likelihoodMat, preLikelihoodMat, allCondProb_t, topologyMatrix, currSites,
											numLineTopMat, currChangedNodes, changedEdges, leftLikelihoodMat, rightLikelihoodMat,
											preLeftLikelihoodMat, preRightLikelihoodMat, num_chars);
				}
			}
			currLogLikelihood = 0.0;
			for (j=0; j<num_chars; j++) {
				currLogLikelihood += likelihoodMat[(numLineTopMat-1)*num_chars+j]*variables.rootNodeFreq[j];
			}
			int numOfSites = alignment.numEachUniqueSites[i];
			if (alignment.isConstantSites[i]==0) {
				currLogLikelihood = log(currLogLikelihood * (variables.alpha[0]))*numOfSites;
			} else {
				var_inv = 0.0;
				for (j=0; j<num_chars; j++) {
					if ((alignment.isConstantSites[i] >> j) & 1)
						var_inv += variables.probXGivenInv[j];
				}
				currLogLikelihood = log(currLogLikelihood * variables.alpha[0] + variables.beta * var_inv)*numOfSites;
			}
			result_Log_Ls[k] += currLogLikelihood;
			
			preLikelihoodMat = &(likelihoodMatArr[k*numInternalNodes*num_chars]);
			preLeftLikelihoodMat = &(leftLikelihoodMatArr[k*numInternalNodes*num_chars]);
			preRightLikelihoodMat = &(rightLikelihoodMatArr[k*numInternalNodes*num_chars]);
		}
	}
	
	free_decimals(likelihoodMatArr);
	free_decimals(leftLikelihoodMatArr);
	free_decimals(rightLikelihoodMatArr);
	free(changedEdges);
	free(changedEdgesAllOne);
	
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_L += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
	
}


// compute the total log-likelihood for a phylogenetic tree
// (for more than one rate category)
double getLogLikelihood(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars) {

	int i, j, k;
	double var_inv;
	int numInternalNodes = alignment.numSeqs-1;
	double*  likelihoodMat = malloc_decimals(numInternalNodes * num_chars);
	double  likelihoodMatCol[4];
	double logLikelihood = 0.0;
	int numSeqs = alignment.numSeqs;
	int numUniqSites = alignment.numUniqueSites;
	
	// variables.compute_rootNodeFreq_alpha();
	variables.compute_invar_constant();
	
	for (i=0; i<numUniqSites; i++) {
		char* currSites = &(alignment.uniqueSites[i*numSeqs]);
		computeLikelihoodMatrix(likelihoodMat, allCondProbSet_t[0], topologyMatrix, currSites, numLineTopMat, num_chars);
		// mul_decimals(likelihoodMatCol, &likelihoodMat[(numLineTopMat-1)*num_chars], &variables.rootNodeFreq_alpha[0], num_chars);
		mul_decimals(likelihoodMatCol, &likelihoodMat[(numLineTopMat-1)*num_chars], &variables.rootNodeFreq[0], variables.alpha[0], num_chars);
		
		for (k=1; k<(int)allCondProbSet_t.size(); k++) {
			computeLikelihoodMatrix(likelihoodMat, allCondProbSet_t[k], topologyMatrix, currSites, numLineTopMat, num_chars);
			// mul_add_decimals(likelihoodMatCol, &likelihoodMat[(numLineTopMat-1)*num_chars], &variables.rootNodeFreq_alpha[k*num_chars], num_chars);
			mul_add_decimals(likelihoodMatCol, &likelihoodMat[(numLineTopMat-1)*num_chars], &variables.rootNodeFreq[k*num_chars], variables.alpha[k], num_chars);
		}
		logLikelihood += log(likelihoodMatCol[0]+likelihoodMatCol[1]+likelihoodMatCol[2]+likelihoodMatCol[3] + variables.invar_constant[alignment.isConstantSites[i]])*alignment.numEachUniqueSites[i];
	}
	free_decimals(likelihoodMat);
	return logLikelihood;
}

// compute the total log-likelihood for a phylogenetic tree
// (for more than one rate category)
// NEWLY ADDED: siteLikelihoods - the likelihoods of each nucleotide of each site for each variable category
//                                dimension: num_chars * numUniqSites * numVarCategories
//
//              updatedCat      - if "updatedCat" >= 0, only update the category "updatedCat"
//                                else update all the categories
//
//              temporary use:  likelihoodMat, leftLikelihoodMat, rightLikelihoodMat

double getLogLikelihood(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars, double* siteLikelihoods, int updatedCat, double* likelihoodMat, double* leftLikelihoodMat, double* rightLikelihoodMat) {
	
	int i, j, k;
	double var_inv;
	int numInternalNodes = alignment.numSeqs-1;
	int numSeqs = alignment.numSeqs;
	int numUniqSites = alignment.numUniqueSites;
	int numVarCategories = (int)allCondProbSet_t.size();
	char* currSites;
	char* currChangedNodes;
	double currLikelihood;
	double currSiteLikelihood;
	int numOfSites;
	
	int startCat = 0;
	int endCat = numVarCategories-1;
	if (updatedCat>=0)
		startCat=endCat=updatedCat;
	
	for (k=startCat; k<=endCat && k<numVarCategories; k++) {
		for (i=0; i<numUniqSites; i++) {
			currSites = &(alignment.uniqueSites[i*numSeqs]);
			currChangedNodes = &(alignment.changedNodes[i*(numInternalNodes+numSeqs+1)]);
			computeLikelihoodMatrix(likelihoodMat, allCondProbSet_t[k], topologyMatrix, currSites, numLineTopMat, currChangedNodes, leftLikelihoodMat, rightLikelihoodMat, num_chars);
			copy_decimals(&siteLikelihoods[(k*numUniqSites+i)*num_chars], &likelihoodMat[(numLineTopMat-1)*num_chars], num_chars);
		}
	}
	
	// variables.compute_rootNodeFreq_alpha();
	variables.compute_invar_constant();
	
	// L(c,i) = sum_r { likelihood(i, r, c) * rootNodeFreq(r, c) * alpha(r) }
	// where c ~ [0,3] i.e. characters
	//       i ~ [0,numUniqSites-1] i.e. column
	//       r ~ [0,numVarCategories-1] i.e. rate category
	double*  L = malloc_decimals(numUniqSites * num_chars);
	
	// for the first rate category
	for (i=0; i<numUniqSites; i++) {
		// mul_decimals(&L[i*4], &siteLikelihoods[i*num_chars], &variables.rootNodeFreq_alpha[0], num_chars);
		mul_decimals(&L[i*4], &siteLikelihoods[i*num_chars], &variables.rootNodeFreq[0], variables.alpha[0], num_chars);
	}
	
	// for the rest rate categories
	for (k=1; k<numVarCategories; k++) {
		for (i=0; i<numUniqSites; i++) {
			// mul_add_decimals(&L[i*4], &siteLikelihoods[(k*numUniqSites+i)*num_chars], &variables.rootNodeFreq_alpha[k*num_chars], num_chars);
			mul_add_decimals(&L[i*4], &siteLikelihoods[(k*numUniqSites+i)*num_chars], &variables.rootNodeFreq[k*num_chars], variables.alpha[k], num_chars);
		}
	}

	// L[i] = L[i*4]+L[i*4+1]+L[i*4+2]+L[i*4+3] where 0 <= i < num_items
	special_op(L, numUniqSites);
	
	for (i=0; i<numUniqSites; i++)
		L[i] = log(L[i] + variables.invar_constant[alignment.isConstantSites[i]]);
	
	double  sum[] = {0,0,0,0};
	
	mul_sum_decimals(sum, L, alignment.numEachUniqueSites, numUniqSites);
	
	free_decimals(L);
	
	return sum[0] + sum[1] + sum[2] + sum[3];
}


void getLogLikelihoodDetails(vector<double*> allCondProbSet_t, int* topologyMatrix, int numLineTopMat, Alignment& alignment, VariableSet& variables, int num_chars) {

	int numInternalNodes = alignment.numSeqs-1;
	double*  likelihoodMat = malloc_decimals(numInternalNodes * num_chars);
	int numSeqs = alignment.numSeqs;
	int numSites = alignment.numSites;
	double currSiteLikelihood;
	for (int i=0; i<numSites; i++) {
		char* currSites = &(alignment.sites[i*numSeqs]);
		for (int k=0; k<(int)allCondProbSet_t.size(); k++) {
			computeLikelihoodMatrix(likelihoodMat, allCondProbSet_t[k], topologyMatrix, currSites, numLineTopMat, num_chars);
			currSiteLikelihood = 0.0;
			for (int j=0; j<num_chars; j++) {
				currSiteLikelihood += likelihoodMat[(numLineTopMat-1)*num_chars+j]*variables.rootNodeFreq[k*num_chars + j];
			}
			cerr << currSiteLikelihood << ",";
		}
		cerr << endl;
	}
	free_decimals(likelihoodMat);
}

// load the rateGrpFile and collect the list of rate groups
// Input: (1) rateGrpFile; (2) topology matrix; (3) leafList
// Output: array of the rate groups (i.e. array of integers) and number of rate group
vector<int>* collectRateGrp(char* rateGrpFile, int* topMatrix, vector<string>* leafList, int& numRateGrp) {
	int numEdges = leafList->size()*2 - 2;
	ifstream fin;
	fin.open(rateGrpFile);
	if (!fin.is_open()) {
		cerr << "Error opening the rate matrix file :" << rateGrpFile << endl;
		exit(1);
	}
	string rateMatrixStr = "";
	// ignore all the lines starting with "[" character
	while (getline(fin, rateMatrixStr)) {
		if (rateMatrixStr.length() > 0 && rateMatrixStr[0]!='[')
			break;
	}
	fin.close(); // close the file

	if (rateMatrixStr.length() == 0) {
		cerr << "Error! The rate matrix file " << rateGrpFile << " is empty" << endl;
		exit(1);
	}

	// get the arrangement of rate matrix
	vector<int>* rateMatArray = NULL;
	if (rateMatrixStr.length() > 0) {
		if (rateMatrixStr[0] == '(')
			rateMatArray = treeFormatToRateMatrix(topMatrix, rateMatrixStr, leafList);
		else {
			rateMatArray = new vector<int>;
			vector<string> token;
			tokenizer(rateMatrixStr, ", ", &token);
			for (int i=0; i<(int)token.size(); i++)
				rateMatArray->push_back(atoi(token[i].c_str()));
		}
	}

	if (rateMatArray == NULL || rateMatArray->size() == 0 || rateMatArray->size() != numEdges) {
		cerr << "Error! The arrangement of rate matrices was loaded unsuccessfully" << endl;
		exit(1);
	}

	numRateGrp = 0;
	for (int i=0; i<(int)rateMatArray->size(); i++) {
		int rateGrp = rateMatArray->at(i);
		if (numRateGrp < rateGrp)
			numRateGrp = rateGrp;
	}

	return rateMatArray;
}



// initialize the thread locks
void initThreadLocks(ThreadLocks* threadLocks) {
	pthread_mutex_init(&(threadLocks->lock_cpuThresRecycle), NULL);
	pthread_mutex_init(&(threadLocks->lock_output), NULL);
	pthread_mutex_init(&(threadLocks->lock_optLK), NULL);
	pthread_mutex_init(&(threadLocks->lock_optInfo), NULL);
	pthread_mutex_init(&(threadLocks->lock_master), NULL);
	pthread_mutex_init(&(threadLocks->lock_maxMasterID), NULL);
	pthread_mutex_init(&(threadLocks->lock_createThread), NULL);
	pthread_mutex_init(&(threadLocks->lock_jobAllocate), NULL);
	pthread_mutex_init(&(threadLocks->lock_checkpoint), NULL);
}

// helping function for displaying the arrangement of rate matrix into a tree format
void subRateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, int line, string& str, vector<string>* leafList) {
	int node1 = topMatrix[line*2];
	int node2 = topMatrix[line*2+1];
	int rate1 = rateMatrix->at(line*2);
	int rate2 = rateMatrix->at(line*2+1);
	str.append("(");
	if (node1 < 0) {
		// a leaf
		if (leafList!=NULL)
			str.append(leafList->at(-node1-1) + "_");
		str.append(intToStr(rate1));
	} else {
		// an internal node
		subRateMatrixToTreeFormat(topMatrix,rateMatrix,node1-1, str, leafList);
		str.append(intToStr(rate1));
	}
	str.append(",");
	if (node2 < 0) {
		// a leaf
		if (leafList!=NULL)
			str.append(leafList->at(-node2-1) + "_");
		str.append(intToStr(rate2));
	} else {
		// an internal node
		subRateMatrixToTreeFormat(topMatrix,rateMatrix,node2-1, str, leafList);
		str.append(intToStr(rate2));
	}
	str.append(")");
}


// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, vector<string>* leafList) {
	int numLineTopMatrix=(int)rateMatrix->size()/2;
	string str = "";

	subRateMatrixToTreeFormat(topMatrix, rateMatrix, numLineTopMatrix-1, str, leafList);
	if (str.length() > 0)
		str.append(";");
	return str;    
}

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, string rateMatrixStr, vector<string>* leafList) {
	vector<int> rateMatrix;
	vector<string> token;
	tokenizer(rateMatrixStr, ", ", &token);
	for (int i=0; i<(int)token.size(); i++)
		rateMatrix.push_back(atoi(token[i].c_str()));
	int numLineTopMatrix=(int)rateMatrix.size()/2;
	string str = "";

	subRateMatrixToTreeFormat(topMatrix, &rateMatrix, numLineTopMatrix-1, str, leafList);
	if (str.length() > 0)
		str.append(";");
	return str;    
}

// helping function for displaying the arrangement of rate matrix into a tree format
void subRateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, int line, string& str, vector<string>* leafList, vector<int>* untouchNodes, int currNode) {
	int node1 = topMatrix[line*2];
	int node2 = topMatrix[line*2+1];
	int rate1 = rateMatrix->at(line*2);
	int rate2 = rateMatrix->at(line*2+1);
	int untouch1 = 0;
	int untouch2 = 0;
	if (node1 >= 0)
		untouch1 = untouchNodes->at(node1);
	if (node2 >= 0)
		untouch2 = untouchNodes->at(node2);
	str.append("(");
	if (node1 < 0) {
		// a leaf
		if (leafList!=NULL)
			str.append(leafList->at(-node1-1) + "_");
		str.append(intToStr(rate1));
	} else {
		// an internal node
		subRateMatrixToTreeFormat(topMatrix,rateMatrix,node1-1, str, leafList, untouchNodes, currNode);
		str.append(intToStr(rate1));
		str.append("[" + intToStr(untouch1) + "]");
		if (node1 == currNode)
			str.append("**");
	}
	str.append(",");
	if (node2 < 0) {
		// a leaf
		if (leafList!=NULL)
			str.append(leafList->at(-node2-1) + "_");
		str.append(intToStr(rate2));
	} else {
		// an internal node
		subRateMatrixToTreeFormat(topMatrix,rateMatrix,node2-1, str, leafList, untouchNodes, currNode);
		str.append(intToStr(rate2));
		str.append("[" + intToStr(untouch2) + "]");
		if (node2 == currNode)
			str.append("**");
	}
	str.append(")");
}

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, vector<int>* rateMatrix, vector<string>* leafList, vector<int>* untouchNodes, int currNode) {
	int numLineTopMatrix=(int)rateMatrix->size()/2;
	string str = "";
	
	subRateMatrixToTreeFormat(topMatrix, rateMatrix, numLineTopMatrix-1, str, leafList, untouchNodes, currNode);
	if (str.length() > 0) {
		if (currNode==0)
			str.append("**");
		str.append(";");
	}
	return str;
}

// display the arrangement of rate matrix into a tree format
string rateMatrixToTreeFormat(int* topMatrix, string rateMatrixStr, vector<string>* leafList, vector<int>* untouchNodes, int currNode) {
	vector<int> rateMatrix;
	vector<string> token;
	tokenizer(rateMatrixStr, ", ", &token);
	for (int i=0; i<(int)token.size(); i++)
		rateMatrix.push_back(atoi(token[i].c_str()));
	int numLineTopMatrix=(int)rateMatrix.size()/2;
	string str = "";
	
	subRateMatrixToTreeFormat(topMatrix, &rateMatrix, numLineTopMatrix-1, str, leafList, untouchNodes, currNode);
	if (str.length() > 0) {
		if (currNode==0)
			str.append("**");
		str.append(";");
	}
	return str;
}


// get the arrangement of rate matrix from the tree format
vector<int>* treeFormatToRateMatrix(int* topMatrix, string rateMatrixTree, vector<string>* leafList) {

	int numSpecies = (int)leafList->size();
	int numEdges = 2 * numSpecies - 2;

	int* theTopMatrix;
	int theLineTopMat;
	vector<string>* theLeafList;
	vector<int>* rateMatrix = new vector<int>;
	
	genTopMatrixFrStr(rateMatrixTree, &theTopMatrix, &theLineTopMat, &theLeafList, rateMatrix, NULL);

	int theNumEdges = theLineTopMat * 2;

	// check whether "topMatrix" and "theTopMatrix" are the same
	bool isSame = true;
	int i;
	if (numEdges != theNumEdges || numEdges==0) {
		isSame = false;
	} else {
		for (i=0; i<numEdges; i++) {
			if (theTopMatrix[i]!=topMatrix[i]) {
				isSame = false;
				break;
			}
		}
	}
	// check whether "leafList" and "theLeafList" are the same
	if (isSame) {
		for (i=0; i<numSpecies; i++) {
			if (leafList->at(i) != theLeafList->at(i)) {
				isSame = false;
				break;
			}
		}
	}
	if (!isSame) {
		cerr << "Error! The tree representing the arrangement of the rate matrix is not the same as the one in the topology file" << endl;
		cerr << rateMatrixTree << endl;
		exit(1);
	}
	return rateMatrix;
}

// return the number of rate groups of a given rate matrix
int getNumRateGrp(vector<int>& rateMatrix) {
	int num = 0;
	for (int i=0; i<(int)rateMatrix.size(); i++) {
		if (rateMatrix.at(i) > num)
			num = rateMatrix.at(i);
	}
	return num;
}

// read the HAS result file
// output: 1. number of site categories; 2. HALHAS model number
// not yet finished
void readHASResultfile(char* file, int& catNum, int& modelNum) {
	ifstream fin;
	fin.open(file);
	string aline;
	while (getline(fin, aline)) {
		if (aline.length() > 18 && aline.substr(0,18)=="HAL-HAS model type") {

		}
	}
}

//===========================================
// for loading the files used for simulation
//===========================================

// Input: topology string (can also include the rate group and the edge lengths)
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length, 4. List of node names
// Return: Topology matrix
int* genTopMatrixFrStr2(string newickStr, int& rowNum, vector<string>& leafList,
		vector<double>& edgeLens, vector<string>& nodeList) {

	trim(newickStr);
	if (newickStr.length()==0) {
		cerr << "Error! The topology string is empty" << endl;
		exit(1);
	}

	// get the number of rows in the matrix
	rowNum = howManyCharInside(newickStr, ',');

	if (rowNum == 0) {
		return NULL;
	}

	// clear all the arrays
	leafList.clear();
	edgeLens.clear();
	nodeList.clear();

	// allocate memory for the topMatrix
	int* topMatrix = new int[rowNum*2];

	// create another array to store the nodes waiting to be processed
	vector<int> waitNodeIDList;
	vector<string> waitNodeNameList;
	vector<double> waitNodeDist;

	// the number of new nodes formed
	int newNodeNum = 0;

	// current node's name and its edge length
	string currNodeName = "";
	string currDist = "";
	int distMode = 0;
	int isLeaf = 1;
	for (int i=0; i<(int)newickStr.length(); i++) {
		char c = newickStr.at(i);
		if (c == ';')
			break;
		if (c == '(' || c == ' ') {
			continue; // do nothing
		} else if (c == ',' || c==')') {
			if (currNodeName=="") {
				// no node name is provided
				cerr << "There is missing node name in the tree file" << endl;
				exit(1);
			}
			if (currDist=="") {
				// no edge length information is provided
				cerr << "No edge length information is provided for the node " << currNodeName << endl;
				exit(1);
			}
			// add the current node name to waitNodeNameList
			waitNodeNameList.push_back(currNodeName);
			// add the current edge length to waitNodeDist
			waitNodeDist.push_back((double) atof(currDist.c_str()));
			if (isLeaf) {
				// add the current leaf to leafList
				leafList.push_back(currNodeName);
				// add the current node ID to waitNodeIDList
				waitNodeIDList.push_back(-1 * (int)(leafList.size()));
			} else {
				// add the current node ID to waitNodeIDList
				waitNodeIDList.push_back(++newNodeNum);
			}
			// the following is applied for ')' only
			if (c == ')') {
				// add a new row into the topMatrix
				int lastID = waitNodeIDList.back();
				waitNodeIDList.pop_back();
				topMatrix[newNodeNum*2] = waitNodeIDList.back();
				waitNodeIDList.pop_back();
				topMatrix[newNodeNum*2+1] = lastID;
				// add the corresponding edge lengths into the edgeLens
				double lastDist = waitNodeDist.back();
				waitNodeDist.pop_back();
				edgeLens.push_back(waitNodeDist.back());
				waitNodeDist.pop_back();
				edgeLens.push_back(lastDist);
				// add the corresponding node name into the nodeList
				string lastName = waitNodeNameList.back();
				waitNodeNameList.pop_back();
				nodeList.push_back(waitNodeNameList.back());
				waitNodeNameList.pop_back();
				nodeList.push_back(lastName);
			}
			// reset the variables
			currNodeName = "";
			currDist = "";
			distMode = 0;
			if (c==')')
				isLeaf = 0;
			else
				isLeaf = 1;
		} else if (c == ':') {
			distMode = 1;
		} else if (distMode==1) {
			currDist.append(1,c);
		} else {
			currNodeName.append(1,c);
		}
	}
	return topMatrix;
}

// for the tree file using for simulation
// Input: Trees file
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length (# of items = # of edges * # of site categories)
//         4. List of node names, 5. List of category
// return: Topology matrix
int* genTopMatrix2(char* treesFile, int& rowNum, vector<string>& leafList, vector<double>& edgeLens,
		vector<string>& nodeList, vector<string>& siteCatNames) {

	// read the trees file
	ifstream fin;
	fin.open(treesFile); // open the file
	if (!fin.is_open()) {
		cerr << "Error opening the topology file :" << treesFile << endl;
		exit(1);
	}
	string aline;
	int* firstTopMatrix = NULL;
	int* currTopMatrix = NULL;
	int currRowNum;
	vector<string> currLeafList;
	vector<double> currEdgeLens;
	vector<string> currNodeList;
	bool sameTree = true;
	int i;
	vector<string> token;

	while (getline(fin, aline) && sameTree) {
		if (aline.length()>0 && aline[0]=='#')
			continue;
		tokenizer(aline, " \t", &token);
		if (token.size() >= 2) {
			if (siteCatNames.size()==0) {
				firstTopMatrix = genTopMatrixFrStr2(token[1], rowNum, leafList, edgeLens, nodeList);                
				if (rowNum>0) {
					siteCatNames.push_back(token[0]);
				}
			} else {
				currTopMatrix = genTopMatrixFrStr2(token[1], currRowNum, currLeafList, currEdgeLens, currNodeList);
				if (currRowNum>0) {
					// check whether all the tree topologies are the same
					if (currRowNum != rowNum) {
						sameTree = false;
						break;
					}
					// topMatrix vs currTopMatrix
					for (i=0; i<currRowNum*2; i++) {
						if (firstTopMatrix[i] != currTopMatrix[i]) {
							sameTree = false;
							break;
						}
					}
					// leafList vs currLeafList
					if (sameTree) {
						for (i=0; i<currRowNum+1; i++) {
							if (leafList[i] != currLeafList[i]) {
								sameTree = false;
								break;
							}
						}
					}
					// nodeList vs currNodeList
					if (sameTree) {
						for (i=0; i<currRowNum*2; i++) {
							if (nodeList[i] != currNodeList[i]) {
								sameTree = false;
								break;
							}
						}
					}
					// if the topology is not the same, then give out error message
					if (!sameTree) {
						cerr << "Error! The topoloies are not the same inside the tree file" << endl;
						exit(1);
					}
					// append the new items in currEdgeLens to edgeLens
					for (i=0; i<currRowNum*2; i++) {
						edgeLens.push_back(currEdgeLens[i]);
					}
					siteCatNames.push_back(token[0]);
					// reset the variables
					delete[] currTopMatrix;
					currTopMatrix = NULL;
				}
			}
		}
	}
	fin.close(); // close the file
	return firstTopMatrix;
}


//============================================
// to check how many changes on the rate matrices (convergent or non-convergent)
// across the tree
// this function is invoked recursively, from the last time of the topMatrix to
// the first line of the topMatrix
//============================================
int numChanges(int* topMatrix, vector<int>* rateMatrix, int line, int parentRate) {
// int numChanges(int* topMatrix, vector<int>* rateMatrix, int line, int parentRate, vector<string>* leafList, string& str) {

	int rate1 = rateMatrix->at(line*2);
	int rate2 = rateMatrix->at(line*2+1);
	int node1 = topMatrix[line*2];
	int node2 = topMatrix[line*2+1];
	int currNumChanges = 0;

	if (rate1==rate2) {
		if (rate1 != parentRate && parentRate != 0)
			currNumChanges++;
	} else {
		if (parentRate==0)
			currNumChanges++;
		else {
			if (rate1 != parentRate)
				currNumChanges++;
			if (rate2 != parentRate)
				currNumChanges++;
		}
	}
	
	if (node1 > 0) {
		// an internal node
		currNumChanges += numChanges(topMatrix,rateMatrix,node1-1, rate1);
	}
	
	if (node2 > 0 ) {
		// an internal node
		currNumChanges += numChanges(topMatrix,rateMatrix,node2-1, rate2);
	}
	
	return currNumChanges;
}
