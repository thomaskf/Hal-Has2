/*
 *
 * RAL_common.cpp
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

#include "RAL_common.h"

// ============================================================================================================================
// FUNCTIONS of rate matrix generation
// ============================================================================================================================

int* getEdgeOrder(ParameterSet* ps, int unRooted) {
	// this function is to compute the order of edges according to the edges length
	// when all the parameters are optimized according to the simplest model
	// (i.e. all edges are assigned to the same rate group)
	// The order of the edges is based on the decreasing order of the corresponding length
	
	// for unrooted tree, the root edges have the lower priority
	
	int i;
	double* curr_w;
	double* curr_pi;
	double base;
	double curr_len;
	multimap<double, int> edge_lens;
	multimap<double, int>::reverse_iterator rit;
	int* edgeOrder = NULL;
	
	if (ps->num_edges==0)
		return edgeOrder;
	
	for (i=0; i<ps->num_edges; i++) {
		if (unRooted && (i==ps->num_edges-2 || i==ps->num_edges-1)) {
			curr_len = -1;
		} else {
			curr_w = &(ps->w[i*6]);
			curr_pi = &(ps->pi[i*4]);
			base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
			curr_len = ps->t[i] * base;
		}
		edge_lens.insert(pair<double,int>(curr_len, i));
	}
	
	edgeOrder = (int*) malloc(sizeof(int)*ps->num_edges);
	
	i=0;
	for (rit=edge_lens.rbegin(); rit!=edge_lens.rend(); rit++) {
		edgeOrder[i] = rit->second;
		i++;
	}
	
	return edgeOrder;
}

int generateRateMat(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp, int unRooted) {

	return generateRateMat(outRateMatrix, inRateMatrix, ps, vs, IC, topMatrix, leafID2line, internalID2line,
						   untouchNodes, numSpecies, numRateGrp, NULL, NULL, unRooted);

}

int generateRateMat(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp,
					vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted) {
	
	// generate a set of rate matrices
	//
	// input: topMatrix - topology matrix
	//        leafID2line - leafID to line number in the parameter file
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        untouchNodes - the list of nodes which cannot be processed (including root)
	//                       rootID = 0
	//        numSpecies - number of species
	//        numRateGrp - number of rate groups
	//        masterRateMatrix - all the rate matrices which have come across before
	// output:
	//      a set of rate matrices together with the corresponding internal node ID
	//      each rate matrix is represented by an integer array (i.e. vector<int>)
	//      Note: the nodeIDs and the groupIDs are 1-based
	// Report:
	//      the number of new rate matrices generated
	
	
	int new_rate_matrices = 0;
	int insertToSet = 1;
	int i;
	
	// for a rooted binary tree, number of internal nodes (including root) = number of species - 1;
	int numInternalNodes = numSpecies-1;
	int numEdges = 2*numSpecies-2;
	
	vector<int> currMatrix;
	vector<int> leftNodes, rightNodes;
	for (int node=0; node<numInternalNodes; node++) {
		
		if (untouchNodes->at(node)==1) {
			// the node is already inside the list of untouchNodes
			continue;
		}
		
		// reset leftNodes and rightNodes
		leftNodes.clear(); rightNodes.clear();
		// get all the nodes inside the subtree rooted at the left child until the node appear in the untouchNodes
		collectNodes(leftNodes, node, topMatrix, numInternalNodes, 1, untouchNodes);
		// get all the nodes inside the subtree rooted at the right child until the node appear in the untouchNodes
		collectNodes(rightNodes, node, topMatrix, numInternalNodes, 2, untouchNodes);
		// set the current node inside "untouchNodes" to 1
		// untouchNodes->at(node)=1;
		
		// -------------------------------------------------------------------
		// assign a new rate group for the nodes inside the left subtree
		// -------------------------------------------------------------------
		// set the content of currMatrix as inRateMatrix
		currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		// assign the new rate group for the nodes inside the left subtree
		assignGrpID(currMatrix, leftNodes, numRateGrp+1, leafID2line, internalID2line);
		// refine the rate matrix
		refineRateMatrix(currMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
		if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(currMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
					new_rate_matrices++;
				} else {
					for (i=0; i<(int) opt_modes->size(); i++) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
		
		// -------------------------------------------------------------------
		// assign a new rate group for the nodes inside the right subtree
		// -------------------------------------------------------------------
		// set the content of currMatrix as inRateMatrix
		currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		// assign the new rate group for the nodes inside the right subtree
		assignGrpID(currMatrix, rightNodes, numRateGrp+1, leafID2line, internalID2line);
		// refine the rate matrix
		refineRateMatrix(currMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
		if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(currMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
		
		// --------------------------------------------------------------------------------------
		// assign a new same rate group for the nodes inside both the left and the right subtree
		// --------------------------------------------------------------------------------------
		// set the content of currMatrix as inRateMatrix
		currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		// assign the new rate group for the nodes inside the left subtree
		assignGrpID(currMatrix, leftNodes, numRateGrp+1, leafID2line, internalID2line);
		// assign the new rate group for the nodes inside the right subtree
		assignGrpID(currMatrix, rightNodes, numRateGrp+1, leafID2line, internalID2line);
		// refine the rate matrix
		refineRateMatrix(currMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
		if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(currMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
		
		// ----------------------------------------------------------------------------------------------
		// assign two new different rate groups for the nodes inside both the left and the right subtree
		// ----------------------------------------------------------------------------------------------
		// set the content of currMatrix as inRateMatrix
		currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		// assign the new rate group for the nodes inside the left subtree
		assignGrpID(currMatrix, leftNodes, numRateGrp+1, leafID2line, internalID2line);
		// assign the new rate group for the nodes inside the right subtree
		assignGrpID(currMatrix, rightNodes, numRateGrp+2, leafID2line, internalID2line);
		// refine the rate matrix
		// refineRateMatrix(currMatrix, numRateGrp+2);
		refineRateMatrix(currMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
		if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(currMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
	}
	
	return new_rate_matrices;
}

int generateRateMatJohn(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
						int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp, int unRooted) {
	
	return generateRateMatJohn(outRateMatrix, inRateMatrix, ps, vs, IC, topMatrix,
						leafID2line, internalID2line, untouchNodes, numSpecies, numRateGrp, NULL, NULL, unRooted);
}

int generateRateMatJohn(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
						int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp,
						vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted) {
	
	// generate a set of rate matrices
	//
	// input: topMatrix - topology matrix
	//        leafID2line - leafID to line number in the parameter file
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        untouchNodes - the list of nodes which cannot be processed (including root)
	//                       rootID = 0
	//        numSpecies - number of species
	//        numRateGrp - number of rate groups
	//        masterRateMatrix - all the rate matrices which have come across before
	// output:
	//      a set of rate matrices together with the corresponding internal node ID
	//      each rate matrix is represented by an integer array (i.e. vector<int>)
	//      Note: the nodeIDs and the groupIDs are 1-based
	// Report:
	//      the number of new rate matrices generated
	
	
	int new_rate_matrices = 0;
	int insertToSet = 1;
	int i;
	
	// for a rooted binary tree, number of internal nodes (including root) = number of species - 1;
	int numInternalNodes = numSpecies-1;
	int numEdges = 2*numSpecies-2;
	int rateGrp, rateGrp2;
	
	vector<int> currMatrix;
	vector<int> leftNodes, rightNodes;
	
	// print out the input rate matrix arrangement
	// cout << "In Rate Matrix:" << endl;
	// cout << rateMatrixToTreeFormat(topMatrix, inRateMatrix, leafList, untouchNodes, -1) << endl;
	
	for (int node=0; node<numInternalNodes; node++) {
		
		if (untouchNodes->at(node)==1) {
			// the node is already inside the list of untouchNodes
			continue;
		}
		
		// reset leftNodes and rightNodes
		leftNodes.clear(); rightNodes.clear();
		// get all the nodes inside the subtree rooted at the left child until the node appear in the untouchNodes
		collectNodes(leftNodes, node, topMatrix, numInternalNodes, 1, untouchNodes);
		// get all the nodes inside the subtree rooted at the right child until the node appear in the untouchNodes
		collectNodes(rightNodes, node, topMatrix, numInternalNodes, 2, untouchNodes);
		// set the current node inside "untouchNodes" to 1
		// untouchNodes->at(node)=1;
		
		// -------------------------------------------------------------------
		// generate a set of new models W2 by
		// assigning a set of rate groups for the nodes inside the left subtree
		// -------------------------------------------------------------------
		
		// cout << "W2:" << endl;
		for (rateGrp=1; rateGrp<=numRateGrp+1; rateGrp++) {
			// set the content of currMatrix as inRateMatrix
			currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
			// assign the new rate group for the nodes inside the left subtree
			assignGrpID(currMatrix, leftNodes, rateGrp, leafID2line, internalID2line);
			// refine the rate matrix
			refineRateMatrix(currMatrix);
			// for unrooted tree, the rate matrices of the root edges have to be the same
			if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
				if (!outRateMatrix->isInside(currMatrix)) {
					if (opt_modes==NULL) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
						// cout << rateMatrixToTreeFormat(topMatrix, &currMatrix, leafList, untouchNodes, node) << endl;
						new_rate_matrices++;
					} else {
						for (i=0; i<(int)opt_modes->size(); i++) {
							outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
							new_rate_matrices++;
						}
					}
				}
			}
		}
		
		// -------------------------------------------------------------------
		// generate a set of new modes W3 by
		// assigning a set of rate groups for the nodes inside the right subtree
		// -------------------------------------------------------------------
		
		// cout << "W3:" << endl;
		for (rateGrp=1; rateGrp<=numRateGrp+1; rateGrp++) {
			// set the content of currMatrix as inRateMatrix
			currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
			// assign the new rate group for the nodes inside the right subtree
			assignGrpID(currMatrix, rightNodes, rateGrp, leafID2line, internalID2line);
			// refine the rate matrix
			refineRateMatrix(currMatrix);
			// for unrooted tree, the rate matrices of the root edges have to be the same
			if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
				if (!outRateMatrix->isInside(currMatrix)) {
					if (opt_modes==NULL) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
						// cout << rateMatrixToTreeFormat(topMatrix, &currMatrix, leafList, untouchNodes, node) << endl;
						new_rate_matrices++;
					} else {
						for (i=0; i<(int)opt_modes->size(); i++) {
							outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
							new_rate_matrices++;
						}
					}
				}
			}
		}
		
		// --------------------------------------------------------------------------------------
		// generate a set of new modes W1 by
		// assigning a set of rate groups for the nodes inside both the left and the right subtree
		// --------------------------------------------------------------------------------------

		// cout << "W1:" << endl;
		for (rateGrp=1; rateGrp<=numRateGrp+1; rateGrp++) {
			// set the content of currMatrix as inRateMatrix
			currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
			// assign the new rate group for the nodes inside the left subtree
			assignGrpID(currMatrix, leftNodes, rateGrp, leafID2line, internalID2line);
			// assign the new rate group for the nodes inside the right subtree
			assignGrpID(currMatrix, rightNodes, rateGrp, leafID2line, internalID2line);
			// refine the rate matrix
			refineRateMatrix(currMatrix);
			if (!outRateMatrix->isInside(currMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
					// cout << rateMatrixToTreeFormat(topMatrix, &currMatrix, leafList, untouchNodes, node) << endl;
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
		
		// ----------------------------------------------------------------------------------------------
		// generate a set of new modes W4 by
		// assigning a set of two different rate groups for the nodes inside both the left and the right subtree
		// ----------------------------------------------------------------------------------------------
		
		// cout << "W4:" << endl;
		for (rateGrp=1; rateGrp<=numRateGrp+1; rateGrp++) {
			for (rateGrp2=1; rateGrp2<=numRateGrp+1; rateGrp2++) {
				if (rateGrp==numRateGrp+1 && rateGrp2==numRateGrp+1)
					rateGrp2 = numRateGrp+2;
				if (rateGrp == rateGrp2)
					continue;
		
				// set the content of currMatrix as inRateMatrix
				currMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
				// assign the new rate group for the nodes inside the left subtree
				assignGrpID(currMatrix, leftNodes, rateGrp, leafID2line, internalID2line);
				// assign the new rate group for the nodes inside the right subtree
				assignGrpID(currMatrix, rightNodes, rateGrp2, leafID2line, internalID2line);
				// refine the rate matrix
				refineRateMatrix(currMatrix);
				// for unrooted tree, the rate matrices of the root edges have to be the same
				if (!unRooted || currMatrix[numEdges-1]==currMatrix[numEdges-2]) {
					if (!outRateMatrix->isInside(currMatrix)) {
						if (opt_modes==NULL) {
							outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet);
							// cout << rateMatrixToTreeFormat(topMatrix, &currMatrix, leafList, untouchNodes, node) << endl;
							new_rate_matrices++;
						} else {
							for (i=0; i<(int)opt_modes->size(); i++) {
								outRateMatrix->insert(currMatrix, *untouchNodes, IC, node, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
								new_rate_matrices++;
							}
						}
					}
				}
			}
		}
	}
	
	return new_rate_matrices;
}

int generateRateMatTD(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					  int* internalID2line, vector<int>* untouchNodes, int numRateGrp, double thresHeight, int unRooted) {
	
	return generateRateMatTD(outRateMatrix, inRateMatrix, ps, vs, IC, topMatrix,
					  internalID2line, untouchNodes, numRateGrp, thresHeight,
					  NULL, NULL, unRooted);
	
}


int generateRateMatTD(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					  int* internalID2line, vector<int>* untouchNodes, int numRateGrp, double thresHeight,
					  vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted) {
	// generate a set of rate matrices using top-down approach
	//
	// input: inRateMatrix - the current rate matrix
	//        ps - a parameter set
	//        topMatrix - topology matrix
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        numRateGrp - number of rate groups
	//        masterRateMatrix - all the rate matrices which have come across before
	//        thresHeight
	// output:
	//      a set of rate matrices together with the corresponding internal node ID
	//      each rate matrix is represented by an integer array (i.e. vector<int>)
	//      Note: the nodeIDs and the groupIDs are 1-based
	// Report:
	//      the number of new rate matrices generated
	
	
	int new_rate_matrices = 0;
	int numEdges = inRateMatrix->size();
	
	int NO_COMBINE_ALGO = 1;
	int insertToSet = 1;
	
	if (numRateGrp<=1)
		return 0;
	
	int i,j;
	ParameterSet currPS(ps->num_chars);
	currPS.copyFrom(*ps);
	currPS.normalize(); // need to normlize the parameters
	
	// =================================================
	// combine the rate group according to Pi value
	// =================================================
	
	int* uniqueGrp = new int[numRateGrp];
	for (i=0; i<numRateGrp; i++)
		uniqueGrp[i] = -1;
	
	double* uniquePi = new double[numRateGrp * 4];
	for (i=0; i<(int)inRateMatrix->size(); i++) {
		int rateGrp = inRateMatrix->at(i) - 1;
		if (uniqueGrp[rateGrp]==-1) {
			uniqueGrp[rateGrp]=i;
			for (j=0; j<4; j++) {
				uniquePi[rateGrp*4+j] = currPS.pi[i*4+j];
			}
		}
	}
	
	// compute the distance matrix for the pi matrix
	double* distMat = new double[numRateGrp*numRateGrp];
	computeDistMatrix(distMat, uniquePi, numRateGrp, 4);
	
	
	int* merge = new int[(numRateGrp-1)*2];
	double* height = new double[numRateGrp-1];
	hierCluster(merge, height, distMat, numRateGrp);
	
	// reduce the number of rate matrices by one
	vector<int> newRateMatrix;
	newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
	int frRateGrp = -merge[0];
	int toRateGrp = -merge[1];
	changeRateMatrix(newRateMatrix, frRateGrp, toRateGrp);
	refineRateMatrix(newRateMatrix);
	
	// for unrooted tree, the rate matrices of the root edges have to be the same
	if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
		if (!outRateMatrix->isInside(newRateMatrix)) {
			if (opt_modes==NULL) {
				outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet);
				// cout << rateMatrixToTreeFormat(topMatrix, &cnewRateMatrix, leafList, untouchNodes, -1) << endl;
				new_rate_matrices++;
			} else {
				for (i=0; i<(int)opt_modes->size(); i++) {
					outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
					new_rate_matrices++;
				}
			}
		}
	}
	
	if (numRateGrp >= 3 && (height[1]-height[0])<thresHeight) {
		// further reduce the number of rate matrices by one
		newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		changeRateMatrix(newRateMatrix, frRateGrp, toRateGrp);
		if (merge[2] < 0)
			frRateGrp = -merge[2];
		if (merge[3] < 0)
			toRateGrp = -merge[3];
		changeRateMatrix(newRateMatrix, frRateGrp, toRateGrp);
		refineRateMatrix(newRateMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
	if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(newRateMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet);
					// cout << rateMatrixToTreeFormat(topMatrix, &cnewRateMatrix, leafList, untouchNodes, -1) << endl;
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
	}
	
	// =====================================================
	// combine the rate group according to the shortest edge
	// =====================================================
	
	
	int smallPGrp = -1;
	int smallCGrp = -1;
	double smallEdgeLen = 1.0;
	
	for (i=0; i<(int)inRateMatrix->size()-2; i++) {
		// only consider the edges NOT connecting to the root
		
		int parent = i/2+1;
		int parentLine = internalID2line[parent];
		int childGrp = inRateMatrix->at(i);
		int parentGrp = inRateMatrix->at(parentLine);
		
		if (childGrp != parentGrp) {
			if (smallEdgeLen > currPS.t[i]) {
				smallEdgeLen = currPS.t[i];
				smallPGrp = parentGrp;
				smallCGrp = childGrp;
			}
		}
	}
	
	if (smallPGrp != -1) {
		newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		changeRateMatrix(newRateMatrix, smallCGrp, smallPGrp);
		refineRateMatrix(newRateMatrix);
		// for unrooted tree, the rate matrices of the root edges have to be the same
		if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
			if (!outRateMatrix->isInside(newRateMatrix)) {
				if (opt_modes==NULL) {
					outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet);
					// cout << rateMatrixToTreeFormat(topMatrix, &cnewRateMatrix, leafList, untouchNodes, -1) << endl;
					new_rate_matrices++;
				} else {
					for (i=0; i<(int)opt_modes->size(); i++) {
						outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
						new_rate_matrices++;
					}
				}
			}
		}
	}
	
	
	// =========================================================
	// combine the rate group according to BOTH S and Pi values
	// =========================================================
	
	if (numRateGrp > 2 && NO_COMBINE_ALGO==0) {
		
		int rowLen = 4 + currPS.num_w;
		double* uniqueSPi = new double[numRateGrp*rowLen];
		for (i=0; i<numRateGrp; i++) {
			int lineNo = uniqueGrp[i];
			// Pi
			for (j=0; j<4; j++) {
				uniqueSPi[i*rowLen+j] = currPS.pi[lineNo*4+j];
			}
			// S
			for (j=0; j<currPS.num_w; j++) {
				uniqueSPi[i*rowLen+(4+j)] = currPS.w[lineNo*currPS.num_w+j];
			}
		}
		
		// compute the distance matrix for the pi matrix
		computeDistMatrix(distMat, uniqueSPi, numRateGrp, rowLen);
		
		// find the two smallest distances
		double smallDists[2];
		int smallRows[2];
		int smallCols[2];
		
		if (distMat[1] < distMat[2]) {
			smallDists[0]=distMat[1]; smallDists[1]=distMat[2];
			smallCols[0]=1; smallCols[1]=2;
		} else {
			smallDists[0]=distMat[2]; smallDists[1]=distMat[1];
			smallCols[0]=2; smallCols[1]=1;
		}
		smallRows[0]=smallRows[1]=0;
		
		for (i=0; i<numRateGrp-1; i++) {
			for (j=i+1; j<numRateGrp; j++) {
				double currDist = distMat[i*numRateGrp+j];
				if (currDist < smallDists[1]) {
					if (currDist < smallDists[0]) {
						smallDists[1] = smallDists[0]; smallRows[1] = smallRows[0]; smallCols[1] = smallCols[0];
						smallDists[0] = currDist; smallRows[0] = i; smallCols[0] = j;
					} else {
						smallDists[1] = currDist; smallRows[1] = i; smallCols[1] = j;
					}
				}
			}
		}
		
		// reduce the number of rate matrices by one
		for (i=0; i<2; i++) {
			newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
			frRateGrp = smallRows[i]+1;
			toRateGrp = smallCols[i]+1;
			changeRateMatrix(newRateMatrix, frRateGrp, toRateGrp);
			refineRateMatrix(newRateMatrix);
			// for unrooted tree, the rate matrices of the root edges have to be the same
			if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
				if (!outRateMatrix->isInside(newRateMatrix)) {
					if (opt_modes==NULL) {
						outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet);
						// cout << rateMatrixToTreeFormat(topMatrix, &cnewRateMatrix, leafList, untouchNodes, -1) << endl;
						new_rate_matrices++;
					} else {
						for (i=0; i<(int)opt_modes->size(); i++) {
							outRateMatrix->insert(newRateMatrix, *untouchNodes, IC, -1, insertToSet, opt_modes->at(i), opt_maxITs->at(i));
							new_rate_matrices++;
						}
					}
				}
			}
		}
		
		delete[] uniqueSPi;
		
	}
	
	
	// release the memory
	delete[] uniqueGrp;
	delete[] uniquePi;
	delete[] distMat;
	delete[] merge;
	delete[] height;
	
	return new_rate_matrices;
}


int nextEdgeConsidered(ParameterSet* ps, int* topMatrix, vector<int>* untouchEdges, int numSpecies) {
	// according to the new bottom-up approach version 1
	// compute the next edges to be considered
	// the returned edge is 0-based
	// if no edge is available, then return -1
	
	int numInterNodes = numSpecies - 1; // number of internal nodes including root
	int numNodes = numInterNodes + numSpecies;
	int numEdges = numNodes - 1;
	
	if (numEdges<=0)
		return -1;
	
	int i;

	/*
	// for testing only
	for (i=0; i<numEdges; i++)
		if (untouchEdges->at(i)==0)
			return i;
	*/
	/*
	 // show the content of untouchEdges
	 cout << "untouchEdges : ";
	 for (i=0; i<numEdges; i++) {
	 if (i>0)
	 cout << ",";
	 cout << untouchEdges->at(i);
	 }
	 cout << endl;
	 */
	
	// check how many available edges
	int numAvailEdges = 0;
	for (i=0; i<numEdges; i++)
		if (untouchEdges->at(i)==0)
			numAvailEdges++;
	
	if (numAvailEdges==0)
		return -1;
	
	else if (numAvailEdges==1)
		for (i=0; i<numEdges; i++)
			if (untouchEdges->at(i)==0)
				return i;
	
	// update the parameters such that
	// Edge length is set to the rate of substitution; OR
	// ps->updateContent(1);
	// ps->showContent();
	
	// =============================================================
	// Computation of Mi
	// =============================================================
	// Given a tree T, let I be the set of nodes inside T
	// Define Mi as the length of the longest path from a
	//     node i \in I to a leaf in a subtree rooted at i
	//
	// For an internal node i, let u,v \in I be its children
	//     Mi = max { Mu + len(i->u) , Mv + len(i->v) }
	//     where len(i->u) and len(i->v) are the lengths
	//     of the edges from i to u, and from i to v, respectively
	//
	// For a leaf node i, Mi = 0
	// =============================================================
	double* M = new double[numNodes];
	// leaves: [0 ... numSpecies-1]; internal nodes: [ntrue.result.txtumSpecies ... numNodes-1]
	memset(M, 0, numNodes*sizeof(double));
	
	int node1, node2, currNode;
	double len1, len2;
	for (i=0; i<numEdges; i+=2) {
		node1 = topMatrix[i];
		if (node1 < 0) {
			node1 = -node1-1; // leaf
		} else {
			node1 = node1 + numSpecies - 1; // internal node
		}
		len1 = M[node1] + ps->t[i];
		node2 = topMatrix[i+1];
		if (node2 < 0) {
			// leaf
			node2 = -node2-1;
		} else {
			node2 = node2 + numSpecies - 1;
		}
		len2 = M[node2] + ps->t[i+1];
		currNode = i/2 + numSpecies;
		if (len2 > len1)
			M[currNode] = len2;
		else
			M[currNode] = len1;
	}
	
	/*
	 // show the values of M
	 cout << "M : ";
	 for (i=0; i<numNodes; i++) {
	 if (i>0)
	 cout << " ";
	 cout << M[i];
	 }
	 cout << endl;
	 */
	
	
	// =============================================================
	// Computation of Wi
	// =============================================================
	// Define Wi as the length of the path from the root of T
	//     to the node i \in I
	//
	// For a node i which is not the root, let p be its parent
	//     Wi = Wp + len(p->i)
	//
	// For a root i, Wi = 0
	// =============================================================
	
	double* W = new double[numNodes];
	// leaves: [0 ... numSpecies-1]; internal nodes: [ntrue.result.txtumSpecies ... numNodes-1]
	memset(W, 0, numNodes*sizeof(double));
	
	int node, parent;
	for (i=numEdges-1; i>=0; i--) {
		parent = i/2 + numSpecies;
		node = topMatrix[i];
		if (node < 0) {
			node = -node-1; // leaf
		} else {
			node = node + numSpecies - 1; // internal node
		}
		W[node] = W[parent] + ps->t[i];
	}
	
	/*
	 // show the values of W
	 cout << "W : ";
	 for (i=0; i<numNodes; i++) {
	 if (i>0)
	 cout << " ";
	 cout << W[i];
	 }
	 cout << endl;
	 */
	
	// =============================================================
	// Computation of Pe
	// =============================================================
	// Define Pe as the length of the longest path between leaves
	//     passing through an edge e
	//
	// Let a,b be the children of the root of T
	// Let Qa and Qb be the set of the nodes inside the subtree
	//     rooted at node a and b, respectively
	//
	// Consider an edge (x->y),
	//     if y \in Qa
	//           Pe = My + Wy + Mb + Wb
	//     else (i.e. y \in Qb)
	//           Pe = My + Wy + Ma + Wa
	// =============================================================
	
	
	double* P = new double[numEdges];
	memset(P, 0, numEdges*sizeof(double));
	
	int* grouping = new int[numNodes];
	memset(grouping, 0, numNodes*sizeof(int));
	
	int a,b;
	// node a
	a = topMatrix[numEdges-1];
	if (a < 0) {
		a = -a-1; // leaf
	} else {
		a = a + numSpecies - 1; // internal node
	}
	grouping[a] = 1;
	// node b
	b = topMatrix[numEdges-2];
	if (b < 0) {
		b = -b-1; // leaf
	} else {
		b = b + numSpecies - 1; // internal node
	}
	grouping[b] = 2;
	P[numEdges-1] = M[a] + W[a] + M[b] + W[b];
	P[numEdges-2] = M[a] + W[a] + M[b] + W[b];
	// other nodes
	for (i=numEdges-3; i>=0; i--) {
		node = topMatrix[i];
		if (node < 0) {
			node = -node-1; // leaf
		} else {
			node = node + numSpecies - 1; // internal node
		}
		parent = i/2 + numSpecies;
		grouping[node] = grouping[parent];
		if (grouping[node]==1) {
			P[i] = M[node] + W[node] + M[b] + W[b];
		} else {
			P[i] = M[node] + W[node] + M[a] + W[a];
		}
	}
	
	/*
	 // show the values of P
	 cout << "P : ";
	 for (i=0; i<numEdges; i++) {
	 if (i>0)
	 cout << " ";
	 cout << P[i];
	 }
	 cout << endl;
	 // show the values of Grouping
	 cout << "grouping : ";
	 for (i=0; i<numNodes; i++) {
	 if (i>0)
	 cout << " ";
	 cout << grouping[i];
	 }
	 cout << endl;
	 */
	
	// return the available one with largest value of Pe
	int resultEdge = -1;
	double resultPe = -1.0;
	for (i=0; i<numEdges; i++) {
		if (untouchEdges->at(i)==0 && isSame(P[i], resultPe)==0 && P[i]>resultPe) {
			resultEdge = i;
			resultPe = P[i];
		}
	}
	
	delete[] M;
	delete[] W;
	delete[] P;
	delete[] grouping;
	
	return resultEdge;
}


// new BU algorithm version 1
int generateRateMatNewBU1(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
						  int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, int* edgeOrder, int unRooted) {
	
	return generateRateMatNewBU1(outRateMatrix, inRateMatrix, ps, vs, logL, parentIC, IC,
								 topMatrix, leafID2line, internalID2line, untouchEdges, numSpecies,
								 numRateGrp, lowPossibleIC, NULL, NULL, edgeOrder, unRooted);
	
}


// new BU algorithm version 1
int generateRateMatNewBU1(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
						  int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int* edgeOrder, int unRooted) {
	// generate a set of rate matrices using NEW bottom-up approach version 1
	//
	// input: topMatrix - topology matrix
	//        leafID2line - leafID to line number in the parameter file
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        untouchNodes - the list of nodes which cannot be processed (including root)
	//                       rootID = 0
	//        numSpecies - number of species
	//        numRateGrp - number of rate groups
	//        masterRateMatrix - all the rate matrices which have come across before
	// output:
	// a set of rate matrices together with the corresponding internal node ID
	// each rate matrix is represented by an integer array (i.e. vector<int>)
	// Note: the nodeIDs and the groupIDs are 1-based
	
	int i,k;
	int new_rate_matrices = 0;
	/*
	cout << "ps's content:" << endl << flush;
	ps->showContent();
	int numEdges = 2*numSpecies-2;
	cout << "numEdges = " << numEdges << endl << flush;
	cout << "untouchEdges: " << flush;
	for (i=0; i<numEdges; i++)
		cout << " " << untouchEdges->at(i) << flush;
	*/
	int edge=-1;
	int numEdges = 2*numSpecies-2;
	
    if (edgeOrder == NULL) {
		edge = nextEdgeConsidered(ps, topMatrix, untouchEdges, numSpecies);
    } else {
		for (i=0; i<numEdges; i++) {
			if (untouchEdges->at(edgeOrder[i])==0) {
				edge = edgeOrder[i];
				break;
			}
		}
	}
	// cout << "Next edge to be considered : " << edge << endl << flush;
	// exit(1);
	
	vector<int> newRateMatrix;
	newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
	
	// insert the same rate matrix to the outRateMatrix
	int insertToSet = 0;
	if (opt_modes==NULL) {
        if (ps == NULL)
            outRateMatrix->insert(newRateMatrix, *untouchEdges, logL, parentIC, IC, edge, insertToSet, lowPossibleIC);
        else
            outRateMatrix->insert(newRateMatrix, *untouchEdges, logL, parentIC, *ps, *vs, IC, edge, insertToSet, lowPossibleIC, *ps, *vs);
		new_rate_matrices++;
	} else {
		for (i=0; i<(int)opt_modes->size(); i++) {
            if (ps == NULL)
                outRateMatrix->insert(newRateMatrix, *untouchEdges, logL, parentIC, IC, edge, insertToSet, lowPossibleIC, opt_modes->at(i), opt_maxITs->at(i));
            else
                outRateMatrix->insert(newRateMatrix, *untouchEdges, logL, parentIC, *ps, *vs, IC, edge, insertToSet, lowPossibleIC, *ps, *vs, opt_modes->at(i), opt_maxITs->at(i));
			new_rate_matrices++;
		}
	}
	insertToSet = 1;
	for (i=1; i<=numRateGrp+1; i++) {
		newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
		if (newRateMatrix[edge]!=i) {
			newRateMatrix[edge] = i;
			refineRateMatrix(newRateMatrix);
			// for unrooted tree, the rate matrices of the root edges have to be the same
			if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
				if (!outRateMatrix->isInside(newRateMatrix)) {
					// not appear before
					if (opt_modes==NULL) {
                        if (ps == NULL)
                            outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet);
                        else
                            outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet, *ps, *vs);
						new_rate_matrices++;
					} else {
						for (k=0; k<(int)opt_modes->size(); k++) {
                            if (ps == NULL)
                                outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet, opt_modes->at(k), opt_maxITs->at(k));
                            else
                                outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet, *ps, *vs, opt_modes->at(k), opt_maxITs->at(k));
							new_rate_matrices++;
						}
					}
				}
			}
		}
	}
	// note that the first one in the outRateMatrix is the same as inRateMatrix
	
	return new_rate_matrices;
}

// new BU algorithm version 2
int generateRateMatNewBU2(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
						  int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, int unRooted) {
	
	return generateRateMatNewBU2(outRateMatrix, inRateMatrix, ps, vs, logL, parentIC, IC, topMatrix, leafID2line, internalID2line, untouchEdges, numSpecies,
						  numRateGrp, lowPossibleIC, NULL, NULL, unRooted);
}

// new BU algorithm version 2
int generateRateMatNewBU2(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
						  int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted) {
	// generate a set of rate matrices using NEW bottom-up approach version 2
	//
	// input: topMatrix - topology matrix
	//        leafID2line - leafID to line number in the parameter file
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        untouchNodes - the list of nodes which cannot be processed (including root)
	//                       rootID = 0
	//        numSpecies - number of species
	//        numRateGrp - number of rate groups
	//        masterRateMatrix - all the rate matrices which have come across before
	// output:
	// a set of rate matrices together with the corresponding internal node ID
	// each rate matrix is represented by an integer array (i.e. vector<int>)
	// Note: the nodeIDs and the groupIDs are 1-based
	
	int numEdges = 2*numSpecies-2;

	
	// print out the input rate matrix arrangement
	// cout << "In Rate Matrix:" << endl;
	// cout << rateMatrixToTreeFormat(topMatrix, inRateMatrix, leafList, untouchEdges, -1) << endl;
	
	vector<int> newRateMatrix;
	newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
	
	int i,k;
	int new_rate_matrices = 0;
	int edge;
	int insertToSet = 1;
	// cout << "size of untouchEdges : " << untouchEdges->size() << endl << flush;
	for (edge=0; edge<numEdges; edge++) {
		// cout << "edge = " << edge << endl << flush;
		if (untouchEdges->at(edge)==1)
			continue;
		for (i=1; i<=numRateGrp+1; i++) {
			newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
			if (newRateMatrix[edge]!=i) {
				newRateMatrix[edge] = i;
				refineRateMatrix(newRateMatrix);
				// for unrooted tree, the rate matrices of the root edges have to be the same
				if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
					if (!outRateMatrix->isInside(newRateMatrix)) {
						if (opt_modes==NULL) {
							outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet, *ps, *vs);
							// cout << rateMatrixToTreeFormat(topMatrix, &cnewRateMatrix, leafList, untouchNodes, -1) << endl;
							new_rate_matrices++;
						} else {
							for (k=0; k<(int)opt_modes->size(); k++) {
								outRateMatrix->insert(newRateMatrix, *untouchEdges, IC, edge, insertToSet,  *ps, *vs, opt_modes->at(k), opt_maxITs->at(k));
								new_rate_matrices++;
							}
						}
					}
				}
			}
		}
	}
	// note that the first one in the outRateMatrix is the same as inRateMatrix
	// cout << "new_rate_matrices: " << new_rate_matrices << endl << flush;
	return new_rate_matrices;
}

// bute-force
int generateRateMatButeForce(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
							 int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
							 int maxNumRateGrp, double lowPossibleIC, int unRooted) {
	
	return generateRateMatButeForce(outRateMatrix, inRateMatrix, ps, vs, logL, parentIC, IC, topMatrix, leafID2line, internalID2line, untouchEdges, numSpecies,
									maxNumRateGrp, lowPossibleIC, NULL, NULL, unRooted);
	
}

// bute-force
int generateRateMatButeForce(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC, double IC,
							 int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
							 int maxNumRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted) {
	// generate a set of rate matrices using bute-force algorithm
	//
	// input: topMatrix - topology matrix
	//        leafID2line - leafID to line number in the parameter file
	//        internalID2line - internalID to line number in the parameter file
	//                          (the leafID2line and the internalID2line are 1-based)
	//        untouchNodes - the list of nodes which cannot be processed (including root)
	//                       rootID = 0
	//        numSpecies - number of species
	//        maxNumRateGrp - maximum number of rate groups allowed
	//        masterRateMatrix - all the rate matrices which have come across before
	// output:
	// a set of rate matrices together with the corresponding internal node ID
	// each rate matrix is represented by an integer array (i.e. vector<int>)
	// Note: the nodeIDs and the groupIDs are 1-based
	
	int new_rate_matrices = 0;
	int numEdges = 2*numSpecies-2;
	
	// int edge = nextEdgeConsidered(ps, topMatrix, untouchEdges, numSpecies);
	// cout << "Next edge to be considered : " << edge << endl;
	// exit(1);
	
	vector<int> newRateMatrix;
	newRateMatrix.assign(inRateMatrix->begin(), inRateMatrix->end());
	
	vector<int> newUntouchEdges(numEdges,1);
	
	// clear the outRateMatrix
	outRateMatrix->clear();
	
	// insert the same rate matrix to the outRateMatrix
	int insertToSet = 1;
	int i,j,k,edge,numItems;
	if (opt_modes==NULL) {
		outRateMatrix->insert(newRateMatrix, newUntouchEdges, logL, parentIC, *ps, *vs, IC, -1, insertToSet, lowPossibleIC);
		new_rate_matrices++;
	} else {
		for (i=0; i<(int)opt_modes->size(); i++) {
			outRateMatrix->insert(newRateMatrix, newUntouchEdges, logL, parentIC, *ps, *vs, IC, -1, insertToSet, lowPossibleIC, opt_modes->at(i), opt_maxITs->at(i));
			new_rate_matrices++;
		}
	}
	
	for (edge=0; edge<numEdges; edge++) {
		if (untouchEdges->at(edge)==1)
			continue;
		numItems = outRateMatrix->size();
		for (j=0; j<numItems; j++) {
			newRateMatrix.assign(outRateMatrix->matrices[j]->begin(), outRateMatrix->matrices[j]->end());
			for (i=1; i<=maxNumRateGrp; i++) {
				if (newRateMatrix[edge]!=i) {
					newRateMatrix[edge] = i;
					refineRateMatrix(newRateMatrix);
					// for unrooted tree, the rate matrices of the root edges have to be the same
					if (!unRooted || newRateMatrix[numEdges-1]==newRateMatrix[numEdges-2]) {
						if (!outRateMatrix->isInside(newRateMatrix)) {
							// not appear before
							if (opt_modes==NULL) {
								outRateMatrix->insert(newRateMatrix, newUntouchEdges, IC, -1, insertToSet);
								new_rate_matrices++;
							} else {
								for (k=0; k<(int)opt_modes->size(); k++) {
									outRateMatrix->insert(newRateMatrix, newUntouchEdges, IC, -1, insertToSet, opt_modes->at(k), opt_maxITs->at(k));
									new_rate_matrices++;
								}
							}
						}
					}
				}
			}
		}
	}
	
	/*
	 // show all items inside outRateMatrix
	 for (j=0; j<outRateMatrix->size(); j++) {
	 cout << "[";
	 for (i=0; i<(int)outRateMatrix->matrices[j]->size(); i++) {
	 if (i>0)
	 cout << ",";
	 cout << outRateMatrix->matrices[j]->at(i);
	 }
	 cout << "]" << endl;
	 }
	 exit(1);
	 */
	
	return new_rate_matrices;
}


// ============================================================================================================================
// OTHER FUNCTIONS
// ============================================================================================================================

void collectNodes(vector<int>& nodes, int currInternalNodeID, int* topMatrix, int numLineTopMat, int type,
				  vector<int>* untouchNodes) {
	// type: 1 - collect all nodes inside the subtree rooted at the left child
	//           until the node (also including that node) appear in the untouchNodes
	//       2 - collect all nodes inside the subtree rooted at the right child
	//           until the node (also including that node) appear in the untouchNodes
	//       3 - collect all nodes inside the subtree rooted at this node
	//           until the node (also including that node) appear in the untouchNodes
	
	if (type==3) {
		// put this node into the list
		nodes.push_back(currInternalNodeID);
	}
	if (currInternalNodeID >= 0 && untouchNodes->at(currInternalNodeID)==0) {
		// it is an internal node
		int row = currInternalNodeID-1;
		if (currInternalNodeID==0)
			row = numLineTopMat-1;
		if (type != 2) {
			int left = topMatrix[row*2];
			collectNodes(nodes,left,topMatrix,numLineTopMat,3,untouchNodes);
		}
		if (type != 1) {
			int right = topMatrix[row*2+1];
			collectNodes(nodes,right,topMatrix,numLineTopMat,3,untouchNodes);
		}
	}
}

void assignGrpID(vector<int>& currMatrix, vector<int>& nodes, int newGrp, int* leafID2line, int* internalID2line) {
	// assign the new rate group for the nodes
	
	// Note: the leafID2line and the internalID2line are 1-based
	//       the currMatrix is 0-based
	
	for (int i=0; i<(int)nodes.size(); i++) {
		if (nodes[i]>0) {
			currMatrix[internalID2line[nodes[i]]]=newGrp;
		} else if (nodes[i] < 0) {
			currMatrix[leafID2line[-nodes[i]]]=newGrp;
		}
	}
}

void changeRateMatrix(vector<int>& rateMatrix, int frRateGrp, int toRateGrp) {
	for (int i=0; i<(int)rateMatrix.size(); i++) {
		if (rateMatrix[i]==frRateGrp) {
			rateMatrix[i] = toRateGrp;
		}
	}
}

void refineRateMatrix(vector<int>& rateMatrix) {
	
	vector<int> mapping;
	
	int numRateGrp = 0;
	for (int i=0; i<(int)rateMatrix.size(); i++) {
		int rateGrp = rateMatrix[i];
		if ((int)mapping.size() <= rateGrp) {
			mapping.resize(rateGrp+1,-1);
		}
		if (mapping[rateGrp]==-1) {
			mapping[rateGrp]=++numRateGrp;
		}
		rateMatrix[i] = mapping[rateGrp];
	}
}

