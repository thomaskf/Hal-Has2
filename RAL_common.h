/*
 *
 * RAL_common.h
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

#ifndef __HAL_HAS__matrix_generation__
#define __HAL_HAS__matrix_generation__

#include <iostream>
#include <map>
#include "rateMatrixResult.h"
#include "parameters.h"

// ============================================================================================================================
// FUNCTIONS for rate matrices generations
// ============================================================================================================================

// generate a set of rate matrices using original bottom-up approach
int generateRateMat(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp, int unRooted);
int generateRateMat(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp,
					vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted);

// generate a set of rate matrices using top-down approach
int generateRateMatTD(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					  int* internalID2line, vector<int>* untouchNodes, int numRateGrp, double thresHeight, int unRooted);
int generateRateMatTD(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
					  int* internalID2line, vector<int>* untouchNodes, int numRateGrp, double thresHeight,
					  vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted);

// generate a set of rate matrices using NEW bottom-up approach version 1
int generateRateMatNewBU1(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
						  double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, int* edgeOrder, int unRooted);
int generateRateMatNewBU1(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
						  double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int* edgeOrder, int unRooted);

// generate a set of rate matrices using NEW bottom-up approach version 2
int generateRateMatNewBU2(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
						  double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, int unRooted);
int generateRateMatNewBU2(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
						  double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
						  int numRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted);

// bute-force
int generateRateMatButeForce(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
							 double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
							 int maxNumRateGrp, double lowPossibleIC, int unRooted);
int generateRateMatButeForce(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double logL, double parentIC,
							 double IC, int* topMatrix, int* leafID2line, int* internalID2line, vector<int>* untouchEdges, int numSpecies,
							 int maxNumRateGrp, double lowPossibleIC, vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted);

// generate a set of rate matrices using John's bottom-up approach
int generateRateMatJohn(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
						int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp, int unRooted);
int generateRateMatJohn(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
						int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp,
						vector<int>* opt_modes, vector<int>* opt_maxITs, int unRooted);


// this function is to compute the order of edges according to the edges length
// when all the parameters are optimized according to the simplest model
// (i.e. all edges are assigned to the same rate group)
// The order of the edges is based on the decreasing order of the corresponding length
int* getEdgeOrder(ParameterSet* ps, int unRooted);


//int generateRateMatJohn(TempRateMats* outRateMatrix, vector<int>* inRateMatrix, ParameterSet* ps, VariableSet* vs, double IC, int* topMatrix,
//                    int* leafID2line, int* internalID2line, vector<int>* untouchNodes, int numSpecies, int numRateGrp, vector<string>* leafList);

// ============================================================================================================================
// OTHER FUNCTIONS
// ============================================================================================================================

void collectNodes(vector<int>& nodes, int currInternalNodeID, int* topMatrix, int numLineTopMat, int type,
				  vector<int>* untouchNodes);
void assignGrpID(vector<int>& currMatrix, vector<int>& nodes, int newGrp, int* leafID2line, int* internalID2line);

void changeRateMatrix(vector<int>& rateMatrix, int frRateGrp, int toRateGrp);
void refineRateMatrix(vector<int>& rateMatrix);

#endif /* defined(__HAL_HAS__matrix_generation__) */
