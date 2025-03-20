/*
 *
 * optimization.cpp
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


#include "optimization.h"

Optim::Optim(VariableSet *vs, ParameterSet *ps, Alignment *algn, int *topMatrix, int isReversible, int num_chars) {
	this->vs = vs;
	this->ps = ps;
	this->algn = algn;
	this->topMatrix = topMatrix;
	this->isReversible = isReversible;
	this->num_chars = num_chars;
	this->modelType = -1;
	this->ps_set = NULL;
	this->isShare = 0;
	this->num_w = ps->num_w;
	this->num_edges = ps->num_edges;
	this->num_rateCat = vs->num_alpha;
	this->num_rateMat = (int) (ps->rateIDs).size();
	edgeScalars = NULL;
	likelihoodMat = NULL;
	leftLikelihoodMat = NULL;
	rightLikelihoodMat = NULL;
	//changedNodes = NULL;
	HASsiteLikelihoods = NULL;
	
}

// constructor
Optim::Optim(Alignment *algn, int *topMatrix, int isReversible, int num_chars) {
	this->vs = NULL;
	this->ps = NULL;
	this->algn = algn;
	this->topMatrix = topMatrix;
	this->isReversible = isReversible;
	this->num_chars = num_chars;
	this->modelType = -1;
	this->ps_set = NULL;
	this->isShare = 0;
	edgeScalars = NULL;
	likelihoodMat = NULL;
	leftLikelihoodMat = NULL;
	rightLikelihoodMat = NULL;
	//changedNodes = NULL;
	HASsiteLikelihoods = NULL;
}

// constructor
Optim::Optim(VariableSet *vs, AllParameterSet *ps_set, Alignment *algn, int *topMatrix, int isReversible, int modelType, int num_chars) {
	this->vs = vs;
	this->ps_set = ps_set;
	this->algn = algn;
	this->topMatrix = topMatrix;
	this->isReversible = isReversible;
	this->modelType = modelType;
	this->num_chars = num_chars;
	this->isShare = 0;

	this->num_w = ps_set->ps[0]->num_w;
	this->num_edges = ps_set->ps[0]->num_edges;
	this->num_rateCat = vs->num_alpha;
	this->num_rateMat = (int) (ps_set->ps[0]->rateIDs).size();

	if (modelType == 5 && num_rateCat > 1) {
		edgeScalars = (double*) malloc (sizeof(double) * (num_rateCat-1));
		for (int i=0; i<num_rateCat-1; i++) {
			edgeScalars[i] = 1.0;
		}
	} else {
		edgeScalars = NULL;
	}

	likelihoodMat = NULL;
	leftLikelihoodMat = NULL;
	rightLikelihoodMat = NULL;
	// changedNodes = NULL;
	HASsiteLikelihoods = NULL;
}

void Optim::setVarParam(VariableSet *vs, ParameterSet *ps) {
	this->vs = vs;
	this->ps = ps;
	this->num_w = ps->num_w;
	this->num_edges = ps->num_edges;
	this->num_rateCat = vs->num_alpha;
	this->num_rateMat = (int) (ps->rateIDs).size();
	this->num_chars = vs->num_chars;
}


void Optim::getOptValues(int type, int isShare, int maxit) {

	// optimizing the values
	// return the log-likelihood value
	// type -- 1 : all edge lengths; 2 : s-matrix along an edge; 3: pi-vector along an endge, 4: edgeScalars
	// type -- 5 : pi-vector as well as the root frequency (for GTR_upsilon model) 
	// type -- 6 : s-matrix and pi-vector along an edge
	// type -- 7 : edge-lengths, s-matrix, and pi-vector along all edges
	// type -- 8 : alpha, beta and nucleotide distribution along the constant sites
	// type -- 9 : nucleotide distribution along the root

	if (type==0)
		return;

	int i,j,k;

	paramType = type;
	this->isShare = isShare;

	int numValues = 0;
	if (type == 1)
		numValues = num_edges;
	else if (type == 2)
		numValues = num_w;
	else if (type == 3 || type == 5)
		numValues = num_chars;
	else if (type == 4) {
		if (num_rateCat <= 1) {
			cerr << "[getOptValues] Error! Number of rate categories <= 1" << endl;
			exit(1);
		}
		numValues = 1;
	} else if (type == 6) {
		numValues = num_w + num_chars;
	} else if (type == 7) {
		numValues = num_edges + num_rateMat * num_w + num_rateMat * num_chars;
	} else if (type == 8 && vs->isGTRUpilson) {
		numValues = num_rateCat;
	} else if (type == 8) {
#ifndef GTRIFO
		numValues = num_rateCat+1+num_chars;
#else
		numValues = num_rateCat+1;
#endif
	} else if (type == 9) {
		if (isShare)
			numValues = num_chars;
		else
			numValues = num_chars*num_rateCat;
	}

	double*  x = malloc_decimals(numValues);
	double*  l = malloc_decimals(numValues);

	ParameterSet* curr_ps;
	if (isShare) {
		curr_ps = ps_set->ps[0]; // refer to the parameters of the first rate category
	} else {
		curr_ps = ps;
	}

	/*
	cout << "[start] getOptValues type : " << type << " isShare : " << isShare << endl << flush;
	switch (type) {
		case 1:
			cout << "Edge" << endl;
			break;
		case 2:
			cout << "S-matrix" << endl;
			break;
		case 3:
			cout << "Pi-vector" << endl;
			break;
		case 4:
			cout << "edge scalar (N/A)" << endl;
			break;
	}
	 */

	if (type==1) {
		/* Initialize the variables */
		for (i=0; i<numValues; i++) {
			x[i] = curr_ps->t[i];
			l[i] = MIN_values;
		}
		lbfgsb_R(numValues, x, l, maxit, this);
	} else if (type==2) {
		for (k=0; k<num_rateMat; k++) {
			// for each rate matrix
			if ((curr_ps->rateIDs)[k]->size() > 0) {
				currRateMatID = k;
				/* Initialize the variables */
				for (i=0; i<numValues; i++) {
					x[i] = curr_ps->w[((curr_ps->rateIDs)[k]->at(0)) * num_w + i];
					l[i] = MIN_values;
				}
				lbfgsb_R(numValues, x, l, maxit, this);
			}
		}
	} else if (type==3 || type==5) {
		for (k=0; k<num_rateMat; k++) {
			// for each rate matrix
			if ((curr_ps->rateIDs)[k]->size() > 0) {
				currRateMatID = k;
				/* Initialize the variables */
				// cout << "x=";
				for (i=0; i<numValues; i++) {
					x[i] = curr_ps->pi_t[((curr_ps->rateIDs)[k]->at(0)) * num_chars + i];
					l[i] = MIN_PI;
					// cout << x[i] << ",";
				}
				// cout << endl << flush;
				lbfgsb_R(numValues, x, l, maxit, this);
			}
		}
	} else if (type==4) {
		for (k=1; k<num_rateCat; k++) {
			currRateCatID = k;
			ps = ps_set->ps[k];
			/* Initialize the variables */
			x[0] = edgeScalars[k-1];
			l[0] = MIN_values;
			lbfgsb_R(numValues, x, l, maxit, this);
		}
	} else if (type==6) {
		for (k=0; k<num_rateMat; k++) {
			// for each rate matrix
			if ((curr_ps->rateIDs)[k]->size() > 0) {
				currRateMatID = k;
				/* Initialize the variables */
				for (i=0; i<num_w; i++) {
					x[i] = curr_ps->w[((curr_ps->rateIDs)[k]->at(0)) * num_w + i];
					l[i] = MIN_values;
				}
				for (i=0; i<num_chars; i++) {
					x[i+num_w] = curr_ps->pi_t[((curr_ps->rateIDs)[k]->at(0)) * num_chars + i];
					l[i+num_w] = MIN_PI;
				}
				lbfgsb_R(numValues, x, l, maxit, this);
			}
		}
	} else if (type==7) {
		// for edge lengths, all s-matrix,and all PIs
		// num_edges + num_rateMat * num_w + num_rateMat * num_chars
		i=0;
		// edge lengths (# : num_edges)
		for (j=0; j<num_edges; j++) {
			x[i] = curr_ps->t[i];
			l[i] = MIN_values;
			i++;
		}
		// all s-matrix (# : num_rateMat * num_w)
		for (j=0; j<num_rateMat; j++) {
			// for each rate matrix
			for (k=0; k<num_w; k++) {
				x[i] = curr_ps->w[((curr_ps->rateIDs)[j]->at(0)) * num_w + k];
				l[i] = MIN_values;
				i++;
			}
		}
		// all PIs (# : num_rateMat * num_chars)
		for (j=0; j<num_rateMat; j++) {
			// for each rate matrix
			for (k=0; k<num_chars; k++) {
				x[i] = curr_ps->pi_t[((curr_ps->rateIDs)[j]->at(0)) * num_chars + k];
				l[i] = MIN_PI;
				i++;
			}
		}
		lbfgsb_R(numValues, x, l, maxit, this);
	} else if (type==8) {
		// for alpha, beta and nucleotide distribution along the constant sites
		for (i=0; i<num_rateCat; i++) {
			// for alpha
			x[i] = vs->alpha_t[i];
			l[i] = MIN_values;
		}
		for (i=num_rateCat; i<num_rateCat+1; i++) {
			// for beta
			x[i] = vs->beta_t;
			l[i] = MIN_values;
		}
#ifndef GTRIFO
		for (i=0; i<num_chars; i++) {
			// for nucleotide distribution along the constant sites
			x[num_rateCat+1+i] = vs->probXGivenInv_t[i];
			l[num_rateCat+1+i] = MIN_values;
		}
#endif
		lbfgsb_R(numValues, x, l, maxit, this);
	} else if (type==9) {
		// the nucleotide frequency along the root
		for (i=0; i<numValues; i++) {
			x[i] = vs->rootNodeFreq_t[i];
			l[i] = MIN_values;
		}
		lbfgsb_R(numValues, x, l, maxit, this);
	}

	free_decimals(x);
	free_decimals(l);
}

void Optim::getOptRootVector(int isShare, int maxit) {

	// optimize the root vector of each rate category
	int type = 9;
	getOptValues(type, isShare, maxit);

	/*
	int numRateCat = vs->num_alpha;
	int numSeqs = algn->numSeqs;
	int numUniqSites = algn->numUniqueSites;
	int numInternalNodes = numSeqs-1;
	int numLineTopMat = numSeqs-1;
	double* likelihoodMat = matrix(numInternalNodes, num_chars);
	double* temp = matrix(numRateCat, num_chars);
	double* updateMatrix = matrix(numRateCat, num_chars);
	resetMat(updateMatrix, numRateCat, num_chars);
	double* categorySum = matrix(1, numRateCat);
	double alphaDemon, updateMatrixSum;


	// optimize the root vector
	for (int i=0; i<numUniqSites; i++) {
		char* currSites = &(algn->uniqueSites[i*numSeqs]);

		for (int rateCat=0; rateCat<numRateCat; rateCat++) {

			if (ps_set != NULL)
				computeLikelihoodMatrix(likelihoodMat, ps_set->ps[rateCat]->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);
			else
				computeLikelihoodMatrix(likelihoodMat, ps->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);

			// multiply by the root vector
			for (int k=0; k<num_chars; k++) {
				temp[rateCat*num_chars + k] = likelihoodMat[(numInternalNodes-1)*num_chars+k] * vs->rootNodeFreq[rateCat*num_chars + k];
			}
		}

		for (int rateCat=0; rateCat<numRateCat; rateCat++) {
			categorySum[rateCat] = 0.0;
			for (int k=0; k<num_chars; k++) {
				categorySum[rateCat] += temp[rateCat*num_chars + k];
			}
		}

		alphaDemon = 0.0;
		for (int rateCat=0; rateCat<numRateCat; rateCat++) {
			alphaDemon += categorySum[rateCat] * vs->alpha[rateCat];
		}

		for (int rateCat=0; rateCat<numRateCat; rateCat++) {
			if (algn->isConstantSites[i]) {
				int base = currSites[0]-1;
				for (int k=0; k<num_chars; k++) {
					double denom = alphaDemon + vs->beta * vs->probXGivenInv[base];
					updateMatrix[rateCat*num_chars+k] += (double) algn->numEachUniqueSites[i] * vs->alpha[rateCat] * temp[rateCat*num_chars+k] / denom;
				}
			} else {
				for (int k=0; k<num_chars; k++) {
					updateMatrix[rateCat*num_chars+k] += (double) algn->numEachUniqueSites[i] * vs->alpha[rateCat] * temp[rateCat*num_chars+k] / alphaDemon;
				}
			}
		}

	}

	updateMatrixSum = 0.0;
	for (int k=0; k<num_chars; k++) {
		for (int rateCat=0; rateCat<numRateCat; rateCat++) {
			updateMatrixSum += updateMatrix[rateCat*num_chars + k];
		}
	}

	if (isShare) {
		for (int k=0; k<num_chars; k++) {
			double k_sum = 0.0;
			for (int rateCat=0; rateCat<numRateCat; rateCat++) {
				k_sum += updateMatrix[rateCat*num_chars + k];
			}
			for (int rateCat=0; rateCat<numRateCat; rateCat++) {
				vs->rootNodeFreq[rateCat*num_chars + k] = k_sum / updateMatrixSum;
			}
		}
	} else {
		for (int rateCat=0; rateCat<numRateCat; rateCat++) {
			double cat_sum = 0.0;
			for (int k=0; k<num_chars; k++) {
				cat_sum += updateMatrix[rateCat*num_chars + k];
			}
			for (int k=0; k<num_chars; k++) {
				vs->rootNodeFreq[rateCat*num_chars + k] = updateMatrix[rateCat*num_chars + k] / cat_sum;
			}
		}
	}

	free(likelihoodMat);
	free(temp);
	free(updateMatrix);
	free(categorySum);
	*/

}


// optimize the alpha and invariant site
void Optim::getOptAlphaInvariantSite(int maxit) {

	// optimize the invariant sites
	int isShare = 0;
	int type = 8;
	getOptValues(type, isShare, maxit);

	/*
	int numSeqs = algn->numSeqs;
	int numUniqSites = algn->numUniqueSites;
	int numInternalNodes = numSeqs-1;
	int numLineTopMat = numSeqs-1;
	double* likelihoodMat = matrix(numInternalNodes, num_chars);
	double* siteLikelihoodInvSite = matrix(num_rateCat, num_chars*num_chars);
	resetMat(siteLikelihoodInvSite, num_rateCat, num_chars);
	double* probXGivenVar = matrix(num_rateCat, num_chars);
	resetMat(probXGivenVar, num_rateCat, num_chars);
	double* probXGivenInvDeriv = matrix(1, num_chars);
	resetMat(probXGivenInvDeriv, 1, num_chars);
	double* alphaDeriv = matrix(1,num_rateCat);
	resetMat(alphaDeriv, 1, num_rateCat);
	double* currLikelihood = matrix(1, num_rateCat);
	resetMat(currLikelihood,1,num_rateCat);
	double betaDeriv = 0.0;
	double probXGivenInvDeriv_sum = 0.0;

	for (int i=0; i<numUniqSites; i++) {
		if (algn->isConstantSites[i]) {
			char* currSites = &(algn->uniqueSites[i*numSeqs]);
			int base = currSites[0]-1;

			for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
				if (ps_set != NULL)
					computeLikelihoodMatrix(likelihoodMat, ps_set->ps[rateCat]->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);
				else
					computeLikelihoodMatrix(likelihoodMat, ps->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);

				for (int k=0; k<num_chars; k++) {
					siteLikelihoodInvSite[rateCat*num_chars*num_chars + base*num_chars + k] = likelihoodMat[(numInternalNodes-1)*num_chars+k];
				}
			}
		}
	}

	for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
		for (int i=0; i<num_chars; i++) {
			for (int k=0; k<num_chars; k++) {
				probXGivenVar[rateCat*num_chars + i] += siteLikelihoodInvSite[rateCat*num_chars*num_chars+i*num_chars+k] * vs->rootNodeFreq[rateCat*num_chars + k];
			}
		}
	}

	for (int i=0; i<numUniqSites; i++) {
		char* currSites = &(algn->uniqueSites[i*numSeqs]);
		if (algn->isConstantSites[i] && vs->beta > 0.0) {
			// constant sites
			int base = currSites[0]-1;

			double denom = 0.0;
			for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
				denom += vs->alpha[rateCat] * probXGivenVar[rateCat*num_chars+base];
			}
			denom += vs->beta * vs->probXGivenInv[base];

			for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
				alphaDeriv[rateCat] += algn->numEachUniqueSites[i] * probXGivenVar[rateCat*num_chars+base]/ denom;
			}

			betaDeriv += algn->numEachUniqueSites[i] * vs->probXGivenInv[base]/ denom;
			probXGivenInvDeriv[base] +=  algn->numEachUniqueSites[i]* vs->beta * vs->probXGivenInv[base]/ denom;

		} else {
			// variable sites
			resetMat(currLikelihood,1,num_rateCat);
			double denom = 0.0;
			for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
				if (ps_set != NULL)
					computeLikelihoodMatrix(likelihoodMat, ps_set->ps[rateCat]->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);
				else
					computeLikelihoodMatrix(likelihoodMat, ps->allCondProb, topMatrix, currSites, numLineTopMat, num_chars);
				for (int k=0; k<num_chars; k++) {
					currLikelihood[rateCat] += likelihoodMat[(numInternalNodes-1)*num_chars+k]* vs->rootNodeFreq[rateCat*num_chars + k];
				}
				denom += currLikelihood[rateCat] * vs->alpha[rateCat];
			}
			for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
				alphaDeriv[rateCat] += algn->numEachUniqueSites[i] * currLikelihood[rateCat] / denom ;
			}
		}
	}


	for (int rateCat=0; rateCat<num_rateCat; rateCat++) {
		vs->alpha[rateCat] = vs->alpha[rateCat] * alphaDeriv[rateCat] / algn->numSites;
	}

	vs->beta = vs->beta * betaDeriv / algn->numSites;

	for (int i=0; i<num_chars; i++)
		probXGivenInvDeriv_sum += probXGivenInvDeriv[i];

	if (vs->beta > 0.0) {
		for (int i=0; i<num_chars; i++)
			vs->probXGivenInv[i] = probXGivenInvDeriv[i]/ probXGivenInvDeriv_sum;
	}


	free(likelihoodMat);
	free(siteLikelihoodInvSite);
	free(probXGivenVar);
	free(probXGivenInvDeriv);
	free(alphaDeriv);
	free(currLikelihood);
	*/

}


void pickOptimizationOrder(int mode, int* optItem, int isRAS) {
	// pick the order of the optimization

	// if the mode is not between 1 and 8, then randomly pick one
	if (mode < 1 || mode > NUMBER_OPT_MODES)
		mode = rand() % NUMBER_OPT_MODES + 1;

	// for RAS only, mode 7 and mode 8 are replaced by some specific modes
	if (isRAS) {
		if (mode==7)
			mode = RAS_OPT_MODE7;
		else if (mode==8)
			mode = RAS_OPT_MODE8;
	}

	switch (mode) {
	case 1:
		optItem[0]=1;optItem[1]=2;optItem[2]=3;
		break;
	case 2:
		optItem[0]=1;optItem[1]=3;optItem[2]=2;
		break;
	case 3:
		optItem[0]=2;optItem[1]=1;optItem[2]=3;
		break;
	case 4:
		optItem[0]=2;optItem[1]=3;optItem[2]=1;
		break;
	case 5:
		optItem[0]=3;optItem[1]=1;optItem[2]=2;
		break;
	case 6:
		optItem[0]=3;optItem[1]=2;optItem[2]=1;
		break;
	case 7:
		optItem[0]=1;optItem[1]=6;optItem[2]=0;
		break;
	case 8:
		optItem[0]=6;optItem[1]=1;optItem[2]=0;
		break;
	}
}

int Optim::optimizeAllParam(int numIterations, double& loglikelihood, double& IC, ThreadLocks* threadLocks, int masterID, int threadID, int maxit, int mode, int ICType, int numRateChanges, int sameRootMat, int precise) {
	// optimize all parameters
	// return the number of iterations processed

	// int minTimesForConverage = 3;
	// int times = 0;
	double thres = 0.1;
	if (precise)
		thres = 0.001;

	// double preThres = 1.0;
	// double estMaxBound = 50.0;

	int numLineTopMat = algn->numSeqs-1;

	//double pre_loglike = loglikelihood;
	double pre_loglike = LARGE_NUMBER;

	int numIter = 0;

	// get the updated value of num_rateMat
	num_rateMat = (int) (ps->rateIDs).size();

	int iterContinue = 1;

	// to reset all values in the parameters
	// ps->reset();
	// vs->resetAllVariables(*algn);
	ps->set_for_opt();
	vs->set_for_opt();

	// show the rate matrix
	// ps->printRateMatID();

	// compute the loglikelihood value
	int theEdge = -1; // for all edges
	ps->computeAllEigenValues_t(theEdge);
	ps->computeAllCondProb(isReversible);
	// cout << "computing the loglikelihood value" << endl << flush;
	// loglikelihood = getLogLikelihood(ps->allCondProb, topMatrix, numLineTopMat, (*algn), (*vs));
	// cout << loglikelihood << endl << flush;

	int  optItem[3]; // 1: edge len; 2: s-matrix; 3: pi

	if (mode >= 1 && mode <= NUMBER_OPT_MODES) {
		pickOptimizationOrder(mode, optItem, 0);
		// get the order of the optimization
	}

	int i;
	int isShareRootVector = 0;

	// VariableSet vs_zero_beta(vs->num_alpha, vs->num_chars);
	// double loglikelihood_zero_beta;

	while (iterContinue) {

		numIter++;

		if (mode < 1 || mode > NUMBER_OPT_MODES) {
			pickOptimizationOrder(mode, optItem, 0);
			// get the order of the optimization
		}

		int isShare = 0; // for RAL, this is always zero

		// optimizing the edge lengths / s-matrix / pi
		for (i=0; i<3; i++) {

			getOptValues(optItem[i], isShare, maxit);
			
			// cout << "optItem[" << i << "]=" << optItem[i] << endl;
			// loglikelihood = getLogLikelihood(ps->allCondProb, topMatrix, numLineTopMat, (*algn), (*vs), num_chars);
			// cout << loglikelihood << endl;

		}

#ifndef GTRIFO
		// optimizing the root vector
		getOptRootVector(isShareRootVector, maxit);
#endif

		// cout << "optimizing the root vector" << endl;
		// loglikelihood = getLogLikelihood(ps->allCondProb, topMatrix, numLineTopMat, (*algn), (*vs), num_chars);
		// cout << loglikelihood << endl;

		getOptAlphaInvariantSite(maxit);
		// cout << "optimizing the alpha and invariant site" << endl;
		loglikelihood = getLogLikelihood(ps->allCondProb_t, topMatrix, numLineTopMat, (*algn), (*vs), num_chars);
		// cout << setprecision(10) << loglikelihood << endl;

		if (numIterations==-1 && numIter > 1 && (loglikelihood-pre_loglike)<thres) {
			iterContinue = 0;
		}

		// break due to the number of iterations
		if (numIterations!=-1 && numIter >= numIterations)
			iterContinue = 0;

		// cout << "# " << loglikelihood << endl;
		pre_loglike = loglikelihood;

	}

	IC = getIC(loglikelihood, ps->num_rate_matrices, algn->numSites, algn->numSeqs, ICType, numRateChanges, sameRootMat);

	return numIter;
}

// optimize all parameters for RAS
// return the number of iterations processed
// Output: 1. loglikelihood; 2. IC; 3. df (degree of freedoms)
int Optim::optimizeAllParamRAS(int numIterations, int numRateMat, double& loglikelihood,
	double& IC, int& df, int isGTRUpsilon, int threadID, ThreadLocks* threadLocks,
	int maxit, int mode, int ICType, int sameRootMat, int numRateChanges, int precise) {

	if (modelType == -1) {

		if (threadLocks!=NULL)
			pthread_mutex_lock(&(threadLocks->lock_output));
		cerr << "[Optim::optimizeAllParamRAS] Error! The model type has not been defined" << endl;
		if (threadLocks!=NULL)
			pthread_mutex_unlock(&(threadLocks->lock_output));

	}

#ifdef LARGE_DATA
	double thres = 0.1;
#else
	double thres = 0.01;
	if (precise)
		thres = 0.001;
#endif

	int numLineTopMat = algn->numSeqs-1;
	int numInternalNodes = numLineTopMat;

	double pre_loglike = loglikelihood;

	int numIter = 0;

	int isShare;

	int i;

	// to reset all values in the parameters
	ps_set->reset();
	vs->resetAllVariables(*algn);
	
	// set up the parameters for optimisation method
	ps_set->set_for_opt();
	vs->set_for_opt();
	
	// set up the temporary arrays for optimisation method
	if (HASsiteLikelihoods == NULL) {
		HASsiteLikelihoods = matrix(algn->numUniqueSites * (int)ps_set->allCondProbSet_t.size(), num_chars);
	}
	if (likelihoodMat==NULL)
		likelihoodMat = matrix(numInternalNodes, num_chars);
	if (leftLikelihoodMat==NULL)
		leftLikelihoodMat = matrix(numInternalNodes,num_chars);
	if (rightLikelihoodMat==NULL)
		rightLikelihoodMat = matrix(numInternalNodes,num_chars);
	
	// show the rate matrix
	// ps_set->ps[0]->printRateMatID();

	// compute the loglikelihood value
	ps_set->computeAllEigenValues_t();
	ps_set->computeAllCondProb(isReversible);
	
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
	printf("Time used before compuring the loglikelihood: %4.2f ms\n", second);
	pre_time = curr_time;
	cout << "computing the loglikelihood value" << endl << flush;
#endif

	loglikelihood = getLogLikelihood(ps_set->allCondProbSet_t, topMatrix, numLineTopMat, (*algn), (*vs), num_chars);
	
#ifdef GET_TIME_STAT_DETAIL
	cout << loglikelihood << endl << flush;
	curr_time = clock();
	second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
	printf("Time used after compuring the loglikelihood: %4.2f ms\n", second);
	pre_time = curr_time;
#endif

	int  optItem[3]; // 1: edge len; 2: s-matrix; 3: pi

	if (mode >= 1 && mode <= NUMBER_OPT_MODES) {
		pickOptimizationOrder(mode, optItem, 1);
		// get the order of the optimization
	}

	while (1) {

		numIter++;

		if (mode < 1 || mode > NUMBER_OPT_MODES) {
			pickOptimizationOrder(mode, optItem, 1);
			// get the order of the optimization
		}

		// optimizing the edge lengths / s-matrix / pi
		for (i=0; i<3; i++) {

			// for RAS, optItem[i] is either 1, 2 or 3

			if (optItem[i]==1) {
				// optimizing the edge lengths
#ifdef GET_TIME_STAT_DETAIL
				cout << "optimizing the edge lengths [numIter : " << numIter << "]" << endl;
#endif
				if (modelType != 5) {
					// not scalar edge length
					isShare = 0;
					for (int i=0; i<ps_set->size(); i++) {
						ps = ps_set->ps[i];
						currRateCatID = i;
						getOptValues(1, isShare, maxit);
					}
				} else {
					// scalar edge length
					isShare = 1;
					if (num_rateCat > 1) {
						getOptValues(4, isShare, maxit); // optimizing the edge scalars
					}
					ps = NULL;
					getOptValues(1, isShare, maxit); // optimizing the edge length of the first rate category
				}
			} else if (optItem[i]==2) {
				// optimizing the s-matrix
#ifdef GET_TIME_STAT_DETAIL
				cout << "optimizing the s-matrix [numIter : " << numIter << "]" << endl;
#endif
				if (modelType == 1) {
					// not sharing the same S
					isShare = 0;
					for (int i=0; i<ps_set->size(); i++) {
						ps = ps_set->ps[i];
						currRateCatID = i;
						getOptValues(2, isShare, maxit);
					}
				} else {
					// sharing the same S
					isShare = 1;
					ps = NULL;
					getOptValues(2, isShare, maxit);
				}
			} else if (optItem[i]==3) {
				// optimizing the pi
#ifdef GET_TIME_STAT_DETAIL
				cout << "optimizing the pi [numIter : " << numIter << "]" << endl;
#endif
				if (modelType <= 2) {
					// not sharing the same Pi
					isShare = 0;
					for (int i=0; i<ps_set->size(); i++) {
						ps = ps_set->ps[i];
						currRateCatID = i;
						if (isGTRUpsilon)
							getOptValues(5, isShare, maxit);
						else
							getOptValues(3, isShare, maxit);
					}
				} else {
					isShare = 1;
					ps = NULL;
					if (isGTRUpsilon)
						getOptValues(5, isShare, maxit);
					else
						getOptValues(3, isShare, maxit);
				}
			}
			
#ifdef GET_TIME_STAT_DETAIL
			curr_time = clock();
			second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
			printf("Time used: %4.2f ms\n", second);
			pre_time = curr_time;
#endif

		}

		// optimizing the root vectors
		if (!isGTRUpsilon) {
			isShare = 0;
			if (modelType > 3) {
				// sharing the same root vectors
				isShare = 1;
			}
			getOptRootVector(isShare, maxit);
		}

		// optimize the alpha and invariant site
		getOptAlphaInvariantSite(maxit);

		loglikelihood = getLogLikelihood(ps_set->allCondProbSet_t, topMatrix, numLineTopMat, (*algn), (*vs), num_chars);

		if (threadLocks!=NULL)
			pthread_mutex_lock(&(threadLocks->lock_output));
		// cout << "[threadID:" << threadID << " - " << numIter << "] " << longDoublToStr(loglikelihood,ANS_DECI) << endl;
		if (threadLocks!=NULL)
			pthread_mutex_unlock(&(threadLocks->lock_output));

		if (numIter > 1 && numIterations==-1 && (loglikelihood-pre_loglike)<thres)
			break;

		if (numIterations!=-1 && numIter >= numIterations)
			break;

		pre_loglike = loglikelihood;

	}

	int nsites = algn->numSites;
	int use_reg = 1;
	IC = getIC(loglikelihood, numRateMat, nsites, algn->numSeqs, num_rateCat, modelType, df, isGTRUpsilon, ICType, numRateChanges, sameRootMat, use_reg);

	return numIter;

}


// destructor
Optim::~Optim() {
	if (edgeScalars!=NULL)
		free(edgeScalars);
	if (likelihoodMat!=NULL)
		free(likelihoodMat);
	if (leftLikelihoodMat!=NULL)
		free(leftLikelihoodMat);
	if (rightLikelihoodMat!=NULL)
		free(rightLikelihoodMat);
	// if (changedNodes!=NULL)
	//    free(changedNodes);
	if (HASsiteLikelihoods!=NULL)
		free(HASsiteLikelihoods);
}

void printEdgeLens(Optim* op) {
	cout << "edge scalars:";
	for (int i=0; i<op->num_rateCat-1; i++) {
		cout << op->edgeScalars[i] << ",";
	}
	cout << endl;
	cout << "edge lengths:" << endl;
	for (int i=0; i<op->num_rateCat; i++) {
		for (int j=0; j<op->num_edges; j++) {
			cout << op->ps_set->ps[i]->t[j] << ",";
		}
		cout << endl;
	}
}

//======================================================================================
// CORE FUNCTIONS FOR OPTIMIZATION
//======================================================================================

// initial update on the parameters and variables
static void initialUpdate(Optim* op, double *values, int num, vector<int>* edges, int num_chars) {

	/*
	// print out the values
	cout << "[1] values: ";
	for (int i=0; i<num; i++)
		cout << values[i] << " ";
	cout << endl;
	 */

	int num_w = 6;
	int num_edges = op->num_edges;
	int num_rateMat = op->num_rateMat;
	ParameterSet *curr_ps = NULL;
	ParameterSet *ps0 = NULL;
	if (op->isShare)
		ps0 = op->ps_set->ps[0];

	double* curr_w;
	double* curr_pi; double* curr_pi_t;
	double sum_values;
	int i,j,k,l;

/*
	// check whether there is any nan value
	for (i=0; i<num; i++) {
		if (values[i]!=values[i]) {
			// it is nan
			// restore the previous value
			if (op->paramType == 1) {
				// edge lens
				if (op->isShare) {
					values[i] = ps0->t[i];
				} else {
					values[i] = op->ps->t[i];
				}
			} else if (op->paramType == 2) {
				// s-matrices
				if (op->isShare) {
					curr_ps = op->ps_set->ps[0];
					curr_w = &(curr_ps->w[(edges->at(0)) * curr_ps->num_w]);
					values[i] = curr_w[i];
				} else {
					curr_w = &(op->ps->w[(edges->at(0)) * op->ps->num_w]);
					values[i] = curr_w[i];
				}
			} else if (op->paramType == 3 || op->paramType == 5) {
				// pi
				if (op->isShare) {
					curr_ps = op->ps_set->ps[0];
					curr_pi = &(curr_ps->pi[(edges->at(0)) * num_chars]);
					values[i] = curr_pi[i];
				} else {
					curr_pi = &(op->ps->pi[(edges->at(0)) * num_chars]);
					values[i] = curr_pi[i];
				}
			} else if (op->paramType == 4) {
				// edge scalars
				values[i] = op->edgeScalars[op->currRateCatID-1];
			} else if (op->paramType == 6) {
				if (i < num_w) {
					// s-matrices
					if (op->isShare) {
						curr_ps = op->ps_set->ps[0];
						curr_w = &(curr_ps->w[(edges->at(0)) * curr_ps->num_w]);
						values[i] = curr_w[i];
					} else {
						curr_w = &(op->ps->w[(edges->at(0)) * op->ps->num_w]);
						values[i] = curr_w[i];
					}
				} else {
					// pi
					if (op->isShare) {
						curr_ps = op->ps_set->ps[0];
						curr_pi = &(curr_ps->pi[(edges->at(0)) * num_chars]);
						values[i] = curr_pi[i-num_w];
					} else {
						curr_pi = &(op->ps->pi[(edges->at(0)) * num_chars]);
						values[i] = curr_pi[i-num_w];
					}
				}
			} else if (op->paramType == 7) {
				// for edge lengths, all s-matrix,and all PIs
				// num_edges + num_rateMat * num_w + num_rateMat * num_chars

				// restore the previous value
				if (i<num_edges) {
					// edge lengths (# : num_edges)
					values[i] = op->ps->t[i];
				} else if (i < num_edges + num_rateMat * num_w) {
					// all s-matrix (# : num_rateMat * num_w)
					j = ( i - num_edges ) / num_w;
					k = ( i - num_edges ) % num_w;
					values[i] = curr_ps->w[((curr_ps->rateIDs)[j]->at(0)) * num_w + k];
				} else {
					// all PIs (# : num_rateMat * num_chars)
					j = ( i - num_edges + num_rateMat * num_w ) / num_chars;
					k = ( i - num_edges + num_rateMat * num_w ) % num_chars;
					values[i] = curr_ps->pi[((curr_ps->rateIDs)[j]->at(0)) * num_chars + k];
				}
			} else if (op->paramType == 8) {
				// alpha values
				if (i < op->num_rateCat)
					values[i] = op->vs->alpha[i];
				else if (i < op->num_rateCat+1)
					values[i] = op->vs->beta;
				else
					values[i] = op->vs->probXGivenInv[i - op->num_rateCat - 1];
			} else if (op->paramType == 9) {
				// root frequency
				values[i] = op->vs->rootNodeFreq[i];
			}
		}
	}
*/
	// reinforce the maximum of "values" to be MAX_values
	if (op->paramType == 3 || op->paramType == 5) {
		max_min_decimals(values, MAX_PI, MIN_PI, num);
	} else if (op->paramType == 6) {
		max_min_decimals(values, MAX_values, MIN_values, num_w);
		max_min_decimals(values, MAX_PI, MIN_PI, num-num_w);
	} else if (op->paramType == 7) {
		max_min_decimals(values, MAX_values, MIN_values, num_edges+num_rateMat*num_w);
		max_min_decimals(values, MAX_PI, MIN_PI, num-(num_edges+num_rateMat*num_w));
	} else {
		max_min_decimals(values, MAX_values, MIN_values, num);
	}

	if (op->paramType == 1) {
		// update the edge lens
		if (op->isShare) {
			// scalar edge length
			for (int i=0; i<num; i++) {
				ps0->t[i] = values[i];
			}
			for (int rateCat=1; rateCat<op->num_rateCat; rateCat++) {
				curr_ps = op->ps_set->ps[rateCat];
				for (int i=0; i<op->num_edges; i++) {
					curr_ps->t[i] = ps0->t[i] * op->edgeScalars[rateCat-1];
				}
			}
		} else {
			for (int i=0; i<num; i++) {
				op->ps->t[i] = values[i];
			}
		}
	} else if (op->paramType == 2) {
		// update the corresponding w values of the s-matrices
		// for all the edges belonging the rate matrix "op->currRateMatID"
		// "num" here equals to the number of w values
		if (edges->size()==0) {
			cerr << "[getUpdateValue] Error! No edge is belonging to rate matrice group :" << op->currRateMatID << endl;
		}
		if (op->isShare) {
			for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
				curr_ps = op->ps_set->ps[rateCat];
				for (int i=0; i<(int)edges->size(); i++) {
					curr_w = &(curr_ps->w[(edges->at(i)) * curr_ps->num_w]);
					for (int k=0; k<num; k++) {
						curr_w[k] = values[k];
					}
				}
			}
		} else {
			for (int i=0; i<(int)edges->size(); i++) {
				curr_w = &(op->ps->w[(edges->at(i)) * op->ps->num_w]);
				for (int k=0; k<num; k++) {
					curr_w[k] = values[k];
				}
			}
		}
	} else if (op->paramType == 3 || op->paramType == 5) {
		// update the Pi values
		// compute the sum of values
		sum_values = 0.0;
		for (int k=0; k<num; k++) {
			sum_values += values[k];
		}

		if (edges->size()==0) {
			cerr << "[getUpdateValue] Error! No edge is belonging to rate matrice group :" << op->currRateMatID << endl;
		}
		
		if (sum_values > 0.0) {
			if (op->isShare) {
				for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
					curr_ps = op->ps_set->ps[rateCat];
					for (int i=0; i<(int)edges->size(); i++) {
						curr_pi = &(curr_ps->pi[(edges->at(i)) * num_chars]);
						curr_pi_t = &(curr_ps->pi_t[(edges->at(i)) * num_chars]);
						for (int k=0; k<num; k++) {
							curr_pi[k] = values[k]/sum_values;
							curr_pi_t[k] = values[k];
						}
					}
				}
#ifdef GTRIFO
				curr_ps = op->ps_set->ps[0];
				curr_pi = &(curr_ps->pi[(edges->at(0)) * num_chars]);
				curr_pi_t = &(curr_ps->pi_t[(edges->at(0)) * num_chars]);
				for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
					for (int k=0; k<num; k++) {
						op->vs->rootNodeFreq[rateCat*num_chars+k] = curr_pi[k];
						op->vs->rootNodeFreq_t[rateCat*num_chars+k] = curr_pi_t[k];
					}
				}
				for (int k=0; k<num; k++) {
					op->vs->probXGivenInv[k] = curr_pi[k];
					op->vs->probXGivenInv_t[k] = curr_pi_t[k];
				}
#endif
			} else {
				for (int i=0; i<(int)edges->size(); i++) {
					curr_pi = &(op->ps->pi[(edges->at(i)) * num_chars]);
					curr_pi_t = &(op->ps->pi_t[(edges->at(i)) * num_chars]);
					for (int k=0; k<num; k++) {
						curr_pi[k] = values[k]/sum_values;
						curr_pi_t[k] = values[k];
					}
				}
			}
#ifdef GTRIFO
			curr_ps = op->ps;
			curr_pi = &(curr_ps->pi[(edges->at(0)) * num_chars]);
			curr_pi_t = &(curr_ps->pi_t[(edges->at(0)) * num_chars]);
			for (int k=0; k<num; k++) {
				op->vs->rootNodeFreq[k] = curr_pi[k];
				op->vs->rootNodeFreq_t[k] = curr_pi_t[k];
			}
			for (int k=0; k<num; k++) {
				op->vs->probXGivenInv[k] = curr_pi[k];
				op->vs->probXGivenInv_t[k] = curr_pi_t[k];
			}
#endif
		}

		// added for upilson model (i.e. op->paramType == 5)
		if (op->paramType == 5) {
			curr_ps = op->ps_set->ps[0];
			curr_pi = &(curr_ps->pi[(edges->at(0)) * num_chars]);
			for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
				for (int k=0; k<num; k++) {
					op->vs->rootNodeFreq[rateCat*num_chars+k] = curr_pi[k];
				}
			}
		}
	} else if (op->paramType == 4) {

		// edge scalars
		op->edgeScalars[op->currRateCatID-1] = fabs(values[0]);
		// cout << "op->edgeScalars[op->"<< op->currRateCatID << "-1] = " << op->edgeScalars[op->currRateCatID-1] << endl << flush;
		curr_ps = op->ps_set->ps[op->currRateCatID];
		for (int i=0; i<op->num_edges; i++) {
			curr_ps->t[i] = ps0->t[i] * op->edgeScalars[op->currRateCatID-1];
		}
	} else if (op->paramType == 6) {
		// (obsoleted)
		// update the corresponding w values of the s-matrices
		// for all the edges belonging the rate matrix "op->currRateMatID"
		// "num" here equals to the number of w values
		if (edges->size()==0) {
			cerr << "[getUpdateValue] Error! No edge is belonging to rate matrice group :" << op->currRateMatID << endl;
		}

		// ==for s-matrix==
		if (op->isShare) {
			for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
				curr_ps = op->ps_set->ps[rateCat];
				for (int i=0; i<(int)edges->size(); i++) {
					curr_w = &(curr_ps->w[(edges->at(i)) * curr_ps->num_w]);
					for (int k=0; k<num_w; k++) {
						curr_w[k] = values[k];
					}
				}
			}
		} else {
			for (int i=0; i<(int)edges->size(); i++) {
				curr_w = &(op->ps->w[(edges->at(i)) * op->ps->num_w]);
				for (int k=0; k<num_w; k++) {
					curr_w[k] = values[k];
				}
			}
		}

		// ==for pi==
		// compute the sum of values
		double sum_values = 0.0;
		for (int k=num_w; k<num_w+num_chars; k++) {
			sum_values += values[k];
		}

		if (op->isShare) {
			for (int rateCat=0; rateCat<op->num_rateCat; rateCat++) {
				curr_ps = op->ps_set->ps[rateCat];
				for (int i=0; i<(int)edges->size(); i++) {
					curr_pi = &(curr_ps->pi[(edges->at(i)) * num_chars]);
					curr_pi_t = &(curr_ps->pi_t[(edges->at(i)) * num_chars]);
					for (int k=0; k<num_chars; k++) {
						curr_pi[k] = values[k+num_w]/sum_values;
						curr_pi_t[k] = values[k+num_w];
					}
				}
			}
		} else {
			for (int i=0; i<(int)edges->size(); i++) {
				curr_pi = &(op->ps->pi[(edges->at(i)) * num_chars]);
				curr_pi_t = &(op->ps->pi_t[(edges->at(i)) * num_chars]);
				for (int k=0; k<num_chars; k++) {
					curr_pi[k] = values[k+num_w]/sum_values;
					curr_pi_t[k] = values[k+num_w];
				}
			}
		}

	} else if (op->paramType == 7) {
		// obsoleted
		for (i=0; i<num; i++) {
			if (i<num_edges) {
				// update the edge lens
				op->ps->t[i] = values[i];
			} else if (i < num_edges + num_rateMat * num_w) {
				// all s-matrix (# : num_rateMat * num_w)
				j = ( i - num_edges ) / num_w;
				k = ( i - num_edges ) % num_w;
				for (l=0; l<(int)(op->ps->rateIDs)[j]->size(); l++)
					op->ps->w[((op->ps->rateIDs)[j]->at(l)) * op->ps->num_w + k] = values[i];
			} else {
				// all PIs (# : num_rateMat * num_chars)
				j = ( i - (num_edges + num_rateMat * num_w) ) / num_chars;
				k = ( i - (num_edges + num_rateMat * num_w) ) % num_chars;
				for (l=0; l<(int)(op->ps->rateIDs)[j]->size(); l++) {
					op->ps->pi[((op->ps->rateIDs)[j]->at(l)) * num_chars + k] = values[i];
					op->ps->pi_t[((op->ps->rateIDs)[j]->at(l)) * num_chars + k] = values[i];
				}
			}
		}
	} else if (op->paramType == 8) {
		// compute the sum for alphas and beta
		sum_values = 0.0;
		for (i=0; i<op->num_rateCat; i++) {
			sum_values += values[i];
		}
		if (!op->vs->isGTRUpilson) {
			sum_values += values[i]; // beta value
		}
		if (sum_values > 0.0) {
			// alpha values
			for (i=0; i<op->num_rateCat; i++) {
				op->vs->alpha[i] = values[i] / sum_values;
				op->vs->alpha_t[i] = values[i];
			}
			if (!op->vs->isGTRUpilson) {
				// beta value
				op->vs->beta = values[i] / sum_values;
				op->vs->beta_t = values[i];
			}
		}
		
#ifndef GTRIFO
			// compute the sum for nucleotide distribution along constant sites
			sum_values = 0.0;
			for (i=0; i<num_chars; i++) {
				sum_values += values[i+op->num_rateCat+1];
			}
			if (sum_values > 0.0) {
				// nucleotide distribution along constant sites
				for (i=0; i<num_chars; i++) {
					op->vs->probXGivenInv[i] = values[i+op->num_rateCat+1]/sum_values;
					op->vs->probXGivenInv_t[i] = values[i+op->num_rateCat+1];
				}
			}
#endif
		}

	} else if (op->paramType == 9) {
		// nucleotide distribution along the root
		sum_values = 0.0;
		for (k=0; k<op->num_rateCat; k++) {
			if (k==0 || !(op->isShare)) {
				sum_values = 0.0;
				for (i=0; i<num_chars; i++)
					sum_values += values[k*num_chars+i];
			}
			for (i=0; i<num_chars; i++) {
				if (op->isShare) {
					op->vs->rootNodeFreq[k*num_chars+i] = values[i]/sum_values;
					op->vs->rootNodeFreq_t[k*num_chars+i] = values[i];
				} else {
					op->vs->rootNodeFreq[k*num_chars+i] = values[k*num_chars+i]/sum_values;
					op->vs->rootNodeFreq_t[k*num_chars+i] = values[k*num_chars+i];
				}
			}
		}
	}

	/*
	// print out the values
	op->ps->showContent();
	op->vs->showContent();
	cout << "[2] values: ";
	for (int i=0; i<num; i++)
		cout << values[i] << " ";
	cout << endl;
	 */
}


// The objective function to be defined
// Input:
// 1. int n          : number of values
// 2. double* x      : values of the parameters
// 3. opt->paramType : a) all edge lengths;
//                     b) s-matrix along an edge;
//                     c) pi-vector along an endge;
//                     d) scalar edge lengths
// Output:
// The corresponding log likelihood value
double fn(int n, double* x, void* opt) {
	Optim* op = (Optim*) opt;

	// for paramType 2 and paramType 3
	vector<int>* edges = NULL;
	if (op->paramType == 2 || op->paramType == 3 || op->paramType == 5 || op->paramType == 6) {
		if (op->isShare) {
			edges = op->ps_set->ps[0]->rateIDs[op->currRateMatID];
		} else if (op->paramType != 7) {
			edges = op->ps->rateIDs[op->currRateMatID];
		}
	}
	int numLineTopMat = op->algn->numSeqs-1;

	// initial update on the parameters and variables
	initialUpdate(op, x, n, edges, op->num_chars);

	// compute the corresponding log likelihood value of this set of variables
	// show the content of ps
	// cout << "content of ps" << endl;
	// op->ps->showContent();
	// cout << endl;
	if (op->ps_set!=NULL && (op->isShare || op->paramType==8 || op->paramType==9))
		op->ps_set->computeAllEigenValues_t();
	else
		op->ps->computeAllEigenValues_t(-1);

	if (op->ps_set!=NULL && (op->isShare || op->paramType==8 || op->paramType==9))
		op->ps_set->computeAllCondProb(op->isReversible);
	else
		op->ps->computeAllCondProb(op->isReversible);

	if (op->ps_set == NULL)
		return -getLogLikelihood(op->ps->allCondProb_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars);
	else
		return -getLogLikelihood(op->ps_set->allCondProbSet_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->HASsiteLikelihoods, -1, op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat);
}


void updateAllEigenAndCondProb(int n, double* x, int updated_item1, int updated_item2, Optim* op) {

	vector<int>* edges = NULL;
	if (op->paramType == 2 || op->paramType == 3 || op->paramType == 5 || op->paramType == 6) {
		if (op->isShare) {
			edges = op->ps_set->ps[0]->rateIDs[op->currRateMatID];
		} else if (op->paramType != 7) {
			edges = op->ps->rateIDs[op->currRateMatID];
		}
	}
	int i,j;
	
	// initial update on the parameters and variables
	initialUpdate(op, x, n, edges, op->num_chars);
	
	// compute the corresponding log likelihood value of this set of variables
	// show the content of ps
	// cout << "content of ps" << endl;
	// op->ps->showContent();
	// cout << endl;
	if (op->ps_set==NULL) {
		if (op->paramType==2 || op->paramType==3)
			op->ps->computeAllEigenValuesForRateGrp_t(op->currRateMatID);
		else if (op->paramType!=1)
			op->ps->computeAllEigenValues_t(-1);
		

		if ((op->paramType==2 || op->paramType==3) && (edges != NULL)) {
			for (i=0; i<(int)edges->size(); i++)
				op->ps->computeCondProbForEdge(edges->at(i));
		} else if (op->paramType==1) {
			if (updated_item1 != -1)
				op->ps->computeCondProbForEdge(updated_item1);
			if (updated_item2 != -1)
				op->ps->computeCondProbForEdge(updated_item2);
		} else
			op->ps->computeAllCondProb(op->isReversible);
	
	} else {
		
		// if (op->isShare || op->paramType==8 || op->paramType==9)
		//	op->ps_set->computeAllEigenValues();
		if (op->paramType==2 || op->paramType==3) {
			if (op->isShare) {
				for (i=0; i<(int)op->ps_set->ps.size(); i++)
					op->ps_set->ps[i]->computeAllEigenValuesForRateGrp_t(op->currRateMatID);
			} else {
				op->ps->computeAllEigenValuesForRateGrp_t(op->currRateMatID);
			}
		} else if (op->paramType==4) {
			op->ps->computeAllEigenValues_t(-1);
		}
		
		
		if (op->paramType==2 || op->paramType==3) {
			if (op->isShare) {
				for (i=0; i<(int)op->ps_set->ps.size(); i++)
					for (j=0; j<(int)edges->size(); j++)
						op->ps_set->ps[i]->computeCondProbForEdge(edges->at(j));
			} else {
				for (j=0; j<(int)edges->size(); j++)
					op->ps->computeCondProbForEdge(edges->at(j));
			}
		} else if (op->paramType==1) {
			if (op->isShare) {
				for (i=0; i<(int)op->ps_set->ps.size(); i++) {
					if (updated_item1 != -1)
						op->ps_set->ps[i]->computeCondProbForEdge(updated_item1);
					if (updated_item2 != -1)
						op->ps_set->ps[i]->computeCondProbForEdge(updated_item2);
				}
			} else {
				if (updated_item1 != -1)
					op->ps->computeCondProbForEdge(updated_item1);
				if (updated_item2 != -1)
					op->ps->computeCondProbForEdge(updated_item2);
			}
		} else if (op->paramType==4) {
			op->ps->computeAllCondProb(op->isReversible);
		}

	}
}


// The objective function used when some variables are updated
// Input:
// 1. int n          : number of values
// 2. double* x      : values of the parameters
// 3. opt->paramType : 1 - all edge lengths;
//                     2 - s-matrix along an edge;
//                     3 - pi-vector along an edge;
//                     4 - scalar edge lengths
// Output:
// The corresponding log likelihood value
double fn_update(int n, double* x, int updated_item1, int updated_item2, void* opt) {

	Optim* op = (Optim*) opt;
	
	updateAllEigenAndCondProb(n, x, updated_item1, updated_item2, op);
	
	int numLineTopMat = op->algn->numSeqs-1;
	
   
	if (op->ps_set == NULL)
		return -getLogLikelihood(op->ps->allCondProb_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat); //, op->align->changedNodes);
	else
		return -getLogLikelihood(op->ps_set->allCondProbSet_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->HASsiteLikelihoods, -1, op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat);
}


// the gradient function
void compute_gr1(int n, double *x, double *g, void *opt) {
	
	
	int i;
	double f1, f2;
	double orig_x;

	Optim* op = (Optim*) opt;

	// cout << "[inside compute_gr] op->paramType: " << op->paramType << endl << flush;
	
	if (op->ps_set == NULL && (op->paramType>=1 && op->paramType<=3)) {
		int num_edges = op->num_edges;
		double* condProbArray = matrix(2*n*num_edges, 16);
		// double* condProb;
		double* x1_array = matrix(n,1);
		double* x2_array = matrix(n,1);
		double* log_Ls = matrix(n,2);
		vector<int>* edges = NULL;
		if (op->paramType != 1) {
			edges = op->ps->rateIDs[op->currRateMatID];
		}
		int numLineTopMat = op->algn->numSeqs-1;

		for (i=0; i<n; i++) {
			orig_x = x[i];
			if (x[i] > DELTA/2) {
				x[i] -= DELTA/2;
			}
			x1_array[i] = x[i];
			// condProb = op->ps->allCondProb;
			// op->ps->allCondProb = &condProbArray[2*i*num_edges*16];
			updateAllEigenAndCondProb(n, x, i-1, i, op);
			memcpy(&condProbArray[2*i*num_edges*16], op->ps->allCondProb_t, num_edges*16*sizeof(double));
			x[i] = orig_x + DELTA/2;
			x2_array[i] = x[i];
			// op->ps->allCondProb = &condProbArray[(2*i+1)*num_edges*16];
			updateAllEigenAndCondProb(n, x, -1, i, op);
			memcpy(&condProbArray[(2*i+1)*num_edges*16], op->ps->allCondProb_t, num_edges*16*sizeof(double));
			// op->ps->allCondProb = condProb;
			x[i] = orig_x;
		}
		
		// compute the total log-likelihood for a phylogenetic tree
		// for many conditional probabilities
		getLogLikelihood(condProbArray, 2*n, edges, op->topMatrix, numLineTopMat, *(op->algn), *(op->vs), op->num_chars, log_Ls, 1);
		
		for (i=0; i<n; i++) {
			// beware that our target is to get the -ve of loglikehood value
			// thus: g[i] = ((-log_Ls[2*i+1]) - (-log_Ls[2*i+1])) / (x2_array[i] - x1_array[i]);
			g[i] = (log_Ls[2*i] - log_Ls[2*i+1]) / (x2_array[i] - x1_array[i]);
		}
		
		free (condProbArray);
		free (x1_array);
		free (x2_array);
		free (log_Ls);        
	} else if (op->ps_set == NULL && (op->paramType==8 || op->paramType==9)) {
		double* x1_array = matrix(n,1);
		double* x2_array = matrix(n,1);
		double* log_Ls = matrix(n,2);
		vector<VariableSet*> vsArray;
		VariableSet* vs;
		int numLineTopMat = op->algn->numSeqs-1;
		for (i=0; i<n; i++) {
			orig_x = x[i];
			if (x[i] > DELTA/2) {
				x[i] -= DELTA/2;
			}
			x1_array[i] = x[i];
			initialUpdate(op, x, n, NULL, op->num_chars);
			vs = new VariableSet(op->num_chars);
			vs->copyFrom(*op->vs);
			vsArray.push_back(vs);
			x[i] = orig_x + DELTA/2;
			x2_array[i] = x[i];
			initialUpdate(op, x, n, NULL, op->num_chars);
			vs = new VariableSet(op->num_chars);
			vs->copyFrom(*op->vs);
			vsArray.push_back(vs);
			x[i] = orig_x;
		}
		
		// op->ps->computeAllEigenValues(-1);
		// op->ps->computeAllCondProb(op->isReversible);
		getLogLikelihood(op->ps->allCondProb_t, op->topMatrix, numLineTopMat,*(op->algn), &vsArray, op->num_chars, log_Ls);
		for (i=0; i<n; i++) {
			// beware that our target is to get the -ve of loglikehood value
			// thus: g[i] = ((-log_Ls[2*i+1]) - (-log_Ls[2*i+1])) / (x2_array[i] - x1_array[i]);
			g[i] = (log_Ls[2*i] - log_Ls[2*i+1]) / (x2_array[i] - x1_array[i]);
		}
		free (x1_array);
		free (x2_array);
		free (log_Ls);
		for (i=0; i<(int)vsArray.size(); i++)
			delete vsArray[i];
	} else {
		double x1, x2;
		for (i=0; i<n; i++) {
			orig_x = x[i];
			if (x[i] > DELTA/2) {
				x[i] -= DELTA/2;
			}
			x1 = x[i];
			f1 = fn_update(n, x, i-1, i, opt);
			x[i] = orig_x + DELTA/2;
			x2 = x[i];
			f2 = fn_update(n, x, -1, i, opt);
			x[i] = orig_x;
			g[i] = ( f2-f1 ) / ( x2 - x1 );
		}
	}
	
	/*
	if (op->paramType==8) {
		for (i=0; i<n; i++) {
			cout << g[i] << ",";
		}
		cout << endl;
	}
	 */
	// cout << "[quit compute_gr] op->paramType: " << op->paramType << endl << flush;
}


// the gradient function
void compute_gr2(int n, double *x, double *g, void *opt) {
	
	
	int i;
	double f0, f1;

	Optim* op = (Optim*) opt;
	int numLineTopMat = op->algn->numSeqs-1;

	// cout << "[inside compute_gr] op->paramType: " << op->paramType << endl << flush;

	if (op->ps_set == NULL && (op->paramType>=1 && op->paramType<=3)) {
		int num_edges = op->num_edges;
		double* condProbArray = matrix(n*num_edges, 16);
		double* log_Ls = matrix(n,1);
		vector<int>* edges = NULL;
		if (op->paramType != 1) {
			edges = op->ps->rateIDs[op->currRateMatID];
		}

		f0 = fn(n,x,opt);
		// printf("f0 : %14.10f\n", f0);
		for (i=0; i<n; i++) {
			x[i] += DELTA;
			updateAllEigenAndCondProb(n, x, i-1, i, op);
			memcpy(&condProbArray[i*num_edges*16], op->ps->allCondProb_t, num_edges*16*sizeof(double));
			x[i] -= DELTA;
		}
		
		// compute the total log-likelihood for a phylogenetic tree
		// for many conditional probabilities
		getLogLikelihood(condProbArray, n, edges, op->topMatrix, numLineTopMat, *(op->algn), *(op->vs), op->num_chars, log_Ls, 2);
		
		for (i=0; i<n; i++) {
			// beware that our target is to get the -ve of loglikehood value
			// printf("log_Ls[%d] : %14.10f\n", i, log_Ls[i]);
			g[i] = (-log_Ls[i] - f0) / DELTA;
		}

		free (condProbArray);
		free (log_Ls);        
	} else if (op->ps_set == NULL && (op->paramType==8 || op->paramType==9)) {
		double* log_Ls = matrix(n,1);
		vector<VariableSet*> vsArray;
		VariableSet* vs;
		int numLineTopMat = op->algn->numSeqs-1;
		f0 = fn(n,x,opt);
		for (i=0; i<n; i++) {
			x[i]+=DELTA;
			initialUpdate(op, x, n, NULL, op->num_chars);
			vs = new VariableSet(op->num_chars);
			vs->copyFrom(*op->vs);
			vsArray.push_back(vs);
			x[i]-=DELTA;
		}
		getLogLikelihood(op->ps->allCondProb_t, op->topMatrix, numLineTopMat,*(op->algn), &vsArray, op->num_chars, log_Ls);
		for (i=0; i<n; i++) {
			// beware that our target is to get the -ve of loglikehood value
			g[i] = (- log_Ls[i] - f0) / DELTA;
		}
		free (log_Ls);
		for (i=0; i<(int)vsArray.size(); i++)
			delete vsArray[i];
	} else if (op->ps_set != NULL && op->isShare==0 && op->paramType>=1 && op->paramType<=3) {
		f0 = fn(n,x,opt);
		for (i=0; i<n; i++) {
			x[i]+=DELTA;
			updateAllEigenAndCondProb(n, x, i-1, i, op);
			f1 = -getLogLikelihood(op->ps_set->allCondProbSet_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->HASsiteLikelihoods, op->currRateCatID, op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat);
			x[i]-=DELTA;
			g[i]=(f1-f0)/DELTA;
		}
	} else if (op->ps_set != NULL && op->paramType==4) {
		f0 = fn(n,x,opt);
		for (i=0; i<n; i++) {
			x[i]+=DELTA;
			updateAllEigenAndCondProb(n, x, i-1, i, op);
			f1 = -getLogLikelihood(op->ps_set->allCondProbSet_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->HASsiteLikelihoods, op->currRateCatID, op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat);
			x[i]-=DELTA;
			g[i]=(f1-f0)/DELTA;
		}
	} else if (op->ps_set != NULL && op->paramType>=8 && op->paramType<=9) {
		f0 = fn(n,x,opt);
		for (i=0; i<n; i++) {
			x[i]+=DELTA;
			initialUpdate(op, x, n, NULL, op->num_chars);
			f1 = -getLogLikelihood(op->ps_set->allCondProbSet_t, op->topMatrix, numLineTopMat,*(op->algn), *(op->vs), op->num_chars, op->HASsiteLikelihoods, (int)op->ps_set->allCondProbSet_t.size(), op->likelihoodMat, op->leftLikelihoodMat, op->rightLikelihoodMat);
			x[i]-=DELTA;
			g[i]=(f1-f0)/DELTA;
		}
	} else {
		f0 = fn(n,x,opt);
		for (i=0; i<n; i++) {
			x[i]+=DELTA;
			f1 = fn_update(n, x, i-1, i, opt);
			x[i]-=DELTA;
			g[i]=(f1-f0)/DELTA;
		}
	}
	
	// cout << "[quit compute_gr] op->paramType: " << op->paramType << endl << flush;
}
