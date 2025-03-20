/*
 *
 * parameters.cpp
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


#include "parameters.h"

VariableSet::VariableSet(int num_chars) {
	// constructor

	this->num_chars = num_chars;
	num_alpha = 0;
	beta = 0.0;
	alpha = NULL;
	probXGivenInv = NULL;
	rootNodeFreq = NULL;
	isGTRUpilson = 0;
	norootfreq = 0;

	// the following variables are for optimization
	alpha_t = NULL;
	probXGivenInv_t = NULL;
	rootNodeFreq_t = NULL;
	
	rootNodeFreq_alpha = NULL;
	invar_constant = NULL;

}


VariableSet::VariableSet(int num_alpha, int num_chars) {
	// constructor

	this->num_chars = num_chars;
	this->num_alpha = num_alpha;
	beta = 0.0;
	alpha = (double*) malloc (sizeof(double)*num_alpha);
	for (int i=0; i<num_alpha; i++)
		alpha[i] = 1.0/num_alpha;
	probXGivenInv = (double*) malloc (num_chars * sizeof(double));
	memset(probXGivenInv,0,num_chars * sizeof(double));
	// rootNodeFreq = (double*) malloc (num_chars*num_alpha*sizeof(double));
	rootNodeFreq = malloc_decimals (num_chars * num_alpha);
	// memset(rootNodeFreq,0,num_chars * sizeof(double)*num_alpha);
	reset_decimals(rootNodeFreq, num_chars*num_alpha);
	isGTRUpilson = 0;
	
	// the following variables are for optimization
	alpha_t = NULL;
	probXGivenInv_t = NULL;
	rootNodeFreq_t = NULL;
	rootNodeFreq_alpha = NULL;
	invar_constant = NULL;
}

VariableSet::~VariableSet() {
	// destructor
	if (alpha != NULL)
		free(alpha);
	if (probXGivenInv != NULL)
		free(probXGivenInv);
	if (rootNodeFreq != NULL)
		free_decimals(rootNodeFreq); //free(rootNodeFreq); // 
	
	if (alpha_t != NULL)
		free(alpha_t);
	if (probXGivenInv_t != NULL)
		free(probXGivenInv_t);
	if (rootNodeFreq_t != NULL)
		free(rootNodeFreq_t);
	if (rootNodeFreq_alpha != NULL)
		free_decimals(rootNodeFreq_alpha);
	if (invar_constant != NULL)
		free(invar_constant);
}

void VariableSet::setGTRUpilson() {
	// set "isGTRUpilson" to TRUE
	isGTRUpilson = 1;
}

void VariableSet::setNoRootFreq() {
	// set "norootfreq" to TRUE
	norootfreq = 1;
}

void VariableSet::readVariableFile(char* fileName) {
	// read the values of the variables from the file

	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the variable file :" << fileName << endl;
		exit(1);
	}
	string aline;
	int itemRead = 0;
	vector<string> token;
	while ((itemRead < 4) && getline(fin,aline)) {
		trim(aline);
		if ((aline.length() > 0) && (aline[0] != '#')) {
			// read the item
			itemRead++;

			switch(itemRead) {
			case 1:
				// beta
				beta = atof(aline.c_str());
				break;
			case 2:
				// alpha
				tokenizer(aline, ", \t", &token);
				num_alpha = (int)token.size();
				for (int i=0; i<num_alpha; i++)
					alpha[i] = atof(token[i].c_str());
				break;
			case 3:
				// probXGivenInv
				tokenizer(aline, ", \t", &token);
				if ((int)token.size() < num_chars) {
					cerr << "Error! The number of values for probXGivenInv is less than " << num_chars << " in the variable file : " << fileName << endl;
				}
				for (int i=0; i < num_chars; i++)
					probXGivenInv[i] = atof(token[i].c_str());
				break;
			case 4:
				// rootNodeFreq
				tokenizer(aline, ", \t", &token);
				if ((int)token.size() < num_chars) {
					cerr << "Error! The number of values for rootNodeFreq is less than " << num_chars << " in the variable file : " << fileName << endl;
				}
				for (int i=0; i < num_chars*num_alpha; i++)
					rootNodeFreq[i] = atof(token[i].c_str());
				break;
			}
		}
	}
	fin.close();
}

void VariableSet::resetAllVariables(Alignment& alignment) {
	// reset the variables according to the constant sites

	// set the value of beta to 0.75 * proportion of invariable sites
	if (isGTRUpilson)
		beta = 0.0;
	else
		beta = 0.75 * (double) alignment.numConstSites / (double) alignment.numSites;

	for (int i=0; i<num_alpha; i++) {
		alpha[i] = (1.0 - beta) / num_alpha;

		for (int k=0; k<num_chars; k++)
			rootNodeFreq[i*num_chars + k] = 1.0 / num_chars;
	}

	for (int i=0; i<num_chars; i++) {
		probXGivenInv[i] = 1.0 / num_chars;
	}

}


// set up the parameters for optimisation method
void VariableSet::set_for_opt() {
	// set up the parameters for optimisation method
	if (alpha_t == NULL) {
		alpha_t = (double*) malloc (sizeof(double) * num_alpha);
	}
	if (probXGivenInv_t == NULL) {
		probXGivenInv_t = (double*) malloc (num_chars * sizeof(double));
	}
	if (rootNodeFreq_t == NULL) {
		rootNodeFreq_t = (double*) malloc (num_chars * sizeof(double) * num_alpha);
	}
	
	beta_t = beta;
	for (int i=0; i<num_alpha; i++) {
		alpha_t[i] = alpha[i];
		for (int k=0; k<num_chars; k++)
			rootNodeFreq_t[i*num_chars + k] = rootNodeFreq[i*num_chars + k];
	}
	for (int i=0; i<num_chars; i++) {
		probXGivenInv_t[i] = probXGivenInv[i];
	}
}


void VariableSet::randomInitAllVariables() {
	// randomly initialize the variables


	// randomly initialize the value of alpha and beta
	int sum_x=0;
	int i,k;
	cout << "[initial values of alpha and beta] ";
	for (i=0; i<num_alpha; i++) {
		int x = genRandInt(1, 100-sum_x-(num_alpha-i));
		sum_x += x;
		alpha[i] = (double)x/100.0;
		cout << alpha[i] << ",";
	}
	beta = 1.0-(double)sum_x/100.0;
	cout << beta << endl;

	// randomly initialize the value of rootNodeFreq
	cout << "[initial values of rootNodeFreq]" << endl;
	for (i=0; i<num_alpha; i++) {
		int sum_y=0;
		for (k=0; k<num_chars-1; k++) {
			int y = genRandInt(1, 100-sum_y-(num_chars-k-1));
			sum_y+=y;
			rootNodeFreq[i*num_chars + k] = (double)y/100.0;
			cout << rootNodeFreq[i*num_chars + k] << ",";
		}
		rootNodeFreq[i*num_chars + num_chars-1] = 1.0 - (double)sum_y/100.0;
		cout << rootNodeFreq[i*num_chars + num_chars-1] << endl;
	}

	// randomly initialize the value of probXGivenInv
	int sum_z=0;
	cout << "[initial values of probXGivenInv] ";
	for (k=0; k<num_chars-1; k++) {
		int z = genRandInt(1, 100-sum_z-(num_chars-k-1));
		sum_z+=z;
		probXGivenInv[k] = (double)z/100.0;
		cout << probXGivenInv[k] << ",";
	}
	probXGivenInv[num_chars-1] = 1.0 - (double)sum_z/100.0;
	cout << probXGivenInv[num_chars-1] << endl;

}

// copy the whole content (except rateIDs) from another instance of ParameterSet
void VariableSet::copyFrom(VariableSet& vs) {
	num_chars = vs.num_chars;
	beta = vs.beta;
	num_alpha = vs.num_alpha;
	if (alpha!=NULL)
		free(alpha);
	alpha = (double*) malloc (sizeof(double)*num_alpha);
	memcpy(alpha,vs.alpha,sizeof(double)*num_alpha);
	if (probXGivenInv==NULL)
		probXGivenInv = (double*) malloc (num_chars * sizeof(double));
	memcpy(probXGivenInv,vs.probXGivenInv,num_chars * sizeof(double));
	if (rootNodeFreq!=NULL)
		free_decimals(rootNodeFreq); // free(rootNodeFreq); // 
	rootNodeFreq = malloc_decimals (num_chars * num_alpha);
	copy_decimals(rootNodeFreq,vs.rootNodeFreq,num_chars*num_alpha);
	// memcpy(rootNodeFreq,vs.rootNodeFreq,num_chars*num_alpha* sizeof(double));
}

void VariableSet::showContent(string &outStr, int format) {
	// output the content to the string
	// format: 0 - output to the checkpoint file
	//         1 - output to the HAS result file

	if (format==0) {
		outStr.append("beta : " + longDoublToStr(beta,PARAM_DECI) + "\n");
		outStr.append("alpha : ");
		for (int i=0; i<num_alpha; i++) {
			if (i > 0)
				outStr.append(",");
			outStr.append(longDoublToStr(alpha[i],PARAM_DECI));
		}
		outStr.append("\n");

		outStr.append("probXGivenInv : ");
		for (int i=0; i<num_chars; i++) {
			if (i > 0)
				outStr.append(",");
			outStr.append(longDoublToStr(probXGivenInv[i],PARAM_DECI));
		}
		outStr.append("\n");

		outStr.append("rootNodeFreq : ");
		for (int i=0; i<num_chars*num_alpha; i++) {
			if (i > 0)
				outStr.append(",");
			outStr.append(longDoublToStr(rootNodeFreq[i],PARAM_DECI));
		}
		outStr.append("\n");
	} else {
		outStr.append("[Nucleotide frequencies in root]\n");
		outStr.append("               \tfreq(A)\tfreq(C)\tfreq(G)\tfreq(T)\n");
		outStr.append("Invariable_site");
		for (int i=0; i<num_chars; i++) {
			outStr.append("\t");
			if (beta==0.0)
				outStr.append("N/A");
			else
				outStr.append(longDoublToStr(probXGivenInv[i],PARAM_DECI));
		}
		for (int i=0; i<num_chars*num_alpha; i++) {
			if (i%4==0)
				outStr.append("\nSite_category_" + intToStr((int)i/4+1));
			outStr.append("\t" + longDoublToStr(rootNodeFreq[i],PARAM_DECI));
		}
		outStr.append("\n\n");
		outStr.append("[Proportion of different sites]\n");
		outStr.append("Invariable_site\t" + longDoublToStr(beta,PARAM_DECI) + "\n");
		for (int i=0; i<num_alpha; i++) {
			outStr.append("Site_category_" + intToStr(i+1) + "\t" + longDoublToStr(alpha[i],PARAM_DECI) + "\n");
		}
	}
}

void VariableSet::showContent() {
	// show the content
	string str = "";
	showContent(str, 0);
	cout << str;
}

void VariableSet::showContent(ofstream* fout) {
	// show the content
	string str="";
	showContent(str, 0);
	(*fout) << str;
}

void VariableSet::roundContent() {
	// roundup all the values
	int i;
	beta = roundNumber(beta, PARAM_DECI);
	if (alpha != NULL) {
		for (i=0; i<num_alpha; i++)
			alpha[i] = roundNumber(alpha[i], PARAM_DECI);
	}
	if (probXGivenInv != NULL) {
		for (i=0; i<num_chars; i++)
			probXGivenInv[i] = roundNumber(probXGivenInv[i], PARAM_DECI);
	}
	if (rootNodeFreq != NULL) {
		for (i=0; i<num_chars; i++)
			rootNodeFreq[i] = roundNumber(rootNodeFreq[i], PARAM_DECI);
	}
}

void VariableSet::compute_rootNodeFreq_alpha() {
	// compute the array of multiplication of rootNodeFreq and alpha
	if (rootNodeFreq_alpha == NULL)
		rootNodeFreq_alpha = malloc_decimals(num_chars*num_alpha);
	for (int i=0; i<num_alpha; i++) {
		mul_decimals(&rootNodeFreq_alpha[i*num_chars], &rootNodeFreq[i*num_chars], alpha[i], num_chars);
	}
}

void VariableSet::compute_invar_constant() {
	// compute the array of invar_constant
	int i,j;
	if (invar_constant == NULL)
		invar_constant = (double*) malloc(num_chars * num_chars * sizeof(double));
	invar_constant[0] = 0.0;
	for (i=1; i<num_chars*num_chars; i++) {
		double var_inv = 0.0;
		for (j=0; j<num_chars; j++) {
			if ((i >> j) & 1)
				var_inv += probXGivenInv[j];
		}
		invar_constant[i] = beta * var_inv;
	}
}

void ParameterSet::readParamFile(char* fileName, int num_w, int num_edges) {
	// read the values of w and pi from the file
	// For reversible, there are six w values (i.e. num_w should be 6)

	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the file :" << fileName << endl;
		exit(1);
	}

	// allocate the space to the arrays
	this->num_w = num_w;
	this->num_edges = num_edges;
	if (w == NULL)
		w = malloc_decimals(num_edges*num_w); // matrix(num_edges, num_w);
	if (pi == NULL)
		pi = malloc_decimals(num_edges*num_chars); // matrix(num_edges, num_chars);
	if (t == NULL)
		t = malloc_decimals(num_edges); // matrix(num_edges, 1);

	// read the file
	vector<string> token;
	string aline;
	int lineNum=0;
	while(getline(fin,aline)) {
		trim(aline);
		if (aline.length() > 0) {
			tokenizer(aline, " \t,", &token);
			if ((int)token.size() < num_w+num_chars) {
				cerr << "Error! The number of items are too few in the following line inside the file " << fileName << " :" << endl;
				cerr << aline << endl;
				exit(1);
			}
			for (int i=0; i<num_w; i++) {
				w[lineNum*num_w+i] = atof(token[i].c_str());
			}
			for (int i=0; i<num_chars; i++) {
				pi[lineNum*num_chars+i] = atof(token[num_w+i].c_str());
			}
			t[lineNum] = atof(token[num_w+4].c_str());
			lineNum++;
		}
	}

	// initalize the rate matrix ID
	initializeRateMatIDFrParam();

	fin.close();
}

void ParameterSet::initialize(int num_w, int num_edges, int isSingleGrp) {
	// initialize the parameters
	// isSingleGrp : 1 - all edges are initialized as the same group
	//               0 - all edges are initizlized as different groups

	// allocate the space to the arrays
	this->num_w = num_w;
	this->num_edges = num_edges;

	if (w == NULL) {
		w = malloc_decimals(num_edges*num_w); // matrix(num_edges, num_w);
		// set all values of w to 1
		assign_decimals(w, num_w*num_edges, (double) 1.0);
		// setMatrixToVal(w, num_w, num_edges, 1.0);
	}
	if (pi == NULL) {
		pi = malloc_decimals(num_edges*num_chars); // matrix(num_edges, num_chars);
		// set all values of pi to 0.25
		assign_decimals(pi, num_chars*num_edges, (double) 0.25);
		// setMatrixToVal(pi, num_chars, num_edges, 0.25);
	}
	if (t == NULL) {
		t = malloc_decimals(num_edges); // matrix(num_edges, 1);
		// set all values of t to 1
		assign_decimals(t, num_edges, (double) 1.0);
		// setMatrixToVal(t, 1, num_edges, 1.0);
	}
	// initalize the rate matrix ID
	initializeRateMatID(isSingleGrp);
}

void ParameterSet::reset() {
	// reset the parameters

	// set all values of w to 1
	assign_decimals(w, num_w*num_edges, (double) 1.0);
	// setMatrixToVal(w, num_w, num_edges, 1.0);
	// set all values of pi to 0.25
	assign_decimals(pi, num_chars*num_edges, (double) 0.25);
	// setMatrixToVal(pi, num_chars, num_edges, 0.25);
	// set all values of t to 1
	assign_decimals(t, num_edges, (double) 1.0);
	// setMatrixToVal(t, 1, num_edges, 1.0);
}

void ParameterSet::set_for_opt() {
	// set up the parameters for optimisation method
	if (pi_t == NULL) {
		pi_t = malloc_decimals(num_chars * num_edges);
		// pi_t = new double[num_chars * num_edges];
	}
	int i;
	copy_decimals(pi_t, pi, num_chars*num_edges);
	// for (i=0; i<num_chars*num_edges; i++)
	//	pi_t[i] = pi[i];
}

void ParameterSet::randomInit(int num_w, int num_edges, int initSeed, int isSingleGrp) {
	// initialize the parameters randomly
	// but still the same set of values for each edge
	// isSingleGrp : 1 - all edges are initialized as the same group
	//               0 - all edges are initizlized as different groups

	this->num_w = num_w;
	this->num_edges = num_edges;

	if (initSeed)
		initialRandSeed(); // initialize the random seed

	int i,j,x,max,sum_x;
	if (w == NULL) {
		w = malloc_decimals(num_edges*num_w); // w = matrix(num_edges, num_w);
		// randomly set the values of w between 0.1 to 1.0
		for (i=0;i<num_w;i++) {
			x = genRandInt(1, 10);
			w[i] = (double)x/10.0;
		}
		for (i=1;i<num_edges;i++) {
			for (j=0;j<num_w;j++) {
				if (isSingleGrp)
					w[i*num_w+j] = w[j];
				else {
					x = genRandInt(1, 10);
					w[i*num_w+j] = (double)x/10.0;
				}
			}
		}
		// print out the values of w
		cout << "[initial values of w]" << endl;
		for (i=0;i<num_edges;i++) {
			for (j=0; j<num_w; j++) {
				if (j>0)
					cout << ",";
				cout << w[i*num_w+j];
			}
			cout << endl;
		}
	}
	if (pi == NULL) {
		pi = malloc_decimals(num_edges*num_chars); // pi = matrix(num_edges, num_chars);
		// randomly set the values of pi
		// (and for each edge, the last one is 1 - the sum of the last three)
		sum_x = 0;
		for (j=0;j<num_chars-1;j++) {
			max = 100-sum_x-(num_chars-j);
			x = genRandInt(1, max);
			sum_x += x;
			pi[j]=(double)x/100.0;
		}
		pi[num_chars-1]=1.0-(double)sum_x/100.0;
		for (i=1;i<num_edges;i++) {
			sum_x=0;
			for (j=0;j<num_chars;j++) {
				if (isSingleGrp) {
					pi[i*num_chars+j] = pi[j];
				} else {
					if (j==num_chars-1) {
						pi[i*num_chars+j]=1.0-(double)sum_x/100.0;
					} else {
						max = 100-sum_x-(num_chars-j);
						x = genRandInt(1, max);
						sum_x += x;
						pi[i*num_chars+j]=(double)x/100.0;
					}
				}
			}
		}
		// print out the value of pi
		cout << "[initial values of pi]" << endl;
		for (i=0;i<num_edges;i++) {
			for (j=0;j<num_chars;j++) {
				if (j>0)
					cout << ",";
				cout << pi[i*num_chars+j];
			}
			cout << endl;
		}
	}
	if (t == NULL) {
		t = malloc_decimals(num_edges); // t = matrix(num_edges, 1);
		cout << "[initial values of t]" << endl;
		// randomly set the values of t between 0.1 to 1.0
		for (i=0;i<num_edges;i++) {
			int x = genRandInt(1, 10);
			t[i] = (double)x/10.0;
			if (i>0)
				cout << ",";
			cout << w[i];
		}
		cout << endl;
	}
	// initalize the rate matrix ID
	initializeRateMatID(isSingleGrp);
}

ParameterSet::ParameterSet(int num_chars) {
	// constructor
	this->num_chars = num_chars;
	w = NULL;
	pi = NULL;
	pi_t = NULL;
	t = NULL;
	allEigenVals = NULL;
	allEigenVects_t = NULL;
	allCondProb_t=NULL;

	// for temporary usage
	eigenMat=NULL;
	s = NULL;
	tmpMat1 = NULL;
	tmpMat2 = NULL;
	sqrtPi = NULL;

	this->num_chars = num_chars;
}

ParameterSet::ParameterSet(int num_w, int num_edges, int num_chars) {
	// constructor
	this->num_w = num_w;
	this->num_edges = num_edges;
	this->num_chars = num_chars;
	w = malloc_decimals(num_edges*num_w);
	pi = malloc_decimals(num_chars*num_edges);
	pi_t = NULL;
	t = malloc_decimals(num_edges);
	allEigenVals = NULL;
	allEigenVects_t = NULL;
	allCondProb_t=NULL;

	// for temporary usage
	eigenMat=NULL;
	s = NULL;
	tmpMat1 = NULL;
	tmpMat2 = NULL;
	sqrtPi = NULL;
}

ParameterSet::~ParameterSet() {
	// destructor
	if (w!=NULL)
		free_decimals(w); // free(w);
	if (pi!=NULL)
		free_decimals(pi); // free(pi);
	if (pi_t!=NULL)
		free_decimals(pi_t); // free(pi_t);
	if (t!=NULL)
		free_decimals(t); // free(t);
	if (allEigenVals!=NULL)
		free_decimals(allEigenVals); // free(allEigenVals);
	if (allEigenVects_t!=NULL)
		free_decimals(allEigenVects_t); // free(allEigenVects);
	if (allCondProb_t!=NULL)
		free_decimals(allCondProb_t); // free(allCondProb);
	for (int i=0; i<(int)rateIDs.size(); i++) {
		delete(rateIDs[i]); // checked
	}

	// for temporary usage
	if (eigenMat!=NULL)
		free_decimals(eigenMat);
	if (s!=NULL)
		free_decimals(s);
	if (tmpMat1!=NULL)
		free_decimals(tmpMat1);
	if (tmpMat2!=NULL)
		free_decimals(tmpMat2);
	if (sqrtPi!=NULL)
		free_decimals(sqrtPi);

}

// update the parameters such that
// outFormat = 1 : Edge length is set to the rate of substitution; OR
//             2 : S6 in the rate matrix is set to 1
void ParameterSet::updateContent(int outFormat) {
	if (w!=NULL && pi!=NULL && t!=NULL) {
		for (int i=0; i<num_edges; i++) {
			if (num_w == 6) {
				double* curr_w = &(w[i*num_w]);
				double* curr_pi = &(pi[i*num_chars]);
				double base = 1.0;
				if (outFormat==1)
					base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
				// this formula is for num_w=6 and num_chars=4
				else if (outFormat==2)
					base = curr_w[num_w-1];
				for (int j=0; j<num_w; j++) {
					curr_w[j] = curr_w[j]/base;
				}
				t[i] = t[i]*base;
			}
		}
	}
}

// roundup the parameters based on the format
// which the edge length is set to the rate of substitution
void ParameterSet::roundContent() {
	if (w!=NULL && pi!=NULL && t!=NULL) {
		for (int i=0; i<num_edges; i++) {
			if (num_w == 6) {
				double* curr_w = &(w[i*num_w]);
				double* curr_pi = &(pi[i*num_chars]);
				double base = 1.0;
				// the following formula is for num_w=6 and num_chars=4
				base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
				for (int j=0; j<num_w; j++) {
					curr_w[j] = roundNumber(curr_w[j]/base, PARAM_DECI);
				}
				t[i] = roundNumber(t[i]*base,PARAM_DECI);
				for (int j=0; j<num_chars; j++) {
					curr_pi[j] = roundNumber(curr_pi[j],PARAM_DECI);
				}
			}
		}
	}
}


void ParameterSet::showContent(string &outStr, int* topMatrix, vector<string>* leafList) {

	if (topMatrix!=NULL) {
		if (w!=NULL) {
			// show the information of the rate matrix
			outStr.append("[GTR rate parameters and nucleotide distribution of edge leading to each node]\n");
			outStr.append("      \t");
			// outStr.append("S1\tS2\tS3\tS4\tS5\tS6");
            outStr.append("S_AC\tS_AG\tS_AT\tS_CG\tS_CT\tS_GT");
			outStr.append("\tPi_A\tPi_C\tPi_G\tPi_T\n");
			for (int i=0; i<num_edges; i++) {
				if (topMatrix!=NULL) {
					if (topMatrix[i] < 0) {
						// a leaf
						outStr.append(leafList->at(-topMatrix[i]-1));
					} else {
						// an internal node
						outStr.append(intToStr(topMatrix[i]));
					}
				}
				for (int j=0; j<num_w; j++) {
					outStr.append("\t" + longDoublToStr(w[i*num_w+j],PARAM_DECI));
				}
				for (int j=0; j<num_chars; j++) {
					outStr.append("\t" + longDoublToStr(pi[i*num_chars+j],PARAM_DECI));
				}
				outStr.append("\n");
			}
		}
		outStr.append("\n");
		/*
		if (pi!=NULL) {
			// show the information of Nucleotide distribution
			outStr.append("[Nucleotide distribution of edge leading to each node]\n");
			if (topMatrix!=NULL)
				outStr.append("      \t");
			outStr.append("Pi_1\tPi_2\tPi_3\tPi_4\n");
			for (int i=0; i<num_edges; i++) {
				if (topMatrix!=NULL) {
					if (topMatrix[i] < 0) {
						// a leaf
						outStr.append(leafList->at(-topMatrix[i]-1));
					} else {
						// an internal node
						outStr.append(intToStr(topMatrix[i]));
					}
				}
				for (int j=0; j<4; j++) {
					outStr.append("\t" + longDoublToStr(pi[i*4+j],PARAM_DECI));
				}
				outStr.append("\n");
			}
		}
		 */
		outStr.append("\n");
	} else {
		if (w!=NULL && pi!=NULL && t!=NULL) {
			for (int i=0; i<num_edges; i++) {
				for (int j=0; j<num_w; j++) {
					outStr.append(longDoublToStr(w[i*num_w+j],PARAM_DECI) + "\t");
				}
				for (int j=0; j<num_chars; j++) {
					outStr.append(longDoublToStr(pi[i*num_chars+j],PARAM_DECI) + "\t");
				}
				outStr.append(longDoublToStr(t[i],PARAM_DECI) + "\n");
			}
		}
	}
}

void ParameterSet::showContent() {
	string str = "";
	showContent(str, NULL, NULL);
	cout << str;
}

void ParameterSet::showContent(ofstream* fout) {
	string str = "";
	showContent(str, NULL, NULL);
	(*fout) << str;
}

struct paramComp {
	double* w;
	double* pi;
	int num_w;
	bool operator() (const paramComp& edge1, const paramComp& edge2) const {
		int i=0;
		while (i<num_w && edge1.w[i]==edge2.w[i])
			i++;
		if (i<num_w)
			return (edge1.w[i]<edge2.w[i]);
		int k=0;
		while (k<4 && edge1.pi[k]==edge2.pi[k])
			k++;
		if (k<4)
			return (edge1.pi[i]<edge2.pi[k]);
		else
			return false;
	}
};

// initialize the rate matrix for each edge
// isSingleGrp : 1 - all edges are initialized as the same group
//               0 - all edges are initizlized as different groups
void ParameterSet::initializeRateMatID(int isSingleGrp) {
	rateIDs.clear();
	if (isSingleGrp==1) {
		rateIDs.push_back(new vector<int>);
		for (int i=0; i<num_edges; i++)
			rateIDs[0]->push_back(i);
	} else {
		for (int i=0; i<num_edges; i++) {
			rateIDs.push_back(new vector<int>);
			rateIDs[i]->push_back(i);
		}
	}
	num_rate_matrices = (int) rateIDs.size();
}

// get the initialize ID of rate matrices for each edge
// before: the values of the parameters (i.e. w and pi) have been read from the file
void ParameterSet::initializeRateMatIDFrParam() {
	map<paramComp,int,paramComp> paramList;
	map<paramComp,int,paramComp>::iterator itr;
	for (int i=0; i<num_edges; i++) {
		paramComp currParam;
		currParam.num_w = num_w;
		currParam.w = &(w[num_w * i]);
		currParam.pi = &(pi[4 * i]);
		itr = paramList.find(currParam);
		if (itr!=paramList.end()) {
			// found
			rateIDs[itr->second]->push_back(i);
		} else {
			// new
			int newID = (int) paramList.size();
			rateIDs.push_back(new vector<int>);
			rateIDs[newID]->push_back(i);
			paramList.insert(pair<paramComp,int>(currParam,newID));
		}
	}
	num_rate_matrices = (int) paramList.size();
}

void ParameterSet::printRateMatID() {
	for (int i=0; i<(int)rateIDs.size(); i++) {
		cout << "Rate matrix " << i << " :";
		for (int j=0; j<(int)rateIDs[i]->size(); j++) {
			cout << " " << rateIDs[i]->at(j);
		}
		cout << endl;
	}
}

// output the rate matrice in single array form
void ParameterSet::outputRateMat(vector<int>& outRateMatArray) {
	outRateMatArray.assign(num_edges,0);
	for (int i=0; i<(int)rateIDs.size(); i++) {
		for (int j=0; j<(int)rateIDs[i]->size(); j++) {
			outRateMatArray[rateIDs[i]->at(j)]=i+1;
		}
	}
}

// update the rate matrix
void ParameterSet::updateRateMat(vector<int>& rateMatArray, int numRateGrp) {
	// clear the rateIDs
	int i,j,k;
	for (i=0; i<(int)rateIDs.size(); i++) {
		rateIDs[i]->clear();
		delete(rateIDs[i]); // checked
	}
	rateIDs.clear();
	// initialize the rateIDs
	for (i=0; i<numRateGrp; i++) {
		vector<int> * new_vector = new vector<int>;
		rateIDs.push_back(new_vector);
	}
	// set the rateIDs
	for (i=0; i<(int)rateMatArray.size(); i++) {
		rateIDs[rateMatArray[i]-1]->push_back(i);
	}
	num_rate_matrices=numRateGrp;
	// update the parameters according to the arrangement of the rate matrices
	int f,g;
	for (i=0; i<numRateGrp; i++) {
		if (rateIDs[i]->size()>1) {
			f = rateIDs[i]->at(0);
			for (j=1; j<(int)rateIDs[i]->size(); j++) {
				g = rateIDs[i]->at(j);
				// w's of edge g are set to be the same as w's of edge f
				for (k=0; k<num_w; k++) {
					w[g*num_w+k] = w[f*num_w+k];
				}
				// pi's of edge g are set to be the same as pi's of edge f
				for (k=0; k<num_chars; k++) {
					pi[g*num_chars+k] = pi[f*num_chars+k];
				}
			}
		}
	}
}

// load the rate matrix from the file
void ParameterSet::loadRateMat(char* rateGrpFile, int* topMatrix, vector<string>* leafList) {

	int numRateGrp;
	vector<int>* rateMatArray = collectRateGrp(rateGrpFile, topMatrix, leafList, numRateGrp);
	updateRateMat(*rateMatArray, numRateGrp);
	delete rateMatArray;
}

// import the rate matrix
// prerequisite: the value of "num_edges" has to be set
void ParameterSet::loadRateMat(vector<int>& rateMat) {
	int maxRateGrp = 0;
	for (int i=0; i<num_edges; i++) {
		if (maxRateGrp < rateMat[i])
			maxRateGrp = rateMat[i];
	}

	updateRateMat(rateMat, maxRateGrp);
}

// import the rate matrix string
// prerequisite: the value of "num_edges" has to be set
void ParameterSet::loadRateMat(string rateMatStr) {
    vector<int> rateMat;
    size_t p1 = 0;
    size_t p2 = rateMatStr.find_first_of(',', p1);
    while (p2 != string::npos && p2 > p1) {
        int r = atoi(rateMatStr.substr(p1,p2-p1).c_str());
        rateMat.push_back(r);
        p1 = p2 + 1;
        p2 = rateMatStr.find_first_of(',', p1);
    }
    if (p1 < rateMatStr.length()) {
        int r = atoi(rateMatStr.substr(p1,rateMatStr.length()-p1).c_str());
        rateMat.push_back(r);
    }
    if (rateMat.size() != num_edges) {
        cerr << "Error! Incorrect size of the input starting rate matrix" << endl;
        cerr << "Input starting rate matrix: " << rateMatStr << endl;
        cerr << "number of edges: " << num_edges << " but number of items inside the input rate matrix: " << rateMat.size() << endl;
        exit(1);
    }
    loadRateMat(rateMat);
}

// compute the eigenvalues and eigenvectors for the edges
// if "theEdge" is -1, then do for all edges
// output: allEigenVals and allEigenVects
void ParameterSet::computeAllEigenValues_t(int theEdge) {

	// allocate memory to the matrices
	if (allEigenVals==NULL)
		allEigenVals = malloc_decimals(num_edges*4); // matrix(num_edges,4);
	if (allEigenVects_t==NULL)
		allEigenVects_t = malloc_decimals(num_edges*16); // matrix(num_edges,16);
	if (eigenMat==NULL)
		eigenMat = malloc_decimals(16); // matrix(4,4);
	if (s==NULL)
		s = malloc_decimals(16); // matrix(4,4);

	vector<int>* edgeList;
	double* curr_w;
	double* curr_pi;
	double* currEigenVals;
	double* currEigenVects_t;
	double* firstEigenVals;
	double* firstEigenVects_t;
	int rateGrp;
	int edge;
	int numRound, rot_num;
	int i;
	if (theEdge == -1) {
		for (rateGrp=0; rateGrp<(int)rateIDs.size(); rateGrp++) {
			edgeList = rateIDs[rateGrp];
			if (edgeList != NULL && edgeList->size() > 0) {
				// for the first item
				edge = edgeList->at(0);
				curr_w = &(w[edge*num_w]);
				curr_pi = &(pi[edge*4]);
				computeSMatrix(s, curr_w, curr_pi);
				computeEigenMatrix(eigenMat, curr_w, curr_pi, s);
				/*
				cout << "eigen-matrix:" << endl;
				for (int m=0; m<4; m++) {
					for (int n=0; n<4; n++) {
						printf("%6.4f ", eigenMat[m*4+n]);
					}
					printf("\n");
				}*/
				firstEigenVals = &(allEigenVals[edge*4]);
				firstEigenVects_t = &(allEigenVects_t[edge*16]);
				// compute the eigenvalues and transpose of eigenvectors
				jacobi_eigenvalue ( 4, eigenMat, MAX_ROUND, firstEigenVects_t, firstEigenVals, &numRound, &rot_num );
				/*
				cout << "firstEigenVals:" << endl;
				for (int m=0; m<4; m++) {
					printf("%6.4f ", firstEigenVals[m]);
				}
				printf("\n");
				cout << "firstEigenVects:" << endl;
				for (int m=0; m<4; m++) {
					for (int n=0; n<4; n++) {
						printf("%6.4f ", firstEigenVects_t[n*4+m]);
					}
					printf("\n");
				}*/
				// computeEigenVect(eigenMat, 4, firstEigenVals, firstEigenVects, numRound);
				// for the rest of the items
				for (i=1; i<(int)edgeList->size(); i++) {
					edge = edgeList->at(i);
					currEigenVals = &(allEigenVals[edge*4]);
					currEigenVects_t = &(allEigenVects_t[edge*16]);
					copy_decimals(currEigenVals, firstEigenVals, 4);
					copy_decimals(currEigenVects_t, firstEigenVects_t, 16);
				}
			}
		}
	} else {
		// just for theEdge
		curr_w = &(w[theEdge*num_w]);
		curr_pi = &(pi[theEdge*4]);
		computeSMatrix(s, curr_w, curr_pi);
		computeEigenMatrix(eigenMat, curr_w, curr_pi, s);
		currEigenVals = &(allEigenVals[theEdge*4]);
		currEigenVects_t = &(allEigenVects_t[theEdge*16]);
		// compute the eigenvalues and transpose of eigenvectors
		jacobi_eigenvalue ( 4, eigenMat, MAX_ROUND, currEigenVects_t, currEigenVals, &numRound, &rot_num );
		// computeEigenVect(eigenMat, 4, currEigenVals, currEigenVects, numRound);
	}
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_eigen += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif

}

// compute the eigenvalues and eigenvectors for the edges of a particular rate group
// output: allEigenVals and allEigenVects
void ParameterSet::computeAllEigenValuesForRateGrp_t(int rateGrp) {
	
	// allocate memory to the matrices
	if (allEigenVals==NULL)
		allEigenVals = malloc_decimals(num_edges*4); // matrix(num_edges,4);
	if (allEigenVects_t==NULL)
		allEigenVects_t = malloc_decimals(num_edges*16); // matrix(num_edges,16);
	if (eigenMat==NULL)
		eigenMat = malloc_decimals(16); // matrix(4,4);
	if (s==NULL)
		s = malloc_decimals(16); // matrix(4,4);
	
	vector<int>* edgeList;
	double* curr_w;
	double* curr_pi;
	double* currEigenVals;
	double* currEigenVects_t;
	double* firstEigenVals;
	double* firstEigenVects_t;
	int edge;
	int numRound, rot_num;
	int i;
	edgeList = rateIDs[rateGrp];
	if (edgeList != NULL && edgeList->size() > 0) {
		// for the first item
		edge = edgeList->at(0);
		curr_w = &(w[edge*num_w]);
		curr_pi = &(pi[edge*4]);
		computeSMatrix(s, curr_w, curr_pi);
		computeEigenMatrix(eigenMat, curr_w, curr_pi, s);
		firstEigenVals = &(allEigenVals[edge*4]);
		firstEigenVects_t = &(allEigenVects_t[edge*16]);
		// compute the eigenvalues and transpose of eigenvectors
		jacobi_eigenvalue ( 4, eigenMat, MAX_ROUND, firstEigenVects_t, firstEigenVals, &numRound, &rot_num );
		// computeEigenVect(eigenMat, 4, firstEigenVals, firstEigenVects, numRound);
		// for the rest of the items
		for (i=1; i<(int)edgeList->size(); i++) {
			edge = edgeList->at(i);
			currEigenVals = &(allEigenVals[edge*4]);
			currEigenVects_t = &(allEigenVects_t[edge*16]);
			copy_decimals(currEigenVals, firstEigenVals, 4);
			// memcpy(currEigenVals,firstEigenVals,4*sizeof(double));
			copy_decimals(currEigenVects_t, firstEigenVects_t, 16);
			// memcpy(currEigenVects,firstEigenVects,16*sizeof(double));
		}
	}
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_eigen += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
}


// OBJECTIVE:             compute the conditional probabilities for all the edges.
// OUTPUT:                allCondProb
// SIMILAR FUNCTION IN R: getCondProb
// PREREQUISITE:          computeAllEigenMatrix()
void ParameterSet::computeAllCondProb(int isReversible) {
	
	if (!isReversible) {
		cerr << "Error! the function ""computeCondProb"" cannot support not isreversible condition yet." << endl;
	}
	if (allCondProb_t==NULL) {
		allCondProb_t=malloc_decimals(num_edges*16); //matrix(num_edges, 16);
	}

	if (tmpMat1==NULL)
		tmpMat1 = malloc_decimals(16); //matrix(4,4);
	if (sqrtPi==NULL)
		sqrtPi = malloc_decimals(4); //matrix(1,4);

	for (int edge=0; edge<num_edges; edge++) {
		double* condProb_t = &(allCondProb_t[edge*16]);
		double* currPi = &(pi[edge*4]);
		double time = t[edge];
		double* eigenVals = &(allEigenVals[edge*4]);
		double* eigenVects_t = &(allEigenVects_t[edge*16]);

		// transpose of condProb = eigen-vector * transpose of {eigen-vector * diag(exp(eigen-values * time))}
		// tmpMat1 = transpose of {eigen-vector * diag(exp(eigen-values * time))}
		for (int i=0; i<4; i++)
			mul_decimals(&tmpMat1[i*4], &eigenVects_t[i*4], (double) exp(eigenVals[i] * time), 4);

		multiply44_t(condProb_t, eigenVects_t, tmpMat1);

		// sqrtPi = sqrt(Pi)
		sqrt_decimals(sqrtPi, currPi, 4);

		// condProb = diag(invserse(sqrtPi)) * condProb * diag(sqrtPi)
		for (int i=0; i<4; i++) {
			div_decimals(&condProb_t[i*4], &condProb_t[i*4], sqrtPi, 4);
			mul_decimals(&condProb_t[i*4], &condProb_t[i*4], sqrtPi[i], 4);
		}
		
		min_decimals(condProb_t, MIN_PROB_VALUE, 16);
		/*
		// PRINT the condProb
		if (edge==0) {
			cout << "condProb (after min):" << endl;
			for (int m=0; m<4; m++) {
				for (int n=0; n<4; n++) {
					printf("%6.4f ", condProb_t[n*4+m]);
				}
				printf("\n");
			}
		}*/
	}
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_prob += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
}

// OBJECTIVE:             compute the conditional probabilities for a particular edge
// OUTPUT:                allCondProb
// SIMILAR FUNCTION IN R: getCondProb
// PREREQUISITE:          computeAllEigenMatrix()
void ParameterSet::computeCondProbForEdge(int edge) {
	if (allCondProb_t==NULL) {
		allCondProb_t=malloc_decimals(num_edges*16); //matrix(num_edges, 16);
	}
	
	if (tmpMat1==NULL)
		tmpMat1 = malloc_decimals(16); //matrix(4,4);
	if (sqrtPi==NULL)
		sqrtPi = malloc_decimals(4); //matrix(1,4);
	
	double* condProb_t = &(allCondProb_t[edge*16]);
	double* currPi = &(pi[edge*4]);
	double time = t[edge];
	double* eigenVals = &(allEigenVals[edge*4]);
	double* eigenVects_t = &(allEigenVects_t[edge*16]);
	
	// transpose of condProb = eigen-vector * transpose of {eigen-vector * diag(exp(eigen-values * time))}
	// tmpMat1 = transpose of {eigen-vector * diag(exp(eigen-values * time))}
	for (int i=0; i<4; i++)
		mul_decimals(&tmpMat1[i*4], &eigenVects_t[i*4], (double) exp(eigenVals[i] * time), 4);
	
	multiply44_t(condProb_t, eigenVects_t, tmpMat1);
	//multiplyQuick44_t(condProb_t, eigenVects_t, tmpMat1);
	
	// sqrtPi = sqrt(Pi)
	sqrt_decimals(sqrtPi, currPi, 4);

	// condProb = diag(invserse(sqrtPi)) * condProb * diag(sqrtPi)
	for (int i=0; i<4; i++) {
		div_decimals(&condProb_t[i*4], &condProb_t[i*4], sqrtPi, 4);
		mul_decimals(&condProb_t[i*4], &condProb_t[i*4], sqrtPi[i], 4);
	}
	
	min_decimals(condProb_t, MIN_PROB_VALUE, 16);
	
#ifdef GET_TIME_STAT_DETAIL
	curr_time = clock();
	compute_prob += (curr_time - pre_optim_time)* 1000.0 / CLOCKS_PER_SEC;
	pre_optim_time = curr_time;
#endif
}

// to normalize S and the edge length
// so that the edge length equals to the number of substitution
// only for number of w = 6
void ParameterSet::normalize() {
	for (int i=0; i<num_edges; i++) {
		double* curr_w = &(w[i*6]);
		double* curr_pi = &(pi[i*4]);
		double base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
		for (int j=0; j<6; j++) {
			curr_w[j] = curr_w[j] / base;
		}
		t[i] = t[i] * base;
	}
}

void ParameterSet::printAllEigens() {
	if (allEigenVects_t != NULL) {
		printf("allEigenVects_t\n");
		for (int i=0; i<num_edges; i++) {
			printf("Edge %d\n", i);
			for (int j=0; j<4; j++) {
				for (int k=0; k<4; k++) {
					printf("%6.4f ", allEigenVects_t[i*16+j*4+k]);
				}
				printf("\n");
			}
		}
		printf("allEigenVals\n");
		for (int i=0; i<num_edges; i++) {
			printf("Edge %d\n", i);
			for (int j=0; j<4; j++) {
				printf("%6.4f ", allEigenVals[i*4+j]);
			}
			printf("\n");
		}
	}
}


void ParameterSet::printAllCondProb() {
	if (allCondProb_t != NULL) {
		// print out in R format
		for (int k=0; k<4; k++) {
			cout << endl << ", , " << k+1 << endl;
			cout << "\t[,1]\t[,2]\t[,3]\t[,4]" << endl;
			for (int i=0; i<num_edges; i++) {
				cout << "[" << i+1 << ",]";
				for (int j=0; j<4; j++) {
					cout << "\t" << allCondProb_t[i*16+k*4+j];
				}
				cout << endl;
			}
		}
	}
}

void ParameterSet::print_t() {
	for (int i=0; i<num_edges; i++) {
		cout << t[i] << " ";
	}
	cout << endl;
}


// copy the whole content (except rateIDs) from another instance of ParameterSet
void ParameterSet::copyFrom(ParameterSet& ps) {

	num_chars = ps.num_chars;
	num_w = ps.num_w;
	num_edges = ps.num_edges;

	if (w == NULL) {
		w = malloc_decimals(num_edges*num_w); //matrix(num_edges, num_w);
	}
	copy_decimals(w, ps.w, num_w*num_edges);
	// memcpy(w,ps.w,num_w*num_edges*sizeof(double));

	if (pi == NULL) {
		pi = malloc_decimals(num_edges*num_chars); //matrix(num_edges, num_chars);
	}
	copy_decimals(pi, ps.pi, num_edges*num_chars);
	// memcpy(pi,ps.pi,num_edges*num_chars*sizeof(double));

	if (t == NULL) {
		t = malloc_decimals(num_edges); //matrix(num_edges, 1);
	}
	copy_decimals(t, ps.t, num_edges);
	// memcpy(t,ps.t,num_edges*sizeof(double));

	if (ps.allEigenVals!=NULL) {
		if (allEigenVals==NULL) {
			allEigenVals = malloc_decimals(num_edges*num_chars); //matrix(num_edges,num_chars);
		}
		copy_decimals(allEigenVals,ps.allEigenVals,num_edges*num_chars);
		// memcpy(allEigenVals,ps.allEigenVals,num_edges*num_chars*sizeof(double));
	}

	if (ps.allEigenVects_t!=NULL) {
		if (allEigenVects_t==NULL) {
			allEigenVects_t = malloc_decimals(num_edges*num_chars*num_chars); //matrix(num_edges,num_chars*num_chars);
		}
		copy_decimals(allEigenVects_t,ps.allEigenVects_t,num_edges*num_chars*num_chars);
		// memcpy(allEigenVects,ps.allEigenVects,num_edges*num_chars*num_chars*sizeof(double));
	}

	if (ps.allCondProb_t!=NULL) {
		if (allCondProb_t==NULL) {
			allCondProb_t=malloc_decimals(num_edges*num_chars*num_chars); //matrix(num_edges, num_chars*num_chars);
		}
		copy_decimals(allCondProb_t,ps.allCondProb_t,num_edges*num_chars*num_chars);
		// memcpy(allCondProb,ps.allCondProb,num_edges*num_chars*num_chars*sizeof(double));
	}
}

// constructor
AllParameterSet::AllParameterSet() {
	this->numRateCat = 0;
}


// constructor
AllParameterSet::AllParameterSet(int numRateCat, int num_chars) {
	this->numRateCat = numRateCat;
	// then create the ParameterSet objects
	for (int i=0; i<numRateCat; i++) {
		ParameterSet* a_ps = new ParameterSet(num_chars);
		ps.push_back(a_ps);
	}
}

// constructor
AllParameterSet::AllParameterSet(int numRateCat, int num_chars, int num_edges) {
	this->numRateCat = numRateCat;
	// then create the ParameterSet objects
	int num_w = 6;
	for (int i=0; i<numRateCat; i++) {
		ParameterSet* a_ps = new ParameterSet(num_w, num_edges, num_chars);
		ps.push_back(a_ps);
	}
}

// destructor
AllParameterSet::~AllParameterSet() {
	for (int i=0; i<numRateCat; i++) {
		delete (ps[i]); // checked
	}
	ps.clear();
}


// read the values of w and pi from all the parameter files
void AllParameterSet::readParamFile(char* paramFile, int num_w, int num_edges) {
    int n = ((int)log10(numRateCat))+strlen(paramFile)+8;
	char* parameterFileName = (char*) malloc (sizeof(char) * n);
	for (int i=0; i<numRateCat; i++) {
		snprintf(parameterFileName,n,"%s%i.txt",paramFile,i+1);
		ps[i]->readParamFile(parameterFileName, num_w, num_edges);
	}
	free(parameterFileName);
}


// load the file of rate groups
// return the number of rate groups
int AllParameterSet::loadRateMat(char* rateGrpFile, int* topMatrix, vector<string>* leafList) {
	if (numRateCat == 0) {
		cerr << "[allParameterSet::loadRateMat] Error! The number of rate category is 0!" << endl;
	}
	for (int i=0; i<numRateCat; i++) {
		ps[i]->loadRateMat(rateGrpFile, topMatrix, leafList);
	}
	return ps[0]->num_rate_matrices;;
}

// update the rate matrix
void AllParameterSet::updateRateMat(vector<int>& rateMatArray, int numRateGrp) {
	for (int i=0; i<numRateCat; i++) {
		ps[i]->updateRateMat(rateMatArray, numRateGrp);
	}
}

void AllParameterSet::showContent() {
	for (int i=0; i<numRateCat; i++) {
		cout << "[parameter set " << i << "]" <<  endl;
		ps[i]->showContent();
	}
}

void AllParameterSet::showContent(ofstream* fout) {
	for (int i=0; i<numRateCat; i++) {
		(*fout) << "[parameter set " << i << "]" <<  endl;
		ps[i]->showContent(fout);
	}
}

void AllParameterSet::showContent(string& outStr, int* topMatrix, vector<string>* leafList) {
	for (int i=0; i<numRateCat; i++) {
		outStr.append("--------------------------------\n");
		outStr.append("Site category " + intToStr(i+1) + "\n");
		outStr.append("--------------------------------\n");
		ps[i]->showContent(outStr, topMatrix,  leafList);
	}
}

void AllParameterSet::updateContent(int outFormat) {
	for (int i=0; i<numRateCat; i++) {
		ps[i]->updateContent(outFormat);
	}
}


void AllParameterSet::computeAllEigenValues_t() {
	for (int i=0; i<numRateCat; i++) {
		ps[i]->computeAllEigenValues_t(-1);
	}
}

void AllParameterSet::computeAllCondProb(int isReversible) {
	allCondProbSet_t.clear();
	for (int i=0; i<numRateCat; i++) {
		ps[i]->computeAllCondProb(isReversible);
		allCondProbSet_t.push_back(ps[i]->allCondProb_t);
	}
}

/*
// to normalize S and the edge length
// so that when the number of substitution equals to the edge length
void AllParameterSet::normalizeAll() {
	for (int i=0; i<numRateCat; i++) {
		ps[i]->normalize();
	}
}
 */

void AllParameterSet::printAllEigens() {
	for (int p=0; p<numRateCat; p++) {
		cout << "Rate Category " << p+1 << endl;
		ps[p]->printAllEigens();
	}
}

void AllParameterSet::printAllCondProb() {
	for (int k=0; k<4; k++) {
		for (int j=0; j<4; j++) {
			cout << endl;
			cout << ", ," << j+1 << "," << k+1 << endl;
			cout << endl;
			for (int s=0; s<2; s++) {
				for (int i=s*7; i<(s+1)*7; i++) {
					cout << "\t[," << i+1 << "]";
				}
				cout << endl;
				for (int p=0; p<numRateCat; p++) {
					cout << "[" << p+1 << ",]";
					for (int i=s*7; i<(s+1)*7; i++) {
						cout << "\t" << ps[p]->allCondProb_t[i*16+k*4+j];
					}
					cout << endl;
				}
			}
		}
	}
}

void AllParameterSet::reset() {
	// reset the parameters
	// set all values of w to 1
	// set all values of pi to 0.25
	// set all values of t to 1
	for (int i=0; i<numRateCat; i++) {
		ps[i]->reset();
	}
}


// set up the parameters for optimisation method
void AllParameterSet::set_for_opt() {
	for (int i=0; i<numRateCat; i++) {
		ps[i]->set_for_opt();
	}
}

void AllParameterSet::initialize(int num_w, int num_edges) {
	// initialize all the parameters
	int isSingleGrp = 1;
	for (int i=0; i<numRateCat; i++) {
		ps[i]->initialize(num_w, num_edges, isSingleGrp);
	}
}

void AllParameterSet::randomInit(int num_w, int num_edges) {
	// initialize all the parameters randomly
	int isSingleGrp = 1;
	for (int i=0; i<numRateCat; i++) {
		ps[i]->randomInit(num_w, num_edges, 0, isSingleGrp);
	}
}

int AllParameterSet::size() {
	return (int) ps.size();
}

// copy the whole content (except rateIDs) from another instance of AllParameterSet
void AllParameterSet::copyFrom(AllParameterSet& all_ps) {
	numRateCat = all_ps.numRateCat;
	ps.resize(numRateCat, NULL);
	int i;
	for (i=0; i<numRateCat; i++) {
		if (ps[i]==NULL) {
			ps[i] = new ParameterSet(all_ps.ps[i]->num_chars);
		}
		ps[i]->copyFrom(*(all_ps.ps[i]));
	}
}

// compute the s-matrix
// output: s_matrix
void computeSMatrix(double* s_matrix, double* w, double * pi) {

#ifdef DEBUG_MODE
	// list the values of w
	cout << "w: ";
	for (int i=0; i<6; i++)
		cout << w[i] << ", ";
	cout << endl;
	// list the value of pi
	cout << "pi: ";
	for (int i=0; i<4; i++)
		cout << pi[i] << ", ";
	cout << endl;
#endif

	// compute the s-matrix (dimension: 4 x 4)
	for (int i=0; i<3; i++)
		s_matrix[0*4+(i+1)]=s_matrix[(i+1)*4+0]=w[i]; //w1,w2,w3
	s_matrix[1*4+2]=s_matrix[2*4+1]=w[3]; //w4
	s_matrix[1*4+3]=s_matrix[3*4+1]=w[4]; //w5
	s_matrix[2*4+3]=s_matrix[3*4+2]=w[5]; //w6
	for (int i=0; i<4; i++) {
		s_matrix[i*4+i]=0.0;
		s_matrix[i*4+i]=(- sumProduct(&(s_matrix[i*4]), pi, 4)) / pi[i];
	}
}

// compute the eigen-matrix
// output: eigenMat
void computeEigenMatrix(double* eigenMat, double* w, double* pi, double* s) {
	int k = 0;
	for (int i=0; i<3; i++) {
		for (int j=i+1; j<4; j++) {
			// (i,j) = sqrt(pi[i] * pi[j]) * w[k]
			eigenMat[i*4+j] = sqrt(pi[i] * pi[j]) * w[k++];
		}
	}
	for (int i=0; i<4; i++)
		// (i,i) = pi[i] * s[i][j]
		eigenMat[i*4+i] = pi[i] * s[i*4+i];
	for (int i=1; i<4; i++)
		for (int j=0; j<i; j++)
			// (i,j) = (j,i)
			eigenMat[i*4+j] = eigenMat[j*4+i];
}

//=================================================================
// for reading the new format of parameter and variable files
//=================================================================

void VariableSet::readSiteInfoFile(char* fileName, vector<string>* siteCatNames) {
	// read the values of the variables from the site info file

	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the variable file :" << fileName << endl;
		exit(1);
	}
	string aline;
	vector<string> token;
	string catName;
	int isVariant;
	int catID=0;
	int* catLoaded = new int[num_alpha+1]; // catLoaded[num_alpha] is for the invariant site
	memset(catLoaded, 0, sizeof(int)*(num_alpha+1));
	while (getline(fin,aline)) {
		if ((aline.length() > 0) && (aline[0] != '#')) {
			tokenizer(aline," \t", &token);
			if (token.size() < 7)
				continue;
			catName = token[0];
			if (token[1] == "invariant") {
				isVariant = 0;
			} else if (token[1] == "variant") {
				isVariant = 1;
				// check for the category name
				for (catID=0; catID<(int)siteCatNames->size(); catID++)
					if (catName == siteCatNames->at(catID))
						break;
				if (catID >= (int)siteCatNames->size()) {
					cerr << "Error! The variant category " << catName << " does not appear in the tree file" << endl;
					exit(1);
				}
			} else {
				cerr << "Error! No indication on whether the category " << catName << " is variant or invariant in the \"site info file\"" << endl;
				exit(1);
			}
			if (isVariant) {
				if (catLoaded[catID]) {
					cerr << "Error! The variant category " << catName << " appears twice in the \"site info file\"" << endl;
					exit(1);
				}
				alpha[catID] = atof(token[2].c_str()); // proportion
				rootNodeFreq[catID*4  ] = atof(token[3].c_str()); // freq(A)
				rootNodeFreq[catID*4+1] = atof(token[4].c_str()); // freq(C)
				rootNodeFreq[catID*4+2] = atof(token[5].c_str()); // freq(G)
				rootNodeFreq[catID*4+3] = atof(token[6].c_str()); // freq(T)
				catLoaded[catID] = 1;
			} else {
				if (catLoaded[num_alpha]) {
					cerr << "Error! There are more than on invariant category in the \"site info file\"" << endl;
					exit(1);
				}
				beta = atof(token[2].c_str()); // proportion
				probXGivenInv[0] = atof(token[3].c_str()); // freq(A)
				probXGivenInv[1] = atof(token[4].c_str()); // freq(C)
				probXGivenInv[2] = atof(token[5].c_str()); // freq(G)
				probXGivenInv[3] = atof(token[6].c_str()); // freq(T)
				catLoaded[num_alpha] = 1;
			}
		}
	}
	fin.close();
	// if there is no constant site
	if (!catLoaded[num_alpha]) {
		beta = 0; // proportion
		probXGivenInv[0] = 0; // freq(A)
		probXGivenInv[1] = 0; // freq(C)
		probXGivenInv[2] = 0; // freq(G)
		probXGivenInv[3] = 0; // freq(T)
	}
}

void ParameterSet::readParamFile(char* fileName, vector<string>& nodeList, vector<double>& edgeLens, 
		int catID, int num_edges, int edgeRepresent) {

	// read the values of w and pi from the parameter file
	// For reversible, there are six w values (i.e. num_w should be 6)

	// edgeRepresent: representation of the edge length
	//                1 - Average number of substitutions per site
	//                2 - Time

	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the file :" << fileName << endl;
		exit(1);
	}

	// this function is designed for num_w = 6
	int num_w = 6;

	// allocate the space to the arrays
	this->num_w = num_w;
	this->num_edges = num_edges;
	if (w == NULL)
		w = malloc_decimals(num_edges*num_w); //matrix(num_edges, num_w);
	if (pi == NULL)
		pi = malloc_decimals(num_edges*4); //matrix(num_edges, 4);
	if (t == NULL)
		t = malloc_decimals(num_edges); //matrix(num_edges, 1);

	// read the file
	vector<string> token;
	string aline;
	string nodeName;
	int nodeID;
	double base;
	double* curr_w;
	double* curr_pi;
	int* accessed = new int[num_edges];
	memset(accessed, 0, num_edges*sizeof(int));
	while(getline(fin,aline)) {
		if (aline.length() > 0 && aline[0]!='#') {
			tokenizer(aline, " \t, ", &token);
			if ((int)token.size() <= num_w+4) {
				cerr << "Error! The number of items are too few in the following line inside the file " << fileName << " :" << endl;
				cerr << aline << endl;
				exit(1);
			}
			// get the node name and the node ID
			nodeName = token[0];
			for (nodeID=0; nodeID<num_edges; nodeID++) {
				if (nodeName == nodeList[nodeID])
					break;
			}
			if (nodeID >= num_edges) {
				cerr << "Error! The node name : " << nodeName << " does not appear in the tree file" << endl;
				exit(1);
			}
			if (accessed[nodeID]) {
				cerr << "Error! The node name : " << nodeName << " appears twice in the file " << fileName << endl;
				exit(1);
			}
			accessed[nodeID] = 1;
			curr_w = &(w[nodeID*num_w]);
			curr_pi = &(pi[nodeID*4]);
			for (int i=0; i<num_w; i++) {
				curr_w[i] = atof(token[i+1].c_str());
			}
			for (int i=0; i<4; i++) {
				curr_pi[i] = atof(token[num_w+i+1].c_str());
			}
			if (edgeRepresent == 1) {
				// edge length represents the average number of substitutions per site
				base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
				t[nodeID] = edgeLens[catID*num_edges + nodeID] / base;
			} else {
				// edge length represents the time
				t[nodeID] = edgeLens[catID*num_edges + nodeID];
			}
			// normalize
			/*
			for (int i=0; i<num_w; i++) {
				curr_w[i] = curr_w[i] / base;
			}
			t[nodeID] = edgeLens[catID*num_edges + nodeID];
			 */
		}
	}
	fin.close();
	// check whether there is missing node name
	for (nodeID=0; nodeID<num_edges; nodeID++) {
		if (accessed[nodeID] == 0) {
			cerr << "Error! The node name : " << nodeList[nodeID] << " does not appear in the file " << fileName << endl;
			exit(1);
		}
	}
	delete[] accessed;
}

// read the values of w and pi from the parameter file list
void AllParameterSet::readParamFileList(char* fileName, vector<string>& categoryList, vector<string>& nodeList, 
		vector<double>& edgeLens, int num_edges, int edgeRepresent) {
	string aline;
	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error! Cannot open the file : " << fileName << endl;
		exit(1);
	}
	vector<string> token;
	int catID;
	int catNum = (int)categoryList.size();
	int* accessed = new int[catNum];
	memset(accessed, 0, catNum*sizeof(int));
	while (getline(fin, aline)) {
		if (aline.length()==0 || aline[0]=='#')
			continue;
		tokenizer(aline, "\t ", &token);
		if (token.size() >= 2) {
			for (catID=0; catID<catNum; catID++) {
				if (token[0] == categoryList[catID])
					break;
			}
			if (catID>=catNum) {
				cerr << "Error! The category name : " << token[0] << " does not appear in the tree file" << endl;
				exit(1);
			}
			if (accessed[catID]) {
				cerr << "Error! The category name : " << token[0] << " appears twice in the file: " << fileName << endl;
				exit(1);
			}
			accessed[catID] = 1;
			// read the values of w and pi from the parameter file list
			ps[catID]->readParamFile((char*) token[1].c_str(), nodeList, edgeLens, catID, num_edges, edgeRepresent);
		}
	}
	// check whether there is category which does not appear in the parameter file list
	for (catID=0; catID<catNum; catID++) {
		if (accessed[catID]==0) {
			cerr << "Error! The cateogry name : " << categoryList[catID] << " does not appear in the file: " << fileName << endl;
			exit(1);
		}
	}
	fin.close();
	delete[] accessed;
}
