/*
 *
 * core_mpi.cpp
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


#include "core_mpi.h"

// =====================
// FOR RAL
// =====================

void sendInfo(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int proceed,
		int num_edges, int num_w, int num_chars, int destID, int tagID, ParameterSet* ps, VariableSet* vs) {

	// proceed:
	// 0 - stop the process
	// 1 - information about the arrangement of rate matrices is about to receive

	int i,k;
	int ps_vs_valid = 0;
	int dim = num_edges+5;
	int* buffer = new int[dim];
	// buffer[0]: proceed
	// buffer[1...num_edges]: rate matrix
	// buffer[num_edges+1]: mode; 
	// buffer[num_edges+2]: maxIT;
	// buffer[num_edges+3]: jobID;
	// buffer[num_edges+4]: 1 - ps and vs are not NULL; 0 - otherwise
	
	buffer[0] = proceed;
	k=1;
	if (proceed) {
		for (i=0; i<num_edges; i++)
			buffer[k++] = rateMatrix->at(i);
		buffer[k++] = currMode;
		buffer[k++] = currMaxIT;
		buffer[k++] = jobID;

		if (ps!=NULL && vs!=NULL) {
			buffer[k++] = 1;
			ps_vs_valid = 1;
		} else {
			buffer[k++] = 0;
		}
	}
	
	MPI_Send(&(buffer[0]), dim, MPI_INT, destID, tagID, MPI_COMM_WORLD);

	delete[] buffer;
	
	if (ps_vs_valid)
		sendPsVs(*ps, *vs, destID, tagID);
}

void receiveInfo(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& proceed, 
		int num_edges, int num_w, int num_chars, int tagID, ParameterSet* ps, VariableSet* vs, int& isPsVsUpdated) {

	// the signal:
	// 0 - stop the process
	// 1 - information about the arrangement of rate matrices is about to receive
	// 2 - initial status of the parameters for tuning is about to receive
	
	MPI_Status status;
	int i,k;
	int dim = num_edges+5;

	int* buffer = new int[dim];
	// buffer[0]: proceed
	// buffer[1...num_edges]: rate matrix
	// buffer[num_edges+1]: mode; 
	// buffer[num_edges+2]: maxIT;
	// buffer[num_edges+3]: jobID;
	// buffer[num_edges+4]: 1 - ps and vs are not NULL; 0 - otherwise

	MPI_Recv(&(buffer[0]), dim, MPI_INT, 0, tagID, MPI_COMM_WORLD, &status);

	proceed = buffer[0];
	k = 1;

	if (proceed) {
		for (i=0; i<num_edges; i++)
			 rateMatrix.at(i) = buffer[k++];
		num_rateMat = getNumRateGrp(rateMatrix);
		currMode = buffer[k++];
		currMaxIT = buffer[k++];
		jobID = buffer[k++];
		isPsVsUpdated = buffer[k++];
		
		if (isPsVsUpdated) {
			receivePsVs(*ps, *vs, tagID);
		}
	}

	delete[] buffer;
}


/*
void sendInfo(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int num_rateMat, int signal,
		int num_edges, int destID, int tagID) {

	// the signal:
	// 0 - stop the process
	// 1 - information about the arrangement of rate matrices is about to receive
	// 2 - initial status of the parameters for tuning is about to receive

	int i;
	int dim = num_edges+4;
	int* buffer = new int[dim];
	// buffer[0...num_edges-1]: rate matrix; buffer[num_edges]: mode; 
	// buffer[num_edges+1]: maxIT; buffer[num_edges+2]: SIGNAL;
	// buffer[num_edges+3]: jobID;
	if (signal == 1) {
		for (i=0; i<num_edges; i++)
			buffer[i] = rateMatrix->at(i);
		buffer[num_edges] = currMode;
		buffer[num_edges+1] = currMaxIT;
		buffer[num_edges+3] = jobID;
	}
	buffer[num_edges+2] = signal;
	// list out the content of buffer
	MPI_Send(&(buffer[0]), dim, MPI_INT, destID, tagID, MPI_COMM_WORLD);

	delete[] buffer;
}
*/

/*
void receiveInfo(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& signal, 
		int num_edges, int tagID) {

	// the signal:
	// 0 - stop the process
	// 1 - information about the arrangement of rate matrices is about to receive
	// 2 - initial status of the parameters for tuning is about to receive
	
	MPI_Status status;
	int i;
	int dim = num_edges+4;
	int* buffer = new int[dim];
	// buffer[0...num_edges-1]: rate matrix; buffer[num_edges]: mode; 
	// buffer[num_edges+1]: maxIT; buffer[num_edges+2]: SIGNAL
	// buffer[num_edges+3]: jobID;

	MPI_Recv(&(buffer[0]), dim, MPI_INT, 0, tagID, MPI_COMM_WORLD, &status);
	// list out the content of buffer

	signal = buffer[num_edges+2];

	if (signal == 1) {
		currMode = buffer[num_edges];
		currMaxIT = buffer[num_edges+1];
		jobID = buffer[num_edges+3];
		for (i=0; i<num_edges; i++)
			rateMatrix[i] = buffer[i];
		num_rateMat = getNumRateGrp(rateMatrix);
	}

	delete[] buffer;
}
*/

void sendPsVs(ParameterSet& ps, VariableSet& vs, int destID, int tagID) {

	int i,k;

	// First "num_w*ps.num_edges" items in buffer : values of w in ps
	// Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	// Next "num_edges" items in buffer : values of t in ps

	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars" items in buffer : values of rootNodeFreq in vs
	
	int dim = ps.num_w*ps.num_edges + ps.num_chars*ps.num_edges + ps.num_edges + 1 + vs.num_alpha + vs.num_chars + vs.num_chars;
	double* buffer = new double[dim];
	for (i=0; i<ps.num_w*ps.num_edges; i++)
		buffer[i] = ps.w[i];
	k=ps.num_w*ps.num_edges;
	for (i=0; i<ps.num_chars*ps.num_edges; i++)
		buffer[i+k] = ps.pi[i];
	k+=ps.num_chars*ps.num_edges;
	for (i=0; i<ps.num_edges; i++)
		buffer[i+k] = ps.t[i];
	k+=ps.num_edges;
	
	buffer[k] = vs.beta;
	k++;
	for (i=0; i<vs.num_alpha; i++)
		buffer[i+k] = vs.alpha[i];
	k+=vs.num_alpha;
	for (i=0; i<vs.num_chars; i++)
		buffer[i+k] = vs.probXGivenInv[i];
	k+=vs.num_chars;
	for (i=0; i<vs.num_chars; i++)
		buffer[i+k] = vs.rootNodeFreq[i];
	k+=vs.num_chars;

	MPI_Send(&(buffer[0]), dim, MPI_DOUBLE, destID, tagID, MPI_COMM_WORLD);

	delete[] buffer;
}

void receivePsVs(ParameterSet& ps, VariableSet& vs, int tagID) {
	
	int dim = ps.num_w*ps.num_edges + ps.num_chars*ps.num_edges + ps.num_edges + 1 + vs.num_alpha + vs.num_chars + vs.num_chars;
	double* buffer = new double[dim];

	MPI_Status status;

	int senderID = 0;
	MPI_Recv(&(buffer[0]), dim, MPI_DOUBLE, senderID, tagID, MPI_COMM_WORLD, &status);
	
	int i,k;

	// First "num_w*ps.num_edges" items in buffer : values of w in ps
	// Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	// Next "num_edges" items in buffer : values of t in ps
	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars" items in buffer : values of rootNodeFreq in vs

	double x;
	for (i=0; i<ps.num_w*ps.num_edges; i++)
		ps.w[i] = buffer[i];
	k=ps.num_w*ps.num_edges;
	for (i=0; i<ps.num_chars*ps.num_edges; i++)
		ps.pi[i] = buffer[i+k];
	k+=ps.num_chars*ps.num_edges;
	for (i=0; i<ps.num_edges; i++)
		ps.t[i] = buffer[i+k];
	k+=ps.num_edges;
	vs.beta = buffer[k];
	k++;
	x=0.0;
	for (i=0; i<vs.num_alpha; i++) {
		vs.alpha[i] = buffer[i+k];
		x += vs.alpha[i];
	}
	// vs.beta = 1.0 - x;
	k+=vs.num_alpha;
	for (i=0; i<vs.num_chars; i++)
		vs.probXGivenInv[i] = buffer[i+k];
	k+=vs.num_chars;
	for (i=0; i<vs.num_chars; i++)
		vs.rootNodeFreq[i] = buffer[i+k];
	k+=vs.num_chars;
	
	delete[] buffer;
}

void sendResult(int jobID, double logL, double ic, ParameterSet& ps, VariableSet& vs, int numIter, int tagID) {

	int i,k;

	// First "num_w*ps.num_edges" items in buffer : values of w in ps
	// Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	// Next "num_edges" items in buffer : values of t in ps
	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars" items in buffer : values of rootNodeFreq in vs
	// Next three items in buffer: logL, ic score and number of iterations
	// Next item in buffer: jobID
	int dim = ps.num_w*ps.num_edges + ps.num_chars*ps.num_edges + ps.num_edges + 1 + vs.num_alpha + vs.num_chars + vs.num_chars + 3 + 1;
	double* buffer = new double[dim];
	for (i=0; i<ps.num_w*ps.num_edges; i++)
		buffer[i] = ps.w[i];
	k=ps.num_w*ps.num_edges;
	for (i=0; i<ps.num_chars*ps.num_edges; i++)
		buffer[i+k] = ps.pi[i];
	k+=ps.num_chars*ps.num_edges;
	for (i=0; i<ps.num_edges; i++)
		buffer[i+k] = ps.t[i];
	k+=ps.num_edges;
	buffer[k] = vs.beta;
	k++;
	for (i=0; i<vs.num_alpha; i++)
		buffer[i+k] = vs.alpha[i];
	k+=vs.num_alpha;
	for (i=0; i<vs.num_chars; i++)
		buffer[i+k] = vs.probXGivenInv[i];
	k+=vs.num_chars;
	for (i=0; i<vs.num_chars; i++)
		buffer[i+k] = vs.rootNodeFreq[i];
	k+=vs.num_chars;
	buffer[k++] = logL;
	buffer[k++] = ic;
	buffer[k++] = numIter;
	buffer[k++] = jobID;

	/*
	// list the content of buffer
	cout << "[sendResult] ";
	for (i=0; i<dim; i++)
		cout << " " << buffer[i];
	cout << endl;
	*/

	MPI_Send(&(buffer[0]), dim, MPI_DOUBLE, 0, tagID, MPI_COMM_WORLD);

	delete[] buffer;
}

void receiveResult(int& jobID, double& logL, double& ic, ParameterSet& ps, VariableSet& vs,
					int& numIter, int tagID, int& sender) {
	
	int dim = ps.num_w*ps.num_edges + ps.num_chars*ps.num_edges + ps.num_edges + 1 + vs.num_alpha + vs.num_chars + vs.num_chars + 3 + 1;
	double* buffer = new double[dim];
	MPI_Status status;

	MPI_Recv(&(buffer[0]), dim, MPI_DOUBLE, MPI_ANY_SOURCE, tagID, MPI_COMM_WORLD, &status);
	
	int i,k;

	/*
	// list the content of buffer
	cout << "[receiveResult] ";
	for (i=0; i<dim; i++)
		cout << " " << buffer[i];
	cout << endl;
	*/
	
	// First "num_w*ps.num_edges" items in buffer : values of w in ps
	// Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	// Next "num_edges" items in buffer : values of t in ps
	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars" items in buffer : values of rootNodeFreq in vs
	// Next three items in buffer: logL, ic score and number of iterations
	// Next item in buffer: jobID

	double x;
	for (i=0; i<ps.num_w*ps.num_edges; i++)
		ps.w[i] = buffer[i];
	k=ps.num_w*ps.num_edges;
	for (i=0; i<ps.num_chars*ps.num_edges; i++)
		ps.pi[i] = buffer[i+k];
	k+=ps.num_chars*ps.num_edges;
	for (i=0; i<ps.num_edges; i++)
		ps.t[i] = buffer[i+k];
	k+=ps.num_edges;
	vs.beta = buffer[k];
	k++;
	x=0.0;
	for (i=0; i<vs.num_alpha; i++) {
		vs.alpha[i] = buffer[i+k];
		x += vs.alpha[i];
	}
	vs.beta = 1.0 - x;
	k+=vs.num_alpha;
	for (i=0; i<vs.num_chars; i++)
		vs.probXGivenInv[i] = buffer[i+k];
	k+=vs.num_chars;
	for (i=0; i<vs.num_chars; i++)
		vs.rootNodeFreq[i] = buffer[i+k];
	k+=vs.num_chars;
	logL = buffer[k++];
	ic = buffer[k++];
	numIter = buffer[k++];
	jobID = buffer[k++];
	sender = status.MPI_SOURCE;
	
	delete[] buffer;
}

// =====================
// For HAS
// =====================

void sendInfoHAS(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int num_rateMat, int proceed,
		int numRateCat, int modelType, int num_edges, int destID, int tagID) {

	// MPI_Status status;
	int i;
	int dim = num_edges+6;
	int* buffer = new int[dim];
	// buffer[0...num_edges-1]: rate matrix;
	// buffer[num_edges]: mode; 
	// buffer[num_edges+1]: maxIT; 
	// buffer[num_edges+2]: numRateCat;
	// buffer[num_edges+3]: modelType;
	// buffer[num_edges+4]: jobID;
	// buffer[num_edges+5]: proceed;
	if (proceed == 1) {
		for (i=0; i<num_edges; i++)
			buffer[i] = rateMatrix->at(i);
		buffer[num_edges] = currMode;
		buffer[num_edges+1] = currMaxIT;
		buffer[num_edges+2] = numRateCat;
		buffer[num_edges+3] = modelType;
		buffer[num_edges+4] = jobID;
	}
	buffer[num_edges+5] = proceed;

	/*
	// list out the content of buffer
	cout << "[sendInfo(RAS)] tagID = " << tagID << " ";
	for (i=0; i<dim; i++)
		cout << buffer[i] << " ";
	cout << endl;
	*/

	MPI_Send(&(buffer[0]), dim, MPI_INT, destID, tagID, MPI_COMM_WORLD);

	delete[] buffer;
}

void receiveInfoHAS(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& proceed,
		int& numRateCat, int& modelType, int num_edges, int tagID) {

	MPI_Status status;
	int i;
	int dim = num_edges+6;
	int* buffer = new int[dim];
	// buffer[0...num_edges-1]: rate matrix; 
	// buffer[num_edges]: mode; 
	// buffer[num_edges+1]: maxIT; 
	// buffer[num_edges+2]: numRateCat;
	// buffer[num_edges+3]: modelType;
	// buffer[num_edges+4]: jobID;
	// buffer[num_edges+5]: proceed;
	
	// cout << "[receiveInfo(RAS)a] tagID = " << tagID << " waiting..." << endl << flush;

	MPI_Recv(&(buffer[0]), dim, MPI_INT, 0, tagID, MPI_COMM_WORLD, &status);

	/*
	// list out the content of buffer
	cout << "[receiveInfo(RAS)b] ";
	for (i=0; i<dim; i++)
		cout << buffer[i] << " ";
	cout << endl;
	*/

	proceed = buffer[num_edges+5];

	if (proceed == 1) {
		currMode   = buffer[num_edges];
		currMaxIT  = buffer[num_edges+1];
		numRateCat = buffer[num_edges+2];
		modelType  = buffer[num_edges+3];
		jobID      = buffer[num_edges+4];
		for (i=0; i<num_edges; i++)
			rateMatrix[i] = buffer[i];
		num_rateMat = getNumRateGrp(rateMatrix);
	}

	delete[] buffer;

}

void sendResultHAS(int jobID, double logL, double ic, int df, AllParameterSet& all_ps, VariableSet& vs, int numIter, int tagID) {
/*
	// show the content of all_ps
	cout << "content of all_ps" << endl;
	cout << "----------------------" << endl;
	all_ps.showContent();
	cout << "----------------------" << endl;
	
	// show the content of vs
	cout << "content of vs" << endl;
	cout << "----------------------" << endl;
	vs.showContent();
	cout << "----------------------" << endl;
*/    
	int i,j,k;

	// for each ps (the number of ps = all_ps.numRateCat) :
	//  First "num_w*ps.num_edges" items in buffer : values of w in ps
	//  Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	//  Next "num_edges" items in buffer : values of t in ps
	
	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars*num_alpha" items in buffer : values of rootNodeFreq in vs
	// Next three items in buffer: logL, ic score and number of iterations
	// Next two items in buffer: jobID and df
	ParameterSet* ps = all_ps.ps[0];
	int dim = (ps->num_w*ps->num_edges + ps->num_chars*ps->num_edges + ps->num_edges)*all_ps.numRateCat + 1 + vs.num_alpha + vs.num_chars + vs.num_chars*vs.num_alpha + 3 + 2;
	
	// first send {dim, num_alpha, num_w, num_edges, num_chars}
	int* bufferInt = new int[5];
	bufferInt[0] = dim; bufferInt[1] = vs.num_alpha; bufferInt[2] = ps->num_w;
	bufferInt[3] = ps->num_edges; bufferInt[4] = ps->num_chars;

	/*
	// list the content of buffer
	cout << "[sendResult(HAS)a] ";
	for (i=0; i<5; i++)
		cout << " " << bufferInt[i];
	cout << endl;
	*/

	MPI_Send(&(bufferInt[0]), 5, MPI_INT, 0, tagID, MPI_COMM_WORLD);
	
	double* buffer = new double[dim];
	k = 0;
	for (j=0; j<all_ps.numRateCat; j++) {
		ps = all_ps.ps[j];
		for (i=0; i<ps->num_w*ps->num_edges; i++)
			buffer[i+k] = ps->w[i];
		k+=ps->num_w*ps->num_edges;
		for (i=0; i<ps->num_chars*ps->num_edges; i++)
			buffer[i+k] = ps->pi[i];
		k+=ps->num_chars*ps->num_edges;
		for (i=0; i<ps->num_edges; i++)
			buffer[i+k] = ps->t[i];
		k+=ps->num_edges;
	}
	buffer[k] = vs.beta;
	k++;
	for (i=0; i<vs.num_alpha; i++)
		buffer[i+k] = vs.alpha[i];
	k+=vs.num_alpha;
	for (i=0; i<vs.num_chars; i++)
		buffer[i+k] = vs.probXGivenInv[i];
	k+=vs.num_chars;
	for (i=0; i<vs.num_chars*vs.num_alpha; i++)
		buffer[i+k] = vs.rootNodeFreq[i];
	k+=vs.num_chars*vs.num_alpha;
	buffer[k++] = logL;
	buffer[k++] = ic;
	buffer[k++] = numIter;
	buffer[k++] = jobID;
	buffer[k++] = df;

	/*
	// list the content of buffer
	cout << "[sendResult(HAS)b] ";
	for (i=0; i<dim; i++)
		cout << " " << buffer[i];
	cout << endl;
	*/

	// send the values of ps and vs
	MPI_Send(&(buffer[0]), dim, MPI_DOUBLE, 0, tagID, MPI_COMM_WORLD);

	delete[] buffer;
}

void receiveResultHAS(int& jobID, double& logL, double& ic, int& df, AllParameterSet& all_ps, VariableSet& vs,
					int& numIter, int tagID, int& sender) {
	
	MPI_Status status;
	int dim, num_alpha, num_w, num_edges, num_chars;
	
	// first receive {dim, num_alpha, num_w, num_edges, num_chars}
	int* bufferInt = new int[5];
	MPI_Recv(&(bufferInt[0]), 5, MPI_INT, MPI_ANY_SOURCE, tagID, MPI_COMM_WORLD, &status);

	/*
	// list the content of buffer
	cout << "[receiveResult(HAS)a] ";
	for (int i=0; i<5; i++)
		cout << " " << bufferInt[i];
	cout << endl;
	*/

	sender = status.MPI_SOURCE;
	dim = bufferInt[0]; num_alpha = bufferInt[1]; num_w = bufferInt[2];
	num_edges = bufferInt[3]; num_chars = bufferInt[4];
	delete[] bufferInt;

	if (dim <= 0)
		return;
		
	double* buffer = new double[dim];

	// get the results
	MPI_Recv(&(buffer[0]), dim, MPI_DOUBLE, sender, tagID, MPI_COMM_WORLD, &status);
	
	int i,j, k;

	/*
	// list the content of buffer
	cout << "[receiveResult(HAS)b] ";
	for (i=0; i<dim; i++)
		cout << " " << buffer[i];
	cout << endl;
	*/

	// for each ps (the number of ps = all_ps.numRateCat) :
	//  First "num_w*ps.num_edges" items in buffer : values of w in ps
	//  Next "num_chars*ps.num_edges" items in buffer : values of pi in ps
	//  Next "num_edges" items in buffer : values of t in ps
	
	// Next item in buffer : value of beta in vs
	// Next "num_alpha" items in buffer : values of alpha in vs
	// Next "num_chars" items in buffer : values of probXGivenInv in vs
	// Next "num_chars*num_alpha" items in buffer : values of rootNodeFreq in vs
	// Next three items in buffer: logL, ic score and number of iterations
	// Next two items in buffer: jobID and df

	all_ps.numRateCat = num_alpha;
	
	ParameterSet* ps;

	k = 0;
	for (j=0; j<num_alpha; j++) {
		
		ps = all_ps.ps[j];
		ps->num_chars = num_chars;
		ps->num_edges = num_edges;
		ps->num_w = num_w;
		for (i=0; i<num_w*num_edges; i++)
			ps->w[i] = buffer[i+k];
		k+=num_w*num_edges;
		for (i=0; i<num_chars*num_edges; i++)
			ps->pi[i] = buffer[i+k];
		k+=num_chars*num_edges;
		for (i=0; i<num_edges; i++)
			ps->t[i] = buffer[i+k];
		k+=num_edges;
	}
	vs.num_alpha = num_alpha;
	vs.num_chars = num_chars;
	// vs.beta = buffer[k];
	k++;
	double x=0.0;
	for (i=0; i<num_alpha; i++) {
		vs.alpha[i] = buffer[i+k];
		x += vs.alpha[i];
	}
	vs.beta = 1.0 - x;
	k+=num_alpha;
	for (i=0; i<num_chars; i++)
		vs.probXGivenInv[i] = buffer[i+k];
	k+=num_chars;
	for (i=0; i<num_chars*num_alpha; i++)
		vs.rootNodeFreq[i] = buffer[i+k];
	k+=num_chars*num_alpha;
	logL = buffer[k++];
	ic = buffer[k++];
	numIter = buffer[k++];
	jobID = round(buffer[k++]);
	df = round(buffer[k++]);
	
	// cout << "leaving receiveResult ..." << endl << flush;
	delete[] buffer;
}

// Terminate all processes
void terminateAllProcs(int num_edges, int num_w, int num_chars, int tagID) {
	// Get the number of processes
	cout << "terminating all processes...." << endl;
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int i;
	for (i=1; i<world_size; i++)
		sendInfo(NULL, 0, 0, 0, 0, num_edges, num_w, num_chars, i, tagID, NULL, NULL);
}

// Terminate all processes (for HAS)
void terminateAllProcsHAS(int num_edges, int tagID) {
	// Get the number of processes
	cout << "terminating all processes...." << endl;
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int i;
	for (i=1; i<world_size; i++)
		sendInfoHAS(NULL , 0, 0, 0, 0, 0, 0, 0, num_edges, i, tagID);
}
