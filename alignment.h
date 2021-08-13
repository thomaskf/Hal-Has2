/*
 *
 * alignment.h
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


#ifndef __RAL_RAS__alignment__
#define __RAL_RAS__alignment__

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include "tool_box.h"
#include "charSet.h"
#include "definitions.h"
#include "simd.h"

using namespace std;

#define MAX_NAME_LEN 64

// #define USE_64_BIT_BINARY

class Alignment {
public:
	char* names;
	char* sites;
	// the i-th sites of all sequences : sites[numSeqs*i ~ numSeqs*(i+1)-1]
	
	// for speeding up the procedure of computing p-values
#ifndef USE_64_BIT_BINARY
	unsigned int* sitesBinary; // dimension: siteB_len * numUniqueSites * 4
	unsigned int* returnSiteB(int siteB_ID, int nucl); // return the binary array of the site
	unsigned int* tmpArray; // for temporary use, length = siteB_len
#else
	unsigned long long* sitesBinary; // dimension: siteB_len * numUniqueSites * 4
	unsigned long long* returnSiteB(int siteB_ID, int nucl); // return the binary array of the site
	unsigned long long* tmpArray; // for temporary use, length = siteB_len
#endif
	int siteB_len;
	int* pdistOrder; // the new order according to the p-distance between the sites
	
	int numSeqs;
	int numSites;
	// int numConstSites;

	char* uniqueSites;
	double*  numEachUniqueSites;     // number of appearance of each unique site
	char* isConstantSites;   // indicate whether the site is a constant site
	int numUniqueSites; // number of different unique sites
	int numConstSites; // number of constant sites

	bool isGTRUpilson; // true if it is GTR Upilson model
	
	// for testing
	// int* prePosWithSmallDiff;
	char* changedNodes; // for computation of log-likelihood

	Alignment(); // constructor
	~Alignment(); // destructor

	void setGTRUpilson(); // set "isGTRUpilson" to TRUE

	// read the sequences file in sequential phylip format
	// the format of the first line: "[number of sequences]    [number of sites]"
	// if one wants to specify the order of the sequences, use name2index
	void readFile(char* fileName, map<string,int>* name2index);

	// read the sequences file in FASTA format
	// the format of the first line: "[number of sequences]    [number of sites]"
	// if one wants to specify the order of the sequences, use name2index
	void readFASTA(char* fileName, map<string,int>* name2index);

	// only keep the sites with predefined set of characters
	void keepSitesWithChars(int includeGap);
	// includeGap: 1 - allow gaps inside the columns; 0 - remove the columns with gap

	// sort the sites (by radix sort),
	// return the sorted site ID
	int* getSortedSitesID();

	void computeNumUniqueSites(int* sortedSiteID);

	// only keep the unique sites
	// input: sortedSiteID returned by getSortedSitesID()
	// output: 1. sites (with only the unique site patterns)
	//         2. numEachPattern (i.e. number of each different pattern)
	//         3. numPattern (i.e. total number of different patterns)
	//         4. rateCats (i.e. the rate categories for each of the unqiue sites)
	void keepUniqueSites(int* sortedSiteID);

	// print out the content according to the site order
	void print(int* siteOrder);

	// print out the statistics of the unique sites
	void printUniqueStat();
	
	// for building the changedNodes array
	void buildChangedNodes(int* topMatrix, int numLineTopMat);
	
	// ==================
	// for testing only
	// ==================
	// across the adjacent sites, compute the total number of differences
	void getNumDiffBWAdjSites();
	// void buildPrePosWithSmallDiff();
	// reorder the unique sites by using greedy approach
	void reorderUniqueSites();
	// build the binary array for computation of p-dist
	void buildSitesArray();
	// get the pdist between two char arrays (newer version)
	int p_dist(int siteID1, int siteID2);
	// original version
	// get the pdist between two char arrays
	int p_dist_orig(int siteID1, int siteID2);

};

// convert the phy format into maze format
void phy2maze(char* file);

// read the alignments in maze format
// remove the columns with gaps
// then output the alignment in maze format
void maze_remove_gaps(char* file);

// convert the maze format into phy format
void maze2phy(char* file);

// convert the phy format into FASTA format
void phy2fasta(char* file);

// check whether all the characters in the character string are the same
int isConstantSite(char* aSite, int len);

// popCount for 32 bits
int popCount(unsigned int bitVector);

// popCount for 64 bits
int popCount64(unsigned long long bitVector);

#endif /* defined(__RAL_RAS__alignment__) */
