/*
 *
 * alignment.cpp
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

#include "alignment.h"

Alignment::Alignment() {
	// constructor
	names = NULL;
	sites = NULL;
	numSeqs=0;
	numSites=0;
	uniqueSites = NULL;
	numEachUniqueSites = NULL;
	isConstantSites = NULL;
	numUniqueSites=0;
	numConstSites=0;
	isGTRUpilson=0;
	// prePosWithSmallDiff = NULL;
	changedNodes = NULL;
	sitesBinary = NULL;
	pdistOrder = NULL;
	tmpArray = NULL;
}

Alignment::~Alignment() {
	// destructor
	if (names != NULL)
		free(names);
	if (sites != NULL)
		free(sites);
	if (uniqueSites!=NULL)
		free(uniqueSites);
	if (numEachUniqueSites != NULL)
		free_decimals(numEachUniqueSites);
	if (isConstantSites != NULL)
		free(isConstantSites);
	// if (prePosWithSmallDiff != NULL)
	//    free(prePosWithSmallDiff);
	if (changedNodes != NULL)
		free(changedNodes);
	if (sitesBinary != NULL)
		free(sitesBinary);
	if (pdistOrder != NULL)
		free(pdistOrder);
	if (tmpArray != NULL)
		free(tmpArray);
}

// set "isGTRUpilson" to TRUE
void Alignment::setGTRUpilson() {
	isGTRUpilson = 1;
}


// read the sequences file
// the format of the first line: "[number of sequences]    [number of sites]"
// if one wants to specify the order of the sequences, use name2index
void Alignment::readFile(char* fileName, map<string,int>* name2index) {

	ifstream fin;
	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the file :" << fileName << endl;
		exit(1);
	}

	string aline;
	vector<string> token;

	// get the number of sequences and the number of sites in the first line
	if (getline(fin,aline) && aline.length()>0) {
		trim(aline); // remove the space or the nonprintable characters at the end
		tokenizer(aline, " \t,", &token);
		if (token.size() == 2) {
			numSeqs = atoi(token[0].c_str());
			numSites = atoi(token[1].c_str());
		} else {
			cerr << "Error! The number of fields in the first line of the sequence file is not 2" << endl;
			exit(1);
		}
	} else {
		cerr << "Error! Either the sequence file is empty or the first line of the sequence file is empty" << endl;
		exit(1);
	}

	// cout << "numSeqs " << numSeqs << endl << flush;
	// cout << "numSites " << numSites << endl << flush;

	// allocate memory to the array sites
	if (sites != NULL)
		free(sites);
	sites = (char*) malloc (sizeof(char) * numSeqs * numSites);
	memset(sites,0,sizeof(char)*numSeqs*numSites);
	if (names != NULL)
		free(names);
	names = (char*) malloc (sizeof(char) * numSeqs * MAX_NAME_LEN);
	memset(names,0,sizeof(char)*numSeqs*MAX_NAME_LEN);

	// get the mapping between the characters and the integers
	char* char2Int = charToIntArray();

	// get the names and the sequences
	int curr_reads = 0;
	map<string,int>::iterator itr;
	while (getline(fin,aline)) {
		trim(aline); // remove the space or the nonprintable characters at the end
		if (aline.length() > 0 && aline[0]!=' ' ) {
			// cannot start with a space

			// read until the space
			int i=0; int j;
			int read_id = curr_reads;
			while (i < (int)aline.length() && aline[i]!=' ' && i<10)
				i++;
			if (i > 0) {
				// name: aline[0...i-1]
				int name_len=i;
				if (name_len > MAX_NAME_LEN) {
					name_len = MAX_NAME_LEN;
					cerr << "Warning! The name has been trimmed" << endl;
				}
				if (name2index != NULL) {
					itr = name2index->find(aline.substr(0,name_len));
					if (itr != name2index->end()) {
						read_id = itr->second;
					} else {
						cerr << "Error! The name : " << aline.substr(0,name_len) << " cannot be found in the topology file" << endl;
						exit(1);
					}
				}
				for (j=0; j<name_len; j++) {
					names[read_id*MAX_NAME_LEN+j] = aline[j];
				}
			} else {
				cerr << "Error! No name is in a line." << endl;
				exit(1);
			}

			// skip all the space
			while (i < (int)aline.length() && aline[i]==' ')
				i++;

			if (i < (int)aline.length()) {
				// seq : aline[i... end]
				int seq_len = (int)aline.length() - i;
				if (seq_len != numSites) {
					cerr << "Error! The " << curr_reads+1 << "-th sequence's length does not match with the number of sites specified in the beginning of the sequence file" << endl;
					exit(1);
				}
				for (j=0; j<seq_len; j++) {
					sites[j*numSeqs+read_id] = (char) char2Int[(int) toupper(aline[i+j])];
				}
			} else {
				cerr << "Error! No sequence is in a line." << endl;
				exit(1);
			}
			curr_reads++;
		}
	}
	// check whether the number of sequences matches
	if (curr_reads != numSeqs) {
		cerr << "Error! The number of sequences does not match with the number specified in the beginning of the sequence file" << endl;
		exit(1);
	}
	
	free(char2Int);
	
	fin.close();
}


// read the sequences file in FASTA format
// the format of the first line: "[number of sequences]    [number of sites]"
// if one wants to specify the order of the sequences, use name2index
void Alignment::readFASTA(char* fileName, map<string,int>* name2index) {

	ifstream fin;
	int i,j;
	string aline;
	vector<string> seqNames; // store all the names of the sequences
	vector<string> seqs; // store all the sequences
	int readID = -1;

	fin.open(fileName);
	if (!fin.is_open()) {
		cerr << "Error opening the file :" << fileName << endl;
		exit(1);
	}

	while (getline(fin,aline)) {
		trim(aline);
		if (aline.length()>0) {
			if (aline[0] == '>') {
				// name (up to the space)
				size_t space_pos = aline.find_first_of(" ");
				if (space_pos!=string::npos && space_pos>1)
					seqNames.push_back(aline.substr(1, space_pos-1));
				else
					seqNames.push_back(aline.substr(1));
				readID++;
				seqs.push_back("");
			} else if (readID >= 0) {
				// sequence
				seqs[readID].append(aline);
			}
		}
	}
	fin.close();

	if (seqs.size() == 0) {
		cerr << "Error! No sequence is loaded from the file " << fileName << endl;
		exit(1);
	}

	// check whether all the sequences are in the same length
	numSites = (int) seqs[0].length();
	numSeqs = (int) seqs.size();
	for (i=1; i<numSeqs; i++) {
		if (numSites != (int) seqs[i].length()) {
			cerr << "Error! The length of " << i+1 << "-th sequence is not the same as those of the prevoius sequences" << endl;
			exit(1);
		}
	}

	// check whether some sequence names are too long
	for (i=0; i<numSeqs; i++) {
		if (seqNames[i].length() > MAX_NAME_LEN) {
			cerr << "Warning! The name of " << i+1 << "-th sequence has been trimmed" << endl;
			seqNames[i] = seqNames[i].substr(0,MAX_NAME_LEN);
		}
	}

	// allocate memory to the array sites
	if (sites != NULL)
		free(sites);
	sites = (char*) malloc (sizeof(char) * numSeqs * numSites);
	memset(sites,0,sizeof(char)*numSeqs*numSites);
	if (names != NULL)
		free(names);
	names = (char*) malloc (sizeof(char) * numSeqs * MAX_NAME_LEN);
	memset(names,0,sizeof(char)*numSeqs*MAX_NAME_LEN);

	// get the mapping between the characters and the integers
	char* char2Int = charToIntArray();

	// save the names and the sequences into the arrays "names" and "sites"
	map<string,int>::iterator itr;

	for (i=0; i<numSeqs; i++) {
		int read_id = i;
		if (name2index != NULL) {
			itr = name2index->find(seqNames[i]);
			if (itr != name2index->end()) {
				read_id = itr->second;
			} else {
				cerr << "Error! The name : " << seqNames[i] << " cannot be found in the topology file" << endl;
			}
		}
		for (j=0; j<(int)seqNames[i].length(); j++) {
			names[read_id*MAX_NAME_LEN+j] = seqNames[i].at(j);
		}
		for (j=0; j<(int)seqs[i].length(); j++) {
			sites[j*numSeqs+read_id] = (char)char2Int[(int)toupper(seqs[i].at(j))];
		}
	}

	// release the memory allocated for "seqNames" and "seqs"
	seqNames.clear();
	seqs.clear();
}

// only keep the sites with predefined set of characters
void Alignment::keepSitesWithChars(int includeGap) {
	// includeGap: 1 - allow gaps and ambiguous chars inside the columns;
	//             0 - remove the columns with gap or ambiguous chars
	// in either case, all the other unrecognized characters will be removed
	int newNumSites = 0;
	int i,j;
	for (i=0; i<numSites; i++) {
		// check whether the current site contains any invalid character
		char* curr_site = &(sites[i*numSeqs]);
		bool isValid = true;
		for (j=0; j<numSeqs; j++) {
			if (curr_site[j]==0) {
				isValid = false;
				break;
			} else if (!includeGap && curr_site[j]!=A_ID && curr_site[j]!=C_ID && curr_site[j]!=G_ID && curr_site[j]!=T_ID) {
				isValid = false;
				break;
			}
		}
		// remove those sites with all gaps
		bool allGaps = true;
		if (isValid) {
			for (j=0; j<numSeqs; j++) {
				if (curr_site[j]!=GAP_ID) {
					allGaps = false;
					break;
				}
			}
			if (allGaps)
				isValid = false;
		}
		if (isValid) {
			if (i>newNumSites) {
				char* to = &(sites[newNumSites*numSeqs]);
				memcpy(to,curr_site,numSeqs*sizeof(char));
			}
			newNumSites++;
		}
	}
	numSites = newNumSites;
}

// sort the sites (by radix sort),
// return the sorted site ID
int* Alignment::getSortedSitesID() {

	int* buffer[2];

	//initialization
	buffer[0] = new int[numSites]; // (int*) malloc (numSites * sizeof(int));
	buffer[1] = new int[numSites]; // (int*) malloc (numSites * sizeof(int));
	for (int i=0; i<numSites; i++)
		(buffer[0])[i] = i;
	int numEachChar[MAXCHAR+1];
	int startPosEachChar[MAXCHAR+1];

	int currBufferID = 0;
	int* currBuffer = buffer[0];
	int* preBuffer = buffer[1];

	for (int i=numSeqs-1; i>=0; i--) {

		// update the pointers of the buffers
		preBuffer = currBuffer;
		currBufferID = (currBufferID+1)%2;
		currBuffer = buffer[currBufferID];

		// compute the starting positions for each character
		for (int k=0; k<=MAXCHAR; k++)
			numEachChar[k] = 0;
		for (int k=0; k<numSites; k++)
			numEachChar[(int)sites[k*numSeqs+i]]++;
		startPosEachChar[0]=0;
		for (int k=1; k<=MAXCHAR; k++) {
			startPosEachChar[k]=startPosEachChar[k-1]+numEachChar[k-1];
		}

		// initialize the array of storing the number of characters coming across
		int numEachCharSoFar[MAXCHAR+1];
		for (int k=0; k<=MAXCHAR; k++) {
			numEachCharSoFar[k]=0;
		}

		// perform classification according to the bases
		for (int k=0; k<numSites; k++) {
			int site = preBuffer[k];
			int base = sites[site*numSeqs+i];
			int pos = startPosEachChar[base]+numEachCharSoFar[base];
			currBuffer[pos] = site;
			numEachCharSoFar[base]++;
		}

	}
	delete[] preBuffer;

	return currBuffer;
}

void Alignment::computeNumUniqueSites(int* sortedSiteID) {
	numUniqueSites = 1;
	if (sortedSiteID == NULL) {
		cerr << "Error! The array ""sortedSiteID"" is NULL" << endl;
	}
	if (numSites==1)
		return;
	char* preSite;
	char* currSite = &(sites[sortedSiteID[0]*numSeqs]);
	for (int i=1; i<numSites; i++) {
		preSite = currSite;
		currSite = &(sites[sortedSiteID[i]*numSeqs]);
		bool isSame = true;
		for (int j=0; j<numSeqs; j++) {
			if (preSite[j]!=currSite[j]) {
				isSame = false;
				break;
			}
		}
		if (!isSame) {
			numUniqueSites++;
		}
	}
}

// check whether all the characters in the character string are the same
/*
int isConstantSite(char* aSite, int len) {
	if (len <= 1)
		return 1;
	char firstBase = aSite[0];
	for (int i=1; i<len; i++) {
		if (firstBase != aSite[i])
			return 0;
	}
	return 1;
}
*/
int isConstantSite(char* aSite, int len) {
	if (len <= 0) {
		return 0;
	}
	int sameChar = 0xf;
	for (int i=0; i<len && sameChar>0; i++) {
		sameChar &= aSite[i];
	}
	return sameChar;
}

// only keep the unique sites
// input: sortedSiteID returned by getSortedSitesID()
// output: 1. sites (with only the unique site patterns)
//         2. numEachPattern (i.e. number of each different pattern)
//         3. numPattern (i.e. total number of different patterns)
void Alignment::keepUniqueSites(int* sortedSiteID) {

	//get the number of unique sites
	computeNumUniqueSites(sortedSiteID);

	// allocate the memory to the array
	if (uniqueSites!=NULL)
		free(uniqueSites);
	uniqueSites = (char*) malloc(numUniqueSites*numSeqs*sizeof(char));
	if (numEachUniqueSites!=NULL)
		free_decimals(numEachUniqueSites);
	numEachUniqueSites = malloc_decimals(numUniqueSites);
	if (isConstantSites!=NULL)
		free(isConstantSites);
	isConstantSites = (char*) malloc(numUniqueSites*sizeof(char));

	// initialize the variables
	numConstSites = 0;
	
	// save the first site to the array
	char* currSite = &(sites[sortedSiteID[0]*numSeqs]);
	char* newSite = &(uniqueSites[0]);
	memcpy(newSite,currSite,numSeqs*sizeof(char));
	numEachUniqueSites[0]=1.0;
	if (isGTRUpilson)
		isConstantSites[0] = 0;
	else
		isConstantSites[0] = isConstantSite(newSite, numSeqs);
	int numUniqueSitesReported = 1;
	if (isConstantSites[0] > 0)
		numConstSites++;

	char* preSite;
	for (int i=1; i<numSites; i++) {
		preSite = currSite;
		currSite = &(sites[sortedSiteID[i]*numSeqs]);
		bool isSame = true;
		for (int j=0; j<numSeqs; j++) {
			if (preSite[j]!=currSite[j]) {
				isSame = false;
				break;
			}
		}
		if (!isSame) {
			newSite = &(uniqueSites[numUniqueSitesReported*numSeqs]);
			memcpy(newSite,currSite,numSeqs*sizeof(char));
			numEachUniqueSites[numUniqueSitesReported]=1.0;
			if (isGTRUpilson)
				isConstantSites[numUniqueSitesReported] = 0;
			else
				isConstantSites[numUniqueSitesReported] = isConstantSite(newSite, numSeqs);
			numUniqueSitesReported++;
		} else {
			numEachUniqueSites[numUniqueSitesReported-1]+=1.0;
		}
		if (isConstantSites[numUniqueSitesReported-1] > 0)
			numConstSites++;
	}
}


// print out the content according to the site order
void Alignment::print(int* siteOrder) {
	cout << numSeqs << "  " << numSites << endl;
	char* int2Char = intToCharArray();
	for (int i=0; i<numSeqs; i++) {
		int j=0;
		while (j<MAX_NAME_LEN && names[i*MAX_NAME_LEN+j]!=0) {
			cout << (char) names[i*MAX_NAME_LEN+j];
			j++;
		}
		cout << "      ";
		j=0;
		while (j<numSites) {
			if (siteOrder==NULL)
				cout << (char) int2Char[(int)sites[j*numSeqs+i]];
			else
				cout << (char) int2Char[(int)sites[siteOrder[j]*numSeqs+i]];
			j++;
		}
		cout << endl;
	}
}

// print out the statistics of the unique sites
void Alignment::printUniqueStat() {
	for (int i=0; i<numUniqueSites; i++) {
		for (int j=0; j<numSeqs; j++) {
			cout << (int) uniqueSites[i*numSeqs+j] << " ";
		}
		cout << "[" << numEachUniqueSites[i] << "]" << endl;
	}
}

//==================================================================================//
// the following functions are not directly related to this software, but are useful
// for converting one format to another format
//==================================================================================//

// convert the phy format into maze format
void phy2maze(char* file) {
	string aline;
	vector<string> seqs;
	vector<string> names;
	ifstream fin;

	// open the file
	fin.open(file);

	// skip the first line
	getline(fin,aline);

	// get the names and the sequences
	while (getline(fin,aline)) {
		trim(aline); // remove the space or the nonprintable characters at the end
		if (aline.length() > 0 && aline[0]!=' ') {
			// cannot start with a space

			// read until the space
			int i=0;
			while (i < (int)aline.length() && aline[i]!=' ')
				i++;

			// name: aline[0...i-1]
			names.push_back(aline.substr(0,i));

			// skip all the space
			while (i < (int)aline.length() && aline[i]==' ')
				i++;

			// seq : aline[i... end]
			seqs.push_back(aline.substr(i));
		}
	}

	// close the file
	fin.close();

	// print out the sequences in maze format
	for (int i=0; i<(int)names.size(); i++) {
		cout << ";no comment" << endl;
		cout << names[i] << endl;
		for (int j=0; j<(int)seqs[i].length(); j++) {
			if (j>0 && (j%60==0))
				cout << endl;
			char c = seqs[i].at(j);
			/*
			 if (c == 'T')
			 c = 'U';
			 */
			cout << c;
		}
		cout << endl;
	}
}


// read the alignments in maze format
// remove the columns with gaps
// then output the alignment in maze format
void maze_remove_gaps(char* file) {
	string aline;
	vector<string> seqs;
	vector<string> names;
	ifstream fin;
	fin.open(file);
	bool new_seq = true;
	while (getline(fin, aline)) {
		trim(aline);
		if (aline.length() == 0)
			continue;
		if (aline[0]==';') {
			new_seq = true;
		} else if (new_seq) {
			// names
			names.push_back(aline);
			seqs.push_back("");
			new_seq = false;
		} else {
			// sequences
			seqs[seqs.size()-1].append(aline);
		}
	}
	int seqNum = (int) names.size();
	if (seqNum != (int)seqs.size()) {
		cerr << "Error! # of names does not match with # of seqs" << endl;
		exit(1);
	}
	int seqLen = 0;
	if (seqNum > 0) {
		seqLen = (int) seqs[0].length();
		for (int i=1; i<(int)seqs.size(); i++) {
			if (seqLen!=(int)seqs[i].length()) {
				cerr << "Error! The " << i+1 << "-th sequence's length is not consistent with the previous ones, which is: " << seqLen << endl;
				exit(1);
			}
		}
	}

	// just keep the valid columns
	vector<int> validCols;
	for (int i=0; i<seqLen; i++) {
		bool isValid = true;
		for (int j=0; j<seqNum; j++) {
			char c = toupper(seqs[j].at(i));
			if (c!='A' && c!='C' && c!='G' && c!='U') {
				isValid = false;
			}
			/*
			if (seqs[j].at(i)=='-')
				isValid = false;
			 */
		}
		if (isValid) {
			validCols.push_back(i);
		}
	}

	cerr << "Number of valid columns is: " << validCols.size() << endl;
	// print out the sequences in maze format
	for (int i=0; i<seqNum; i++) {
		cout << ";no comment" << endl;
		cout << names[i] << endl;
		for (int j=0; j<(int)validCols.size(); j++) {
			if (j>0 && (j%60==0))
				cout << endl;
			char c = seqs[i].at(validCols[j]);
			/*
			if (c == 'U')
				c = 'T';
			 */
			cout << c;
		}
		cout << endl;
	}
}

// convert the maze format into phy format
void maze2phy(char* file) {
	string aline;
	vector<string> seqs;
	vector<string> names;
	ifstream fin;
	fin.open(file);
	bool new_seq = true;
	while (getline(fin, aline)) {
		trim(aline);
		if (aline.length() == 0)
			continue;
		if (aline[0]==';') {
			new_seq = true;
		} else if (new_seq) {
			// names
			names.push_back(aline);
			seqs.push_back("");
			new_seq = false;
		} else {
			// sequences
			seqs[seqs.size()-1].append(aline);
		}
	}
	int seqNum = (int) names.size();
	if (seqNum != (int)seqs.size()) {
		cerr << "Error! # of names does not match with # of seqs" << endl;
		exit(1);
	}
	int seqLen = 0;
	if (seqNum > 0) {
		seqLen = (int) seqs[0].length();
		for (int i=1; i<(int)seqs.size(); i++) {
			if (seqLen!=(int)seqs[i].length()) {
				cerr << "Error! The " << i+1 << "-th sequence's length is not consistent with the previous ones, which is: " << seqLen << endl;
				exit(1);
			}
		}
	}
	// print out the alignment in phy format
	cout << seqNum << "  " << seqLen << endl;
	for (int i=0; i<seqNum; i++) {
		cout << names[i] << "      " << seqs[i] << endl;
	}
}

// convert the phy format into FASTA format
void phy2fasta(char* file) {
	string aline;
	ifstream fin;
	vector<string> token;
	int i;

	fin.open(file);
	// skip the first line
	getline(fin,aline);
	while (getline(fin,aline)) {
		trim(aline);
		if (aline.length() > 0) {
			tokenizer(aline, " \t", &token);
			// token[0] - name
			// token[1] - seq
			cout << ">" << token[0] << endl;
			for (i=0; i<(int)token[1].length(); i++) {
				if (i>0 && i%60==0)
					cout << endl;
				cout << token[1].at(i);
			}
			cout << endl;
		}
	}
	fin.close();
}

//=============================================================
// For testing only
//=============================================================

/*
// get the pdist between two char arrays
int p_dist(char* site1, char* site2, int length) {
	int d = 0;
	int i;
	for (i=0; i<length; i++) {
		if (site1[i] != site2[i])
			d++;
	}
	return d;
}
*/

// return the binary array of the site
// nucl : [0,3]
#ifndef USE_64_BIT_BINARY
unsigned int* Alignment::returnSiteB(int siteB_ID, int nucl) {
	// return &(sitesBinary[siteB_len*numUniqueSites*nucl + siteB_ID*siteB_len]);
	return &(sitesBinary[siteB_len*4*siteB_ID + nucl*siteB_len]);
}
#else
unsigned long long* Alignment::returnSiteB(int siteB_ID, int nucl) {
	// return &(sitesBinary[siteB_len*numUniqueSites*nucl + siteB_ID*siteB_len]);
	return &(sitesBinary[siteB_len*4*siteB_ID + nucl*siteB_len]);
}
#endif

// build the binary array for computation of p-dist
void Alignment::buildSitesArray() {
	// cout << "[enter buildSitesArray]" << endl << flush;
	if (sitesBinary == NULL) {
#ifndef USE_64_BIT_BINARY
		siteB_len = (numSeqs + 31) / 32;
		sitesBinary = (unsigned int*) malloc(siteB_len*numUniqueSites*4*sizeof(unsigned int));
		memset(sitesBinary, 0, siteB_len*numUniqueSites*4*sizeof(unsigned int));
		unsigned int* currBinSites;
		// array for temporary use
		tmpArray = (unsigned int*) malloc(siteB_len * sizeof(unsigned int));
		unsigned int t;
#else
		siteB_len = (numSeqs + 63) / 64;
		sitesBinary = (unsigned long long*) malloc(siteB_len*numUniqueSites*4*sizeof(unsigned long long));
		memset(sitesBinary, 0, siteB_len*numUniqueSites*4*sizeof(unsigned long long));
		unsigned long long* currBinSites;
		// array for temporary use
		tmpArray = (unsigned long long*) malloc(siteB_len * sizeof(unsigned long long));
		unsigned long long t;
#endif
		int i,j,k;
		int p,q;
		char* currSites;
		char mask;
		for (i=0; i<4; i++) {
			mask = (1 << i);
			for (k=0; k<numUniqueSites; k++) {
				currBinSites = returnSiteB(k, i);
				currSites = &(uniqueSites[k*numSeqs]);
				for (j=0; j<numSeqs; j++) {
#ifndef USE_64_BIT_BINARY
					p=j/32;
					q=j%32;
#else
					p=j/64;
					q=j%64;
#endif
					t = ((currSites[j] & mask) >> i);
					currBinSites[p] |= (t << q);
				}
			}
		}
	}
	// cout << "[leave buildSitesArray]" << endl << flush;
}

// get the pdist between two char arrays (newer version)
// length of tmpArray = siteB_len
int Alignment::p_dist(int siteID1, int siteID2) {
	// cout << "[enter p_dist siteID1=" << siteID1 << " siteID2=" << siteID2 << "]" << endl << flush;
	int d=0;
	int i,j;
#ifndef USE_64_BIT_BINARY
	unsigned int *siteBin1[4];
	unsigned int *siteBin2[4];
#else
	unsigned long long *siteBin1[4];
	unsigned long long *siteBin2[4];
#endif
	for (i=0; i<4; i++) {
		siteBin1[i] = returnSiteB(siteID1, i);
		siteBin2[i] = returnSiteB(siteID2, i);
	}
	for (j=0; j<siteB_len; j++) {
		tmpArray[j] = (siteBin1[0][j] ^ siteBin2[0][j]) | (siteBin1[1][j] ^ siteBin2[1][j]) | (siteBin1[2][j] ^ siteBin2[2][j]) | (siteBin1[3][j] ^ siteBin2[3][j]);
	}
	for (j=0; j<siteB_len; j++) {
#ifndef USE_64_BIT_BINARY
		d+= popCount(tmpArray[j]);
		// d += __builtin_popcount(tmpArray[j]);
#else
		d += popCount64(tmpArray[j]);
		//d += __builtin_popcountll(tmpArray[j]);
#endif
	}
	// cout << "[quit p_dist]" << endl << flush;
	// cerr << siteID1 << " " << siteID2 << " " << d << endl;
	return d;
}

// original version
// get the pdist between two char arrays
int Alignment::p_dist_orig(int siteID1, int siteID2) {
	int d = 0;
	int i;
	char* site1 = &(uniqueSites[siteID1 * numSeqs]);
	char* site2 = &(uniqueSites[siteID2 * numSeqs]);
	for (i=0; i<numSeqs; i++) {
		if (site1[i] != site2[i])
			d++;
	}
	// cerr << siteID1 << " " << siteID2 << " " << d << endl;
	return d;
}


// across the adjacent sites, compute the total number of differences
// this is not completed
void Alignment::getNumDiffBWAdjSites() {
	if (pdistOrder==NULL)
		return;
	int i;
	int diff = 0;
	for (i=1; i<numUniqueSites; i++) {
		diff += p_dist(i-1, i);
	}
	cout << "Total number of differences: " << diff << endl;
	cout << "Total number of uniqueSites: " << numUniqueSites << endl;
}

// reorder the unique sites by using greedy approach
void Alignment::reorderUniqueSites() {
	// cout << "[enter reorderUniqueSites]" << endl << flush;
	 int i,j;
	 int tmp_i;
	 int min_dist, curr_dist;
	 int min_j;
	
	 if (numUniqueSites > 0) {
		 buildSitesArray(); // build the binary sites
		#ifdef GET_TIME_STAT
				 curr_time = clock();
				 second = (curr_time - pre_time) * 1000.0 / CLOCKS_PER_SEC;
				 printf("Time used after building the binary site arrays: %4.2f ms\n", second);
				 pre_time = curr_time;
		#endif
		 if (pdistOrder != NULL)
			 free (pdistOrder);
		 pdistOrder = (int*) malloc (numUniqueSites * sizeof(int));
		 for (i=0; i<numUniqueSites; i++)
			 pdistOrder[i] = i;
		 for (i=0; i<numUniqueSites-2; i++) {
			 // cout << "i=" << i << endl << flush;
			 min_dist = numSeqs * 4;
			 min_j = -1;
			 for (j=i+1; j<numUniqueSites && min_dist > 1; j++) {
			 // for (j=i+1; j<numUniqueSites; j++) {
				 curr_dist = p_dist(pdistOrder[i], pdistOrder[j]);
				 // curr_dist = p_dist_orig(pdistOrder[i], pdistOrder[j]);
				 if (curr_dist < min_dist) {
					 min_dist = curr_dist;
					 min_j = j;
				 }
			 }
			 if (min_j > i+1) {
				 // swapping between pos i+1 and pos min_j
				 tmp_i = pdistOrder[i+1];
				 pdistOrder[i+1] = pdistOrder[min_j];
				 pdistOrder[min_j] = tmp_i;
			 }
		 }
	 }
	
	// update the uniqueSites, numEachUniqueSites, and isConstantSites
	char *siteFr, *siteTo;
	char* uniqueSitesNew = (char*) malloc(numUniqueSites*numSeqs*sizeof(char));
	double*  numEachUniqueSitesNew = malloc_decimals(numUniqueSites);
	char* isConstantSitesNew = (char*) malloc(numUniqueSites*sizeof(char));
	for (i=0; i<numUniqueSites; i++) {
		siteFr = &(uniqueSites[pdistOrder[i]*numSeqs]);
		siteTo = &(uniqueSitesNew[i*numSeqs]);
		memcpy(siteTo, siteFr, numSeqs*sizeof(char));
		numEachUniqueSitesNew[i] = numEachUniqueSites[pdistOrder[i]];
		isConstantSitesNew[i] = isConstantSites[pdistOrder[i]];
	}
	free(uniqueSites);
	free_decimals(numEachUniqueSites);
	free(isConstantSites);
	uniqueSites = uniqueSitesNew;
	numEachUniqueSites = numEachUniqueSitesNew;
	isConstantSites = isConstantSitesNew;
}

/*
// build the array buildPrePosWithSmallDiff()
void Alignment::buildPrePosWithSmallDiff() {
	if (numUniqueSites > 0) {
		prePosWithSmallDiff = (int *) malloc (numUniqueSites * sizeof(int));
		int i,j;
		prePosWithSmallDiff[0] = 0;
		if (numUniqueSites > 1)
			prePosWithSmallDiff[1] = 0;
		int min_diff; int min_j; int curr_diff;
		char* sitei;
		char* sitej;
		for (i=2; i<numUniqueSites; i++) {
			min_j = -1;
			min_diff = numSeqs;
			sitei = &(uniqueSites[i*numSeqs]);
			for (j=0; j<i; j++) {
				sitej = &(uniqueSites[j*numSeqs]);
				curr_diff = p_dist(sitei, sitej, numSeqs);
				if (curr_diff < min_diff) {
					min_diff = curr_diff;
					min_j = j;
					if (curr_diff == 1)
						break;
				}
			}
			prePosWithSmallDiff[i] = min_j;
		}
		
		// get the statistics (for testing only)
		curr_diff = 0;
		for (i=1; i<numUniqueSites; i++) {
			sitei = &(uniqueSites[i*numSeqs]);
			sitej = &(uniqueSites[prePosWithSmallDiff[i]*numSeqs]);
			curr_diff += p_dist(sitei, sitej, numSeqs);
		}
		cout << "By using prePosWithSmallDiff, total number of differences: " << curr_diff << endl;
	}
}
*/

// for building the changedNodes array
// after building (and sorting) the uniqueSites
void Alignment::buildChangedNodes(int* topMatrix, int numLineTopMat) {
	int numInternalNodes = numSeqs-1;
	changedNodes = (char*) malloc (sizeof(char)*(numInternalNodes+numSeqs+1)*numUniqueSites);
	memset(changedNodes, 0, sizeof(char)*(numInternalNodes+numSeqs+1)*numUniqueSites);
	int i,j;
	char* currChangedNodes;
	char* preSites;
	char* currSites;
	int leftNode, rightNode, l, r;
	currChangedNodes = &(changedNodes[0]);
	for (j=0; j<numSeqs+numInternalNodes+1; j++) {
		 currChangedNodes[j] = 1;
	}
	currSites = &(uniqueSites[0]);
	for (i=1; i<numUniqueSites; i++) {
		 currChangedNodes = &(changedNodes[i*(numInternalNodes+numSeqs+1)]);
		 preSites = currSites;
		 currSites = &(uniqueSites[i*numSeqs]);
		 for (j=0; j<numSeqs; j++) {
			  if (preSites[j] != currSites[j]) {
					currChangedNodes[numSeqs-j-1] = 1;
			  }
		 }
		for (j=0; j<numLineTopMat; j++) {
			leftNode = topMatrix[j*2];
			rightNode = topMatrix[j*2+1];
			// l = (leftNode < 0)?(-leftNode-1):(leftNode-1+numSeqs);
			l = leftNode+numSeqs;
			// r = (rightNode < 0)?(-rightNode-1):(rightNode-1+numSeqs);
			r = rightNode+numSeqs;
			if (currChangedNodes[l] | currChangedNodes[r]) {
				currChangedNodes[j+numSeqs+1] = 1;
			}
		}
	}
}

int popCount(unsigned int x) {
	
	x -= ((x >> 1) & 0x55555555);
	x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
	x = (((x >> 4) + x) & 0x0f0f0f0f);
	x += (x >> 8);
	x += (x >> 16);
	return(x & 0x0000003f);
	
}

int popCount64(unsigned long long x) {

	x -= ((x >> 1) & 0x5555555555555555);
	x = (((x >> 2) & 0x3333333333333333) + (x & 0x3333333333333333));
	x = (((x >> 4) + x) & 0x0f0f0f0f0f0f0f0f);
	x += (x >>  8);
	x += (x >> 16);
	x += (x >> 32);
	return (x & 0x000000000000007f);
}

/*
//This is better when most bits in x are 0
//It uses 3 arithmetic operations and one comparison/branch per "1" bit in x.
int popCount64(unsigned long long x) {
	int count;
	for (count=0; x; count++)
		x &= x-1;
	return count;
}
*/