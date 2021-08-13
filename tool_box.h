/*
 *
 * tool_box.h
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

#ifndef __sim_tree__tool_box__
#define __sim_tree__tool_box__

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cstdio>
#include "definitions.h"

using namespace std;

#define SMALLVARIATION 0.00001

// split the seq into many parts by "separators"
// the vector<string> *result cannot be NULL
void tokenizer(string seq, string separators, vector<string>* result);

// return true if c is space or non-printable character
bool isSpaceOrNonprintableChar(char c);

// remove the non-printable characters at the end of the string
void trim(string& str);

// sum of the product of two vectors
double sumProduct(double* vectorA, double* vectorB, int len);

map<string,int>* genReverseMap(vector<string>* strArray);

void initialRandSeed();

// generate a random integer from the integer "frInt" to the integer "toInt"
// suppose the the random seed has been initialized
// (i.e. the function "initialRandSeed()" has been executed)
int genRandInt(int frInt, int toInt);

// check whether str is a non-negative integer
bool isNonNegativeInt(char* str);

string removeExtension(string str);

// check the number of bits of the OS
// if the number of bits is not 64,
// then output the warning message
void checkOSNumOfBits();

// integer to string
string intToStr(int i);

// integer to string
// minLen : minimum length of the integer.
// If the number of digits < minimum length, then the integer will be displayed with leading zeros, unless the integer is zero
string intToStr(int i, int minLen);

// double to string
string longDoublToStr(double d, double decimalPlace);

// round the double to a certain number of decimal place
double roundNumber(double d, int decimalPlace);

// compare two double numbers
// if their difference <= SMALLVARIATION
// then return 1, otherwise 0
int isSame(double d1, double d2);

// copy file
void filecopy(char * dest, char * src);

#endif /* defined(__sim_tree__tool_box__) */
