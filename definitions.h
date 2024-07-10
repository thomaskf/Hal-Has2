/*
 *
 * definitions.h
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


#ifndef RAL_RAS_definitions_h
#define RAL_RAS_definitions_h

#include <config.h>

// ==================================================================================================
// The version of the program
// ==================================================================================================

// The version number is defined in the CMakeLists.txt file
//#ifndef VERSION
//    #define VERSION "2.7"
//#endif

// ==================================================================================================
// The contact person of the program
// ==================================================================================================

// #define CONTACTPERSON "Lars Jermiin <Lars.Jermiin@anu.edu.au>, Thomas Wong <Thomas.Wong@anu.edu.au>"
// LSJermiin 2024-Mar-10: Replaced line above with line below
#define CONTACTPERSON "Lars Jermiin <Lars.Jermiin@universityofgalway.ie>, Thomas Wong <Thomas.Wong@anu.edu.au>"

// ==================================================================================================
// Uncomment the following line if the dataset is very large, in order to speed up the program
// #define LARGE_DATA
// ==================================================================================================

// ==================================================================================================
// Number of decimal places
// ==================================================================================================

#define PARAM_DECI 5.0 // the number of decimal places for the resulting parameter values
// #define ANS_DECI 2.0 // the number of decimal places for the resulting IC and log-likelihood values
// LSJermiin 2024-Mar-10: Replaced line above with line below
#define ANS_DECI 5.0 // the number of decimal places for the resulting IC and log-likelihood values


// ==================================================================================================
// Uncomment the following line if this is a release vesion
#define RELEASE_VERSION
// ==================================================================================================


// ==================================================================================================
// Define the default configuration of optimization method
// ==================================================================================================

// uncomment the following line to use old set of parameters
// #define USE_OLD_PARAM

#ifdef USE_OLD_PARAM
	#define DEFAULT_OPT_MODE "1"  // the mode of the optimization method
#else
	#define DEFAULT_OPT_MODE "5"  // the mode of the optimization method
#endif

#ifdef USE_OLD_PARAM
	#define DEFAULT_OPT_MAXIT "1000" // the maxIT of the optimization method
#else
	#define DEFAULT_OPT_MAXIT "10" // the maxIT of the optimization method
#endif

#define PRECISE_MODE "1,3,6"
#define PRECISE_MAXIT "1000,10,1"

// ==================================================================================================
// Define the default threshold for the log-likelihood value
// ==================================================================================================

#define DEFAULT_LOG_L_THRES 0.1  // the threshold for the log-likelihood value
#define DEFAULT_LOG_L_THRES_FOR_CONSENSUS 10.0

// ==================================================================================================
// Define the default topK
// ==================================================================================================

#define DEFAULT_TOP_K_OLD 3  // the old BU algorithm
#define DEFAULT_TOP_K 10 // for the new BU algorithm

// ==================================================================================================
// Define the default Information Criteria
// ==================================================================================================

#define DEFAULT_IC 3  // BIC
#define TOTAL_NUM_IC 4 // total number of different information criteria available to choose

// ==================================================================================================
// Define the default tag ID for MPI communication
// ==================================================================================================

#define DEFAULT_RAL_TAGID 777
#define DEFAULT_RAS_TAGID 888

// ==================================================================================================
// Define the default starting state for the optimization used in NEW-BU algorithm
// ==================================================================================================

#define DEFAULT_OPT_START 0
#define DEFAULT_OPT_START_OLD 0
// 0: from the original state;
// 1: from the tunned state of the parent configuration
// 2: from the tunned state of involving one rate group

// ==================================================================================================
// Define the value of very large number
// ==================================================================================================

#define LARGE_NUMBER 1e77

// ==================================================================================================
// For the regularization of BIC
// Define the value of Creg
// BICreg = -2 (log L) + (K + C) (log N) where L is the likelihood,
//          K is the number of free parameters, N is the number of sites
//          C is defined in the following, which is between [ 0 , K'-K ]
//          K' is the number of free parameters when an additional new rate matrix is used
//          For HAL model, K'-K is always 8
// ==================================================================================================

#ifndef C_REG
	#define C_REG 0 // has to be between 0 and 8
#endif

// ==================================================================================================
// Uncomment the following line to fix the model as GTR+I+FO during the HAL analysis
// #define GTRIFO
// ==================================================================================================

// ==================================================================================================
// Uncomment the following line to only consider the model type 5 during the HAS analysis
// #define MODEL5ONLY
// ==================================================================================================

// ==================================================================================================
// The locks designed for multi-threads
// ==================================================================================================

// a set of locks used in each CPU thread
typedef struct ThreadLocks
{
	pthread_mutex_t lock_cpuThresRecycle;
	pthread_mutex_t lock_output;
	pthread_mutex_t lock_optLK;
	pthread_mutex_t lock_optInfo;
	pthread_mutex_t lock_master;
	pthread_mutex_t lock_maxMasterID;
	pthread_mutex_t lock_createThread;
	pthread_mutex_t lock_jobAllocate;
	pthread_mutex_t lock_checkpoint;
} ThreadLocks;

// for SSE functions
// define the data type for SSE
// #define ALIGN16 __attribute__((aligned(16)))

// uncomment the following two lines if want to get the time statistics
// #define GET_TIME_STAT
// #define GET_TIME_STAT_DETAIL

#ifdef GET_TIME_STAT
#include "global.h"
#endif

static double  base_arrays[64] = {
	0, 0, 0, 0, // unclassified
	1, 0, 0, 0, // A
	0, 1, 0, 0, // C
	1, 1, 0, 0, // A or C
	0, 0, 1, 0, // G
	1, 0, 1, 0, // A or G
	0, 1, 1, 0, // C or G
	1, 1, 1, 0, // A, C or G
	0, 0, 0, 1, // T
	1, 0, 0, 1, // A or T
	0, 1, 0, 1, // C or T
	1, 1, 0, 1, // A, C or T
	0, 0, 1, 1, // G or T
	1, 0, 1, 1, // A, G or T
	0, 1, 1, 1, // C, G or T
	1, 1, 1, 1  // A, C, G or T
};

#endif
