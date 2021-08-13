#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "parameters.h"
#include "jacobi_eigenvalue.h"
#include "matrix.h"

using namespace std;

// compute the eigen value for symmetric 4x4 matrix
// if success, return 1, otherwise return 0
// input: in[4x4], 4x4 symmetric matrix
// output: eigenVectors[4x4], the matrix of eigenvectors.
// output: eigenValues[4], the eigenvalues, in descending order.
int eigenvalue_4_4(double* in, double* eigenVectors, double* eigenValues);

// To solve a quartic equation: a x^4 + b x^3 + c x^2 + d x + e = 0
// Return 1 if success, 0 otherwise
// Input: in[0] = a; in[1] = b; in[2] = c; in[3] = d; in[4] = e; 
// Output: The values of the roots
int quartic_solver(double* in, double* out);

// To get the real answers for a cubic equation
// Return 1 if success, 0 otherwise
// Input: in[0] = a; in[1] = b; in[2] = c; in[3] = d
//         i.e. a x^3 + b x^2 + c x + d = 0
// Output: out[0] = real ans 1; out[1] = real ans 2; out[2] = real ans 3;
//         num_real_ans = the number of real answers
int cubic_solver(double* in, double* out, int& num_real_ans);

// for testing the function eigenvalue_4_4 for num_times times
void test_eigenvalue_4_4(int num_times);

// for testing the quartic solver for num_times times
void test_quartic_solver(int num_times);

// for testing the cubic solver for num_times times
void test_cubic_solver(int num_times);

// For solving one equation
// If there is more than one answer, return one of the possible answers.
// Return 1 when success; return 0 when failure
// Input: numTerms - number of terms
//        A[0] X0 + A[1] X1 + ... + A[numTerms-2] X(numTerms-2) + A[numTerms-1] = 0
// Output: out[0] = X0; out[1] = X1; ....  ; out[numTerms-2] = X(numTerms-2)
// accept_all_zeros : 1 - accept the results that X,Y = 0; 0 - otherwise
int solve_1_equation(double* A, double* out, int num_terms, int accept_all_zeros);

// For solving equations with two variables
// If there is more than one answer, return one of the possible answers.
// Return 1 when success; return 0 when failure
// Input: numEq - number of equations
//        (1) A[0] X + A[1] Y + A[2] = 0
//        (2) A[3] X + A[4] Y + A[5] = 0
//                       ....
//        (N) A[(N-1)*3] X + A[(N-1)*3+1] Y A[(N-1)*3+2] = 0
// Output: out[0] = X; out[1] = Y 
// accept_all_zeros : 1 - accept the results that X,Y = 0; 0 - otherwise
int solve_equations_2_variables(int numEq, double* A, double* out, int accept_all_zeros);

// For testing the function "solve_equations_2_variables"
void testing_solve_equations_2_variables(int num_times);

// For testing the function "solve_N_equations"
void testing_solve_N_equations(int num_times, int num_equations);

// For solving N equations with N variables (for N >= 2)
// If there is more than one answer, return one of the possible answers.
// Return 1 when success; return 0 when failure
// Input: 
//        (1)   A[0] X1 +   A[1] X2 +   A[2] X3 + ... +  A[N-1] XN +    A[N] = 0
//        (2) A[N+1] X1 + A[N+2] X2 + A[N+3] X3 + ... +   A[2N] XN + A[2N+1] = 0
//        (3) ......
//        (N) A[(N-1)(N+1)] X1 + A[(N-1)(N+1)+1] X2 + ...+ A[N(N+1)-2] XN + A[N(N+1)-1] = 0
// Output: out[0] = X1; out[1] = X2; out[2] = X3; ...; out[N-1] = XN
// accept_all_zeros : 1 - accept the results that X1,X2,...,XN = 0; 0 - otherwise
int solve_N_equations(double* A, double* out, int num_equations, int num_terms, int accept_all_zeros);

// ===========================================================
// Other supporting functions
// ===========================================================
int isTheSameLoosely(double x1, double x2);

int isTheSame(double x1, double x2);

int isTheSame2(double x1, double x2);

void roundToZero(double& x);

void roundToZero2(double& x);

// void roundToZero3(double& x);

// Return the status of two equations, where
// (1) not parallel
// (2) parallel, and no intersection point
// (3) parallel and overlap each other
int check_two_equations(double* A, double* B, int num_terms);

// Return the status of two equations, where
// (1) not parallel
// (2) parallel, and no intersection point
// (3) parallel and overlap each other
// int check_two_equations_loosely(double* A, double* B, int num_terms);

// Check the equation
// Return 1 if all terms are zeros except the last term (i.e. invalid equation)
// Return 2 if all terms are zeros (i.e. the equation has no effect)
// Return 0 otherwise
int check_equation(double* A, int num_terms);

// reduce the variables in the equation
void reduce_variable(double* A, double* new_A, int num_terms, int remove_term);
