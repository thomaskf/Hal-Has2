#include "eigenvalue_4_by_4.h"

#define PI          3.141592653589793238462643383279502884L /* pi */
#define QUITE_SMALL_NUMBER 1.0E-3
#define SMALL_NUMBER 1.0E-5
#define VERY_SMALL_NUMBER 1.0E-10
#define SUPER_SMALL_NUMBER 1.0E-30

// #define DEBUG_MODE

// compute the eigen value for symmetric 4x4 matrix
// if success, return 1, otherwise return 0
// input: in[4x4], 4x4 symmetric matrix
// output: eigenVectors[4x4], the matrix of eigenvectors.
// output: eigenValues[4], the eigenvalues, in descending order.
int eigenvalue_4_4(double* in, double* eigenVectors, double* eigenValues) {

	int x,y;

	// assume the matrix is symmetric
	double a,b,c,d,e,f,g,h,i,j;
	a=in[0];b=in[1];c=in[2];d=in[3];e=in[5];f=in[6];g=in[7];h=in[10];i=in[11];j=in[15];

	// compute the determinant of the matrix: (IN - \lambda I)
	// it is a quartic equation in terms of \lambda
	double quartic[5];
	quartic[0]=1.0;
	quartic[1]=-(a+e+h+j);
	quartic[2]=a*e + e*h + a*h + a*j + e*j + h*j - i*i - f*f - g*g - b*b - c*c - d*d;
	quartic[3]=-(a*e*h + a*e*j + a*h*j + e*h*j) - 2.0*f*g*i + e*i*i + a*i*i + a*f*f + j*f*f + a*g*g + h*g*g + h*b*b + j*b*b - 2.0*b*f*c - 2.0*b*d*g - 2.0*c*d*i + e*c*c + j*c*c + e*d*d + h*d*d;
	quartic[4]=h*a*e*j + 2.0*a*f*g*i - a*e*i*i - a*f*f*j - a*g*g*h - b*b*h*j - 2.0*b*f*d*i - 2.0*b*c*g*i + b*b*i*i + 2.0*b*f*c*j + 2.0*b*d*g*h + 2.0*c*d*i*e + c*c*g*g - 2.0*c*d*f*g - e*j*c*c - e*h*d*d + d*d*f*f;

	for (x=0; x<5; x++)
		roundToZero2(quartic[x]);

#ifdef DEBUG_MODE
	cout << "Quartic equation:" << endl;
	cout << "(" << quartic[0] << ") X^4 +"; 
	cout << "(" << quartic[1] << ") X^3 +"; 
	cout << "(" << quartic[2] << ") X^2 +"; 
	cout << "(" << quartic[3] << ") X^1 +"; 
	cout << "(" << quartic[4] << ")" << endl;; 
#endif

	// solve the quartic equation and get the values of lambdas
	double lambda[4];
	double t;
	int success;
	success = quartic_solver(quartic, lambda);
	if (!success)
		return 0;

	// the smallest value of lambda cannot be smaller than "SUPER_SMALL_NUMBER"
	for (x=0; x<4; x++) {
		if (fabsl(lambda[x]) < SUPER_SMALL_NUMBER) {
			lambda[x] = -SUPER_SMALL_NUMBER;
		}
	}

	// sort lambdas in ascending order
	for (x=0; x<3; x++) {
		for (y=0; y<3; y++) {
			if (lambda[y] > lambda[y+1]) {
				t = lambda[y];
				lambda[y] = lambda[y+1];
				lambda[y+1] = t;
			}
		}
	}

#ifdef DEBUG_MODE
	// show the value of lambdas
	cout << "lambda are: ";
	for (x=0; x<4; x++) {
		if (x>0)
			cout << ", ";
		cout << lambda[x];
	}
	cout << endl;
#endif

	for (x=0; x<4; x++)
		eigenValues[x] = lambda[x];

	// for each lambda value, compute the corresponding eigenVector
	// the equations with four variables
	// (1) (a-l)W +     bX +     cY +     dZ + 0 = 0
	// (2)     bW + (e-l)X +     fY +     gZ + 0 = 0
	// (3)     cW +     fX + (h-l)Y +     iZ + 0 = 0
	// (4)     dW +     gX +     iY + (j-l)Z + 0 = 0

	// (5)  (W0)W +  (X0)X +  (Y0)Y +  (Z0)Z + 0 = 0  where W0,X0,Y0,Z0 are the answers for lambda[0]
	// (6)  (W1)W +  (X1)X +  (Y1)Y +  (Z1)Z + 0 = 0  where W1,X1,Y1,Z1 are the answers for lambda[1]
	// (7)  (W2)W +  (X2)X +  (Y2)Y +  (Z2)Z + 0 = 0  where W2,X2,Y2,Z2 are the answers for lambda[2]

	// for computation of W,X,Y,Z for lambda[0], equations (1)-(4) are needed
	// for computation of W,X,Y,Z for lambda[1], equations (1)-(5) are needed
	// for computation of W,X,Y,Z for lambda[2], equations (1)-(6) are needed
	// for computation of W,X,Y,Z for lambda[3], equations (1)-(7) are needed

	double A[35];
	double R[4];
	double R_magnitude;
	double l;
	double eigenVect[16];
	double* ev;
	int num_equations;
	int num_terms = 5;
	for (x=0; x<4; x++) {
		l = lambda[x];
		// standard equations
		A[ 0]=a-l; A[ 1]=b;   A[ 2]=c;   A[ 3]=d;
		A[ 5]=b;   A[ 6]=e-l; A[ 7]=f;   A[ 8]=g;
		A[10]=c;   A[11]=f;   A[12]=h-l; A[13]=i;
		A[15]=d;   A[16]=g;   A[17]=i;   A[18]=j-l;
		// all constant terms are zeros
		A[4]=A[9]=A[14]=A[19]=A[24]=A[29]=A[34]=0.0;
		// for computation of W,X,Y,Z for lambda[x>0] if lambda[x] = lambda[x-1]
		// if (x>0 && isTheSame(lambda[x],lambda[x-1])) {
		if (x>0) {
			num_equations++;
			for (y=0; y<4; y++)
				A[15+5*(num_equations-4)+y] = ev[y]; // previous values of eigen vectors
		} else {
			num_equations = 4;
		}
		ev = &(eigenVect[x*4]);

		/* roundToZero2(A[0]);
		roundToZero2(A[6]);
		roundToZero2(A[12]);
		roundToZero2(A[18]);*/
#ifdef DEBUG_MODE
		// list out the equations
		int u,w;
		double *curr_A;
		for (u=0; u<num_equations; u++) {
			curr_A = &(A[u*num_terms]);
			for (w=0; w<num_terms-1; w++) {
				if (w>0) {
					cout << " + ";
				}
				cout << "(" << curr_A[w] << ") X" << w+1;
			}
			cout << " + (" << curr_A[num_terms-1] << ") = 0" << endl;
		}
#endif
		if (solve_N_equations(A, R, num_equations, num_terms, 0)) {
#ifdef DEBUG_MODE
			// list out the answer for the equations
			cout << "Result of the equations: ";
			for (y=0; y<4; y++) {
				if (y>0) {
					cout << ", ";
				}
				cout << R[y];
			}
			cout << endl;
#endif
			// normalize the R vector
			R_magnitude = 0.0;
			for (y=0; y<4; y++) {
				R_magnitude += R[y]*R[y];
			}
			R_magnitude = sqrt(R_magnitude);
			if (R_magnitude > 0.0) {
				for (y=0; y<4; y++) {
					R[y] = R[y] / R_magnitude;
				}
			}

			for (y=0; y<4; y++)
				ev[y] = R[y];
		} else {
			return 0;
		}
	}
	for (x=0; x<16; x++)
		eigenVectors[x] = eigenVect[x];
	return 1;
}


// Converting a quartic equation to a depressed quartic
int toDepressedQuartic(double* in, double* out) {
	// Return 1 if success, 0 otherwise
	// Input: in[0] = a4; in[1] = a3; in[2] = a2; in[3] = a1; in[4] = a0;
	//         i.e. a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0 = 0;
	// Output: out[0] = p; out[1] = q; out[2] = r
	//         i.e. y^4 + p y^2 + q y + r = 0; where x = y - a3 / (4 a4)

	double a0,a1,a2,a3,a4;
	a4=in[0]; a3=in[1]; a2=in[2]; a1=in[3]; a0=in[4];

	if (a4 == 0)
		return 0;

	out[0] = (8.0 * a2 * a4 - 3.0 * a3 * a3) / (8.0 * a4 * a4);
	out[1] = (a3 * a3 * a3 - 4.0 * a2 * a3 * a4 + 8.0 * a1 * a4 * a4) / (8.0 * a4 * a4 * a4);
	out[2] = ((-3.0) * a3 * a3 * a3 * a3 + 256.0 * a0 * a4 * a4 * a4 - 64.0 * a1 * a3 * a4 * a4 + 16.0 * a2 * a3 * a3 * a4) / (256.0 * a4 * a4 * a4 * a4);
	return 1;
}

// To get the real answers for a cubic equation
int cubic_solver(double* in, double* out, int& num_real_ans) {
	// Return 1 if success, 0 otherwise
	// Input: in[0] = a; in[1] = b; in[2] = c; in[3] = d
	//         i.e. a x^3 + b x^2 + c x + d = 0
	// Output: out[0] = real ans 1; out[1] = real ans 2; out[2] = real ans 3;
	double a,b,c,d;
	a=in[0];b=in[1];c=in[2];d=in[3];
	if (a==0)
		return 0;
	double delta0 = b*b - 3.0*a*c;
	double delta1 = 2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d;
	double t = delta1*delta1 - 4.0*delta0*delta0*delta0;
	roundToZero2(t);
#ifdef DEBUG_MODE
	cout << "delta0 = " << delta0 << endl;
	cout << "delta1 = " << delta1 << endl;
	cout << "t = " << t << endl;
#endif
	if (t==0.0 && delta0==0) {
		// special case 1
		out[0]=out[1]=out[2]=b/(-3.0*a);
		num_real_ans=3;
		return 1;
	}
	if (t==0.0 && delta0!=0) {
		// special case 2
		out[0]=out[1]=(9.0*a*d-b*c)/(2.0*delta0);
		out[2]=(4.0*a*b*c-9.0*a*a*d-b*b*b)/(a*delta0);
		num_real_ans=3;
		return 1;
	}
	if (t<0.0) {
		// there will be three real answers
		double p = (3.0*a*c - b*b)/(3*a*a);
		double q = (2.0*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
#ifdef DEBUG_MODE
		cout << "p=" << p << endl;
		cout << "q=" << q << endl;
#endif
		if (p<0) {
			num_real_ans=3;
			t = 3.0*q*sqrt(-3.0/p)/(2.0*p);
			if (t >= 1.0)
				t = 0.0;
			else if (t <= -1.0)
				t = PI/3.0;
			else
				t = acos(3.0*q*sqrt(-3.0/p)/(2.0*p))/3.0;
			double t1 = 2.0*sqrt(p/(-3.0));
			double t2 = b/(3.0*a);
			out[0] = t1*cos(t)-t2;
			out[1] = t1*cos(t-2.0*PI/3.0)-t2;
			out[2] = t1*cos(t-4.0*PI/3.0)-t2;
			return 1;
		} else {
			num_real_ans=0;
			return 0;
		}
	}

	// t > 0.0
	if (delta0==0.0)
		t = (delta1 + delta1) / 2.0;
	else
		t = (delta1 + sqrt(t)) / 2.0;
	if (t<0.0) {
		t = -pow(-t,1.0/3.0);
	} else {
		t = pow(t,1.0/3.0);
	}
	out[0] = (b + t + delta0/t)/(-3.0*a);
	num_real_ans=1;
	return 1;
}


// To solve a depressed quartic
int depressed_quartic_solver(double* in, double* out) {
	// Return 1 if success, 0 otherwise
	// Input: in[0] = a; out[1] = b; out[2] = c
	//         i.e. y^4 + a y^2 + b y + c = 0
	// Output: The values of the roots
	double a,b,c;
	a = in[0]; b = in[1]; c = in[2];
	double cubic_in[4];
	cubic_in[0] = 1.0;
	cubic_in[1] = 2.5*a;
	cubic_in[2] = 2.0*a*a - c;
	cubic_in[3] = a*a*a/2.0 - a*c/2.0 - b*b/8.0;
	int success, num_real_ans;
	success = cubic_solver(cubic_in, out, num_real_ans);
	if (!success || num_real_ans < 1)
		return 0;
	double y = out[0];
	double t1, t2, t3, t4, t5;
	t1 = a + 2.0*y;
	if (t1 < 0.0) {
#ifdef DEBUG_MODE
		cout << "Error A" << endl;
#endif
		return 0;
	}
	t1 = sqrt(t1);
	t2 = 3.0*a + 2.0*y;
	t3 = a + 2.0*y;
	if (t3 <= 0.0) {
#ifdef DEBUG_MODE
		cout << "Error B" << endl;
#endif
		return 0;
	}
	t3 = 2.0*b/(sqrt(t3));

	t4 = -(t2+t3);
	t5 = -(t2-t3);
	roundToZero2(t4);
	if (t4<0.0) {
#ifdef DEBUG_MODE
		cout << "t4 = " << t4 << endl;
		cout << "Error C" << endl;
#endif
		return 0;
	}

	roundToZero2(t5);
	if (t5<0.0) {
#ifdef DEBUG_MODE
		cout << "t5=" << t5 << endl;
		cout << "Error D" << endl;
#endif
		return 0;
	}
	t4 = sqrt(t4);
	t5 = sqrt(t5);
	out[0] = (t1+t4)/2.0;
	out[1] = (t1-t4)/2.0;
	out[2] = (-t1+t5)/2.0;
	out[3] = (-t1-t5)/2.0;
	return 1;
}

// To solve a quartic equation: a x^4 + b x^3 + c x^2 + d x + e = 0
// Only consider the situation when all roots are real
// Return 1 if success, 0 otherwise
// Input: in[0] = a; in[1] = b; in[2] = c; in[3] = d; in[4] = e; 
// Output: The values of the real roots
int quartic_solver(double* in, double* out) {

	double a,b,c,d,e;
	double t1,t2,t3,t4;
	double delta0,delta1,delta,D;
	double p,q,S,Q;
	int i;
	int use_depressed_quartic_method = 0;

	double divisor = 0.0;
	for (i=0; i<5; i++) {
		if (divisor < fabsl(in[i]))
			divisor = fabsl(in[i]);
	}

	a = in[0]/divisor;
	b = in[1]/divisor;
	c = in[2]/divisor;
	d = in[3]/divisor;
	e = in[4]/divisor;

#ifdef DEBUG_MODE
	cout << "a = " << a << "; b = " << b << "; c = " << c << "; d = " << d << "; e = " << e << endl;
#endif

	if (a==0.0) {
#ifdef DEBUG_MODE
		cout << "error 1" << endl;
#endif
		return 0;
	}


	// based on the formula in http://en.wikipedia.org/wiki/Quartic_function
	delta0 = c*c - 3.0*b*d + 12.0*a*e;
	delta1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*a*d*d + 27.0*b*b*e - 72.0*a*c*e;
	D = 64.0*a*a*a*e - 16.0*a*a*c*c + 16.0*a*b*b*c - 16.0*a*a*b*d - 3.0*b*b*b*b;
	t1 = delta1*delta1 - 4.0*delta0*delta0*delta0;
	delta = t1/(-27.0);
	p = (8.0*a*c - 3.0*b*b) / (8.0*a*a);

#ifdef DEBUG_MODE
	cout << "delta = " << delta << endl;
#endif

	// for some figures which are too closed to zero, round them to zero
	roundToZero2(D);
	roundToZero2(delta);
	roundToZero2(delta0);

#ifdef DEBUG_MODE
	cout << "delta = " << delta << endl;
	cout << "delta0 = " << delta0 << endl;
	cout << "delta1 = " << delta1 << endl;
	cout << "p = " << p << endl;
	cout << "D = " << D << endl;
	cout << "t1 = " << t1 << endl;
#endif

	// consider the special cases
	/*
	if (delta < 0.0) {
	#ifdef DEBUG_MODE
		cout << "The equation has two real roots and two complex conjugate roots" << endl;
	#endif
		return 0;
	} else */ 
	if (delta > 0.0) {
		if (p > 0.0 || D > 0.0) {
#ifdef DEBUG_MODE
			cout << "The equation has two pairs of complex conjugate roots" << endl;
#endif
			return 0;
		} // else: all four roots are real
	} else if (delta==0.0) {
		// delta = 0
		if (D==0.0) {
			if (delta0 == 0.0) {
				// all four roots are the same
				out[0]=out[1]=out[2]=out[3]=b/(-4.0 * a);
				return 1;
			}
			if (p > 0.0) {
#ifdef DEBUG_MODE
				cout << "The equation has two complex conjugate double roots" << endl;
#endif
				return 0;
			}
		} else {
			if (delta0 == 0.0) {
				// there is a triple root and a simple root, all real.
				t1 = 9.0*b*b - 24.0*a*c;
				roundToZero2(t1);
				if (t1 < 0.0) {
#ifdef DEBUG_MODE
					cout << "Error in finding triple root answer!" << endl;
#endif
					return 0;
				}
				t2=(-3.0*b + sqrt(t1))/(12.0*a);
				t3=b/(-a)-3.0*t2;
				if (!isTheSame(t2*t2*t2*t3, e/a)) {
					t2=(-3.0*b - sqrt(t1))/(12.0*a);
					t3=b/(-a)-3.0*t2;
				}
				out[0]=out[1]=out[2]=t2;
				out[3]=t3;
				return 1;
			} else if (p>0.0 || D>0.0) {
#ifdef DEBUG_MODE
				cout << "The equation has a real double root and two complex conjugate roots" << endl;
#endif
				return 0;
			}
		}
	}

	if (t1 < 0.0) {
#ifdef DEBUG_MODE
		cout << "here 1" << endl;
#endif
		if (delta0<=0.0) {
#ifdef DEBUG_MODE
			cout << "error 2a" << endl;
#endif
			return 0;
		}
		t2 = delta1/(2.0*sqrt(delta0*delta0*delta0));
		if (t2 <= -1.0)
			t2 = PI;
		else if (t2 >= 1.0)
			t2 = 0.0;
		else
			t2 = acos(t2);
		t2 = (-2.0/3.0)*p + (2.0*sqrt(delta0)*cos(t2/3.0))/(3.0*a);
		roundToZero(t2);
		if (t2 < 0.0) {
#ifdef DEBUG_MODE
			cout << "t2 = " << t2 << endl;
			cout << "error 2b" << endl;
#endif
			return 0;
		}
		S = 0.5*sqrt(t2);
	} else {
		// t1 >= 0.0
		double f1;
		if (delta0 == 0.0)
			f1 = delta1;
		else
			f1 = (delta1 + sqrt(t1))/2.0;
		double f2 = 1.0/3.0;
#ifdef DEBUG_MODE
		cout << "here 2" << endl;
		cout << "sqrt(t1) = " << sqrt(t1) << endl;
		cout << "f1 = " << f1 << endl;
		cout << "f2 = " << f2 << endl;
#endif
		if (f1 < 0.0)
			Q = -pow(-f1,f2);
		else
			Q = pow(f1,f2);
		if (Q == 0.0) {
			use_depressed_quartic_method = 1;
		} else {
			t2 = (-2.0/3.0)*p + (Q + delta0/Q) / (3.0*a);
#ifdef DEBUG_MODE
			cout << "Q = " << Q << endl;
			cout << "t2 = " << t2 << endl;
#endif
			roundToZero(t2);
			if (t2 < 0.0) {
#ifdef DEBUG_MODE
				cout << "error 3" << endl;
#endif
				return 0;
			}
			S = 0.5*sqrt(t2);
		}
	}

	if (use_depressed_quartic_method || (S <= SMALL_NUMBER && S >= -SMALL_NUMBER)) {
		// use depressed quartic and Ferrari's solution
#ifdef DEBUG_MODE
		cout << "using depressed quartic and Ferrari's solution" << endl;
#endif

		int success = toDepressedQuartic(in, out);
		if (!success) {
#ifdef DEBUG_MODE
			cout << "error 3a" << endl;
#endif
			return 0;
		}
		success = depressed_quartic_solver(out, out);
		if (!success) {
#ifdef DEBUG_MODE
			cout << "error 3b" << endl;
#endif
			return 0;
		}
		for (i=0; i<4; i++) {
			out[i] -= b/(4.0*a);
		}

		return 1;
	}

	q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d) / (8.0*a*a*a);
	t3 = (-4.0)*S*S - 2.0*p + q/S;
	t4 = (-4.0)*S*S - 2.0*p - q/S;
#ifdef DEBUG_MODE
	cout << "S = " << S << endl;
	cout << "p = " << p << endl;
	cout << "q = " << q << endl;
	cout << "t3 = " << t3 << endl;
	cout << "t4 = " << t4 << endl;
#endif
	if (t3 < 0 && t3 > -SMALL_NUMBER)
		t3 = 0.0;
	if (t4 < 0 && t4 > -SMALL_NUMBER)
		t4 = 0.0;
	if (t3 < 0.0 && t4 < 0.0) {
#ifdef DEBUG_MODE
		cout << "error 4" << endl;
#endif
		return 0;
	} else if (t3 >= 0.0 && t4 < 0.0) {
		out[0] = (-b)/(4.0*a) - S + 0.5*sqrt(t3);
		out[1] = (-b)/(4.0*a) - S - 0.5*sqrt(t3);
		out[2] = out[0];
		out[3] = out[1];
	} else if (t4 >= 0.0 && t3 < 0.0) {
		out[0] = (-b)/(4.0*a) + S + 0.5*sqrt(t4);
		out[1] = (-b)/(4.0*a) + S - 0.5*sqrt(t4);
		out[2] = out[0];
		out[3] = out[1];
	} else {
		out[0] = (-b)/(4.0*a) - S + 0.5*sqrt(t3);
		out[1] = (-b)/(4.0*a) - S - 0.5*sqrt(t3);
		out[2] = (-b)/(4.0*a) + S + 0.5*sqrt(t4);
		out[3] = (-b)/(4.0*a) + S - 0.5*sqrt(t4);
	}
	return 1;
}

double myRand() {
	int largest_number = 100;
	if (rand()%2==0)
		return rand()%largest_number+1.0;
	else
		return -(rand()%largest_number+1.0);
}

double myRandPos() {
	int largest_number = 100;
	return rand()%largest_number+1.0;
}

int myRandPos(int largest_number) {
	return rand()%largest_number+1.0;
}
// for testing the quartic solver for num_times times
void test_quartic_solver(int num_times) {

	srand(time(NULL));

	double r1,r2,r3,r4,c;
	int i,j,k;

	double* r = new double[4];
	double* in = new double[5];
	double* out = new double[4];
	double* rand1 = new double[5];
	double* rand2 = new double[5];
	int* corr = new int[4];
	int success;

	for (i=0; i<num_times; i++) {

		for (j=0; j<4; j++) {
			rand1[j] = myRand();
			rand2[j] = myRand();
		}
		rand1[4] = myRand();
		rand2[4] = myRand();

#ifdef DEBUG_MODE
		int z1[] = {0,0,2,2,1};
		int z2[] = {1,1,1,1,1};
		for (j=0; j<5; j++) {
			rand1[j] = z1[j];
			rand2[j] = z2[j];
		}
#endif

		for (j=0; j<4; j++)
			r[j] = rand1[j]/rand2[j];
		c = rand1[4]/rand2[4];

		r1 = r[0];r2 = r[1];r3 = r[2];r4 = r[3];

		in[0] = c;
		in[1] = -c*(r1+r2+r3+r4);
		in[2] = c*(r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4);
		in[3] = -c*(r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4);
		in[4] = c*(r1*r2*r3*r4);

		// cout << "A=" << in[0] << " B=" << in[1] << " C=" << in[2] << " D=" << in[3] << " E=" << in[4] << endl;
		success = 0;
		success = quartic_solver(in, out);
		// cout << "num_ans = " << num_ans << endl;

		cout << "R  : ";
		for (j=0; j<4; j++) {
			if (j>0)
				cout << ", ";
			cout << r[j];
		}
		cout << "; C : " << c;
		cout << endl;

		cout << "Rand1: ";
		for (j=0; j<=4; j++) {
			if (j>0)
				cout << ", ";
			cout << rand1[j];
		}
		cout << endl;

		cout << "Rand2: ";
		for (j=0; j<=4; j++) {
			if (j>0)
				cout << ", ";
			cout << rand2[j];
		}
		cout << endl;

		if (success) {

			// show the results
			cout << "Out: ";
			for (j=0; j<4; j++) {
				if (j>0)
					cout << ", ";
				cout << out[j];
			}
			cout << endl;

			// check the correctness
			corr[0]=corr[1]=corr[2]=corr[3]=0;
			for (j=0; j<4; j++) {
				// see whether the j-th answer matches with any r1, r2, r3, or r4
				for (k=0; k<4; k++) {
					if (corr[k]==0 && isTheSame(out[j], r[k])) {
						corr[k] = 1;
						break;
					}
				}
			}
			if (corr[0] && corr[1] && corr[2] && corr[3]) {
				cout << "All correct!" << endl;
			} else {
				cout << "NOT all correct!" << endl;
			}
		} else {
			cout << "No success" << endl;
		}


	}

	delete[] r;
	delete[] in;
	delete[] out;
	delete[] corr;
}

// for testing the cubic solver for num_times times
void test_cubic_solver(int num_times) {

	srand(time(NULL));

	double r1,r2,r3,c;
	int i,j,k;

	double* r = new double[3];
	double* in = new double[4];
	double* out = new double[3]; // only the real answers of the cubic equation
	int* corr = new int[3];
	int num_real_ans;
	int num_correct;
	int success;

	for (i=0; i<num_times; i++) {

		for (j=0; j<3; j++)
			r[j] = myRand()/myRand();
		c = myRand()/myRand();

		/*
		r[0] = 0.13;
		r[1] = 0.13;
		r[2] = 2.65625;
		c = -1.94118;
		 */

		r1 = r[0];r2 = r[1];r3 = r[2];

		in[0] = c;
		in[1] = -c*(r1+r2+r3);
		in[2] = c*(r1*r2 + r1*r3 + r2*r3);
		in[3] = -c*(r1*r2*r3);

		// cout << "A=" << in[0] << " B=" << in[1] << " C=" << in[2] << " D=" << in[3] << endl;
		success = 0; corr[0]=corr[1]=corr[2]= 0; num_correct=0;

		success = cubic_solver(in, out, num_real_ans);

		cout << "    R  : ";
		for (j=0; j<3; j++) {
			if (j>0)
				cout << ", ";
			cout << r[j];
		}
		cout << "; C : " << c;
		cout << endl;

		if (success) {

			// show the real result
			cout << "[" << num_real_ans << "] Out: ";
			for (j=0; j<num_real_ans; j++) {
				if (j > 0)
					cout << ", ";
				cout << out[j];
			}
			cout << endl;

			// check the correctness
			for (j=0; j<num_real_ans; j++) {
				for (k=0; k<3; k++) {
					if (corr[k]==0 && isTheSame(out[j], r[k])) {
						num_correct++;
						corr[k]=1;
						break;
					}
				}
			}

			if (num_correct == num_real_ans)
				cout << "All correct!" << endl;
			else
				cout << "NOT all correct!" << endl;

		} else {
			cout << "No success" << endl;
		}
	}
	delete[] r;
	delete[] in;
	delete[] out;
	delete[] corr;
}

// for testing the function eigenvalue_4_4 for num_times times
void test_eigenvalue_4_4(int num_times) {

#define MULTI_TIME 10

	srand(time(NULL));

	int i,j,k;

	double* in = new double[16];
	double* eigenVectors = new double[16];
	double* eigenValues = new double[4];
	double* s_matrix = new double[16];
	double* pi = new double[4];
	double* w = new double[6];

	double* eigenVectors2 = new double[16];
	double* transposeEigenVect2 = new double[16];
	double* in2 = new double[16];
	double* in_multi = new double[16];

	for (i=0; i<num_times; i++) {

		// simulate pi vector
		int orig_tot = 100;
		int total = 100;
		int t;
		for (j=0; j<3; j++) {
			t = rand()%(total-3+j)+1;
			// t = rand()%2+1;
			total-=t;
			pi[j] = (double)t/orig_tot;
		}
		pi[3]=(double)total/orig_tot;

#ifdef DEBUG_MODE
		double tpi[] = {0.25, 0.25, 0.25, 0.25};
		for (j=0; j<4; j++)
			pi[j] = tpi[j];
#endif

		// show the pi vector
		cout << "pi: ";
		for (j=0; j<4; j++) {
			if (j>0)
				cout << ", ";
			cout << pi[j];
		}
		cout << endl;

		/*
		for (j=0; j<4; j++)
			pi[j]=1.0;
		 */
		// simulate w matrix
		for (j=0; j<6; j++) {
			// w[j] = (double)(rand()%3+1)/3;
			// w[j] = myRandPos()/myRandPos();
			w[j] = (double) myRandPos(100) / myRandPos(10000);
		}

#ifdef DEBUG_MODE
		double tw[] = {1.02237, 1.00147, 1.02449, 0.999993, 0.999993, 0.999993};
		for (j=0; j<6; j++)
			w[j] = tw[j];
#endif

		// show the w matrix
		cout << "w: ";
		for (j=0; j<6; j++) {
			if (j>0)
				cout << ", ";
			cout << w[j];
		}
		cout << endl;

		// compute the S matrix
		computeSMatrix(s_matrix, w, pi);

		// compute the eigen-matrix
		computeEigenMatrix(in, w, pi, s_matrix);

		/*
	   #ifdef DEBUG_MODE
	   double sm[] = {-0.75,  0.25,  0.25,  0.25,
							0.25, -0.75,  0.25,  0.25,
							0.25,  0.25, -0.75,  0.25,
							0.25,  0.25,  0.25, -0.75};
	   for (j=0; j<16; j++)
		   in[j] = sm[j];
	   #endif*/

		// show the symmetric eigen-matrix
		cout << "Q:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << in[j*4+k];
			}
			cout << endl;
		}

		int success = 1;
		// double* ev;

		// compute the eigen values for symmetric 4x4 matrix
		// OUR METHOD
		// success = eigenvalue_4_4(in, eigenVectors, eigenValues);

		// JACOBI METHOD
		int it_max = 100;
		int n = 4;
		int it_num, rot_num;
		// cout << "Use Jacobi method" << endl;
		jacobi_eigenvalue ( n, in, it_max, eigenVectors, eigenValues, &it_num, &rot_num );
		if (!success) {
			cout << "No success" << endl;
		}/* else {
			// show the eigenVectors and eigenValues
			cout << "eigenVectors:" << endl;
			for (j=0; j<4; j++) {
				for (k=0; k<4; k++) {
					if (k>0)
						cout << ", ";
					cout << eigenVectors[j*4+k];
				}
				cout << endl;
			}
			cout << "eigenValues: " << endl;
			for (j=0; j<4; j++) {
				if (j>0)
					cout << ", ";
				cout << eigenValues[j];
			}
			cout << endl;
		}*/

		// verify the result

		// U
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				eigenVectors2[j*4+k] = eigenVectors[k*4+j];
			}
		}

		cout << "U:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << eigenVectors2[j*4+k];
			}
			cout << endl;
		}

		// U-1
		transpose(transposeEigenVect2, eigenVectors2, 4, 4);

		cout << "U-1:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << transposeEigenVect2[j*4+k];
			}
			cout << endl;
		}

		// in^MULTI_TIME
		memset(in_multi, 0, 16 * sizeof(double));
		for (j=0; j<4; j++)
			in_multi[j*4+j] = 1.0;
		for (j=0; j<MULTI_TIME; j++) {
			multiply(in_multi, in_multi, in, 4, 4, 4);
		}

		// the multiplication of U A^MULTI_TIME U-1
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				eigenVectors2[j*4+k] *= pow(eigenValues[k],MULTI_TIME);
			}
		}
		multiplyQuick(in2, eigenVectors2, transposeEigenVect2, 4, 4, 4);

		// show the matrix IN_MULTI
		cout << "Q^" << MULTI_TIME << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << in_multi[j*4+k];
			}
			cout << endl;
		}

		// show the matrix U A^MULTI_TIME U-1
		cout << "U A^" << MULTI_TIME << " U-1:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << in2[j*4+k];
			}
			cout << endl;
		}

		// verify the difference between eigen-matrix^MULTI_TIME and U A^MULTI_TIME U-1
		double diff = 0.0;
		// int isAllSame = 1;
		for (j=0; j<16; j++) {
			diff = fabsl(in_multi[j] - in2[j]) / in_multi[j];
			cout << "diff=" << diff << endl;
		}
		/*
			if (!isTheSameLoosely(in_multi[j], in2[j])) {
				cout << in_multi[j] << " , " << in2[j] << " diff: " << in_multi[j]-in2[j] << endl;
				isAllSame = 0;
				break;
			}
		 */
		/*
		if (isAllSame) {
			cout << "Yes, it is correct" << endl;
		} else {
			cout << "No, it is incorrect." << endl;
		}*/




		/*
		cout << "using the other method:" << endl;

		// by using the other method
		int it_max = 100;
		int n = 4;
		int it_num, rot_num;
		jacobi_eigenvalue ( n, in, it_max, eigenVectors, eigenValues, &it_num, &rot_num );

		// show the eigenVectors and eigenValues
		cout << "eigenVectors:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << eigenVectors[j*4+k];
			}
			cout << endl;
		}
		cout << "eigenValues: " << endl;
		for (j=0; j<4; j++) {
			if (j>0)
				cout << ", ";
			cout << eigenValues[j];
		}
		cout << endl;
		// verify the result
		allcorrect = 1;
		for (j=0; j<4; j++) {
			ev = &(eigenVectors[j*4]);
			for (k=0; k<4; k++) {
				v1[k] = 0.0;
				for (l=0; l<4; l++) {
					v1[k] += in[k*4+l]*ev[l];
				}
				v2[k] = eigenValues[j]*ev[k];
			}
			cout << "v1: ";
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << v1[k];
			}
			cout << endl;
			cout << "v2: ";
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << v2[k];
			}
			cout << endl;
			for (k=0; k<4; k++) {
				if (!isTheSame(v1[k],v2[k]))
					allcorrect = 0;
			}
		}
		if (allcorrect) {
			cout << "Yes, all correct" << endl;
		} else {
			cout << "No, not all correct!" << endl;
		}

		// U
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				eigenVectors2[j*4+k] = eigenVectors[k*4+j];
			}
		}

		// U-1
		transpose(transposeEigenVect2, eigenVectors2, 4, 4);

		// the multiplication of U U-1
		multiplyQuick(in2, eigenVectors2, transposeEigenVect2, 4, 4, 4);

		// show the matrix U U-1
		cout << "U U-1:" << endl;
		for (j=0; j<4; j++) {
			for (k=0; k<4; k++) {
				if (k>0)
					cout << ", ";
				cout << in2[j*4+k];
			}
			cout << endl;
		}
		 */

	}
}

int isTheSameLoosely(double x1, double x2) {
	if (x1 != 0.0 && x2 != 0.0)
		if (fabsl(x1) >= fabsl(x2))
			return (fabsl((x1 - x2)/x1) <= QUITE_SMALL_NUMBER);
		else
			return (fabsl((x1 - x2)/x2) <= QUITE_SMALL_NUMBER);
	else
		return (fabsl(x1 - x2) <= QUITE_SMALL_NUMBER);
}

int isTheSame(double x1, double x2) {
	if (x1 != 0.0 && x2!=0.0)
		if (fabsl(x1) >= fabsl(x2))
			return (fabsl((x1 - x2)/x1) <= SMALL_NUMBER);
		else
			return (fabsl((x1 - x2)/x2) <= SMALL_NUMBER);
	else
		return (fabsl(x1 - x2) <= SMALL_NUMBER);
}

void roundToZero(double& x) {
	if (x > -SMALL_NUMBER && x < SMALL_NUMBER) {
		x = 0.0;
	}
}

int isTheSame2(double x1, double x2) {
	if (x1 != 0.0 && x2 != 0.0)
		if (fabsl(x1) >= fabsl(x2))
			return (fabsl((x1 - x2)/x1) <= VERY_SMALL_NUMBER);
		else
			return (fabsl((x1 - x2)/x2) <= VERY_SMALL_NUMBER);
	else
		return (fabsl(x1 - x2) <= VERY_SMALL_NUMBER);
}

void roundToZero2(double& x) {
	if (x > -VERY_SMALL_NUMBER && x < VERY_SMALL_NUMBER) {
		x = 0.0;
	}
}


// print out the equations
void printEq(double* A, int num_terms, int num_eqs) {
	int i,j;
	for (i=0; i<num_eqs; i++) {
		for (j=0; j<num_terms; j++) {
			cout << "(" << A[i*num_terms+j] << ") ";
			if (j < num_terms-1)
				cout << "X" << j+1 << " + ";
		}
		cout << " = 0" << endl;
	}

}

// check whether the answers fit the equation
// return 1 if yes, 0 otherwise
int check_ans(double* A, int num_terms, double* ans) {
	double sum = 0.0;
	int i;
	for (i=0; i<num_terms-1; i++)
		sum += ans[i]*A[i];
	sum += A[num_terms-1];
	return isTheSame(sum, 0.0);
}

// Check the equation
// Return 1 if all terms are zeros except the last term (i.e. invalid equation)
// Return 2 if all terms are zeros (i.e. the equation has no effect)
// Return 0 otherwise
int check_equation(double* A, int num_terms) {
	int i;
	for (i=0; i<num_terms-1; i++) {
		if (!isTheSame2(A[i],0.0))
			return 0;
	}
	if (isTheSame2(A[num_terms-1], 0.0))
		return 2;
	else
		return 1;
}

// Return the status of two equations, where
// (1) not parallel
// (2) parallel, and no intersection point
// (3) parallel and overlap each other
// assuming two equations are normalized already
int check_two_equations(double* A, double* B, int num_terms) {
	int i, j;
	int isParallel = 1;
	for (i=0; i<num_terms-1; i++) {
		if (!isTheSame(A[i],0.0) || !isTheSame(B[i],0.0)) {
			if (!isTheSame(A[i],B[i])) {
				// cout << "Not the same: A[" << i << "]=" << A[i] << " B[" << i << "]=" << B[i] << endl;
				isParallel = 0;
				break;
			}
		}
	}
	if (isParallel) {
		if (isTheSame(A[num_terms-1],0.0) && isTheSame(B[num_terms-1],0.0))
			return 3; // parallel and overlap each other
		else if (isTheSame(A[num_terms-1],B[num_terms-1]))
			return 3; // parallel and overlap each other
		else
			return 2; // parallel, and no intersection point
	}
	return 1;
}

// normalize the equations (by dividing each term with the max term)
void normalizeEquation(double* eq, int num_terms) {

	int v;
	double max_coeff = 0.0;
	for (v=0; v<num_terms; v++) {
		if (max_coeff < fabsl(eq[v])) {
			max_coeff = fabsl(eq[v]);
		}
	}

	if (max_coeff > 0.0) {
		for (v=0; v<num_terms; v++) {
			eq[v] = eq[v] / max_coeff;
		}
	} else {
		return; // all terms are zeros
	}


	// round those terms to zero again if their new values are very closed to zero
	for (v=0; v<num_terms; v++)
		roundToZero(eq[v]);

	// consider the sign
	int sign=0;
	for (v=0; v<num_terms; v++) {
		if (eq[v] > 0.0) {
			sign = 1;
			break;
		} else if (eq[v] < 0.0) {
			sign = -1;
			break;
		}
	}
	if (sign==0)
		return;

	if (sign == -1) {
		for (v=0; v<num_terms; v++)
			eq[v] = -eq[v];
	}

}


/*
// Return the status of two equations, where
// (1) not parallel
// (2) parallel, and no intersection point
// (3) parallel and overlap each other
// assuming two equations are normalized already
int check_two_equations_loosely(double* A, double* B, int num_terms) {
	int i, j;
	int isParallel;
	for (i=0; i<num_terms-1; i++) {
		if (!isTheSame2(A[i],0.0) && !isTheSame2(B[i],0.0)) {
			isParallel = 1;
			for (j=0; j<num_terms-1; j++) {
				if (j!=i) {
					if (!isTheSame(A[j]/A[i],B[j]/B[i]))
						isParallel = 0;
				}
			}
			if (isParallel) {
				if (isTheSame(A[num_terms-1]/A[i],B[num_terms-1]/B[i]))
					return 3; // parallel and overlap each other
				else
					return 2; // parallel, and no intersection point
			}
		}
	}
	return 1;
}*/

// For solving one equation
// If there is more than one answer, return one of the possible answers.
// Return 1 when success; return 0 when failure
// Input: numTerms - number of terms
//        A[0] X0 + A[1] X1 + ... + A[numTerms-2] X(numTerms-2) + A[numTerms-1] = 0
// Output: out[0] = X0; out[1] = X1; ....  ; out[numTerms-2] = X(numTerms-2)
// accept_all_zeros : 1 - accept the results that X,Y = 0; 0 - otherwise
int solve_1_equation(double* A, double* out, int num_terms, int accept_all_zeros) {

	int i;
	int single_status = check_equation(A, num_terms);

	if (single_status==1) {
		// invalid equation
		return 0;
	} else if (single_status==2) {
		// all terms are zeros
		// just set all variables to be 1.0
		for (i=0; i<num_terms-1; i++)
			out[i] = 1.0;
		return 1;
	}

	// check how many non-zero variable terms
	int numNonZeroVarTerm = 0;
	for (i=0; i<num_terms-1; i++) {
		if (!isTheSame2(A[i],0.0))
			numNonZeroVarTerm++;
	}

	// if NumNonZeroVarTerm is 0, then single_status should be 1 or 2
	// and it would not reach here
	// thus, NumNonZeroVarTerm >= 1
	double sum = 0.0;
	int all_zeros = 1;
	for (i=0; i<num_terms-1; i++) {
		if (isTheSame2(A[i],0.0)) {
			out[i]=1.0; // can be any value
			all_zeros = 0;
		} else if (numNonZeroVarTerm > 1) {
			out[i]=1.0; // can be any value
			numNonZeroVarTerm--;
			sum += A[i]*out[i];
			all_zeros = 0;
		} else {
			out[i] = (sum + A[num_terms-1]) / (-A[i]);
			if (out[i] != 0.0)
				all_zeros = 0;
		}
	}

#ifdef DEBUG_MODE
	// list out the equations
	printEq(A, num_terms, 1);
	// list out the results
	cout << "temp result: ";
	for (i=0; i<num_terms-1; i++) {
		if (i>0)
			cout << ", ";
		cout << out[i];
	}
	cout << endl;
#endif

	if (all_zeros)
		return accept_all_zeros;
	else
		return 1;
}


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
int solve_equations_2_variables(int numEq, double* A, double* out, int accept_all_zeros) {

#ifdef DEBUG_MODE
	cout << "[solve_equations_2_variables] accept_all_zeros = " << accept_all_zeros << endl;
#endif

	int i,j;
	int num_terms = 3;

	/*
	if (loosely) {
		cout << "LOOSELY" << endl;
		for (i=0; i<numEq; i++) {
			for (j=0; j<3; j++) {
				if (fabsl(A[i*3+j]) < QUITE_SMALL_NUMBER)
					A[i*3+j] = 0.0;
			}
		}
		// print out the equations
		cout << "=====================" << endl;
		printEq(A, num_terms, numEq);
		cout << "=====================" << endl;
	}
	 */

	double* eq1;
	double* eq2;
	int *pair_status = new int[numEq*numEq];
	int *single_status = new int[numEq];
	// check whether there is no solution
	// (1) any of the formula is not valid
	for (i=0; i<numEq; i++) {
		eq1 = &(A[i*3]);
		single_status[i] = check_equation(eq1, num_terms);
		if (single_status[i]==1) { // the formula is not valid
			delete[] pair_status;
			delete[] single_status;
			return 0;
		}
	}
	// (2) see whether they are parallel and no intersection point
	for (i=0; i<numEq; i++) {
		eq1 = &(A[i*3]);
		for (j=i+1; j<numEq; j++) {
			eq2 = &(A[j*3]);
			pair_status[i*numEq+j] = check_two_equations(eq1, eq2, num_terms);
			if (pair_status[i*numEq+j]==2) { // pair_status=2 : parallel, and no intersection point
				delete[] pair_status;
				delete[] single_status;
				return 0;
			}
		}
	}

	// pick two non-zero equations which are not parallel to each other
	int found = 0;
	for (i=0; i<numEq && (!found); i++) {
		if (single_status[i]==0) { // non-zero equation
			eq1 = &(A[i*3]);
			for (j=i+1; j<numEq && (!found); j++) {
				if (single_status[j]==0) { // non-zero equation
					eq2 = &(A[j*3]);
					if (pair_status[i*numEq+j]==1) { // pair_status=1 : not parallel
						found = 1;
					}
				}
			}
		}
	}
	double a,b,c,d,e,f;
	if (found) {
		// there is a unique solution (or no solution)
		a=eq1[0];b=eq1[1];c=eq1[2];
		d=eq2[0];e=eq2[1];f=eq2[2];
		out[1]=(a*f-c*d)/(b*d-a*e);
		if (!isTheSame2(a,0.0)) {
			out[0]=(-c-b*(out[1]))/a;
		} else { // d!=0.0
			out[0]=(-f-e*(out[1]))/d;
		}
		if (!accept_all_zeros && out[0]==0 && out[1]==0) {
			delete[] pair_status;
			delete[] single_status;
			// out[0] = SUPER_SMALL_NUMBER;
			// out[1] = SUPER_SMALL_NUMBER;
			return 0;
		} else {
			// check whether the answer fits to the other equations
			for (i=0;i<numEq; i++) {
				eq1 = &(A[i*3]);
				if (!check_ans(eq1, 3, out)) {
					delete[] pair_status;
					delete[] single_status;
					printEq(A, 3, numEq);
					cout << "answer: x1=" << out[0] << "; x2=" << out[1] << endl;
					cout << "didn't pass the " << i << "-th equation" << endl;
					return 0;
				}
			}
			delete[] pair_status;
			delete[] single_status;
			return 1; // the answer fits all the equations!
		}
	}

	// pick any non-zero equation
	found = 0;
	for (i=0; i<numEq && (!found); i++) {
		if (single_status[i]==0) { // non-zero equation
			eq1 = &(A[i*3]);
			found = 1;
		}
	}
	if (found) {
		a = eq1[0]; b = eq1[1]; c=eq1[2];
		if (!isTheSame2(a,0.0)) {
			out[1]=1.0;out[0]=(-c-b)/a;
			delete[] pair_status;
			delete[] single_status;
			return 1;
		} else { // b!=0.0
			out[0]=1.0;out[1]=(-c-a)/b;
			delete[] pair_status;
			delete[] single_status;
			return 1;
		}
	}

	// All are zero equations
	out[0]=out[1]=1.0;
	delete[] pair_status;
	delete[] single_status;
	return 1;
}

// For testing the function "solve_equations_2_variables"
void testing_solve_equations_2_variables(int num_times) {

	srand(time(NULL));
	int i;
	int maxNumEq = 10;
	double* A = new double[maxNumEq*3];
	for (i=0; i<num_times; i++) {

		double x[2];
		x[0] = rand()%2; // myRand()/myRand();
		x[1] = rand()%2; // myRand()/myRand();
		int numEq = rand()%maxNumEq+1; // between 1 and maxNumEq
		int j,k;
		double constant_term;
		for (j=0; j<numEq; j++) {
			constant_term = 0.0;
			for (k=0; k<2; k++) {
				A[j*3+k] = rand()%2; // myRand()/myRand();
				constant_term += A[j*3+k]*x[k];
			}
			A[j*3+2] = -constant_term;
		}

		double out[2];
		// show the value of x[0] and x[1]
		cout << "Real values: x1 = " << x[0] << "; x2 = " << x[1] << endl;
		// solve the equations
		int success = solve_equations_2_variables(numEq, A, out, 1);
		if (!success) {
			cout << "No success" << endl;
		} else if (isTheSame2(x[0],out[0]) && isTheSame2(x[1],out[1])) {
			cout << "Yes correct!" << endl;
		} else {
			cout << "Result: x1 = " << out[0] << "; x2 = " << out[1] << endl;
			// checking whether the result still suit the equations
			int correct = 1;
			for (j=0; j<numEq; j++) {
				double* eq;
				eq = &(A[3*j]);
				if (!check_ans(eq, 3, out)) {
					correct = 0;
					break;
				}
			}
			if (correct)
				cout << "Yes, correct, although not the same as the original one" << endl;
			else
				cout << "No, incorrect!" << endl;
		}
	}
}


// For testing the function "solve_N_equations"
void testing_solve_N_equations(int num_times, int num_equations) {

	if (num_equations == 0)
		return;

	srand(time(NULL));
	int u,v,w;
	int num_terms = rand()%10+3; // at least 3
#ifdef DEBUG_MODE
	num_terms = 5;
#endif
	double* X = new double[num_terms-1];
	double* A = new double[num_equations * num_terms];
	double* out = new double[num_terms-1];
	double* curr_A;
	double sum;
	int success;

	for (v=0; v<num_times; v++) {

		// simulate the values of Xs
		for (u=0; u<num_terms-1; u++) {
			X[u] = rand()%10;
			// X[u] = myRand()/myRand();
		}

		// simulate the values of coefficients
		for (u=0; u<num_equations; u++) {
			curr_A = &(A[u*num_terms]);
			sum = 0.0;
			for (w=0; w<num_terms-1; w++) {
				curr_A[w] = rand()%10;
				// curr_A[w] = myRand()/myRand();
				sum += curr_A[w]*X[w];
			}
			curr_A[num_terms-1] = -sum;
		}
		// assign the values of coefficients if necessary

#ifdef DEBUG_MODE
		double values[] = {1,1,1,1,0,
				1,1,1,1,0,
				1,1,1,1,0,
				1,1,1,1,0,
				1.0/3.0,1.0/3.0,1.0/3.0,-1,0  };

		for (u=0; u<num_equations*num_terms; u++)
			A[u] = values[u];
		double Xvalues[] = {0,0,0,0,0};
		for (u=0; u<num_terms-1; u++)
			X[u] = Xvalues[u];
#endif

		// list out the equations
		for (u=0; u<num_equations; u++) {
			curr_A = &(A[u*num_terms]);
			for (w=0; w<num_terms-1; w++) {
				if (w>0) {
					cout << " + ";
				}
				cout << "(" << curr_A[w] << ") X" << w+1;
			}
			cout << " + (" << curr_A[num_terms-1] << ") = 0" << endl;
		}
		// list out the results
		for (u=0; u<num_terms-1; u++) {
			if (u>0)
				cout << "; ";
			cout << "X" << u+1 << " = " << X[u];
		}
		cout << endl;

		success = solve_N_equations(A,out,num_equations,num_terms,0);

		if (success) {
			cout << "Result: ";
			for (u=0; u<num_terms-1; u++) {
				if (u>0)
					cout << "; ";
				cout << "X" << u+1 << " = " << out[u];
			}
			cout << endl;
		}
		if (!success) {
			cout << "No success" << endl;
		} else {
			// check whether the results are the same as the X's
			int allCorrect = 1;
			int allZeros = 1;
			for (u=0; u<num_terms-1; u++) {
				if (!isTheSame2(out[u], X[u])) {
					allCorrect = 0;
					break;
				}
			}
			if (allCorrect) {
				cout << "Yes, all are correct!" << endl;
			} else {
				// first check whether the results are all zeros
				for (u=0; u<num_terms-1; u++) {
					if (out[u]!=0.0)
						allZeros = 0;
				}
				if (allZeros) {
					cout << "No, they are all zeros!" << endl;
				} else {
					allCorrect = 1;
					for (u=0; u<num_equations; u++) {
						curr_A = &(A[u*num_terms]);
						sum = 0.0;
						for (w=0; w<num_terms-1; w++) {
							sum += curr_A[w]*out[w];
						}
						if (!isTheSame2(curr_A[num_terms-1], -sum)) {
							allCorrect = 0;
							break;
						}
					}
					if (allCorrect) {
						cout << "Yes, all are correct, although some are not the same as the original values." << endl;
					} else {
						cout << "No, incorrect!" << endl;
						// list out the value of A
						cout << "A: ";
						for (u=0; u<num_equations*num_terms; u++) {
							if (u>0)
								cout << ", ";
							cout << A[u];
						}
						cout << endl;
						// list out the value of X
						cout << "X: ";
						for (u=0; u<num_terms-1; u++) {
							if (u>0)
								cout << ", ";
							cout << X[u];
						}
						cout << endl;
					}
				}
			}
		}
	}
}

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
int solve_N_equations(double* A, double* out, int num_equations, int num_terms, int accept_all_zeros) {

	if (num_terms < 3)
		return 0;

	int u,v;
	int status;
	double *eq1, *eq2;

	// normalize the equations (by dividing each term with the max term)
	for (u=0; u<num_equations; u++) {
		eq1 = &(A[num_terms*u]);
		normalizeEquation(eq1, num_terms);
	}

#ifdef DEBUG_MODE
	// print out the equations
	cout << "=====================" << endl;
	printEq(A, num_terms, num_equations);
	cout << "=====================" << endl;
#endif

	// if number of equations == 1, redirect to another function
	if (num_equations==1)
		return solve_1_equation(A, out, num_terms, accept_all_zeros);

	// if number of terms == 3, redirect to another function
	if (num_terms==3) {
		int status = solve_equations_2_variables(num_equations, A, out, accept_all_zeros);
#ifdef DEBUG_MODE
		cout << "Result of 2 equations: " << out[0] << ", " << out[1] << " status=" << status << endl;
#endif
		return status;
	}


	// check whether there is no solution
	// (1) check whether any of the equations is not valid
	int* single_status = new int[num_equations];
	for (u=0; u<num_equations; u++) {
		eq1 = &(A[num_terms*u]);
		single_status[u] = check_equation(eq1, num_terms);
		if (single_status[u]==1) {
			delete[] single_status;
			return 0;
		}
	}
	// (2) see whether they are parallel and no intersection point
	int* pair_status = new int[num_equations*num_equations];
#ifdef DEBUG_MODE
	cout << "pair_status: ";
#endif
	for (u=0; u<num_equations; u++) {
		eq1 = &(A[num_terms*u]);
		for (v=u+1; v<num_equations; v++) {
			eq2 = &(A[num_terms*v]);
			pair_status[u*num_equations+v] = check_two_equations(eq1, eq2, num_terms);
			if (pair_status[u*num_equations+v]==2) {
				delete[] single_status;
				delete[] pair_status;
				return 0; // two equations are parallel with no intersection point
			}
#ifdef DEBUG_MODE
			cout << pair_status[u*num_equations+v] << ", ";
#endif
		}
	}
#ifdef DEBUG_MODE
	cout << endl;
#endif
	// delete any redundant equation
	int* toDelete = new int[num_equations];
	memset(toDelete, 0, num_equations*sizeof(int));
	int new_num_equations = num_equations;
	for (u=0; u<num_equations; u++) {
		if (single_status[u]==2) {
			// all terms are zeros (i.e. the equation has no effect)
			toDelete[u] = 1;
			new_num_equations--;
		}
	}
	for (u=0; u<num_equations; u++) {
		if (!toDelete[u]) {
			for (v=u+1; v<num_equations; v++) {
				if (!toDelete[v] && pair_status[u*num_equations+v]==3) {
					// parallel and overlap each other
					toDelete[v] = 1;
					new_num_equations--;
				}
			}
		}
	}
	delete[] single_status;
	delete[] pair_status;

	if (new_num_equations==0) {
		// just set all answers as 1.0
		for (u=0; u<num_terms; u++)
			out[u] = 1.0;
		return 1;
	} else if (new_num_equations==1) {
		for (u=0; u<num_equations; u++) {
			if (toDelete[u]==0) {
				double* new_A = new double[num_terms];
				for (v=0; v<num_terms; v++) {
					new_A[v] = A[u*num_terms+v];
				}
				status = solve_1_equation(new_A, out, num_terms, accept_all_zeros);
				delete[] new_A;
				return status;
			}
		}
	}

	// new_num_equations >= 2
	double* new_A = new double[(new_num_equations-1)*(num_terms-1)];

	// we are going to remove the first undeleted equation
	// see which term should we remove too
	int rm_tm;
	int rm_eq;
	int found = 0;
	for (u=0; u<num_equations&&(!found); u++) {
		if (toDelete[u] == 0) {
			for (v=0; v<num_terms&&(!found); v++) {
				if (!isTheSame2(A[u*num_terms + v],0.0)) {
					rm_tm = v;
					rm_eq = u;
					found = 1;
				}
			}
		}
	}
	toDelete[rm_eq] = 1;
	int eq=0;
	for (u=0; u<num_equations; u++) {
		if (toDelete[u]==0) {
			int tm=0;
			for (v=0; v<num_terms; v++) {
				if (v!=rm_tm) {
					new_A[eq*(num_terms-1) + tm] = A[u*num_terms+v]*A[rm_eq*num_terms+rm_tm]-A[rm_eq*num_terms+v]*A[u*num_terms+rm_tm];
					tm++;
				}
			}
			eq++;
		}
	}
	if (accept_all_zeros || !isTheSame2(A[rm_eq*num_terms+(num_terms-1)], 0.0))
		status = solve_N_equations(new_A, out, new_num_equations-1, num_terms-1, 1);
	else
		status = solve_N_equations(new_A, out, new_num_equations-1, num_terms-1, accept_all_zeros);

	double sum = A[rm_eq*num_terms + (num_terms-1)];
	if (status) {
		for (u=num_terms-2; u>rm_tm; u--) {
			out[u] = out[u-1];
			sum += out[u] * A[rm_eq*num_terms + u];
		}
		out[rm_tm] = sum/(-A[rm_eq*num_terms+rm_tm]);
	}

	delete[] new_A;
	return status;
}


// reduce the variables in the equation
void reduce_variable(double* A, double* new_A, int num_terms, int remove_term) {
	int i,j;
	j=0;
	for (i=0; i<num_terms-1; i++) {
		if (remove_term!=i) {
			new_A[j++] = A[i];
		}
	}
	// for constant term
	new_A[j++] = A[num_terms-1]+A[remove_term];
}
