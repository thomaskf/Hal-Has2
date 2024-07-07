#ifndef __RAL_RAS__simd__
#define __RAL_RAS__simd__

#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <smmintrin.h>
#include <nmmintrin.h>
#endif

#ifdef __linux__     // Linux
#include <malloc.h>
#endif

inline double * malloc_decimals(int num_items)
{
	num_items = ((num_items+1)/2)*2;
	
#ifdef WIN32           // WIN32
	return (double *) _aligned_malloc(num_items*sizeof(double), 16);
#elif defined __linux__     // Linux
	return (double *) memalign(16, num_items*sizeof(double)); // aligned_alloc(16, num_items*sizeof(double));
#elif defined __MACH__      // Mac OS X
	return (double *) malloc(num_items*sizeof(double));
#else                       // other (use valloc for page-aligned memory)
	return (double *) valloc(num_items*sizeof(double));
#endif
}

inline void reset_decimals(double * pt, int num_items) {
	num_items = ((num_items+1)/2)*2;
	memset(pt, 0, sizeof(double)*num_items);
}

inline void free_decimals(double * pt)
{
#ifdef WIN32           // WIN32
	_aligned_free(pt);
#else                       // other (use valloc for page-aligned memory)
	free(pt);
#endif
}

inline void assign_decimals(double * a, int num_items, double value)
{
	__m128d* a_arr = (__m128d *) a;
	for (int i=0; i<(num_items+1)/2; i++) {
		a_arr[i] = _mm_set1_pd(value);
	}
}

inline void copy_decimals(double * target, double * source, int num_items) {
	memcpy(target, source, num_items*sizeof(double));
}

inline void multiply41(double *r, const double* a, const double* b){
	// multiple of a[4x4] and b[4x1]
	
	__m128d *r_arr;
	r_arr = (__m128d *) r;
	
	__m128d a_128, b_128;
	
	b_128 = _mm_set1_pd(b[0]);
	a_128 = _mm_set_pd(a[4],a[0]);
	r_arr[0] = _mm_mul_pd(a_128, b_128);
	a_128 = _mm_set_pd(a[12],a[8]);
	r_arr[1] = _mm_mul_pd(a_128, b_128);
	
	for (int i=1; i<4; i++) {
		b_128 = _mm_set1_pd(b[i]);
		a_128 = _mm_set_pd(a[4+i],a[i]);
		r_arr[0] = _mm_add_pd(_mm_mul_pd(a_128, b_128), r_arr[0]);
		a_128 = _mm_set_pd(a[12+i],a[8+i]);
		r_arr[1] = _mm_add_pd(_mm_mul_pd(a_128, b_128), r_arr[1]);
	}
	
}



inline void multiply41_t(double *r, const double* a, const double* b){
	// multiple of transpose of a[4x4] and b[4x1]
	
	__m128d *a_arr, *r_arr;
	a_arr = (__m128d *) a;
	r_arr = (__m128d *) r;
	
	__m128d b_128;
	
	b_128 = _mm_set1_pd(b[0]);
	r_arr[0] = _mm_mul_pd(a_arr[0], b_128);
	r_arr[1] = _mm_mul_pd(a_arr[1], b_128);
	
	for (int i=1; i<4; i++) {
		b_128 = _mm_set1_pd(b[i]);
		r_arr[0] = _mm_add_pd(_mm_mul_pd(a_arr[2*i], b_128), r_arr[0]);
		r_arr[1] = _mm_add_pd(_mm_mul_pd(a_arr[2*i+1], b_128), r_arr[1]);
	}
}

inline void multiply44_t(double *r, const double* a, const double* b){
	// multiple of transpose of a[4x4] and b[4x4]
	
	__m128d *b_arr, *r_arr;
	b_arr = (__m128d *) b;
	r_arr = (__m128d *) r;
	
	__m128d a_128;
	
	for (int i=0; i<4; i++) {
		a_128 = _mm_set1_pd(a[i]);
		r_arr[i*2] = _mm_mul_pd(a_128, b_arr[0]);
		r_arr[i*2+1] = _mm_mul_pd(a_128, b_arr[1]);
	
		for (int j=1; j<4; j++) {
			a_128 = _mm_set1_pd(a[j*4+i]);
			r_arr[i*2] = _mm_add_pd(_mm_mul_pd(a_128, b_arr[j*2]), r_arr[i*2]);
			r_arr[i*2+1] = _mm_add_pd(_mm_mul_pd(a_128, b_arr[j*2+1]), r_arr[i*2+1]);
		}
	
	}
}

inline void multiply_t(double *out, const double* a, const double* b, const int e, const int f) {
	// out[e x 1] = multiple of transpose of a[f x e] and b[f x 1]
	// dimensions e has to be divisible by 4
	
	__m128d *a_arr, *out_arr;
	a_arr = (__m128d *) a;
	out_arr = (__m128d *) out;
	
	__m128d b_128;
	
	b_128 = _mm_set1_pd(b[0]);
	int num_words = e/2;
	for (int i=0; i<num_words; i++)
		out_arr[i] = _mm_mul_pd(a_arr[i],b_128);
	
	for (int j=1; j<f; j++) {
		b_128 = _mm_set1_pd(b[j]);
		for (int i=0; i<num_words; i++) {
			out_arr[i] = _mm_add_pd(_mm_mul_pd(a_arr[num_words*j+i],b_128), out_arr[i]);
		}
	}
}


inline void round_decimals(double *a, int num_items, int decimalPlace) {

	__m128d *a_arr;
	__m128d x1, x2;
	
	a_arr = (__m128d *) a;
	x1 = _mm_set1_pd(pow(10.0,decimalPlace));
	x2 = _mm_set1_pd(1.0/pow(10.0,decimalPlace));
	for (int i=0; i<(num_items+1)/2; i++) {
		a_arr[i] = _mm_div_pd(_mm_round_pd(_mm_mul_pd(a_arr[i], x1),0),x1);
	}
	
}

inline void sqrt_decimals(double *target, double *src, int num_items) {
	
	__m128d *src_arr = (__m128d *) src;
	__m128d *target_arr = (__m128d *) target;
	for (int i=0; i<(num_items+1)/2; i++) {
		target_arr[i] = _mm_sqrt_pd(src_arr[i]);
	}
}

// target += src1 * src2 * v
inline void mul_add_decimals(double *target, double *src1, double *src2, double v, int num_items) {
	
	__m128d *src1_arr = (__m128d *) src1;
	__m128d *src2_arr = (__m128d *) src2;
	__m128d v_128 = _mm_set1_pd(v);
	__m128d *target_arr = (__m128d *) target;
	for (int i=0; i<(num_items+1)/2; i++) {
		target_arr[i] = _mm_add_pd(_mm_mul_pd(_mm_mul_pd(src1_arr[i], src2_arr[i]), v_128), target_arr[i]);
	}
}

// target = src1 * src2 * v
inline void mul_decimals(double *target, double *src1, double *src2, double v, int num_items) {
	
	__m128d *src1_arr = (__m128d *) src1;
	__m128d *src2_arr = (__m128d *) src2;
	__m128d v_128 = _mm_set1_pd(v);
	__m128d *target_arr = (__m128d *) target;
	for (int i=0; i<(num_items+1)/2; i++) {
		target_arr[i] = _mm_mul_pd(_mm_mul_pd(src1_arr[i], src2_arr[i]), v_128);
	}
	
}


inline void mul_decimals(double *a, int num_items, double factor) {

	__m128d *a_arr = (__m128d *) a;
	__m128d b = _mm_set1_pd(factor);
	for (int i=0; i<(num_items+1)/2; i++) {
		a_arr[i] = _mm_mul_pd(a_arr[i],b);
	}
	
}

// out = a * b
inline void mul_decimals(double* out, double *a, double b, int num_items) {

	__m128d *out_arr = (__m128d *) out;
	__m128d *a_arr = (__m128d *) a;
	__m128d b_128 = _mm_set1_pd(b);
	for (int i=0; i<(num_items+1)/2; i++) {
		out_arr[i] = _mm_mul_pd(a_arr[i],b_128);
	}
	
}

// out = a * b
inline void mul_decimals(double* out, double *a, double *b, int num_items) {

	__m128d *out_arr = (__m128d *) out;
	__m128d *a_arr = (__m128d *) a;
	__m128d *b_arr = (__m128d *) b;
	for (int i=0; i<(num_items+1)/2; i++) {
		out_arr[i] = _mm_mul_pd(a_arr[i],b_arr[i]);
	}
	
}

// out = a / b
inline void div_decimals(double* out, double *a, double *b, int num_items) {

	__m128d *out_arr = (__m128d *) out;
	__m128d *a_arr = (__m128d *) a;
	__m128d *b_arr = (__m128d *) b;
	for (int i=0; i<(num_items+1)/2; i++) {
		out_arr[i] = _mm_div_pd(a_arr[i],b_arr[i]);
	}
	
}

// out = sum(a * b)
inline void mul_sum_decimals(double* out, double *a, double *b, int num_items) {

	__m128d *out_arr = (__m128d *) out;
	__m128d *a_arr = (__m128d *) a;
	__m128d *b_arr = (__m128d *) b;
	int i;
	out_arr[0] = _mm_mul_pd(a_arr[0],b_arr[0]);
	for (i=1; i<num_items/2; i++) {
		out_arr[0] = _mm_add_pd(_mm_mul_pd(a_arr[i],b_arr[i]),out_arr[0]);
	}
	for (i=(num_items/2)*2; i<num_items; i++)
		out[0] += (a[i] * b[i]);
	
}


// out += a * b
inline void mul_add_decimals(double* out, double *a, double *b, int num_items) {

	__m128d *out_arr = (__m128d *) out;
	__m128d *a_arr = (__m128d *) a;
	__m128d *b_arr = (__m128d *) b;
	for (int i=0; i<(num_items+1)/2; i++) {
		out_arr[i] = _mm_add_pd(_mm_mul_pd(a_arr[i],b_arr[i]), out_arr[i]);
	}
	
}

// a[i] = a[i*4]+a[i*4+1]+a[i*4+2]+a[i*4+3] where 0 <= i < num_items
inline void special_op(double* a, int num_items) {
	
	__m128d * out_arr = (__m128d *) a;
	__m128d a_128;
	__m128d b_128;
	__m128d c_128;
	__m128d d_128;
	int i;
	for (i=0; i<num_items/2; i++) {
		a_128 = _mm_set_pd(a[i*8+4], a[i*8]);
		b_128 = _mm_set_pd(a[i*8+5], a[i*8+1]);
		c_128 = _mm_set_pd(a[i*8+6], a[i*8+2]);
		d_128 = _mm_set_pd(a[i*8+7], a[i*8+3]);
		out_arr[i] = _mm_add_pd(d_128,_mm_add_pd(c_128,_mm_add_pd(b_128,a_128)));
	}
	for (i=(num_items/2)*2; i<num_items; i++) {
		a[i] = a[i*4]+a[i*4+1]+a[i*4+2]+a[i*4+3];
	}

}

// if (a[i] > mx) a[i] = mx
// if (a[i] < mn) a[i] = mn
inline void max_min_decimals(double* a, double mx, double mn, int num_items) {

	__m128d * a_arr = (__m128d *) a;
	__m128d mx_128 = _mm_set1_pd(mx);
	__m128d mn_128 = _mm_set1_pd(mn);
	
	for (int i=0; i<(num_items+1)/2; i++) {
		a_arr[i] = _mm_min_pd(a_arr[i], mx_128);
		a_arr[i] = _mm_max_pd(a_arr[i], mn_128);
	}
}

inline void min_decimals(double* a, double mn, int num_items) {

	__m128d * a_arr = (__m128d *) a;
	__m128d mn_128 = _mm_set1_pd(mn);
	
	for (int i=0; i<(num_items+1)/2; i++) {
		a_arr[i] = _mm_max_pd(a_arr[i], mn_128);
	}
}

#endif // __RAL_RAS__simd__
