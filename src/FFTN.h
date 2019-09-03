/**
 * @file FFTN.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Definition of 2D and 1D Fast Fourier Fransformation (FFT) based on FFTW (see http://www.fftw.org)
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "fftw3.h"
#include "vector"
#include "math.h"
using namespace std;
#ifndef FFTN
#define FFTN

void FFT2d(vector< vector<double> > invector_r, fftw_complex* out_c);

void FFT2d(int rows, int cols, double* realin, fftw_complex* out_c);

void FFT1d(vector<double> invector_r, fftw_complex* out_c);

void FFT1d(int N, double* in, fftw_complex* out_c);

void IFFT1d(int N, fftw_complex* fftout, double* out);

void IFFT1d(int N, fftw_complex* fftout, vector<double>& reoutvector);

void IFFT2d(int rows, int cols, fftw_complex* fftout, vector<double>& reoutvector);

void IFFT2d(int rows, int cols, fftw_complex* fftout, double* out);

#endif
