

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

void GetSpectrum1d(int N, double dx, fftw_complex* fftout, vector<double>& fvector, vector<double>& SpectrumVector);

#endif
