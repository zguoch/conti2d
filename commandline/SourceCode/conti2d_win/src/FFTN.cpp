
#include"FFTN.h"
void FFT2d(vector< vector<double> > invector_r, fftw_complex* out_c)
{

	fftw_complex *in;

	fftw_plan p;//,q;

	int rows = invector_r.size();
	int cols = invector_r[0].size();
	int N = rows*cols;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			in[j + i*cols][0] = invector_r[i][j];
			in[j + i*cols][1] = 0.0;
		}

	}

	p = fftw_plan_dft_2d(rows, cols, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/


	fftw_destroy_plan(p);
	//fftw_destroy_plan(q);
	fftw_free(in);
}
void FFT2d(int rows, int cols, double* realin, fftw_complex* out_c)
{

	fftw_complex *in;

	fftw_plan p;//,q;

	int N = rows*cols;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			in[j + i*cols][0] = realin[i*cols+j];
			in[j + i*cols][1] = 0.0;
		}

	}

	p = fftw_plan_dft_2d(rows, cols, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/


	fftw_destroy_plan(p);

	fftw_free(in);
}
void FFT1d(vector<double> invector_r, fftw_complex* out_c)
{

	fftw_complex *in;

	fftw_plan p;//,q;

	int N = invector_r.size();

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);

	for (int i = 0; i < N; i++)
	{
		in[i][0] = invector_r[i];
		in[i][1] = 0.0;

	}

	p = fftw_plan_dft_1d(N, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/

	fftw_destroy_plan(p);
	//fftw_destroy_plan(q);
	fftw_free(in);
}
void FFT1d(int N, double* in, fftw_complex* out_c)
{

	fftw_plan p;//,q;

	p = fftw_plan_dft_r2c_1d(N, in, out_c, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/

	fftw_destroy_plan(p);
	fftw_free(in);
}
void IFFT1d(int N, fftw_complex* fftout, double* out)
{
	fftw_plan p;//,q;
	p = fftw_plan_dft_c2r_1d(N, fftout, out, FFTW_ESTIMATE);
	fftw_execute(p);

	fftw_destroy_plan(p);
}
void IFFT1d(int N, fftw_complex* fftout, vector<double>& reoutvector)
{
	reoutvector.clear();

	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_1d(N, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (int i = 0; i < N; i++)
	{
		reoutvector.push_back(reout[i][0] / N);
	}

	fftw_destroy_plan(p);
	fftw_free(reout);
}
void IFFT2d(int rows, int cols, fftw_complex* fftout, vector<double>& reoutvector)
{
	reoutvector.clear();

	int N = rows*cols;
	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_2d(rows, cols, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			reoutvector.push_back(reout[j + i*cols][0] / N);
		}
	}

	fftw_destroy_plan(p);
	fftw_free(reout);
}
void IFFT2d(int rows, int cols, fftw_complex* fftout, double* out)
{
	int N = rows*cols;
	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_2d(rows, cols, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out[i*cols+j]=(reout[j + i*cols][0] / N);
		}
	}

	fftw_destroy_plan(p);
	fftw_free(reout);
}
void GetSpectrum1d(int N, double dx, fftw_complex* fftout, vector<double>& fvector, vector<double>& SpectrumVector)
{
	fvector.clear();
	SpectrumVector.clear();

	double FS = 1.0 / dx;
	double df = FS / (N - 1);
	int Nhalf = (int)floor(N / 2.0);
	double f = 0, spectrum = 0;
	for (int i = 0; i < Nhalf; i++)
	{
		f = i*df;
		spectrum = sqrt(pow(fftout[i][0], 2.0) + pow(fftout[i][1], 2.0));
		fvector.push_back(f);
		SpectrumVector.push_back(spectrum);
	}
	return;
}