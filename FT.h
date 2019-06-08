#pragma once
#include <iostream>
#include<complex>
#include<algorithm>
#include<vector>
#define M_PI 3.1415926
class FT
{
private:
	
public:
	FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void FFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);
	void IORFFT1D_forwrow(double** pFreqReal, double** pFreqImag, int columnsize,const int row, bool INVERSE);
	void IORFFT1D_forcolumn(double** pFreqReal, double** pFreqImag, int rowsize,const int column, bool INVERSE);
	void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);
	void bitreversal(std::complex<double>*x, int size);
	void LowpassFilter(double** Real, double** Img, double** filter);
	void HighpassFilter(double** Real, double** Img, double** filter);

private:

};



