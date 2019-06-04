#pragma once
#include <iostream>
#define pi 3.1415926
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
	void FFT1D(double* pFreqReal, double* pFreqImag, int h, int w);
	void iFFT1D(double* pFreqReal, double* pFreqImag, int h, int w);
	void bitreversal(double* pFreqReal, double* pFreqImag, int h, int w);
	void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void LowpassFilter(double** Real, double** Img, double** filter);
	void HighpassFilter(double** Real, double** Img, double** filter);

private:

};



