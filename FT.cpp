#include "FT.h"

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{

	double** pFreq = new double*[h];
	double** temReal = new double*[h];
	double** temImag = new double*[h];
	for (int newcnt = 0; newcnt < h; newcnt++)
	{
		pFreq[newcnt] = new double[w]; // 傅立葉頻率陣列
		temReal[newcnt] = new double[w];
		temImag[newcnt] = new double[w];
	}

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			pFreq[i][j] = 0;
		}
	}

	//row
	double *record_r = new double[w];
	double *record_i = new double[w];
	for (int i = 0; i < h; i++)
	{

		for (int j = 0; j < w; j++)
		{
			record_r[j] = InputImage[i][j];
			record_i[j] = 0.0;
		}
		FFT1D(record_r, record_i, w,false);
		for (int j = 0; j < w; j++)
		{
			temReal[i][j] = record_r[j];
			temImag[i][j] = record_i[j];
		}
	}
	delete[]record_i;
	delete[]record_r;

	//
	double *srecord_r = new double[w];
	double *srecord_i = new double[w];
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			srecord_r[j] = temReal[j][i];
			srecord_i[j] = temImag[j][i];
		}
		FFT1D(srecord_r, srecord_i, h, false);
		for (int j = 0; j < h; j++)
		{
			FreqReal[i][j] = srecord_r[j];
			FreqImag[i][j] = srecord_i[j];
		}
	}
	delete[]srecord_i;
	delete[]srecord_r;
	//
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], 2.0) + pow(FreqImag[i][j], 2.0));
			OutputImage[i][j] = pFreq[i][j];
		}
	}


	



}

void FT::FFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{

}

void FT::FFT1D(double * real, double * imag, int size, bool inverse)
{
	std::complex<double>*x = new std::complex<double>[size];
	for (int i = 0; i < size; i++)
	{
		std::complex<double>tem(real[i], imag[i]);
		x[i] = tem;
	}
	bitreversal(x, size);

	for (int k = 2; k <= size; k *= 2)
	{
		double w;
		if (inverse == true)
			w = 2.0 * 3.1415926 / k;
		else
			w = -2.0 * 3.1415926 / k;
		std::complex<double> dtheata(cos(w), sin(w));
		for (int j = 0; j < size; j += k)
		{
			std::complex<double> theata(1.0, 0.0);
			for (int i = j; i < (j + k / 2); i++)
			{
				std::complex<double> a = x[i],
					b = x[i + k / 2] * theata;
				x[i] = a + b;
				x[i + k / 2] = a - b;
				theata *= dtheata;
			}
		}
	}

	if (inverse == true)
		for (int i = 0; i < size; i++)
		{

			real[i] = x[i].real() / std::sqrt(size);
			imag[i] = x[i].imag() / std::sqrt(size);
		}
	else
		for (int i = 0; i < size; i++)
		{

			real[i] = x[i].real() / std::sqrt(size);
			imag[i] = x[i].imag() / std::sqrt(size);

		}



}

void FT::bitreversal(std::complex<double>* x, int size)
{
	for (int i = 1, j = 0; i < size; i++)
	{
		for (int k = size / 2; !((j ^= k)&k); k /= 2);
		if (i > j)
		{
			std::complex<double> temp = x[i];
			x[i] = x[j];
			x[j] = temp;
		}
	}
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	double** InverseReal = new double*[h];
	double** InverseImag = new double*[h];
	double** pFreq = new double*[w];
	double** temReal = new double*[h];
	double** temImag = new double*[h];

	for (int i = 0; i < h; i++)
	{
		InverseReal[i] = new double[w];
		InverseImag[i] = new double[w];
		pFreq[i] = new double[w]; // 傅立葉頻率陣列
		temReal[i] = new double[w];
		temImag[i] = new double[w];
	}

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//row
	double *record_r = new double[w];
	double *record_i = new double[w];
	for (int i = 0; i < h; i++)
	{

		for (int j = 0; j < w; j++)
		{
			record_r[j] = FreqReal[i][j];
			record_i[j] = FreqImag[i][j];
		}
		FFT1D(record_r, record_i, w, true);
		for (int j = 0; j < w; j++)
		{
			temReal[i][j] = record_r[j];
			temImag[i][j] = record_i[j];
		}
	}
	delete[]record_i;
	delete[]record_r;

	//
	double *srecord_r = new double[w];
	double *srecord_i = new double[w];
	for (int i = 0; i < w; i++)
	{
		for (int j = 0; j < h; j++)
		{
			srecord_r[j] = temReal[j][i];
			srecord_i[j] = temImag[j][i];
		}
		FFT1D(srecord_r, srecord_i, h, true);
		for (int j = 0; j < h; j++)
		{
			FreqReal[i][j] = srecord_r[j];
			FreqImag[i][j] = srecord_i[j];
		}
	}
	delete[]srecord_i;
	delete[]srecord_r;

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], 2.0) + pow(FreqImag[i][j], 2.0));
			OutputImage[i][j] = pFreq[i][j];
		}
	}

}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
}


void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{
}
