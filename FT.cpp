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
	//M&N必須是方陣
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			FFT(FreqReal, FreqImag, InputImage, M, N, j, i);
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

void FT::FFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	int N = h;
	int M = w;
	//N=M 且 必須為 2的幂次方
	double *real = new double[w];
	double *imag = new double[w];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			real[j] = pFreqReal[i][j];
			imag[j] = pFreqImag[i][j];
		}
		FFT1D(real,imag,h, w);
		for (int j = 0; j < M; j++)
		{
			 pFreqReal[i][j]= real[j];
			 pFreqImag[i][j]= imag[j];
		}
	}

	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < M; j++)
	//	{
	//		real[j] = pFreqReal[j][i];
	//		imag[j] = pFreqImag[j][i];
	//	}
	//	FFT1D(real, imag, h, w);
	//	for (int j = 0; j < M; j++)
	//	{
	//		pFreqReal[j][i] = real[j];
	//		pFreqImag[j][i] = imag[j];
	//	}
	//}

	delete[]real;
	delete[]imag;
	
}

void FT::FFT1D(double * pFreqReal, double * pFreqImag, int h, int w)
{
	//reverse bit
	bitreversal(pFreqReal, pFreqImag, h, w);
	//FFT 
	int m = log2(h);
	int c1 = -1, c2 = 0,step=1;
	
	for (int i = 0; i < m; i++)
	{
		int l = step;
		step <<= 1;
		int u1 = 1, u2 = 0;
		for (int j = 0; j <l; j++)
		{
			for (int k = j; k < h; k += step)
			{
				int place = k + l;
				int t1 = u1 * pFreqReal[place] - u2 * pFreqImag[place];
				int t2 = u1 * pFreqImag[place] + u2 * pFreqReal[place];
				pFreqReal[place] = pFreqReal[k] - t1;
				pFreqImag[place] = pFreqImag[k] - t2;
				pFreqReal[k] += t1;
				pFreqImag[k] += t2;

			}
			int z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = std::sqrt((1.0-c1)/2.0)*-1;
		c1 = std::sqrt((1.0 + c1) / 2.0);
	}

	for (int i = 0; i <h; i++)
	{
		pFreqReal[i] /= (double)h;
		pFreqImag[i] /= (double)h;
	}

}

void FT::iFFT1D(double * pFreqReal, double * pFreqImag, int h, int w)
{
	//reverse bit
	bitreversal(pFreqReal, pFreqImag, h, w);
	//FFT 
	int m = log2(h);
	int c1 = -1, c2 = 0, step = 1;

	for (int i = 0; i < m; i++)
	{
		int l = step;
		step <<= 1;
		int u1 = 1, u2 = 0;
		for (int j = 0; j < l; j++)
		{
			for (int k = j; k < h; k += step)
			{
				int place = k + l;
				int t1 = u1 * pFreqReal[place] - u2 * pFreqImag[place];
				int t2 = u1 * pFreqImag[place] + u2 * pFreqReal[place];
				pFreqReal[place] = pFreqReal[k] - t1;
				pFreqImag[place] = pFreqImag[k] - t2;
				pFreqReal[k] += t1;
				pFreqImag[k] += t2;

			}
			int z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = std::sqrt((1.0 - c1) / 2.0);
		c1= std::sqrt((1.0 + c1) / 2.0);
	}

	
}

void FT::bitreversal(double * pFreqReal, double * pFreqImag, int h, int w)
{
	for (int i = 1,j=0; i <h; i++)
	{
		for (int k = h >> 1; !((j ^= k)&k); k >>= 1)
		{
			if (i > j) {
				std::swap(pFreqReal[i],pFreqReal[j]);
				std::swap(pFreqImag[i], pFreqImag[j]);
			}
		}
	}
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i < M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
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
			InverseFFT(InverseReal, InverseImag, FreqReal, FreqImag, M, N, j, i);
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

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	int N = h;
	int M = w;
	//N=M 且 必須為 2的幂次方
	double *real = new double[w];
	double *imag = new double[w];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			real[j] = pFreqReal[i][j];
			imag[j] = pFreqImag[i][j];
		}
		iFFT1D(real, imag, h, w);
		for (int j = 0; j < M; j++)
		{
			pFreqReal[i][j] = real[j];
			pFreqImag[i][j] = imag[j];
		}
	}
	InverseReal = pFreqReal;
	InverseImag = pFreqImag;
	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < M; j++)
	//	{
	//		real[j] = pFreqReal[j][i];
	//		imag[j] = pFreqImag[j][i];
	//	}
	//	FFT1D(real, imag, h, w);
	//	for (int j = 0; j < M; j++)
	//	{
	//		pFreqReal[j][i] = real[j];
	//		pFreqImag[j][i] = imag[j];
	//	}
	//}

	delete[]real;
	delete[]imag;

}


void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{
}
