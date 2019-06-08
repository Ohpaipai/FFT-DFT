#include "FT.h"

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
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
			DFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
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
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
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

	for (int i = 0; i < M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // �ť߸��W�v�}�C
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
			InverseDFT(InverseReal, InverseImag, FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
			//�s�U�ϳť߸���ƻP��Ƴ���
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
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
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
	for (int newcnt = 0; newcnt < h; newcnt++)
	{
		pFreq[newcnt] = new double[w]; // �ť߸��W�v�}�C
	}

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			FreqReal[i][j] = InputImage[i][j];
			FreqImag[i][j] = 0.0;
			pFreq[i][j] = 0;
		}
	}

	//row
	for (int i = 0; i < h; i++)
	{

		IORFFT1D_forwrow(FreqReal, FreqImag, w, i,false);
	}

	// column
	for (int i = 0; i < w; i++)
	{
		IORFFT1D_forcolumn(FreqReal, FreqImag, h, i,false);

	}

	//
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			pFreq[i][j] = sqrt(pow(FreqReal[j][i], 2.0) + pow(FreqImag[j][i], 2.0)) * w;
			OutputImage[i][j] = pFreq[i][j];
		}
	}
}

void FT::FFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
}

void FT::IORFFT1D_forwrow(double ** pFreqReal, double ** pFreqImag, int columnsize, const int row,bool INVERSE)
{

	std::complex<double>*x = new std::complex<double>[columnsize];
	for (int i = 0; i < columnsize; i++)
	{
		std::complex<double>tem(pFreqReal[row][i], pFreqImag[row][i]);
		x[i] = tem;
	}
	bitreversal(x, columnsize);

	for (int k = 2; k <= columnsize; k *= 2)
	{
		double w;
		if(INVERSE==true)
			 w= 2.0 * 3.1415926 / k;
		else
			 w = -2.0 * 3.1415926 / k;
		std::complex<double> dtheata(cos(w), sin(w));
		for (int j = 0; j < columnsize; j += k)
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

	if (INVERSE==true)
		for (int i = 0; i < columnsize; i++)
		{

			pFreqReal[row][i] = x[i].real() ;
			pFreqImag[row][i] = x[i].imag() ;
		}
	else
		for (int i = 0; i < columnsize; i++)
		{

			pFreqReal[row][i] = x[i].real() / columnsize;
			pFreqImag[row][i] = x[i].imag() / columnsize;

		}

}

void FT::IORFFT1D_forcolumn(double ** pFreqReal, double ** pFreqImag, int rowsize, const int column, bool INVERSE)
{

	std::complex<double>*x = new std::complex<double>[rowsize];
	for (int i = 0; i < rowsize; i++)
	{
		std::complex<double>tem(pFreqReal[i][column], pFreqImag[i][column]);
		x[i] = tem;
	}
	bitreversal(x, rowsize);

	for (int k = 2; k <= rowsize; k *= 2)
	{
		double w;
		if (INVERSE == true)
			w = 2.0 * 3.1415926 / k;
		else
			w = -2.0 * 3.1415926 / k;
		std::complex<double> dtheata(cos(w), sin(w));
		for (int j = 0; j < rowsize; j += k)
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
	if(INVERSE==true)
		for (int i = 0; i < rowsize; i++)
		{

			pFreqReal[i][column] = x[i].real() ;
			pFreqImag[i][column] = x[i].imag() ;

		}
	else
		for (int i = 0; i < rowsize; i++)
		{

			pFreqReal[i][column] = x[i].real() / rowsize;
			pFreqImag[i][column] = x[i].imag() / rowsize;

		}


}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	
	double** InverseReal = new double*[h];
	double** InverseImag = new double*[h];
	double** pFreq = new double*[w];

	for (int i = 0; i < h; i++)
	{
		InverseReal[i] = new double[w];
		InverseImag[i] = new double[w];
		pFreq[i] = new double[w]; // �ť߸��W�v�}�C
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
	for (int i = 0; i <h; i++)
	{
		
		IORFFT1D_forwrow(FreqReal, FreqImag, h, i, true);
	}

	//column
	for (int i = 0; i < w; i++)
	{
		
		IORFFT1D_forcolumn(FreqReal, FreqImag, w, i, true);
	}
	
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

void FT::bitreversal(std::complex<double>*x, int size)
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


void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{
}
