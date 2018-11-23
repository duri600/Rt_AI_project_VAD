#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "RTIVA.h"
#include "header.h"
#include "sigproc.h"



realtime_IVA::realtime_IVA()
{
	int sample, i, j, re, im;
	int channel, ch, freq;
	Nfft = Nwin; //1024
	nOverlap = Nwin - NBufferSize;
	Nhalf = Nfft + 2; //1026
	Nfreq = Nfft / 2 + 1; // 513
	epsi = 0.000001;
	beta = 0.3;
	eta = 0.3;
	gamma = 0.001;

	W = new double**[Nch];
	dd = new double**[Nch];
	dW = new double**[Nch];

	for (channel = 0; channel < Nch; channel++)
	{
		W[channel] = new double*[Nch];
		dd[channel] = new double*[Nch];
		dW[channel] = new double*[Nch];
		for (ch = 0; ch < Nch; ch++)
		{
			W[channel][ch] = new double[Nfreq * 2 + 2];
			dd[channel][ch] = new double[Nfreq * 2 + 2];
			dW[channel][ch] = new double[Nfreq * 2 + 2];
		}
	}

	ksi = new double[Nfreq];
	for (i = 0; i < Nfreq; i++)
	{
		ksi[i] = 0.0;
	}


	iksi = new double[Nfreq];
	Sum_Y = new double[Nch];

	hanningWin = new double[Nwin];
	x1 = new double[Nwin * 2 + 2];
	x2 = new double[Nwin * 2 + 2];

	conjY = new double *[Nch];
	Y = new double *[Nch];
	Y_buff = new double *[Nch];
	phi = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		conjY[i] = new double[Nfreq * 2 + 2];
		Y[i] = new double[Nfreq * 2 + 2];
		phi[i] = new double[Nfreq * 2 + 2];
		Y_buff[i] = new double[Nwin];
	}
	for (i = 0; i < Nch; i++)
	{
		for (sample = 0; sample < Nfreq * 2 + 2; sample++)
		{
			Y[i][sample] = 0.0;
		}
		for (sample = 0; sample < Nwin; sample++)
		{
			Y_buff[i][sample] = 0.0;
		}
	}

	for (sample = 0; sample < Nwin; sample++)
	{
		//hanningWin[sample] = 0.5 * (1.0 - cos(2.0 * (double)M_PI*(double)(sample) / (Nwin))) * ((double)NBufferSize / (Nwin / 2));
		hanningWin[sample] = 0.5 * (1.0 - cos(2.0 * (double)M_PI*(double)(sample) / (Nwin)));
	}

	//W 정의
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < Nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				if (i == j)
				{
					W[i][j][re] = 1.0;
					W[i][j][im] = 0.0;
				}
				else
				{
					W[i][j][re] = 0.0;
					W[i][j][im] = 0.0;
				}
			}
		}
	}



}

realtime_IVA::~realtime_IVA()
{
	int ch, channel;
	for (channel = 0; channel < Nch; channel++)
	{
		for (ch = 0; ch < Nch; ch++)
		{
			delete[] W[channel][ch];
			delete[] dd[channel][ch];
			delete[] dW[channel][ch];
		}
		delete[] W[channel];
		delete[] dd[channel];
		delete[] dW[channel];
	}
	delete[] W;
	delete[] dd;
	delete[] dW;

	for (ch = 0; ch < Nch; ch++)
	{
		delete[] Y[ch];
		delete[] conjY[ch];
		delete[] phi[ch];
		delete[] Y_buff[ch];
	}
	delete[] Y;
	delete[] Y_buff;
	delete[] conjY;
	delete[] phi;
	delete[] ksi;
	delete[] iksi;
	delete[] Sum_Y;
	delete[] hanningWin;
	delete[] x1;
	delete[] x2;
}

void realtime_IVA::IVAprocessing(double **input, double **output)
{
	int i, j, sample, freq, ch;
	int re, im;

	//ksi 정의
	for (freq = 0; freq < Nfreq; freq++)
	{
		ksi[freq] += ksi[freq] + epsi;
	}

	//STFT
	for (sample = 0; sample < Nwin; sample++)
	{
		x1[sample] = hanningWin[sample] * input[0][sample];
		x2[sample] = hanningWin[sample] * input[1][sample];
	}

	hfft3(x1, Nfft, 1);
	hfft3(x2, Nfft, 1);


	for (i = 0; i < Nfreq; i++)
	{
		re = i + i;
		im = re + 1;
		Y[0][re] = W[0][0][re] * x1[re] - W[0][0][im] * x1[im] + W[0][1][re] * x2[re] - W[0][1][im] * x2[im];
		Y[0][im] = W[0][0][re] * x1[im] + W[0][0][im] * x1[re] + W[0][1][re] * x2[im] + W[0][1][im] * x2[re];
		Y[1][re] = W[1][0][re] * x1[re] - W[1][0][im] * x1[im] + W[1][1][re] * x2[re] - W[1][1][im] * x2[im];
		Y[1][im] = W[1][0][re] * x1[im] + W[1][0][im] * x1[re] + W[1][1][re] * x2[im] + W[1][1][im] * x2[re];
	}

	for (i = 0; i < Nch; i++)
	{
		Sum_Y[i] = 0.0;
		for (freq = 0; freq < Nfreq; freq++)
		{
			re = freq + freq;
			im = re + 1;
			Sum_Y[i] += pow(Y[i][re], 2) + pow(Y[i][im], 2);
		}
		Sum_Y[i] += epsi;
		Sum_Y[i] = sqrt(Sum_Y[i]); 
	}

	for (i = 0; i < Nch; i++)
	{
		for (freq = 0; freq < Nfreq; freq++)
		{
			re = freq + freq;
			im = re + 1;
			phi[i][re] = Y[i][re] / Sum_Y[i]; 
			phi[i][im] = Y[i][im] / Sum_Y[i];
		}
	}

	for (freq = 0; freq < Nfreq; freq++)
	{
		re = freq + freq;
		im = re + 1;
		ksi[freq] = ksi[freq] * beta + (1 - beta) * ((pow(x1[re], 2) + pow(x1[im], 2) + pow(x2[re], 2) + pow(x2[im], 2)) / 2);
		iksi[freq] = sqrt(1 / ksi[freq]);
	}

	for (i = 0; i < Nch; i++)
	{
		for (freq = 0; freq < Nfreq; freq++)
		{
			re = freq + freq;
			im = re + 1;
			conjY[i][re] = Y[i][re];
			conjY[i][im] = -Y[i][im];
		}
	}

	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < Nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				if (i == 0 && j == 0)
				{
					dd[i][j][re] = eta * iksi[freq] * (-gamma * (pow(Y[i][re], 2) - pow(Y[i][im], 2) + x1[im] * Y[i][im] - x1[re] * Y[i][re]));
					dd[i][j][im] = eta * iksi[freq] * (-gamma * (2 * Y[i][re] * Y[i][im] - x1[re] * Y[i][im] - x1[im] * Y[i][re]));
				}
				else if (i == 0 && j == 1) //i = 0, j = 1
				{
					dd[i][j][re] = eta * iksi[freq] * (-gamma * (Y[i][re] * Y[j][re] - x1[re] * Y[j][re] - Y[i][im] * Y[j][im] + x1[im] * Y[j][im]) - (phi[i][re] * conjY[j][re] - phi[i][im] * conjY[j][im]));
					dd[i][j][im] = eta * iksi[freq] * (-gamma * (Y[i][im] * Y[j][re] - x1[im] * Y[j][re] + Y[i][re] * Y[j][im] - x1[re] * Y[j][im]) - (phi[i][re] * conjY[j][im] + phi[i][im] * conjY[j][re]));
				}
				else if (i == 1 && j == 0) //i = 1, j = 0
				{
					dd[i][j][re] = eta * iksi[freq] * (-gamma * (Y[i][re] * Y[j][re] - Y[j][re] * x2[re] - Y[j][im] * Y[i][im] + Y[j][im] * x2[im]) - (phi[i][re] * conjY[j][re] - phi[i][im] * conjY[j][im]));
					dd[i][j][im] = eta * iksi[freq] * (-gamma * (Y[i][im] * Y[j][re] - Y[j][re] * x2[im] + Y[j][im] * Y[i][re] - Y[j][im] * x2[re]) - (phi[i][re] * conjY[j][im] + phi[i][im] * conjY[j][re]));
				}
				else
				{
					dd[i][j][re] = eta * iksi[freq] * (-gamma * (pow(Y[i][re], 2) - pow(Y[i][im], 2) + x2[im] * Y[i][im] - x2[re] * Y[i][re]));
					dd[i][j][im] = eta * iksi[freq] * (-gamma * (2 * Y[i][re] * Y[i][im] - x2[re] * Y[i][im] - x2[im] * Y[i][re]));
				}
			}
		}
	}

	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < Nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				dW[i][j][re] = (dd[i][0][re] * W[0][j][re] + dd[i][1][re] * W[1][j][re]) - (dd[i][0][im] * W[0][j][im] + dd[i][1][im] * W[1][j][im]);
				dW[i][j][im] = (dd[i][0][re] * W[0][j][im] + dd[i][1][re] * W[1][j][im]) + (dd[i][0][im] * W[0][j][re] + dd[i][1][im] * W[1][j][re]);
			}
		}
	}

	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < Nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				W[i][j][re] = W[i][j][re] + dW[i][j][re];
				W[i][j][im] = W[i][j][im] + dW[i][j][im];
			}
		}
	}

	for (i = 0; i < Nfreq; i++)
	{
		re = i + i;
		im = re + 1;
		Y[0][re] = W[0][0][re] * x1[re] - W[0][0][im] * x1[im] + W[0][1][re] * x2[re] - W[0][1][im] * x2[im];
		Y[0][im] = W[0][0][re] * x1[im] + W[0][0][im] * x1[re] + W[0][1][re] * x2[im] + W[0][1][im] * x2[re];
		Y[1][re] = W[1][0][re] * x1[re] - W[1][0][im] * x1[im] + W[1][1][re] * x2[re] - W[1][1][im] * x2[im];
		Y[1][im] = W[1][0][re] * x1[im] + W[1][0][im] * x1[re] + W[1][1][re] * x2[im] + W[1][1][im] * x2[re];
	}

	for (ch = 0; ch < Nch; ch++)
	{
		hfft3(Y[ch], Nfft, -1);
		for (i = 0; i < Nwin - NBufferSize; i++)
		{
			Y_buff[ch][i] = Y_buff[ch][NBufferSize + i];
			Y_buff[ch][i] += Y[ch][i];
		}
		for (; i < Nwin; i++)
		{
			Y_buff[ch][i] = Y[ch][i];
		}
	}

	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < NBufferSize; i++)
		{
			output[ch][i] = Y_buff[ch][i];
		}
	}
}