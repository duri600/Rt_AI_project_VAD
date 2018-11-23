#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define SWAP(a,b)	tempr=(a);(a)=(b);(b)=tempr
#define Nch			2
#define Nwin		1024
#define NBufferSize	256
#define NchASIOinput 2
#define NchASIOoutput 0
#define Mic 2

/* tools.cpp */
double *CreateVectorD(int size);
int VectorSizeD(double *v);
void ZeroVectorD(double *v);
double **Create2dVectorD(int row, int col);
void Vector2dSizeD(double **v, int *row, int *col);

/* Short-Time Fourier Transform */
double **stft(double *x, int Nfft, double *win, int Noverlap);
double *istft(double **X, int N, double *window, int Noverlap);