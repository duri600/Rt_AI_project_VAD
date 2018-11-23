double *zeros(int len);
int length(double *v);

void fft(double *data, const int n, const int isign);
void hfft1(double *data, const int isign);
void hfft1_a(double data[], const int n, const int isign);
void hfft3(double data[], const int n, const int isign);
double *hanning(int nfft);
void free2(double **v);
void size(double **v, int *rows, int *cols);
double **zeros2(int rows, int cols);
/* File I/O */

void pcm2wav(char *pcm_name, char *wav_name, long Fs);
double *wavread(char *filename);

typedef struct /* WAV File Header */
{
	char	riff_id[4];
	long	riff_size;
	char	wave_id[4];
	char	format_id[4];
	long	format_size;
	short	format_type;
	short	ChannelNo;
	long	SamplePerSec;
	long	AvgBytesPerSec;
	short	BytesPerSample;
	short	BitsPerSample;
	char	data_id[4];
	long	data_size;
}	WAVEHEADER;