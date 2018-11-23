#define VAD_THRESHOLD  200   
#define FramesHangUp   2
#define FramesHangDn   15

class oldVAD {

private:

	double **X, **N;
	double *hanning_win;
	int Buffer_Size;
	int Nol;
	int count;
	double Nvar_mean;
	int record;
	int isSpeech;
	int EPDstate;
	int CNThangup;
	int CNThangdown;
	double PSY;
	int Mic_num;
	int Freq;
	short *input_short;
	int EPDcase;
	int length = 0;

public:
	oldVAD(int BuffSize, int micnum, int SamFreq);
	~oldVAD();
	void Estimate_Noise(double *input);
	int Calculate_VAD(double *input);
	int EPD_print();
	int power_test(double *input);


};