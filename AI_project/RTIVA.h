class realtime_IVA {

private:
	int Nfft;
	int nOverlap;
	int Nhalf;
	int Nfreq;
	double epsi;
	double beta;
	double eta;
	double gamma;

	double *x1;
	double *x2;
	double *ksi;
	double *iksi;
	double *win;
	double **Y;
	double **Y_buff;
	double **conjY;
	double **phi;
	double ***dd;
	double ***dW;
	double ***W;
	double *Sum_Y, *Avg_Y;
	double *hanningWin;

public:
	realtime_IVA();
	~realtime_IVA();
	void realtime_IVA::IVAprocessing(double **input, double **output);
};