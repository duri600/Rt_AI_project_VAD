#define n_file      1
#define sf_file      1
#define sm_file      1
#define str_len      50
#define SamplingFreq 16000
#define init_N_len 10
#define frame_len 1024
#define frame_shift 256
#define size_basis 513  //nFreq = 1024 / 2 + 1
#define blk_len_sep 1
#define R_x 100
#define R_d 100
#define R_a 20
#define m_a 100
#define P_len_l 10
#define P_len_k 60
#define preemph 0.9200
#define fftlen 1024
#define DCbin 5
#define nonzerofloor 1.0e-09
#define n_sparsity 5
#define max_itera 25
#define p_beta 1
#define p_beta_max 1000
#define p_alpha_d 0.85
#define p_alpha_eta 0.30
#define p_overlap_m_a 0.0100
#define p_overlapscale 0.5
#define tr_min 1
#define p_conv_eps 1e-3
#define tr_sec 384
#define P_blk_gap 7
#define P_alpha_p 0.4
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define cluster_buff 1
#define train_Exemplar 0 //training with SNMF
#define p_splice 0
#define p_conv_eps 1e-3
#define RANDOM_OFF 0      //RANDOM 0   ONES 1 


class SNMF {

private:
	int i, j, re, im;
	int fftlen2;
	double **h_w;
	double *in_mat, *in_w, *in_h;

	char sf_strs[sf_file][str_len] = {"female_train.wav"};
	
	char sm_strs[sm_file][str_len] = {"male_train.wav"};

	char n_strs[n_file][str_len] = {"noise_train.wav"};

	double** A;
	
	//buff_init
	double** g_Xm_hat;
	double** g_Dm_hat;
	double* g_Xm_tilde;
	double* g_x_hat;
	double* g_d_hat;
	double* g_x_tilde;
	int g_blk_cnt;
	double* g_lambda_dav;
	int g_update_switch;
	int g_l_mod_lswitch;
	double** g_r_blk;

	//hanning
	double *win_STFT;// win_STFT
	 //SNMF_test
	int N_length_bound;
	int g_W;
	double* output;
	double* d_frame;
	double* x_tilde;

	//bnmf_sep_event_RT
	double *y, *Ym, *Yp;

	double **Ym_mat;
	////////////////////////////////////////////////공유후 추가
	double **A_R_x;
	double **A_R_d;

	//1st array
	double  *wn, *div, *cost, *h_ind, *w_ind;
	double *sum1, *sum2;


	//2nd array
	double **matrix_mul_out1; //w*h -> m x n
	double **matrix_mul_out3; //v_lam * h_w -> m x rank
	double **lambda, **v, **v_lam, **div_sum;
	double **h, **dph, **dmh;//rank x n
	//n x rank
	double **h_w_a, **h_w_t;
	//m x rank
	double **w, **dpw, **dmw;
	//rank x m
	double **w_trans;
	////////////////////////////////////////////////
	//Training
	double **TF_mag;
	//Test
	double** test_init_w;
	double** test_init_h;
	double* eta, *G;
	double *Q;
	double *tmp_Dm_hat, *tmp_Xm_hat;


public:
	SNMF();
	~SNMF();
	void basis_stack(int total_len, double* f_stack, double* m_stack, double* n_stack);
	void SNMF_training(double*input, int input_len, int R, double **tr_out);
	void SNMF_test(double*input, double*output, double**basis_x, double**basis_n);
	void bnmf_sep_event_RT(double*input, double*d_frame, int l, double**basis_x, double**basis_n);
	void synth_ifft_buff(double *TF_mag, double *TF_phase, double *s_buff);
	void blk_sparse(double **X, double **D, int l, double *Q);
	void sparse_nmf(double** v, double** init_w, double** init_h, int w_update_ind, int h_update_ind, int row_v, int col_v, int rank);
	void max_matrix(double**input, int row, int col, double floor, double**output);
	void matrix_mul(double **a, double** b, int row_a, int col_a, int col_b, double** out);
	void rand_func(int row, int col, double** output);
	void stft_fft(double*s, int input_length, double** S_mag, int sz, int shift, int lenfft, int DC_bin, double* win, int ppreemph);
	void frame_splice(double**Feat, int row, int col, int splice);
};