#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "SNMF.h"
#include "sigproc.h"
#include "header.h"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

SNMF::SNMF()
{
	/* Signal Stack*/
	fftlen2 = fftlen / 2 + 1;

	g_Xm_hat = new double*[size_basis];
	g_Dm_hat = new double*[size_basis];
	g_Xm_tilde = new double[size_basis];

	g_lambda_dav = new double[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		g_Xm_hat[i] = new double[blk_len_sep];
		g_Dm_hat[i] = new double[blk_len_sep];
		g_Xm_tilde[i] = 0;
		g_lambda_dav[i] = 0;
		for (j = 0; j < blk_len_sep; j++)
		{
			g_Xm_hat[i][j] = 0;
			g_Dm_hat[i][j] = 0;
		}
	}

	g_x_hat = new double[frame_len];
	g_d_hat = new double[frame_len];
	g_x_tilde = new double[frame_len];
	for (i = 0; i < frame_len; i++)
	{
		g_x_hat[i] = 0;
		g_d_hat[i] = 0;
		g_x_tilde[i] = 0;
	}

	g_blk_cnt = 1;
	g_update_switch = 1;
	g_l_mod_lswitch = 0;

	g_r_blk = new double*[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		g_r_blk[i] = new double[P_len_l];
		for (j = 0; j < P_len_l; j++)
		{
			g_r_blk[i][j] = 0.0;
		}
	}

	win_STFT = new double[frame_len];
	for (i = 0; i < frame_len; i++)
	{
		win_STFT[i] = (1.0 - cos(2.0 * 3.14159265358979323846*(double)(i) / (frame_len))) * ((double)frame_shift / (frame_len / 2));
		win_STFT[i] = sqrt(win_STFT[i]);
	}

	//SNMF_test
	N_length_bound = 1;

	output = new double[frame_shift];//최종출력
	d_frame = new double[frame_len];//처리후 결과
	x_tilde = new double[frame_len];//처리결과 쌓는 버퍼
	for (i = 0; i < frame_len; i++)
	{
		x_tilde[i] = 0;
	}


	//bnmf_sep_event_RT
	y = new double[fftlen + 2];
	Ym = new double[fftlen2];
	Yp = new double[fftlen2];

	Ym_mat = new double*[fftlen2];
	for (i = 0; i < fftlen2; i++)
	{
		Ym_mat[i] = new double[1];
	}
	A_R_x = new double*[R_x];
	for (i = 0; i < R_x; i++)
	{
		A_R_x[i] = new double[1];
	}
	A_R_d = new double*[R_d];
	for (i = 0; i < R_d; i++)
	{
		A_R_d[i] = new double[1];
	}

	in_mat = new double[(tr_sec * SamplingFreq - frame_len) / frame_shift*fftlen2];
	in_w = new double[(R_x + R_d)*fftlen2];
	in_h = new double[(tr_sec * SamplingFreq - frame_len) / frame_shift*(R_x + R_d)];

	//Training
	TF_mag = new double *[fftlen2];
	for (i = 0; i < fftlen2; i++)
	{
		TF_mag[i] = new double[(tr_sec * SamplingFreq - frame_len) / frame_shift];
	}

	//Test
	test_init_w = new double*[size_basis];
	test_init_h = new double*[R_x + R_d];
	for (i = 0; i < size_basis; i++)
	{
		test_init_w[i] = new double[R_x + R_d];
	}
	for (i = 0; i < R_x + R_d; i++)
	{
		test_init_h[i] = new double[1];
	}
	Q = new double[size_basis];
	//eta
	eta = new double[size_basis];
	G = new double[size_basis];
	tmp_Dm_hat = new double[fftlen2];
	tmp_Xm_hat = new double[fftlen2];

}

SNMF::~SNMF()
{
	int i, j;

	for (i = 0; i < size_basis; i++)
	{

		delete[] g_Xm_hat[i];
		delete[] g_Dm_hat[i];
		delete[] g_r_blk[i];
	}

	delete[] g_Xm_hat;
	delete[] g_Dm_hat;
	delete[] g_Xm_tilde;
	delete[] g_lambda_dav;
	delete[] g_r_blk;
	delete[] win_STFT;
	delete[] g_x_hat;
	delete[] g_d_hat;
	delete[] g_x_tilde;
	delete[] output;
	delete[] d_frame;
	delete[] x_tilde;
	delete[] y;
	delete[] Ym;
	delete[] Yp;

	for (i = 0; i < fftlen2; i++)
	{
		delete[] Ym_mat[i];
	}
	delete[] Ym_mat;

	for (i = 0; i < R_x; i++)
	{
		delete[] A_R_x[i];
	}
	delete[] A_R_x;

	for (i = 0; i < R_d; i++)
	{
		delete[] A_R_d[i];
	}
	delete[] A_R_d;

	delete[] in_mat;
	delete[] in_w;
	delete[] in_h;

	//training
	for (i = 0; i < fftlen2; i++)
	{
		delete[] TF_mag[i];
	}
	delete[] TF_mag;

	//test

	delete[] eta;
	delete[] G;
	delete[] Q;
	delete[] tmp_Dm_hat;
	delete[] tmp_Xm_hat;
	for (i = 0; i < size_basis; i++)
	{
		delete[] test_init_w[i];
	}
	delete[] test_init_w;

	for (i = 0; i < R_x + R_d; i++)
	{
		delete[] test_init_h[i];
	}
	delete[] test_init_h;
}

void SNMF::basis_stack(int total_len, double* f_stack, double* m_stack, double* n_stack)
{
	int i, j, sample, random, rand_fin;
	int buffer_len = 3 * SamplingFreq;
	int signal_len_f;
	int signal_len_m;
	int signal_len_n;
	int iteration_f = total_len / (buffer_len*sf_file);
	int iteration_m = total_len / (buffer_len*sm_file);
	int iteration_n = total_len / (buffer_len*n_file);
	double *signal_m[sm_file];
	double *signal_f[sf_file];
	double *signal_n[n_file];
	char dir[10] = "./tr_db/";
	char data_dir[500];

    //저장된 6분 24초 데이터로 basis 만들기 위한 사전 과정
	std::cout << "reading " << sf_strs[0] << std::endl;
	strcpy(data_dir, dir);
	strcat(data_dir, sf_strs[0]);
	signal_f[0] = wavread(data_dir);
	signal_len_f = length(signal_f[0]);
	for ( sample = 0; sample < signal_len_f; sample++)
	{
		f_stack[sample] = signal_f[0][sample + 1] * 32768.0;
	}
	free(signal_f[0]);

	printf("Finish!!\n");

	std::cout << "reading " << sm_strs[0] << std::endl;
	strcpy(data_dir, dir);
	strcat(data_dir, sm_strs[0]);
	signal_m[0] = wavread(data_dir);
	signal_len_m = length(signal_m[0]);
	for (sample = 0; sample < signal_len_m; sample++)
	{
		m_stack[sample] = signal_m[0][sample + 1] * 32768.0;
	}
	free(signal_m[0]);

	printf("Finish!!\n");

	std::cout << "reading " << n_strs[0] << std::endl;
	strcpy(data_dir, dir);
	strcat(data_dir, n_strs[0]);
	signal_n[0] = wavread(data_dir);
	signal_len_n = length(signal_n[0]);
	for (sample = 0; sample < signal_len_n; sample++)
	{
		n_stack[sample] = signal_n[0][sample + 1] * 32768.0;
	}
	free(signal_n[0]);

	printf("Finish!!\n");

}


void SNMF::SNMF_training(double*input, int in_len, int R, double **tr_out)
{
	int *non_zero_id; //non zero index
	srand(time(NULL));
	int i, j, k, random;
	double **TF_mag_new;

	int frame_num = (in_len - frame_len) / frame_shift;

	double **B_DFT_init, **p_init_h, **wn_sq;
	B_DFT_init = new double *[fftlen / 2 + 1];
	wn_sq = new double *[fftlen / 2 + 1];
	for (i = 0; i < fftlen / 2 + 1; i++)
	{
		B_DFT_init[i] = new double[cluster_buff*R];
		wn_sq[i] = new double[cluster_buff*R];
	}

	stft_fft(input, in_len, TF_mag, frame_len, frame_shift, fftlen, DCbin, win_STFT, preemph);

	non_zero_id = new int[frame_num];
	double col_sum = 0.0;
	double *wn;
	int p_w_update_ind = 1, p_h_update_ind = 1;

	wn = new double[cluster_buff * R];

	int sum_nz = 0;

	for (i = 0; i < frame_num; i++)
	{
		col_sum = 0.0;
		for (k = 0; k < fftlen / 2 + 1; k++)
		{
			col_sum += TF_mag[k][i];
		}
		if (col_sum = !0.0)
		{
			non_zero_id[i] = 1;
			sum_nz++;
		}
		else
		{
			non_zero_id[i] = 0;
		}
	}

	int cnt = 0;
	TF_mag_new = new double*[fftlen / 2 + 1];
	for (i = 0; i < fftlen / 2 + 1; i++)
	{
		TF_mag_new[i] = new double[sum_nz];
	}

	for (i = 0; i < frame_num; i++)
	{
		if (non_zero_id[i] == 1)
		{
			for (j = 0; j < fftlen / 2 + 1; j++)
			{
				TF_mag_new[j][cnt] = TF_mag[j][i];
			}
			cnt++;
		}
	}
	delete[] non_zero_id;

	//frame splice 시작
	frame_splice(TF_mag_new, fftlen / 2 + 1, sum_nz, p_splice);

	for (i = 0; i < fftlen / 2 + 1; i++)
	{
		for (j = 0; j < sum_nz; j++)
		{
			TF_mag_new[i][j] = pow(TF_mag_new[i][j], 2) + nonzerofloor;
		}
	}
	for (j = 0; j < cluster_buff * R; j++)
	{

		for (i = 0; i < fftlen / 2 + 1; i++)
		{
			random = rand() % sum_nz;
			B_DFT_init[i][j] = TF_mag_new[i][random];
		}
	}

	//sparse nmf시작
	if (train_Exemplar == 0)
	{
		p_init_h = new double *[cluster_buff*R];
		for (i = 0; i < cluster_buff*R; i++)
		{
			p_init_h[i] = new double[sum_nz];
		}

		rand_func(R, sum_nz, p_init_h);
		sparse_nmf(TF_mag_new, B_DFT_init, p_init_h, p_w_update_ind, p_h_update_ind, fftlen / 2 + 1, sum_nz, R);

		for (j = 0; j < cluster_buff * R; j++)
		{
			for (i = 0; i < fftlen / 2 + 1; i++)
			{
				wn_sq[i][j] = sqrt(pow(B_DFT_init[i][j], 2));
			}
		}

		for (j = 0; j < cluster_buff * R; j++)
		{
			wn[j] = 0.0;
			for (i = 0; i < fftlen / 2 + 1; i++)
			{
				wn[j] += wn_sq[i][j];
			}
		}
		double flr = 1e-9;
		//normalization
		for (j = 0; j < R; j++)
		{
			for (i = 0; i < fftlen / 2 + 1; i++)
			{
				B_DFT_init[i][j] = B_DFT_init[i][j] / wn[j] + flr;
			}
		}
		for (i = 0; i < fftlen / 2 + 1; i++)
		{
			for (j = 0; j < R; j++)
			{
				tr_out[i][j] = B_DFT_init[i][j];
			}
		}
	}

	// delete!!
	for (i = 0; i < fftlen / 2 + 1; i++)
	{
		delete[] B_DFT_init[i];
		delete[] TF_mag_new[i];
		delete[] wn_sq[i];
	}
	delete[] B_DFT_init;
	delete[] TF_mag_new;
	delete[] wn_sq;
	delete[] wn;

	for (i = 0; i < cluster_buff*R; i++)
	{
		delete[] p_init_h[i];
	}
	delete[] p_init_h;

}


//input - frame_len, output - frame_len
void SNMF::SNMF_test(double*output_tmp, double*output_test, double**basis_x, double**basis_n)
{

	if (N_length_bound <= init_N_len)
	{
		g_W = 1;
	}
	else
	{
		g_W = 0;
	}

	bnmf_sep_event_RT(output_tmp, d_frame, N_length_bound, basis_x, basis_n);


	for (i = 0; i < frame_len - frame_shift; i++)
	{
		x_tilde[i] = x_tilde[frame_shift + i];
		x_tilde[i] += d_frame[i];
	}
	for (; i < frame_len; i++)
	{
		x_tilde[i] = d_frame[i];
	}

	for (i = 0; i < frame_shift; i++)
	{
		output_test[i] = x_tilde[i];
	}
	N_length_bound++;
}

void SNMF::bnmf_sep_event_RT(double*input, double*d_frame, int l, double**basis_x, double**basis_n)
{
	srand(time(NULL));

	int test_w_update_ind = 0;
	int test_h_update_ind = 1;

	//input buffer 처리
	for (i = 0; i < 1; i++)
	{
		y[i] = input[i] * win_STFT[i];
	}
	for (; i < frame_len; i++)
	{
		y[i] = (input[i] - preemph*input[i - 1]) * win_STFT[i];
	}
	for (; i < fftlen; i++)
	{
		y[i] = 0;
	}
	//hfft
	hfft3(y, fftlen, 1);

	//mag, phase
	for (i = 0; i < fftlen2; i++)
	{
		re = 2 * i;
		im = re + 1;
		Ym[i] = y[re] * y[re] + y[im] * y[im] + nonzerofloor;
		Yp[i] = atan2(y[im], y[re]);
	}
	for (i = 0; i < DCbin; i++)
	{
		Ym[i] = nonzerofloor;
	}

	//basis합치기
	for (i = 0; i < size_basis; i++)
	{
		for (j = 0; j < R_x; j++)
		{
			test_init_w[i][j] = basis_x[i][j];
		}
		for (; j < R_x + R_d; j++)
		{
			test_init_w[i][j] = basis_n[i][j - R_x];
		}
	}

	rand_func(R_x + R_d, 1, test_init_h);

	//for test
	for (i = 0; i < R_x + R_d; i++)
	{
		test_init_h[i][0] = 1;
	}

	for (i = 0; i < fftlen2; i++)
	{
		Ym_mat[i][0] = Ym[i];
	}

	//sparse_nmf
	sparse_nmf(Ym_mat, test_init_w, test_init_h, test_w_update_ind, test_h_update_ind, fftlen2, 1, R_x + R_d);

	for (i = 0; i < R_x; i++)
	{
		//A_R_x[i][0] = A[i][0];
		A_R_x[i][0] = test_init_h[i][0];
	}
	for (; i < R_x + R_d; i++)
	{
		//A_R_d[i - R_x][0] = A[i][0];
		A_R_d[i - R_x][0] = test_init_h[i][0];
	}

	matrix_mul(basis_x, A_R_x, size_basis, R_x, 1, g_Xm_hat);
	matrix_mul(basis_n, A_R_d, size_basis, R_d, 1, g_Dm_hat);

	//Calculate Block Sparsity

	blk_sparse(g_Xm_hat, g_Dm_hat, l, Q);

	//Enhancement Filter Construction
	double A_d_mag = 0.0, A_x_mag = 0.0, beta;

	if (l == 1)
	{
		for (i = 0; i < size_basis; i++)
		{
			g_lambda_dav[i] = Ym[i];
		}
	}

	for (i = 0; i < R_x; i++)
	{
		A_x_mag += (A_R_x[i][0] / R_x);
	}
	for (i = 0; i < R_d; i++)
	{
		A_d_mag += (A_R_d[i][0] / R_d);
	}
	beta = 20 * log10(A_d_mag / A_x_mag);
	beta *= p_beta;
	if (beta < p_beta)
	{
		beta = p_beta;
	}
	else if (beta >= p_beta_max)
	{
		beta = p_beta_max;
	}

	for (i = 0; i < size_basis; i++)
	{
		g_lambda_dav[i] = p_alpha_d*g_lambda_dav[i] + (1 - p_alpha_d)*g_Dm_hat[i][0] * beta;
	}

	for (i = 0; i < size_basis; i++)
	{
		eta[i] = MAX(0.0031, (p_alpha_eta*g_Xm_tilde[i] + (1 - p_alpha_eta)*g_Xm_hat[i][0] * Q[i]) / MAX(g_lambda_dav[i], nonzerofloor));
		G[i] = MIN(1, (eta[i] / (eta[i] + 1.0)));
	}

	if (l <= init_N_len)
	{
		for (i = 0; i < size_basis; i++)
		{
			G[i] = nonzerofloor;
			A_x_mag = nonzerofloor;
		}
	}
	//Xm_tilde
	for (i = 0; i < size_basis; i++)
	{
		g_Xm_tilde[i] = G[i] * Ym[i];
	}

	//--------block - wise inverse STFT-------------

	for (i = 0; i < fftlen2; i++)
	{
		tmp_Dm_hat[i] = g_Dm_hat[i][0];
		tmp_Xm_hat[i] = g_Xm_hat[i][0];
	}

	synth_ifft_buff(tmp_Xm_hat, Yp, g_x_hat);
	synth_ifft_buff(tmp_Dm_hat, Yp, g_d_hat);
	synth_ifft_buff(g_Xm_tilde, Yp, g_x_tilde);

	for (i = 0; i < frame_len; i++)
	{
		g_x_hat[i] = g_x_hat[i] * p_overlapscale;
		g_d_hat[i] = g_d_hat[i] * p_overlapscale;
		g_x_tilde[i] = g_x_tilde[i] * p_overlapscale;
		d_frame[i] = g_x_tilde[i];
	}
}

void SNMF::synth_ifft_buff(double *TF_input, double *TF_phase, double *s_buff)
{
	double* TF;
	TF = new double[fftlen + 2];

	for (i = 0; i < DCbin; i++)
	{
		TF_input[i] = 0;
	}
	for (; i < fftlen2; i++)
	{
		TF_input[i] = sqrt(TF_input[i]);
	}
	for (i = 0; i < fftlen2; i++)
	{
		re = i * 2;
		im = re + 1;
		TF[re] = TF_input[i] * cos(TF_phase[i]);
		TF[im] = TF_input[i] * sin(TF_phase[i]);
	}

	hfft3(TF, fftlen, -1);

	for (i = 0; i < 1; i++)
	{
		s_buff[i] = TF[i] * win_STFT[i];
	}
	for (; i < frame_len; i++)
	{
		s_buff[i] = TF[i] * win_STFT[i] + preemph*(s_buff[i - 1]);
	}

	delete[] TF;
}

void SNMF::blk_sparse(double **X, double **D, int l, double *Q)
{
	double *SNR_local;
	double MAX_SNR;
	double SNR_l1 = 0, SNR_l2 = 0, P_tmp, P_val;
	int n = P_len_k*P_len_l, P_len_k2 = floor(0.5 * P_len_k);
	SNR_local = new double[size_basis];

	for (i = 0; i < size_basis; i++)
	{
		SNR_local[i] = X[i][0] / MAX(D[i][0], nonzerofloor);
		if (i == 0)
		{
			MAX_SNR = SNR_local[0];
		}
		else if (MAX_SNR < SNR_local[i])
		{
			MAX_SNR = SNR_local[i];
		}
	}

	for (i = 0; i < size_basis; i++)
	{
		for (j = 0; j < P_len_l - 1; j++)
		{
			g_r_blk[i][j] = g_r_blk[i][j + 1];
		}
		SNR_local[i] = SNR_local[i] / MAX_SNR;
		g_r_blk[i][P_len_l - 1] = SNR_local[i];
		if (i < DCbin)
		{
			Q[i] = 0;
		}
		else
		{
			Q[i] = 0.1;
		}
	}

	int blk_gapN2 = (P_blk_gap - 1) / 2;
	if (l > P_len_l)
	{
		for (i = P_len_k2 + DCbin; i < size_basis - P_len_k2; i = i + P_blk_gap)
		{
			for (j = i - P_len_k2; j < i + P_len_k2; j++)
			{
				for (int k = 0; k < P_len_l; k++)
				{
					SNR_l1 = +g_r_blk[j][k];
					SNR_l2 = +g_r_blk[j][k] * g_r_blk[j][k];
				}
			}
			SNR_l2 = sqrt(SNR_l2);
			P_tmp = (sqrt(n) - SNR_l1 / SNR_l2) / (sqrt(n) - 1);
			P_val = P_alpha_p*Q[i - 2] + (1 - P_alpha_p) * P_tmp;
			for (int j = i - blk_gapN2 - 1; j < i; j++)
			{
				Q[j] = P_val;
			}
			for (int j = i - 1; j < i + blk_gapN2; j++)
			{
				Q[j] = P_val;
			}
		}
		for (int i = 0; i < P_len_k - 1; i++)
		{
			Q[i] = Q[P_len_k + DCbin - 1];
		}
	}

	for (int i = 0; i < DCbin; i++)
	{
		Q[i] = 0;
	}
	delete[] SNR_local;
}

//들어오는 프레임에 대해서 w와 h에 대한 정보를 update하는 함수이다.
//연산량이 많아 연산시간을 줄이기 위해 Armadillo library를 사용하였다.
void SNMF::sparse_nmf(double** input_mat, double** init_w, double** init_h, int w_update_ind, int h_update_ind, int m, int n, int rank)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			in_mat[i*m + j] = input_mat[j][i];
		}
	}
	for (i = 0; i < rank; i++)
	{
		for (j = 0; j < m; j++)
		{
			in_w[i*m + j] = init_w[j][i];
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < rank; j++)
		{
			in_h[i*rank + j] = init_h[j][i];
		}
	}

	arma::mat v = arma::mat(in_mat, m, n);
	arma::mat w = arma::mat(in_w, m, rank);
	arma::mat h = arma::mat(in_h, rank, n);
	arma::mat wn;
	arma::mat sum1, sum1_t;
	arma::mat sum3;
	arma::mat wnn_n;
	arma::mat w_h_mul;
	arma::mat lambda;
	arma::mat dph, dmh, dpw, dmw;
	arma::mat w_s;
	arma::mat div = arma::zeros<arma::mat>(1, max_itera);
	arma::mat cost = arma::zeros<arma::mat>(1, max_itera);
	arma::mat p_sparsity = 5 * arma::ones<arma::mat>(rank, n);
	int iter;
	double flr = 1e-9;
	double div_db, cost_db, last_cost = 1e+100;

	arma::mat pow_w = pow(w, 2);
	arma::mat sum_pow_w = sum(pow_w);
	wn = sqrt(sum_pow_w);
	w = w.each_row() / wn;
	h = h.each_col() % wn.t();
	w_h_mul = w*h;
	lambda = clamp(w_h_mul, flr, w_h_mul.max());
	v = clamp(v, flr, v.max());

	for (iter = 0; iter < max_itera; iter++)
	{
		// h에 대한 matrix가 update되는 부분
		if (h_update_ind > 0)
		{
			w_s = sum(w).t();
			dph = w_s + p_sparsity.each_col();
			dph = clamp(dph, flr, dph.max());
			dmh = w.t() * (v / lambda);
			h = h % dmh / dph;
			w_h_mul = w*h;
			lambda = clamp(w_h_mul, flr, w_h_mul.max());
		}
		// w에 대한 matrix가 update되는 부분
		if (w_update_ind > 0)
		{
			sum1 = sum(((v / lambda) * h.t()) % w); // 1* rank
			sum1_t = sum(h, 1).t();
			dpw = sum1_t + (sum1 % w.each_row()).each_row();
			dpw = clamp(dpw, flr, dpw.max());
			sum3 = sum(sum1_t % w.each_row());
			dmw = v / lambda * h.t() + sum3 % w.each_row();
			w = w % dmw / dpw;
			wnn_n = sqrt(sum(pow(w, 2)));
			w = w.each_row() / wnn_n;
			w_h_mul = w*h;
			lambda = clamp(w_h_mul, flr, w_h_mul.max());
		}

		div_db = sum(sum(v % log(v / lambda) - v + lambda));
		cost_db = div_db + sum(sum(p_sparsity % h));

		div(0, iter) = div_db;
		cost(0, iter) = cost_db;

		if (iter > 0 && p_conv_eps > 0)
		{
			double e;
			e = abs(cost_db - last_cost) / last_cost;
			if (e < p_conv_eps)
			{
				goto fin;
			}
		}
		last_cost = cost_db;
	}
fin:

	for (i = 0; i < rank; i++)
	{
		for (j = 0; j < m; j++)
		{
			init_w[j][i] = w(j, i);
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < rank; j++)
		{
			init_h[j][i] = h(j, i);
		}
	}

}

void SNMF::matrix_mul(double **a, double** b, int row_a, int col_a, int col_b, double** out)
{
	double sum = 0;
	int i, j, k;

	for (i = 0; i < row_a; i++)
	{
		for (k = 0; k < col_b; k++)
		{
			for (j = 0; j < col_a; j++)
			{
				sum = sum + a[i][j] * b[j][k];
			}
			out[i][k] = sum;
			sum = 0;
		}
	}
}

void SNMF::max_matrix(double**input, int row, int col, double flr, double**output)
{
	int i, j;
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			output[i][j] = MAX(flr, input[i][j]);
		}
	}
}

void SNMF::rand_func(int row, int col, double** output)
{

	int i, j, random;

	srand(time(NULL));
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
#if RANDOM_OFF==0
			output[i][j] = (double)rand() / (double)RAND_MAX;
#else
			output[i][j] = 1;
#endif
		}
	}
}

void SNMF::stft_fft(double*s, int input_length, double** S_mag, int sz, int shift, int lenfft, int DC_bin, double* win, int ppreemph)
{
	int size_crnt = 0;
	int i = 0;
	int j, k, re, im;

	double *s_frame_pad, *S_frame_mag;
	s_frame_pad = new double[lenfft + 2];
	S_frame_mag = new double[lenfft / 2 + 1];

	while (size_crnt < input_length - lenfft) // row랑 column 혼동 주의
	{
		for (j = 0; j < 1; j++)
		{
			s_frame_pad[j] = s[size_crnt + j] * win[j];
		}
		for (; j < sz; j++)
		{
			s_frame_pad[j] = (s[size_crnt + j] - ppreemph * s[size_crnt + j - 1]) * win[j];
		}
		for (; j < lenfft; j++)
		{
			s_frame_pad[j] = 0.0;
		}

		hfft3(s_frame_pad, lenfft, 1);

		for (j = 0; j < lenfft / 2 + 1; j++)
		{
			re = j + j;
			im = re + 1;
			S_frame_mag[j] = sqrt(pow(s_frame_pad[re], 2) + pow(s_frame_pad[im], 2));
		}

		for (k = 0; k < DC_bin; k++)
		{
			S_frame_mag[k] = 0.000001;
		}


		for (k = 0; k < lenfft / 2 + 1; k++)
		{
			S_mag[k][i] = S_frame_mag[k];
		}

		size_crnt = size_crnt + shift;
		i = i + 1;
	}
	delete[] s_frame_pad;
	delete[] S_frame_mag;
}

void SNMF::frame_splice(double**Feat, int row, int col, int splice)
{
	int i, j;
	int K_s;
	int s = 0;
	double *s_tmp;

	K_s = (2 * splice + 1)*row;
	s_tmp = new double[K_s];

	for (i = 0; i < col; i++)
	{
		if (i < 0)
		{
			for (j = 0; j < row; j++)
			{
				s_tmp[row*(splice + s) + j] = Feat[j][i + s];
				s_tmp[row*(splice - s) + j] = 0.0;
			}

		}
		else if (i > col - 1)
		{
			for (j = 0; j < row; j++)
			{
				s_tmp[row*(splice + s) + j] = 0.0;
				s_tmp[row*(splice - s) + j] = Feat[j][i - s];
			}
		}
		else
		{
			for (j = 0; j < row; j++)
			{
				s_tmp[row*(splice + s) + j] = Feat[j][i + s];
				s_tmp[row*(splice - s) + j] = Feat[j][i - s];
			}
		}
		for (j = 0; j < row; j++)
		{
			Feat[j][i] = s_tmp[j];
		}
	}
	delete[] s_tmp;
}