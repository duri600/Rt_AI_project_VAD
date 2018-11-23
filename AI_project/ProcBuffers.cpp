#define _CRT_SECURE_NO_WARNINGS
#define MAKE_FILE	1		//1 : make

#include <stdio.h>
#include "ProcBuffers.h"
#include "VAD.h"
#include "RTIVA.h"
#include "SNMF.h"
#include "sigproc.h"
#include <io.h>
#include <time.h>
#include <iostream>
#include "header.h"

oldVAD *my_VAD;
realtime_IVA *rtIVA;
SNMF *m_snmf;
SNMF *f_snmf;
#if MAKE_FILE == 1
FILE **nmf, **iva, **rec, **upsample;
double **nmf_test_buff;
double **out_buff;
short **rec_out;
short **iva_out;
short **nmf_out;
#endif
double **Input;
double **basis;
double **basis_f, **basis_m, **basis_fn, **basis_mn;
double **output, **output_temp, **input_temp;
double **InputBuffers, **xx_lp, **InitCond, *XX_LP, *XX, **x, *LPF;
double **vad_input_tmp;
double **test2_pad;
short **test2_48;//

				 //ASIO output
HWAVEOUT hWaveOut;
WAVEFORMATEX wfx = { WAVE_FORMAT_PCM, 1, 16000, 32000, 2, 16, 0 }; //{ WAVE_FORMAT_PCM, 1, 48000, 96000, 2, 16, 0 };
char *buffer1;
WAVEHDR header;


double lpf[256] = { 2.08798576826208e-06, 2.29954344986799e-05, -8.35216069314825e-06, -1.01421599408546e-05, 1.75825071131748e-05, 6.68506278162495e-07, -2.37466843099458e-05, 1.65292262827601e-05, 1.86603000155727e-05, -3.47276454619664e-05, 1.39791928852911e-06, 4.27828841291473e-05,	-3.24679407522093e-05, -2.99725427486013e-05, 6.15260248626211e-05, -6.94391496384368e-06, -7.03040741114150e-05, 5.82408009886648e-05, 4.39850560826189e-05, -0.000101231819760140, 1.84188639093049e-05,
0.000108128664125297, -9.77414224812796e-05, -5.99649810253215e-05, 0.000157527792910078, -3.91650080197067e-05, -0.000157914699984545, 0.000155735031049291, 7.64078744251690e-05, -0.000234514782604306, 7.36122297632366e-05, 0.000221009284717929, -0.000237857827130511, -9.08819603229543e-05, 0.000336524435543129, -0.000127344522646467, -0.000298213623651684, 0.000350655612780098, 9.97613501424823e-05, -0.000468059888805500, 0.000207292215334781, 0.000389557076800700, -0.000501506101607508,
-9.80403878092228e-05, 0.000633560527987403, -0.000321828720343297, -0.000494095184236539, 0.000698643273967565, 7.90873763046442e-05, -0.000837285503953346, 0.000480844420057879, 0.000609619959544147, -0.000951134221010459, -3.44598790385002e-05, 0.00108314551209364, -0.000696005652956184, -0.000732455819102594, 0.00126890531102301, -4.64362037536257e-05, -0.00137463837825748, 0.000980942207594894, 0.000857217253013655, -0.00166293374478135, 0.000176679007271611, 0.00171482336592332, -0.00135171031885412,
-0.000976555872744570, 0.00214562278599450, -0.000372365974724178, -0.00210647594676115, 0.00182754320491227, 0.00108084917623704, -0.00273156581117595, 0.000653525908926587, 0.00255245356503915, -0.00243222943518740, -0.00115775411275296, 0.00343894766908327, -0.00104561120424174, -0.00305643069554936, 0.00319651199968515, 0.00119146660897295, -0.00429207443129812, 0.00158219963360274, 0.00362422367640721, -0.00416242631985649, -0.00116137980775673, 0.00532598247979084, -0.00230995159742169,
-0.00426620546541984, 0.00539142505709083, 0.00103945666930688, -0.00659509904844149, 0.00329821462490733, 0.00500185874499057, -0.00698024820349760, -0.000784864830481446, 0.00819058918837528, -0.00465870264975467, -0.00586892292500419, 0.00909441684396274, 0.000332170044372207, -0.0102780677249779, 0.00658974395002571, 0.00694399559445460, -0.0120466066251864, 0.000437482118368211, 0.0131904224264900, -0.00948907025662504, -0.00839595749341352, 0.0165095446428109, -0.00177826835482988,
-0.0176999110017040, 0.0143004638937376, 0.0106573513166513, -0.0242397828308403, 0.00437616097914038, 0.0260668422700029, -0.0239555296472238, -0.0151991935796053, 0.0416959813296214, -0.0110150817261730, -0.0488064516064210, 0.0544275764049411, 0.0319437385822245, -0.127073765350224, 0.0615707177784359, 0.517443329375340, 0.517443329375340, 0.0615707177784359, -0.127073765350224, 0.0319437385822245, 0.0544275764049411, -0.0488064516064210, -0.0110150817261730, 0.0416959813296214,
-0.0151991935796053, -0.0239555296472238, 0.0260668422700029, 0.00437616097914038, -0.0242397828308403, 0.0106573513166513, 0.0143004638937376, -0.0176999110017040, -0.00177826835482988, 0.0165095446428109, -0.00839595749341352, -0.00948907025662504, 0.0131904224264900, 0.000437482118368211, -0.0120466066251864, 0.00694399559445460, 0.00658974395002571, -0.0102780677249779, 0.000332170044372207, 0.00909441684396274, -0.00586892292500419, -0.00465870264975467, 0.00819058918837528,
-0.000784864830481446, -0.00698024820349760, 0.00500185874499057, 0.00329821462490733, -0.00659509904844149, 0.00103945666930688, 0.00539142505709083, -0.00426620546541984, -0.00230995159742169, 0.00532598247979084, -0.00116137980775673, -0.00416242631985649, 0.00362422367640721, 0.00158219963360274, -0.00429207443129812, 0.00119146660897295, 0.00319651199968515, -0.00305643069554936, -0.00104561120424174, 0.00343894766908327, -0.00115775411275296, -0.00243222943518740, 0.00255245356503915,
0.000653525908926587, -0.00273156581117595, 0.00108084917623704, 0.00182754320491227, -0.00210647594676115, -0.000372365974724178, 0.00214562278599450, -0.000976555872744570, -0.00135171031885412, 0.00171482336592332, 0.000176679007271611, -0.00166293374478135, 0.000857217253013655, 0.000980942207594894, -0.00137463837825748, -4.64362037536257e-05, 0.00126890531102301, -0.000732455819102594, -0.000696005652956184, 0.00108314551209364, -3.44598790385002e-05, -0.000951134221010459, 0.000609619959544147,
0.000480844420057879, -0.000837285503953346, 7.90873763046442e-05, 0.000698643273967565, -0.000494095184236539, -0.000321828720343297, 0.000633560527987403, -9.80403878092228e-05, -0.000501506101607508, 0.000389557076800700, 0.000207292215334781, -0.000468059888805500, 9.97613501424823e-05, 0.000350655612780098, -0.000298213623651684, -0.000127344522646467, 0.000336524435543129, -9.08819603229543e-05, -0.000237857827130511, 0.000221009284717929, 7.36122297632366e-05, -0.000234514782604306, 7.64078744251690e-05,
0.000155735031049291, -0.000157914699984545, -3.91650080197067e-05, 0.000157527792910078, -5.99649810253215e-05, -9.77414224812796e-05, 0.000108128664125297, 1.84188639093049e-05, -0.000101231819760140, 4.39850560826189e-05, 5.82408009886648e-05, -7.03040741114150e-05, -6.94391496384368e-06, 6.15260248626211e-05, -2.99725427486013e-05, -3.24679407522093e-05, 4.27828841291473e-05, 1.39791928852911e-06, -3.47276454619664e-05, 1.86603000155727e-05, 1.65292262827601e-05, -2.37466843099458e-05, 6.68506278162495e-07,
1.75825071131748e-05, -1.01421599408546e-05, -8.35216069314825e-06, 2.29954344986799e-05, 2.08798576826208e-06 };

double up_filt[61] = { -1.43192112886976e-18, -0.00152272743066009, -0.00215148452403510, 3.70892844228060e-18, 0.00382749278365261, 0.00490897726088287, -6.98891002933101e-18, -0.00765507062444613, -0.00936349772317343, 1.12138964279313e-17, 0.0135792681040757, 0.0161482412402103, -1.61834608401010e-17, -0.0224018657807913, -0.0261873394322584, 2.15641756761223e-17, 0.0354229093905987, 0.0410720109334154, -2.69217414371986e-17, -0.0551954652449262, -0.0641574889289079, 3.17720179305444e-17, 0.0880231152256258,
0.104524016396168, -3.56435739509582e-17, -0.155484275560936, -0.198789716325378, 3.81419749332889e-17, 0.409656645427193, 0.825443168036801, 1.00060617355378, 0.825443168036801, 0.409656645427193, 3.81419749332889e-17, -0.198789716325378, -0.155484275560936, -3.56435739509582e-17, 0.104524016396168, 0.0880231152256258, 3.17720179305444e-17, -0.0641574889289079, -0.0551954652449262, -2.69217414371986e-17, 0.0410720109334154, 0.0354229093905987, 2.15641756761223e-17, -0.0261873394322584, -0.0224018657807913,
-1.61834608401010e-17, 0.0161482412402103, 0.0135792681040757, 1.12138964279313e-17, -0.00936349772317343, -0.00765507062444613, -6.98891002933101e-18, 0.00490897726088287, 0.00382749278365261, 3.70892844228060e-18, -0.00215148452403510, -0.00152272743066009, -1.43192112886976e-18
};


ProcBuffers::ProcBuffers()
{
	long NumDevice;
	long Selected;
	int ch, i, j;
	my_VAD = new oldVAD(4 * NBufferSize, Nch, SamplingFreq);
	rtIVA = new realtime_IVA();
	m_snmf = new SNMF();
	f_snmf = new SNMF();
	int total_stack = tr_sec * SamplingFreq;

	Input = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		Input[ch] = new double[4 * NBufferSize];
	}
	basis = new double*[size_basis];
	basis_f = new double*[size_basis];
	basis_m = new double*[size_basis];
	basis_fn = new double*[size_basis];
	basis_mn = new double*[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		basis[i] = new double[2 * R_x + 3 * R_d / 2];
		basis_f[i] = new double[R_x];
		basis_m[i] = new double[R_x];
		basis_fn[i] = new double[R_d];
		basis_mn[i] = new double[R_d];
	}
	output = new double*[Nch];
	output_temp = new double*[Nch];
	input_temp = new double*[Nch];
#if MAKE_FILE == 1
	nmf_test_buff = new double*[Nch];
	out_buff = new double*[Nch];
	rec_out = new short*[Nch];
	iva_out = new short*[Nch];
	nmf_out = new short*[Nch];
#endif
	for (ch = 0; ch < Nch; ch++)
	{
		output[ch] = new double[NBufferSize];
		output_temp[ch] = new double[Nwin];
		input_temp[ch] = new double[Nwin];
#if MAKE_FILE == 1
		nmf_test_buff[ch] = new double[frame_shift];
		out_buff[ch] = new double[NBufferSize];
		rec_out[ch] = new short[frame_shift];
		iva_out[ch] = new short[frame_shift];
		nmf_out[ch] = new short[frame_shift];
#endif
	}
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < Nwin; i++)
		{
			output_temp[ch][i] = 0.0;
			input_temp[ch][i] = 0.0;
		}
	}

#if MAKE_FILE == 1
	char file_name[2][500];
	nmf = new FILE*[Nch];
	iva = new FILE*[Nch];
	rec = new FILE*[Nch];
	upsample = new FILE*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
		nmf[ch] = fopen(file_name[0], "wb");
		sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
		iva[ch] = fopen(file_name[0], "wb");
		sprintf(file_name[0], ".\\output\\REC_ch%d.pcm", ch + 1);
		rec[ch] = fopen(file_name[0], "wb");
		sprintf(file_name[0], ".\\output\\upsample_ch%d.pcm", ch + 1);
		upsample[ch] = fopen(file_name[0], "wb");
	}
#endif
	// ASIO
	InputBuffers = new double *[NchASIOinput];
	InitCond = new double *[NchASIOinput];
	x = new double *[NchASIOinput];
	for (i = 0; i < NchASIOinput; i++)
	{
		InputBuffers[i] = new double[NBufferSize];
		InitCond[i] = new double[NBufferSize];
		x[i] = new double[NBufferSize];
		for (j = 0; j < NBufferSize; j++)
		{
			InputBuffers[i][j] = 0;
			InitCond[i][j] = 0;
			x[i][j] = 0;
		}
	}
	xx_lp = new double *[NchASIOinput];
	for (i = 0; i < NchASIOinput; i++)
	{
		xx_lp[i] = new double[3 * NBufferSize];
		for (j = 0; j < 3 * NBufferSize; j++)
			xx_lp[i][j] = 0;
	}

	XX = new double[NBufferSize * 2 + 2];
	XX_LP = new double[NBufferSize * 2 + 2];
	LPF = new double[NBufferSize * 2 + 2];
	for (i = 0; i < NBufferSize * 2 + 2; i++)
	{
		XX[i] = 0;
		XX_LP[i] = 0;
		LPF[i] = 0;
	}
	for (i = 0; i < NBufferSize; i++)
		LPF[i] = lpf[i];
	hfft1_a(LPF, NBufferSize * 2, 1);

	vad_input_tmp = zeros2(Nch, Nwin);
	buffer1 = new char[NBufferSize * 2];

	test2_48 = new short *[Nch];//
	test2_pad = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		test2_48[ch] = new short[NBufferSize * 3];//
		test2_pad[ch] = new double[3 * NBufferSize + 60];
	}
	for (ch = 0; ch < Nch; ch++) //NMF input
	{
		for (i = 0; i < 3 * NBufferSize + 60; i++)
		{
			test2_pad[ch][i] = 0.0;
		}
	}
}

ProcBuffers::~ProcBuffers()
{
	int ch, i;
	delete my_VAD;
	delete rtIVA;
	delete m_snmf;
	delete f_snmf;

	for (i = 0; i < NchASIOinput; i++)
	{
		delete[] xx_lp[i];
		delete[] x[i];
		delete[] InitCond[i];
	}
	delete[] xx_lp;
	delete[] x;
	delete[] InitCond;
	delete[] LPF;
	delete[] XX_LP;
	delete[] XX;

	for (ch = 0; ch < Nch; ch++)
	{
		delete[] Input[ch];
		delete[] output[ch];
		delete[] output_temp[ch];
		delete[] input_temp[ch];
#if MAKE_FILE == 1
		delete[] nmf_test_buff[ch];
		delete[] out_buff[ch];
		delete[] rec_out[ch];
		delete[] iva_out[ch];
		delete[] nmf_out[ch];
#endif
	}
#if MAKE_FILE == 1
	delete[] nmf_test_buff;
	delete[] out_buff;
	delete[] rec_out;
	delete[] iva_out;
	delete[] nmf_out;
#endif
	delete[] Input;
	delete[] output;
	delete[] output_temp;
	delete[] input_temp;

#if MAKE_FILE == 1
	char file_name[2][500];
	for (ch = 0; ch < Nch; ch++)
	{
		fclose(rec[ch]);
		sprintf(file_name[0], ".\\output\\REC_ch%d.pcm", ch + 1);
		sprintf(file_name[1], ".\\output\\REC_ch%d.wav", ch + 1);
		pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
		remove(file_name[0]);
		char iva1_pcm[500] = ".\\output\\IVA_ch1.pcm";
		char iva2_pcm[500] = ".\\output\\IVA_ch2.pcm";
		if ((_access(iva1_pcm, 0) == 0) || (_access(iva2_pcm,0)==0))
		{
			fclose(iva[ch]);
			sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
			sprintf(file_name[1], ".\\output\\IVA_ch%d_fin.wav", ch + 1);
			pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
			remove(file_name[0]);
			fclose(nmf[ch]);
			sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
			sprintf(file_name[1], ".\\output\\NMF_ch%d_fin.wav", ch + 1);
			pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
			remove(file_name[0]);
			fclose(upsample[ch]);
			sprintf(file_name[0], ".\\output\\upsample_ch%d.pcm", ch + 1);
			sprintf(file_name[1], ".\\output\\upsample_ch%d_fin.wav", ch + 1);
			pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
			remove(file_name[0]);
		}
	}
#endif
	for (i = 0; i < size_basis; i++)
	{
		delete[] basis[i];
		delete[] basis_f[i];
		delete[] basis_m[i];
		delete[] basis_fn[i];
		delete[] basis_mn[i];
	}
	delete[] basis;
	delete[] basis_f;
	delete[] basis_m;
	delete[] basis_fn;
	delete[] basis_mn;
#if MAKE_FILE == 1
	delete[] nmf;
	delete[] iva;
	delete[] rec;
#endif
	free2(vad_input_tmp);
	delete[] buffer1;

	for (ch = 0; ch < Nch; ch++)
	{
		delete[] test2_48[ch];
		delete[] test2_pad[ch];
	}
	delete[] test2_48;
	delete[] test2_pad;
}

int ProcBuffers::Process(double**input, double**test2_48k)
{
	FILE *nFile;

	static int Process_count = 0;
	static int VAD_count = 0;
	static int rec_num = 0;
	static int BuffCnt = 0, isNew16k = 0;
	int EPDcase, ch, i, j;
	int total_stack = tr_sec * SamplingFreq;

	//처음 process를 시작하면서 train후 저장된 basis를 불러온다.
	if (Process_count == 0)
	{
		nFile = fopen("basis_mat.txt", "rt");
		if (nFile != NULL)
		{
			for (i = 0; i < size_basis; i++)
			{
				for (j = 0; j < 2 * R_x + 3 * R_d / 2; j++)
				{
					fscanf(nFile, "%lf", &basis[i][j]);
				}
			}
		}

		for (i = 0; i < fftlen / 2 + 1; i++)
		{
			for (j = 0; j < R_x; j++)
			{
				basis_f[i][j] = basis[i][j];
			}
			for (; j < 2 * R_x; j++)
			{
				basis_m[i][j - R_x] = basis[i][j];
			}
			for (; j < 2 * R_x + R_d; j++) //female noise = noise 50rank female_noise 50rank
			{
				basis_fn[i][j - 2 * R_x] = basis[i][j];
			}
			//male noise = noise 50rank male_noise 50rank
			for (j = 2 * R_x; j < 2 * R_x + R_d / 2; j++)//male noise
			{
				basis_mn[i][j - 2 * R_x] = basis[i][j];
			}
			for (j = 2 * R_x + R_d; j < 2 * R_x + 3 * R_d / 2; j++) //male noise
			{
				basis_mn[i][j - (2 * R_x + R_d / 2)] = basis[i][j];
			}

		}
	}

	//OCTA-Capture로 들어오는 input 48k이므로 process를 진행하기 위해 16k로 Down Sampling을 해야한다.
	isNew16k = (BuffCnt == 2);
	for (ch = 0; ch < NchASIOinput; ch++)
	{
		for (i = 0; i < NBufferSize; i++)
		{
			XX[i + NBufferSize] = 0;
			XX[i] = input[ch][i];
		}

		for (i = 0; i < NBufferSize; i++)
		{
			xx_lp[ch][BuffCnt*NBufferSize + i] = InitCond[ch][i] + XX[i];
			InitCond[ch][i] = XX_LP[NBufferSize + i];
		}
		if (isNew16k == 1)
		{
			for (i = 0; i < NBufferSize; i++)
				x[ch][i] = xx_lp[ch][i + i + i];
		}
	}
	BuffCnt = (BuffCnt + 1) % 3;


	if (isNew16k == 1)
	{
		for (ch = 0; ch < Nch; ch++) //NMF input
		{
			for (i = 0; i < 3 * NBufferSize; i++)
			{
				input_temp[ch][i] = input_temp[ch][NBufferSize + i];
			}
			for (i = 0; i < NBufferSize; i++)
			{
				input_temp[ch][3 * NBufferSize + i] = x[ch][i];
			}
		}
#if MAKE_FILE ==1
		for (i = 0; i < Nch; i++)
		{
			for (j = 0; j < NBufferSize; j++)
			{
				out_buff[i][j] = x[i][j] * 32768.0;
				rec_out[i][j] = (short)(out_buff[i][j]);
			}
			fwrite(rec_out[i], sizeof(short), NBufferSize, rec[i]);
		}
#endif

		for (ch = 0; ch < Nch; ch++)
		{
			for (i = 0; i < Nwin; i++)
			{
				vad_input_tmp[ch + 1][i + 1] = input_temp[ch][i];
			}
		}

		Process_count++;
		if (Process_count >= 1 && Process_count <= 10)//10번째 프레임까지 노이즈 추측
		{
			my_VAD->Estimate_Noise(vad_input_tmp[1]);
		}
		else
		{
			EPDcase = my_VAD->Calculate_VAD(vad_input_tmp[1]);

			//발화구간일 때 처리한다.
			if (my_VAD->EPD_print() == 2 || my_VAD->EPD_print() == 1)
			{
				//IVA와 NMF결과를 wav로 저장하도록 생성
				VAD_count++;
				char iva1_pcm[500] = ".\\output\\IVA_ch1.pcm";
				char iva2_pcm[500] = ".\\output\\IVA_ch2.pcm";
				if ((_access(iva1_pcm, 0) != 0) && (_access(iva2_pcm, 0) != 0))
				{
#if MAKE_FILE == 1
					char file_name[2][500];
					nmf = new FILE*[Nch];
					iva = new FILE*[Nch];
					upsample = new FILE*[Nch];
					for (ch = 0; ch < Nch; ch++)
					{
						sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
						nmf[ch] = fopen(file_name[0], "wb");
						sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
						iva[ch] = fopen(file_name[0], "wb");
						sprintf(file_name[0], ".\\output\\upsample_ch%d.pcm", ch + 1);
						upsample[ch] = fopen(file_name[0], "wb");
					}
#endif
				}
				//화자 분리 (IVA)
				rtIVA->IVAprocessing(input_temp, output);

#if MAKE_FILE == 1
				for (i = 0; i < Nch; i++)
				{
					for (j = 0; j < NBufferSize; j++)
					{
						out_buff[i][j] = output[i][j] * 32768.0;
						iva_out[i][j] = (short)(out_buff[i][j]);
					}
					fwrite(iva_out[i], sizeof(short), NBufferSize, iva[i]);
				}
#endif
				for (ch = 0; ch < Nch; ch++) //NMF input
				{
					for (i = 0; i < 3 * NBufferSize; i++)
					{
						output_temp[ch][i] = output_temp[ch][NBufferSize + i];
					}
					for (i = 0; i < NBufferSize; i++)
					{
						output_temp[ch][3 * NBufferSize + i] = output[ch][i] * 32768.0;
					}
				}

				//SNMF Test - training된 basis를 이용하여 SNMF를 한다.
				//여성화자에 대해 enhancing을 할 경우 잡음 basis로 남성화자와 잡음에 대한 basis를 사용한다.
				f_snmf->SNMF_test(output_temp[0], nmf_test_buff[0], basis_f, basis_mn);
				//남성화자에 대해 enhancing을 할 경우 잡음 basis로 여성화자와 잡음에 대한 basis를 사용한다.
				m_snmf->SNMF_test(output_temp[1], nmf_test_buff[1], basis_m, basis_fn);
				for (ch = 0; ch < Nch; ch++)
				{
					for (i = 0; i < frame_shift; i++)
					{
						nmf_out[ch][i] = (short)(nmf_test_buff[ch][i]);
					}
#if MAKE_FILE == 1
					fwrite(nmf_out[ch], sizeof(short), frame_shift, nmf[ch]);
#endif
				}

				//Upsampling 16k -> 48k
				for (ch = 0; ch < Nch; ch++) //NMF input
				{
					for (i = 0; i < NBufferSize; i++)
					{
						test2_pad[ch][30 + i * 3] = (double)nmf_out[ch][i];
					}
				}

				for (ch = 0; ch < Nch; ch++)
				{
					for (i = 0; i < NBufferSize * 3; i++)
					{
						test2_48k[ch][i] = 0.0;
						for (j = 0; j < 61; j++)
						{
							test2_48k[ch][i] += test2_pad[ch][i + j] * up_filt[61 - j];
						}
					}
				}
#if MAKE_FILE == 1
				for (i = 0; i < Nch; i++)
				{
					for (j = 0; j < 3 * NBufferSize; j++)
					{
						test2_48[i][j] = (short)test2_48k[i][j];
					}
					fwrite(test2_48[i], sizeof(short), NBufferSize * 3, upsample[i]);
				}
#endif
				return 1;
			}
			else
			{
				//발화구간이 끝나는 대로 wav파일 저장
				if (VAD_count!=0)
				{
					rec_num++;
					VAD_count = 0;
#if MAKE_FILE == 1
					char file_name[2][500];
					for (ch = 0; ch < Nch; ch++)
					{
						fclose(nmf[ch]);
						sprintf(file_name[0], ".\\output\\NMF_ch%d.pcm", ch + 1);
						sprintf(file_name[1], ".\\output\\NMF_ch%d_%d.wav", ch + 1, rec_num);
						pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
						remove(file_name[0]);
						fclose(iva[ch]);
						sprintf(file_name[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
						sprintf(file_name[1], ".\\output\\IVA_ch%d_%d.wav", ch + 1, rec_num);
						pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq));
						remove(file_name[0]);
						fclose(upsample[ch]);
						sprintf(file_name[0], ".\\output\\upsample_ch%d.pcm", ch + 1);
						sprintf(file_name[1], ".\\output\\upsample_ch%d_%d.wav", ch + 1, rec_num);
						pcm2wav(file_name[0], file_name[1], (long)(SamplingFreq * 3));
						remove(file_name[0]);
					}
#endif
				}
			}

		}
	}
	return 0;

}

