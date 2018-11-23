/*
duplex.cpp
by Gary P. Scavone, 2006-2007.

This program opens a duplex stream and passes
input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <io.h>
#include <iostream>
#include "ProcBuffers.h"
#include "sigproc.h"
#include "SNMF.h"
#include "header.h"
/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/

SNMF *n_snmf;

typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16

#define CHANNEL 2
#define BUFFERFRAME 256

int record_num = 0;
int copyend = 0;
/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64
*/

// Platform-dependent sleep routines.
#if defined( __WINDOWS_ASIO__ ) || defined( __WINDOWS_DS__ ) || defined( __WINDOWS_WASAPI__ )
#include <windows.h>
#define SLEEP( milliseconds ) Sleep( (DWORD) milliseconds ) 
#else // Unix variants
#include <unistd.h>
#define SLEEP( milliseconds ) usleep( (unsigned long) (milliseconds * 1000.0) )
#endif

void usage(void) {
	// Error function in case of incorrect command-line
	// argument specifications
	std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
	std::cout << "    where N = number of channels,\n";
	std::cout << "    fs = the sample rate,\n";
	std::cout << "    iDevice = optional input device to use (default = 0),\n";
	std::cout << "    oDevice = optional output device to use (default = 0),\n";
	std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
	std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
	exit(0);
}

struct InputData {
	MY_TYPE* in_buffer;
	MY_TYPE* out_buffer;
	MY_TYPE* z_buffer;
	unsigned long bufferBytes;
	unsigned long totalFrames;
	unsigned long iframeCounter;
	unsigned long oframeCounter;
	unsigned int channels;
};


// 실시간으로 MIC로 들어온 입력을 버퍼에 저장하고 저장한 데이터를 processing에서 처리한 후 실시간으로 출력한다.
int inout(void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	double /*streamTime*/, RtAudioStreamStatus status, void *data)
{
	InputData *iData = (InputData *)data;

	////실시간으로 입력을 받는 부분
	unsigned int frames = nBufferFrames;
	//최대 버퍼사이즈를 넘을 경우 크기를 넘는 데이터는 저장하지 않는다.
	if (iData->iframeCounter + nBufferFrames > iData->totalFrames) 
	{
		frames = iData->totalFrames - iData->iframeCounter;
		iData->bufferBytes = frames * iData->channels * sizeof(MY_TYPE);
	}
	unsigned long in_offset = iData->iframeCounter * iData->channels;
	//저장된 버퍼가 존재하면 이를 process하기 위해 복사하며 처리할 데이터가 존재한다고 record_num으로 확인한다.
	memcpy(iData->in_buffer + in_offset, inputBuffer, iData->bufferBytes);
	iData->iframeCounter += frames;
	record_num++;

	//저장된 데이터가 최대 버퍼사이즈를 넘기면 다시 offset을 0으로 하여 overwriting한다.
	//이미 앞의 데이터는 process를 위해 복사되었다.
	if (iData->iframeCounter >= iData->totalFrames)
	{
		return 2; //현재 시간(32초)이 다 되면 프로그램을 멈추도록 데모하였다.
		iData->iframeCounter = 0;
	}

	//process에서 처리된 데이터가 실시간 출력을 위해 처리 후 발생하는 copyend라는 신호가 뜨면 callback함수의 output으로 복사한다.
	if (copyend)
	{
		if (iData->oframeCounter + nBufferFrames > iData->totalFrames)
		{
			frames = iData->totalFrames - iData->oframeCounter;
			iData->bufferBytes = frames * iData->channels * sizeof(MY_TYPE);
		}
		unsigned long out_offset = iData->oframeCounter * iData->channels;
		memcpy(outputBuffer, iData->out_buffer + out_offset, iData->bufferBytes);
		iData->oframeCounter += frames;
		//처리된 데이터가 출력되면 copyend값을 줄인다.
		copyend--;
		if (iData->oframeCounter >= iData->totalFrames)
		{
			iData->oframeCounter = 0;
		}
	}
	//처리된 데이터가 없으면 (발화구간이 없는경우) 0을 출력한다.
	else if (copyend == 0)
	{
		memcpy(outputBuffer, iData->z_buffer, iData->bufferBytes);
	}
	return 0;
}

void Training()
{
	n_snmf = new SNMF();
	int i, j;
	double *f_stack, *m_stack, *n_stack;
	double **f_tr_out, **m_tr_out, **n_tr_out, **mn_tr_out, **fn_tr_out;
	double **total_basis_x;
	int total_stack = tr_sec * SamplingFreq;
	f_stack = new double[total_stack];
	m_stack = new double[total_stack];
	n_stack = new double[total_stack];

	f_tr_out = new double*[size_basis];
	m_tr_out = new double*[size_basis];
	n_tr_out = new double*[size_basis];
	mn_tr_out = new double*[size_basis];
	fn_tr_out = new double*[size_basis];
	total_basis_x = new double*[size_basis];
	for (i = 0; i < size_basis; i++)
	{
		f_tr_out[i] = new double[R_x];
		m_tr_out[i] = new double[R_x];
		n_tr_out[i] = new double[R_d / 2];
		mn_tr_out[i] = new double[R_d / 2];
		fn_tr_out[i] = new double[R_d / 2];
		total_basis_x[i] = new double[2 * R_x + 3 * R_d / 2];
	}
	//Training은 프로그램 시작 시 최초 한번 동작한다. 처리시간은 평균 1분에서 2분정도 소요된다.
	//만약 기존에 Training하여 basis_mat.txt파일이 존재하면 바로 training함수를 나간다.
	if (_access("basis_mat.txt", 0) != 0)
	{
		//basis training을 하기 위한 6분 24초 데이터 stack을 쌓아 가져온다.
		n_snmf->basis_stack(total_stack, f_stack, m_stack, n_stack);

		//총 basis는 350개로 여성화자 basis 100rank, 남성화자 basis 100rank, 잡음 basis 50rank
		//남성화자 + 잡음 basis 50 + 50 ranks, 여성화자 + 잡음 basis 50 + 50 ranks로 training한다.

		n_snmf->SNMF_training(f_stack, total_stack, R_x, f_tr_out); //female data training

		std::cout << "f_tr_out fin!" << std::endl;

		n_snmf->SNMF_training(m_stack, total_stack, R_x, m_tr_out); //male data training

		std::cout << "m_tr_out fin!" << std::endl;

		n_snmf->SNMF_training(n_stack, total_stack, R_d / 2, n_tr_out);

		std::cout << "n_tr_out fin!" << std::endl;

		n_snmf->SNMF_training(m_stack, total_stack, R_d / 2, mn_tr_out);

		std::cout << "mn_tr_out fin!" << std::endl;

		n_snmf->SNMF_training(f_stack, total_stack, R_d / 2, fn_tr_out);

		std::cout << "fn_tr_out fin!" << std::endl;

		for (i = 0; i < size_basis; i++)
		{
			for (j = 0; j < R_x; j++)//fmale_basis 100rank
			{
				total_basis_x[i][j] = f_tr_out[i][j];
			}
			for (; j < 2 * R_x; j++)//male_basis 100rank
			{
				total_basis_x[i][j] = m_tr_out[i][j - R_x];
			}
			for (; j < 2 * R_x + R_d / 2; j++)//noise_basis 50rank
			{
				total_basis_x[i][j] = n_tr_out[i][j - 2 * R_x];
			}
			for (; j < 2 * R_x + R_d; j++)//female_noise_basis 50rank
			{
				total_basis_x[i][j] = fn_tr_out[i][j - (2 * R_x + R_d / 2)];
			}
			for (; j < 2 * R_x + 3 * R_d / 2; j++)//male_noise_basis 50rank
			{
				total_basis_x[i][j] = mn_tr_out[i][j - (2 * R_x + R_d)];
			}
		}
		/*basis 저장하기*/
		FILE *pFile;
		pFile = fopen("basis_mat.txt", "wt");
		if (pFile != NULL)
		{
			for (i = 0; i < size_basis; i++)
			{
				for (j = 0; j < 2 * R_x + 3 * R_d / 2; j++)
				{
					fprintf(pFile, "%E\t", total_basis_x[i][j]);
				}
			}
		}
		fclose(pFile);
	}
	for (i = 0; i < size_basis; i++)
	{
		delete[] f_tr_out[i];
		delete[] m_tr_out[i];
		delete[] n_tr_out[i];
		delete[] mn_tr_out[i];
		delete[] fn_tr_out[i];
		delete[] total_basis_x[i];
	}
	delete[] f_tr_out;
	delete[] m_tr_out;
	delete[] n_tr_out;
	delete[] mn_tr_out;
	delete[] fn_tr_out;
	delete[] total_basis_x;
	delete[] f_stack;
	delete[] m_stack;
	delete[] n_stack;
	delete n_snmf;
}

int main(void)
{
	//Device정보 0 : PC, 1 : OCTA-CAPTURE
	unsigned int channels, fs, bufferBytes, oDevice = 1, iDevice = 1, iOffset = 0, oOffset = 0;
	double time = 16.0;
	int in_buffer_cnt = 0;
	int out_buffer_cnt = 0;
	int i, j, ch;
	int proc_end = 0;

	double **input, **proc_output;

	if (_access("basis_mat.txt", 0) != 0)
	{
		Training();
	}

	input = new double *[CHANNEL];
	proc_output = new double *[CHANNEL];
	for (i = 0; i < CHANNEL; i++)
	{
		input[i] = new double[BUFFERFRAME];
		proc_output[i] = new double[3 * BUFFERFRAME];
		for (j = 0; j < BUFFERFRAME; j++)
		{
			input[i][j] = 0.0;
		}
	}

	ProcBuffers *proc;
	proc = new ProcBuffers();

	RtAudio adac;
	if (adac.getDeviceCount() < 1) {
		std::cout << "\nNo audio devices found!\n";
		exit(1);
	}
	channels = 2;
	fs = 48000;

	adac.showWarnings(true);

	// Set the same number of channels for both input and output.
	unsigned int bufferFrames = 512;
	RtAudio::StreamParameters iParams, oParams;
	iParams.deviceId = iDevice;
	iParams.nChannels = channels;
	iParams.firstChannel = iOffset;
	oParams.deviceId = oDevice;
	oParams.nChannels = channels;
	oParams.firstChannel = oOffset;

	if (iDevice == 0)
		iParams.deviceId = adac.getDefaultInputDevice();
	if (oDevice == 0)
		oParams.deviceId = adac.getDefaultOutputDevice();

	RtAudio::StreamOptions options;
	//options.flags |= RTAUDIO_NONINTERLEAVED;

	InputData data;
	data.in_buffer = 0;
	data.out_buffer = 0;
	data.z_buffer = 0;

	//위의 inout이라는 callback함수를 스트리밍 하기 위해 여러 argument를 입력하고 open한다.
	try {
		adac.openStream(&oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)&data, &options);
	}
	catch (RtAudioError& e) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		exit(1);
	}

	data.bufferBytes = bufferFrames * channels * sizeof(MY_TYPE);
	data.totalFrames = (unsigned long)(fs * time);
	data.iframeCounter = 0;
	data.oframeCounter = 0;
	data.channels = channels;
	unsigned long totalBytes;
	totalBytes = data.totalFrames * channels * sizeof(MY_TYPE);

	// Allocate the entire data buffer before starting stream.
	data.in_buffer = (MY_TYPE*)malloc(totalBytes);
	data.out_buffer = (MY_TYPE*)malloc(totalBytes);
	data.z_buffer = new MY_TYPE[NBufferSize * channels];
	for ( i = 0; i <NBufferSize * channels; i++)
	{
		data.z_buffer[i] = 0.0;
	}

	if (data.in_buffer == 0 || data.out_buffer == 0) {
		std::cout << "Memory allocation error ... quitting!\n";
		goto cleanup;
	}

	// Test RtAudio functionality for reporting latency.
	std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

	//준비된 callback함수를 streaming 시작한다.
	try {
		adac.startStream();

	}
	catch (RtAudioError& e) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		goto cleanup;
	}

	//streaming이 돌면서 앞의 callback(inout)함수에서 
	while (adac.isStreamRunning())
	{
		if (record_num)
		{
			if (in_buffer_cnt >= fs * time * channels)
			{
				in_buffer_cnt = 0;
			}
			for (ch = 0; ch < channels; ch++)
			{
				for (i = 0; i < bufferFrames; i++)
				{
					input[ch][i] = (data.in_buffer[channels*i + ch + in_buffer_cnt]) / 32768.0; //input자료형에 맞게 변환하여 저장
				}
			}
			in_buffer_cnt += bufferFrames * channels;
			
			//Process에서 발화구간에서는 proc_end에 1을 return한다. 발화구간이 아니면 0을 return한다.
			proc_end = proc->Process(input, proc_output);
			//발화구간에서 처리된 데이터에 대해 출력을 위해 callback함수의 outbuffer에 넣어준다.
			if (proc_end == 1)
			{
				if (out_buffer_cnt >= fs * time * channels)
				{
					out_buffer_cnt = 0;
				}
				for (ch = 0; ch < channels; ch++)
				{
					for (i = 0; i < 3 * bufferFrames; i++)
					{
						data.out_buffer[channels*i + ch + out_buffer_cnt] = (MY_TYPE)proc_output[ch][i]; //input자료형에 맞게 변환하여 저장
					}
				}
				out_buffer_cnt += 3 * bufferFrames * channels;
				copyend += 3;
			}
			record_num--;
		}
		else
		{
			SLEEP(16);
		}
	}



	delete proc;
	for (i = 0; i < CHANNEL; i++)
	{
		delete[] input[i];
	}
	delete[] input;

cleanup:
	if (adac.isStreamOpen()) adac.closeStream();

	return 0;
}
