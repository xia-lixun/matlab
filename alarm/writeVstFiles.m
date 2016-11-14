%write files for use in VST
clc
clear all
close all

vstfolder = 'D:\VST\VADVST\ReferenceCode\';
fs = 16000;
fs2 = 14700;
frameLen   = round(32e-3*fs);
frameShift = round(16e-3*fs);
frameOvlp  = frameLen-frameShift;
nFFT       = max(512,frameLen);
halfNFFT   = nFFT/2;
blockTime  = frameShift/fs;

%MFCC 
tc          = 0.5;
order       = 8;
numBanks    = 19;
maxFreq     = fs/2;
MDCT_SIZE   = numBanks*4;

%bandEng
bandFreq         = [120 200 400 800 1600 3200 8000];
nBands           = length(bandFreq)-1;

%cls
modelfname = 'D:\\VST\\VADVST\\fileLabel_Train_1000.model';

preEmph = preEmphVec(2048);
fid = fopen([vstfolder 'preEmph.dat'], 'w');
fprintf(fid, '%1.8ff,\n', preEmph);
fclose(fid);
clear preEmph

fid = fopen([vstfolder 'dec3filt.dat'], 'w');
tx = [0:fs];           % Time vector for original signal
x48 = sin(2*pi*tx*100/fs);   % Define a sinusoid 
[x16, filt] = resample(x48,1,3);    % Change sampling rate
fprintf(fid, '%1.16ff,\n', filt);
fclose(fid);
clear tx x48 x16 filt

mdct_coeff = cos(2*pi*((0:numBanks-1)'*2+1)*(1:order)/MDCT_SIZE);
fid = fopen([vstfolder 'mdct_coeff.dat'], 'w');
for i=1:numBanks
fprintf(fid, '%1.8ff,\n', mdct_coeff(i,:));
end
fclose(fid);
clear mdct_coeff

foldername = [vstfolder '480\256\'];
mkdir(foldername);
writeBandEngFiles(256,fs,foldername)
foldername = [vstfolder '480\512\'];
mkdir(foldername);
writeBandEngFiles(512,fs,foldername)
foldername = [vstfolder '480/1024/'];
mkdir(foldername);
writeBandEngFiles(1024,fs,foldername)

foldername = [vstfolder '441\256\'];
mkdir(foldername);
writeBandEngFiles(256,fs2,foldername)
foldername = [vstfolder '441\512\'];
mkdir(foldername);
writeBandEngFiles(512,fs2,foldername)
foldername = [vstfolder '441\1024\'];
mkdir(foldername);
writeBandEngFiles(1024,fs2,foldername)

load bandFreq2
nBands2 = length(bandFreq2);
clear bandFreq2

fid = fopen([vstfolder 'vad_defines.h'], 'w');
fprintf(fid, '#define NINPUTS 2\n');
fprintf(fid, '//32 ms \n');
fprintf(fid, '#define VAD_FRAMELEN %d\n', frameLen);
fprintf(fid, '//16 ms \n');
fprintf(fid, '#define VAD_FRAMESHIFT %d\n', frameShift);
fprintf(fid, '#define MAX_VAD_FRAMELEN 1024\n');
fprintf(fid, '#define MAX_VAD_FRAMESHIFT 1024\n');
fprintf(fid, '#define VAD_NFFT %d\n', nFFT);
fprintf(fid, '#define VAD_HALFNFFT %d\n', nFFT/2);
fprintf(fid, '#define MFCC_ORDER %d\n', order);
fprintf(fid, '#define MFCC_NBANKS %d\n', numBanks);
fprintf(fid, '#define MDCT_SIZE (MFCC_NBANKS*4)\n');
fprintf(fid, '#define SHR_MINF0 75\n');
fprintf(fid, '#define SHR_MAXF0 450\n');
fprintf(fid, '#define BANDENG_NBANDS %d\n', nBands);
fprintf(fid, '#define BANDENG_NBANDS2 %d\n', nBands2);
fprintf(fid, '#define MAX_BANDFLUX_NBANDS 64\n');
fprintf(fid, '#define BANDFLUX_NBANDS 28\n');
fprintf(fid, '#define BANDFLUX_MINBAND 4\n');
fprintf(fid, '#define BANDFLUX_SMOOTHTIME 0.0315\n');
fprintf(fid, '#define BANDFLUX_FLOORUPTIME 22.4944\n');
fprintf(fid, '#define BANDFLUX_FLOORDOWNTIME 0.155\n');
fprintf(fid, '#define BANDFLUX_LOWERCUTOFF 3\n');
fprintf(fid, '#define BANDFLUX_UPPERCUTOFF 15\n');
fprintf(fid, '#define SMCONF_VoiceAttTime 0.01\n');
fprintf(fid, '#define SMCONF_VoiceRelTime 0.5\n');
fprintf(fid, '#define SMCONF_HOLDTIME 0.35\n');
fprintf(fid, '#define SMCONF_ContinuousVoiceTime 0.7\n');
fprintf(fid, '#define SMCONF_DecDelta 0.9\n');

fprintf(fid, '#define MAT_EPS 2.220446049250313e-16\n');
fprintf(fid, '#define VAD_EPS 7.6294e-006\n');
fprintf(fid, '#define VAD_EOK    1\n');
fprintf(fid, '#define VAD_EFAIL -1\n');
fprintf(fid, '#define MAXVADFEATURES 200\n');
fprintf(fid, '#define MODELFILENAME "%s"\n', modelfname);
fprintf(fid, '#define MIN(a,b)  (((a)<(b))? (a):(b))\n');
fprintf(fid, '#define MAX(a,b)  (((a)>(b))? (a):(b))\n');
fprintf(fid, '#define ABS(a)  (((a)>0)? (a):(-a))\n');
fprintf(fid, 'void vad_compStatistics(float *data, float w, float *mean, float *std, int len);\n');
fprintf(fid, 'extern float binFreqs[MAX_VAD_FRAMELEN/2 + 1];\n');
 
fclose(fid);

