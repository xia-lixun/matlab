%compare vst mat features
clear all
close all
vstfoldername = 'D:\Branch\VADVST\Standalone\';
% vstfoldername = 'D:\AudioMulch 2.2.4\';

matpow = dlmread('pow.txt');
mattim = dlmread('tim.txt');
vstpow = dlmread([vstfoldername 'vstfea_pow.txt']);
vsttim = dlmread([vstfoldername 'vstfea_tim.txt']);
len = length(vsttim);
figure('Name','Frame Time Pow','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
subplot(211), plot(mattim(1:len-1)'-vsttim(1:len-1)), ylabel('time');
% subplot(211), plot(1:len-1, mattim(1:len-1), 1:len-1, vsttim(1:len-1)), ylabel('time');
subplot(212), plot(matpow(1:len-1)-vstpow(1:len-1)), ylabel('pow');

clear matpow mattim vstpow vsttim

mat = dlmread('fea.txt');
vst = dlmread([vstfoldername 'vstfea.txt']);
numfeats = 76;
len = length(vst);
idx = 0;
figure('Name','bandEng frameval','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:9
    matdat = mat(j:numfeats:len);
    vstdat = vst(j:numfeats:len);
    subplot(3,3,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','bandEng winval0','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:9
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(3,3,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','bandEng winval1','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:9
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(3,3,j), plot(matdat-vstdat), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','bandFlux','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
matdat = mat(j+idx:numfeats:len);
vstdat = vst(j+idx:numfeats:len);
subplot(1,3,1), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('frameVal'); 

idx = idx+1;
matdat = mat(j+idx:numfeats:len);
vstdat = vst(j+idx:numfeats:len);
subplot(1,3,2), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('winVal0'); 

idx = idx+1;
matdat = mat(j+idx:numfeats:len);
vstdat = vst(j+idx:numfeats:len);
subplot(1,3,3), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('winVal1'); 

idx = idx+1;
figure('Name','mfcc frameval','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:8
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,4,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','mfcc winval0','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:8
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,4,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','mfcc winval1','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:8
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,4,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(int2str(j)); 
end

idx = idx+j;
figure('Name','SHR','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:2
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,3,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(['frameval' int2str(j)]); 
end

idx = idx+j;
for j=1:4
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,3,j+2), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(['winval' int2str(j)]); 
end

idx = idx+j;
figure('Name','harmonicity ','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
for j=1:2
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,3,j), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(['frameval' int2str(j)]); 
end

idx = idx+j;
for j=1:4
    matdat = mat(j+idx:numfeats:len);
    vstdat = vst(j+idx:numfeats:len);
    subplot(2,3,j+2), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel(['winval' int2str(j)]); 
end

figure('Name','Others','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
idx = idx+j+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(221), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('SpecFlux FrameVal');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(222), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('Specflux winval');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(223), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('PowSpecSkew frameval');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(224), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('PowSpecSkew winval');

figure('Name','Others','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(321), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('pausect winval');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(322), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('zcr frameval');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(323), plot(matdat-vstdat), ylabel('zcr winval 1');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(324), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('zcr winval 2');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(325), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('zcr winval 3');
idx = idx+1;
matdat = mat(idx:numfeats:len);
vstdat = vst(idx:numfeats:len);
subplot(326), plot(matdat(1:end-1)-vstdat(1:end-1)), ylabel('zcr winval 4');

% figure('Name','CLS','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
% idx = idx+1;
% matdat = mat(idx:numfeats:len);
% vstdat = vst(idx:numfeats:len);
% subplot(1,3,1), plot(matdat-vstdat), ylabel('CLS');
% idx = idx+1;
% matdat = mat(idx:numfeats:len);
% vstdat = vst(idx:numfeats:len);
% subplot(1,3,2), plot(matdat-vstdat), ylabel('Score');
% idx = idx+1;
% matdat = mat(idx:numfeats:len);
% vstdat = vst(idx:numfeats:len);
% subplot(1,3,3), plot(matdat-vstdat), ylabel('Confidence');

figure('Name','classification','NumberTitle','off', 'units','normalized','outerposition',[0 0 1 1]);
matcla = dlmread('cla.txt');
matvad = dlmread('smvad.txt');
vstcla = dlmread([vstfoldername 'vstcla.txt']);
vstvad = dlmread([vstfoldername 'vstcla_smvad.txt']);
len = length(vstcla);
subplot(211), plot(matcla(1:len-1)-vstcla(1:len-1)), ylabel('confidence');
% subplot(211), plot(1:len-1, matcla(1:len-1), 1:len-1, vstcla(1:len-1)), ylabel('confidence');
subplot(212), plot(matvad(1:len-1)'-vstvad(1:len-1)), ylabel('vad');

clear matcla matvad vstcla vstvad

