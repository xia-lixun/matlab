%
% demo of car horn/siren detection feature based on coherence of log-x
% power spectral density estimation
% Lixun.Xia2@harman.com
% 2016-06-06
% 
close all; clear all; clc;


% resample the wav file to 16kHz, only mono channel is used
%[x, fs] = audioread('D:/Workspace/BenchMark/ParaExt/x64/Release/mix/horn.wav');
[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/LeaveOut/Honking/CarHorn_04r.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Whistle/whistle_01-100_0178.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Speech/Speech_2_china.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Engine/EngineNoise_01r.wav');

x = resample(x, 16000, fs);
fs = 16000;


% psd estimate parameters
frame_size = 512;
step_size = frame_size / 2;
m = floor((numel(x) - frame_size) / step_size);
x = x(1 : frame_size + m * step_size);
m = m+1; % we have m frames of spectra

% use hamming window, calculate coherence gain for signal RMS estimate
win = hamming(frame_size);
cg = sum(win)/frame_size;
ng = sum(win.^2)/frame_size;

% calculate the periodogram, half-frame step size
u = buffer(x, frame_size, frame_size - step_size, 'nodelay');
u = u - ones(frame_size,1) * mean(u); %centering the time series
u_rms = rms(u);                      %standard deviation of each frame
u = u ./ ( ones(frame_size,1) * u_rms ); %normalize to 0 dBV, w.r.t. 1 Vrms.

% power spectral density estimate
h = u .* repmat(win, 1, m);
h = abs(fft(h)) / frame_size;

puu(1,:) = h(1,:).^2 / ng;
puu(2:frame_size/2,:) = 2 * (h(2:frame_size/2,:).^2) / ng; 
plot(log10(1:frame_size/2), 10*log10(puu)); grid on;

%% find all peaks and sort them from 1 to 8 descending

idx = 120:160;
y = median(puu(32:256,idx),2);
plot(y); hold on; 
plot(1:256, ones(256,1)*rms(y), 'r--');
grid on;
%%

ratio = zeros(m,1);
for k = 1:m
%k = 100; 
    
s = puu(:,k);
ps = sum(s(13:256));

[value, i] = max(s(13:256));
i = i + 12;
ratio(k) =  sum(s(i-3:i+3)) / ps;

end

ratio( find(isnan(ratio)==1) ) = 0;
figure; plot(ratio);


tc = 0.032;
alpha = 1-exp(-step_size/fs/tc);
feature_sm = zeros(size(ratio));
for i = 2:length(ratio)
    feature_sm(i) = ratio(i-1) * alpha + (1 - alpha) * feature_sm(i-1);
end

figure; plot(feature_sm);
