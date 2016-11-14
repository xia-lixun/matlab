% feature verification
close all; clear all; clc;

[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Siren/wail_and_yelp.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Honk/atv_horn_part_66065.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Speech/Speech_2_china.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Nox/Nox_4_street_noise_potrero_1.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Engine/EngineNoise_01r.wav');

x = resample(x, 16000, fs);
fs = 16000;


frame_size = 512;
step_size = frame_size / 2;
m = floor((numel(x) - frame_size) / step_size);
x = x(1 : frame_size + m * step_size);
m = m+1;

u = buffer(x, frame_size, frame_size - step_size, 'nodelay');
u = u .* repmat(hamming(frame_size), 1, m);
u = fft(u);
u = u(1:frame_size/2, :);

%%
plot(angle(u(:,300)), 'b'); hold on;
plot(angle(u(:,400)), 'r--');




