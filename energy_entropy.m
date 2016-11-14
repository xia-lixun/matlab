% feature verification
close all; clear all; clc;

[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Siren/315_03.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Honk/atv_horn_part_66065.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Speech/Speech_2_china.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Nox/Nox_4_windnoise.wav');

x = resample(x, 16000, fs);
fs = 16000;




frame_size = 512;
step_size = frame_size / 2;
m = floor((numel(x) - frame_size) / step_size);
x = x(1 : frame_size + m * step_size);
m = m+1;

u = buffer(x, frame_size, frame_size - step_size, 'nodelay');
u = u - ones(frame_size,1) * mean(u);

ste = sum(u.^2);


i = 1;
j = 0;
while j + 64 <= 512
    s2(i,:) = sum(u(j+1:j+64,:).^2) ./ ste;
    j = j + 64;
    i = i + 1;
end

ee = -sum(s2 .* log2(s2));

figure; plot(ee); grid on;
