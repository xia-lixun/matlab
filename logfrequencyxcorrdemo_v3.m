% feature verification
close all; clear all; clc;

[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Siren/wail.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Honk/atv_horn_part_66065.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Speech/Speech_2_china.wav');
%[x, fs] = audioread('D:/Workspace/Checkout/Alarm/audioFeaExt/data/Pool/Nox/Nox_4_windnoise.wav');

x = resample(x, 16000, fs);
fs = 16000;
x = x + (rand(size(x)) - 0.5);



%%
frame_size = 512;
step_size = frame_size / 2;
m = floor((numel(x) - frame_size) / step_size);
x = x(1 : frame_size + m * step_size);
m = m+1;

u = buffer(x, frame_size, frame_size - step_size, 'nodelay');
u = u .* repmat(hamming(frame_size), 1, m);
u = abs(fft(u));
u = u(1:frame_size/2, :);



% preempsis the spectra
premp = load('preEmph.dat');
premp = premp(1:2048/(frame_size/2):end);
u = u .* repmat(premp, 1, m);



% convert x axis to log and reshape the spectrum for correlation
axis_real = log2(1:256);
tau = min(diff(axis_real));

j = 1;
axis_discrete(j) = 1; j = j + 1;
for i = 14:256  %14:256
    width = round((axis_real(i) - axis_real(i-1)) / tau);
    axis_discrete(j) = axis_discrete(j-1) + width;
    j = j + 1;
end

%%
%plot(axis_discrete, u(13:256,100), 'b'); hold on;
%plot(axis_discrete, u(13:256,100+32), 'r--');

%%
%1. spectral autocorrelation
tc = 0.008;
alpha = 1-exp(-step_size/fs/tc);
for j = 2:256
    u(:,j) = u(:,j) * alpha + (1 - alpha) * u(:,j-1);
end


feature = [];
for j = 1:m-32


a = u(13:256,j);     %13:256
b = u(13:256,j+32);


a_reshaped = [];
b_reshaped = [];
for i = 1:256-13
    a_reshaped = [a_reshaped; a(i) * ones(axis_discrete(i+1) - axis_discrete(i),1)];
    b_reshaped = [b_reshaped; b(i) * ones(axis_discrete(i+1) - axis_discrete(i),1)];
end

    
%figure;
%plot(a_reshaped, 'b'); hold on;
%plot(b_reshaped, 'r--'); grid on;



a_reshaped = a_reshaped - mean(a_reshaped);
b_reshaped = b_reshaped - mean(b_reshaped);
len = numel(a_reshaped);

Rtl = zeros(floor(len/4), 1);
Rtr = zeros(floor(len/4), 1);
for tau = 1:floor(len/4)
    xmm0 = a_reshaped(1:len-tau);
    xmm1 = b_reshaped(1+tau:len);
    Rtl(tau) = sum( xmm0 .* xmm1 ) / ( sqrt( sum(xmm0.*xmm0) * sum(xmm1.*xmm1) ) + eps);
end
for tau = 1:floor(len/4)
    xmm0 = b_reshaped(1:len-tau);
    xmm1 = a_reshaped(1+tau:len);
    Rtr(tau) = sum( xmm0 .* xmm1 ) / ( sqrt( sum(xmm0.*xmm0) * sum(xmm1.*xmm1) ) + eps);
end

%figure;
%subplot(2,1,1); plot(Rtl);
%subplot(2,1,2); plot(Rtr);
feature = [feature max([Rtl; Rtr])];

end

subplot(2,1,1); plot(feature);

tc = 0.064;
alpha = 1-exp(-step_size/fs/tc);
feature_sm = zeros(size(feature));
for i = 2:length(feature)
    feature_sm(i) = feature(i-1) * alpha + (1 - alpha) * feature_sm(i-1);
end
subplot(2,1,2); plot(feature_sm); xlim([1 numel(feature)]); grid on;