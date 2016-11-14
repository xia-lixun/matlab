close all; clear all; clc;

a = load('feature1.dat');
b = load('feature-1.dat');

%%
bin = 63;
la = size(a(:,bin));
lb = size(b(:,bin));


figure;
plot([a(:,bin); zeros(lb)], 'r'); hold on;
plot([zeros(la); b(:,bin)], 'b');
grid on;

%%
figure;
k = 1;
start = 62;
stop = start + 7;
for bin = start:1:stop
    subplot(stop-start+1,1,k);
    histogram(a(:,bin)); hold on;
    histogram(b(:,bin));
    k = k + 1;
end


%%
l2_siren = load('L2-siren2horn.dat');
l2_horn = load('L2-horn2horn.dat');
l2_horn_n = load('L2-horn2horn-n.dat');

%l2_siren_II = load('L2-siren2horn-II.dat');
%l2_horn_II = load('L2-horn2horn-II.dat');


figure; hold on;

plot(l2_siren, 'r'); 
plot(l2_horn, 'b');
plot(l2_horn_n, 'k--');
grid on;

figure; hold on;
histogram(l2_horn_n);
histogram(l2_horn);
histogram(l2_siren);

%%
mu = load('mu-mat.dat');
mustd = load('mu-std-mat.dat');

idx = 131;
x = mu(idx,:);
xstd = mustd(idx, :);

plot(x, 'r'); hold on;
plot(xstd, 'b--');
grid on;


%% spectrum optimization for template matching
close all; clear all; clc;

[x, fs] = audioread('D:\Workspace\Checkout\Alarm\audioFeaExt\data\Negative\Honk1\audio_2911223_safety_sport.wav');

frame_sps = 4096;
shift_sps = 256;
ovlp_sps = frame_sps - shift_sps;
frame_n = floor((length(x)-frame_sps)/shift_sps)+1;

x = single(x(1:(frame_n-1)*shift_sps+frame_sps));
y = buffer(x, frame_sps, ovlp_sps, 'nodelay');

%[apply hd window function]
y = y.*repmat(blackmanharris(frame_sps), 1, frame_n);

% [cavg] 
c = zeros(512, frame_n);
for i = 1:8
    c = c + y((i-1)*512+1:i*512, :);
end
c = c / 8;

% [zero padding]
c = [c; zeros(frame_sps-512, frame_n)];

% [dft]
H = abs(fft(c)) ./ frame_sps; % here we use nFFT to scale instead of sqrt(nfft) because we want to make sure the total energy is less than 1 when the signal amplitude is less than 1.
H = H(1:frame_sps/2, :);

% [de-emphsis]
deempsis = -60:0.625:-0.625;
deempsis = 10.^(deempsis/20);
H(1:96,:) = H(1:96,:) .* repmat(deempsis', 1, frame_n);

% [spectrum as weight ratio]
for i = 1:frame_n
    H(:,i) = H(:,i) ./ sum(H(:,i));
end
LogH = single(20*log10(H+eps));

plot(LogH(:, 200)); grid on;

%%

u = x(10001:10512);
u = u .* hamming(512);
u = abs(fft(u)) / 512;
u = u(1:256);
u = u / sum(u);
u = 20*log10(u);
