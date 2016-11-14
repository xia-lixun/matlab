% demostration of time aliasing issue in stft processors
% lixun, 2016

close all; clear all; clc;

fs = 16000;
t = 2;
nsps = fs * t;


nfft = 512;
nhop = nfft / 4;
%% show that hamming window with %75 overlap for COLA condition
x = zeros(nsps, 1);
for k = 0:20
    x = x + [zeros(nhop*k, 1); hamming(nfft); zeros(nsps-nhop*k-nfft,1)];
end
plot(x); grid on;


%% show perfect decompose/reconstruction of white noise
cola = sum(hamming(nfft)) / nhop;

x = rand(nsps, 1) - 0.5;
nu = floor((nsps - nfft) / nhop) * nhop + nfft;
u = buffer(x(1:nu), nfft, nfft - nhop, 'nodelay');

m = size(u, 2);
u = u .* repmat(hamming(nfft), 1, m);

y = zeros(nu,1);
for k = 1:m
    uk = (1/cola) * ifft(fft(u(:,k)));
    y((k-1)*nhop + 1 : (k-1)*nhop + nfft) = y((k-1)*nhop + 1 : (k-1)*nhop + nfft) + uk;
end

%% show low pass filtering of white noise
h = fir1(nfft - 1, 0.25); %freqz(h, 1, 512)
H = fft(h.');

cola = sum(hamming(nfft)) / nhop;
cola = 1.6;

x = rand(nsps, 1) - 0.5;
nu = floor((nsps - nfft) / nhop) * nhop + nfft;
u = buffer(x(1:nu), nfft, nfft - nhop, 'nodelay');

m = size(u, 2);
u = u .* repmat(hamming(nfft), 1, m);

y = zeros(nu,1);
for k = 1:m
    uk = (1/cola) * ifft(H .* fft(u(:,k)));
    y((k-1)*nhop + 1 : (k-1)*nhop + nfft) = y((k-1)*nhop + 1 : (k-1)*nhop + nfft) + uk;
end

pxx = pwelch(x, hamming(nfft));
pyy = pwelch(y, hamming(nfft));
figure; plot(pxx, 'r'); hold on; plot(pyy, 'b--'); grid on;

%% show time-aliasing stft processor in perfect decompose/reconstruct of white nosie
ka = 4;
cola = sum(hamming(nfft*ka)) / nfft;

x = rand(nsps, 1) - 0.5;
nu = floor((nsps - nfft*ka) / nfft) * nfft + nfft*ka;
u = buffer(x(1:nu), nfft*ka, nfft*ka - nfft, 'nodelay');

m = size(u, 2);
u = u .* repmat(hamming(nfft*ka), 1, m);

y = zeros(nu,1);
for k = 1:m
    % add time aliasing here
    ua = reshape(u(:,k), nfft, ka);
    ua = sum(ua, 2) / ka;
    uk = (1/cola) * ifft(fft(ua));
    
    uk = repmat(uk, ka, 1);
    y((k-1)*nfft + 1 : (k-1)*nfft + nfft*ka) = y((k-1)*nfft + 1 : (k-1)*nfft + nfft*ka) + uk;
end

figure; plot(x, 'r'); hold on; plot(y, 'b--'); grid on;