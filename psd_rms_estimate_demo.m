close all; clear all; clc;


fs = 10e3;
n = 0:fs-1;
x = sqrt(2) * sin(2*pi*(100/fs)*n'); % amplitude / sqrt(2) = rms value of sine wave
y = ( rand(fs, 1) - 0.5 ) * 2 * 1e-3;
z = x + y;


win = hamming(fs);
cg = sum(win) / fs;
ng = sum(win.*win) / fs;

h = abs(fft(z .* win));
h = h / fs / cg;

phh(1) = h(1) .* h(1);
phh(2:fs/2) = 2 * h(2:fs/2) .* h(2:fs/2);

plot(log10(1:numel(phh)), 10*log10(phh));
hold on; grid on;


% noise floor
pyy = sum(y .* y) / fs / (fs/2); % true noise spectral density
plot(log10(1:numel(phh)), 10*log10(pyy) * ones(1, numel(phh)), 'r--');  % real noise power
plot(log10(1:numel(phh)), 10*log10(pyy*ng/cg/cg) * ones(1, numel(phh)), 'r');  % real noise power




%% normalization for reading noise values
close all; clear all; clc;

fs = 10e3;
n = 0:fs-1;
x = sqrt(2) * sin(2*pi*(100/fs)*n'); % amplitude / sqrt(2) = rms value of sine wave
y = ( rand(fs, 1) - 0.5 ) * 2 * 1e-3;
z = x + y;


win = hamming(fs);
cg = sum(win) / fs;
ng = sum(win.*win) / fs;

h = abs(fft(z .* win));
h = h / fs;

phh(1) = h(1) .* h(1) / ng;
phh(2:fs/2) = 2 * h(2:fs/2) .* h(2:fs/2) / ng;

plot(log10(1:numel(phh)), 10*log10(phh));
hold on; grid on;


% signal compensation
psig = 10*log10(phh(101) * ng / cg / cg);  % compensation: + 10log10(ng/cg^2), for hamming win: +1.345dB
plot(log10(101), psig, 's');               % real signal power
%plot(log10(1:numel(phh)), 10*log10(pyy*ng/cg/cg) * ones(1, numel(phh)), 'r');  % real noise power    




%%
energy_t = sum(x.^2);

h = abs(fft(x));
energy_f = sum(h.^2) / 512;

% add window
win = hamming(512);
cpg = sum(win) / 512;


