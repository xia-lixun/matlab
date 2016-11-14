% correlation demo
close all; clear all; clc;

fs = 16000;
f = 1000;
t = 0.016;
nsps = fs * t;

% 1000Hz to 1512Hz, 16 samples to 10.58 samples
phi = rand(1,1) * 2 * pi
amp = 0.5/512;
for i = 1:nsps
    %x(i) = (1.0 - amp * i) * cos(2*pi*(f/fs)*(i-1)+phi);
    x(i) = cos(2*pi*(f/fs)*(i-1)+phi);
    f = f + 0.2;
end

% f1 = 1000;  % 16 samples
% f2 = 2000-100;  % 8 samples
% %f3 = 2000;  % 8 samples
% x1 = sin(2*pi*(f1/fs)*(0:nsps-1) + phi);
% x2 = sin(2*pi*(f2/fs)*(0:nsps-1) + phi);
% %x3 = sin(2*pi*(f3/fs)*(0:nsps-1));
% x = (x1 + x2) / 2;

x = awgn(x, 100, 'measured');
figure;
plot(x); grid on;


rxx = xcorr(x, x,'coeff');
rxx = rxx(nsps:end);
[ymax, xmax] = quadmaxlocate(rxx, min(rxx));

figure;
plot(rxx); grid on;

figure; grid on;
diffxmax = diff(xmax);
plot(diffxmax(1:end-7), '-s', 'MarkerSize', 4);
mean(diffxmax(1:end-7))
std(diffxmax(1:end-7))
std(diffxmax(1:end-7))/mean(diffxmax(1:end-7))



%% Conclusion:
%  (a)
%  512 samples @16000Hz: min frequency for alarm signal detection need two
%  periods (one local maxima), so Fmin = 16000/(512/2) = 62.5Hz. For peaks
%  of two, we need three periods...For peaks of 10, we need 11 periods, so
%  we have 16000/(512/11) = 343.75Hz. The lower bound of ISO7731 alarm
%  specific frequency region is found to be around 400~500Hz. So window
%  length of 512 is a good choice to start with.
%
%  nframe = 256. 16000/(256/11) = 687.5

%  (b)
%  Stationary single tones have the lowest feature value. "Stationary" here
%  we mean constant frequency and amplitude. If two frequency components
%  appears within the correlation subband, the feature we are working on
%  only works when the frequencies of the two can be divided without
%  residue. Of cause jitters can be tolerated to a certain extent --- our
%  case it is 15Hz. (16000/512/2)?



