close all; clear all; clc;


N = 512;
s = rand(N,1) - 0.5;
st = s - mean(s);

%1. temporal autocorrelation
for tau = 0:127
    xmm0 = st(1:N-tau);
    xmm1 = st(1+tau:N);
    Rt(tau+1) = sum( xmm0 .* xmm1 ) / sqrt( sum(xmm0.*xmm0) * sum(xmm1.*xmm1) );
end

%2. spectral autocorrelation
st = st .* hamming(512);
sw = abs(fft(st));
sw = sw(1:N/2);
sw = sw - mean(sw);

for tau = 0:127
    xmm0 = sw(1:N/2-tau);
    xmm1 = sw(1+tau:N/2);
    Rs(tau+1) = sum( xmm0 .* xmm1 ) / sqrt( sum(xmm0.*xmm0) * sum(xmm1.*xmm1) );
end

STA = 0.5 * Rt + 0.5 * Rs;