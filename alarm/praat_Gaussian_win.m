close all
clear all

fs = 16000;
effectiveAnalysisWidth = 0.01; % in second

physicalAnalysisWidth = 2*effectiveAnalysisWidth;
nsamp_window = floor (physicalAnalysisWidth*fs);
halfnsamp_window = nsamp_window / 2 - 1;
nsamp_window = halfnsamp_window * 2;

for k = 1:nsamp_window
    nSamplesPerWindow_f = physicalAnalysisWidth*fs;
    phase = k / nSamplesPerWindow_f;       
    imid = 0.5 * (nsamp_window + 1);
    edge = exp (-12.0);
    phase = (k - imid) / nSamplesPerWindow_f;   % -0.5 .. +0.5 
    value = (exp (-48.0 * phase * phase) - edge) / (1.0 - edge);
    window(k) = value;
end


x = sin(2*pi*1000/fs*(0:16000-1));

% time = (0:0.01:length(x)/fs)*fs;
% leftSample = time*fs;
% rightSample = leftSample+1;
% startSample = rightSample-halfnsamp_window;
% endSample = leftSample + halfnsamp_window;

figure;spectrogram(x,window,halfnsamp_window,2^(ceil(log2(nsamp_window))),fs,'yaxis');colorbar
figure;spectrogram(x,hamming(nsamp_window),halfnsamp_window,2^(ceil(log2(nsamp_window))),fs,'yaxis');colorbar
