function pens = penalty(x, fs, nframe)

% [x, fs] = audioread(clip);
% x = resample(x(1:end,1), 16000, fs);
% x = awgn(x, 90, 'measured');
% fs = 16000;
% spectrogram(x, 512, 256, 512, fs, 'yaxis');
% hold on;


order = 6;
[b,a] = butter(order, [430  7000]/(fs/2), 'bandpass');
x = filter(b,a,x);


m = 2;
%nframe = 512;
nshift = nframe / m;
x = x(1:floor((length(x) - nframe) / nshift) * nshift + nframe); 
x = buffer(x, nframe, nframe-nshift, 'nodelay');
[~, ncol] = size(x);
x = fft(x.*repmat(bh7(nframe), 1, ncol));
x = abs(x(1:nframe/2,:));

pen = zeros(ncol,1);
for i = 1:ncol
    th = (sum(x(:,i).^2)/(nframe/2)).^0.5;
    [~, rmax] = quadmaxlocate(20*log10(x(:,i)'), 20*log10(th));
    
    if numel(rmax) == 1
        pid = 0;
    else 
        %pid = numel(rmax) / max(diff(rmax)); 
        pid = (numel(rmax).^2) / (max(rmax) - min(rmax));  % spread out
    end
    
    pen(i) = pid;
    
    %tt = i * nshift / fs;           % default specgram step is NFFT/2 i.e. 128
    %F = rmax * fs / nframe;                  % Convert R from bins to Hz
    %plot(tt*ones(length(F),1), F', '+');              % the tracks follow the specgram peaks
end

% exponential filtering the feature
tc = 0.016*64;
alpha = 1-exp(-nshift/fs/tc);
pens = zeros(size(pen));
for i = 2:length(pen)
    pens(i) = pen(i-1) * alpha + (1 - alpha) * pens(i-1);
end


end