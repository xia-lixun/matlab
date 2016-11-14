%% sub-band autocorr feature extraction
close all; clear all; clc;



dirpath = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\annotated\';
category = 'CarHorn_*';
filetype = '.wav';
wavfiles = dir(strcat(dirpath,category,filetype));
filetype = '.TextGrid2Lab';
gridfiles = dir(strcat(dirpath,category,filetype));

assert(numel(gridfiles) == numel(wavfiles));
nfiles = numel(gridfiles);




%[+ scan all files and extract all clips of true and false +]
x_true = [];
x_false = [];

for ifile = 1:nfiles
    
    % read audio samples
    [x, fs] = audioread(strcat(dirpath, wavfiles(ifile).name));
    if size(x, 2) == 2
        x = sum(x, 2)./2;
    end
    x = resample(x, 16000, fs);
    fs = 16000;
    
    %parse clip contents from textgrid file
    annotation_mat = load(strcat(dirpath,gridfiles(ifile).name));
    annotation_mat(:,1) = floor(annotation_mat(:,1) * fs)+1;
    annotation_mat(:,2) = floor(annotation_mat(:,2) * fs);
    
    for itie = 1:size(annotation_mat,1)
        if annotation_mat(itie,3) == 1
            x_true = [x_true; x(annotation_mat(itie,1):annotation_mat(itie,2))];
        elseif annotation_mat(itie,3) == 0
            x_false = [x_false; x(annotation_mat(itie,1):annotation_mat(itie,2))];
        end
    end
end

figure;
subplot(2,1,1);
spectrogram(x_true, 512, 512 - 256, 512, fs, 'yaxis');

subplot(2,1,2);
spectrogram(x_false, 512, 512 - 256, 512, fs, 'yaxis');





%% spectrum enhance ment
x = x_true;

PolyB = [1 -0.99];
PolyA = 1;
x = filter(PolyB, PolyA, x); 


m = 4;
nframe = 512;
nframe_cyc = nframe * m;
nshift = nframe/2;

figure; hold;
spectrogram(x, nframe, nframe - nshift, nframe, fs, 'yaxis');


x = x(1:floor((length(x) - nframe) / nshift) * nshift + nframe); 
x = buffer(x, nframe_cyc, nframe_cyc-nshift, 'nodelay');
[~, ncol] = size(x);

x = x.*repmat(bh7(nframe_cyc), 1, ncol);
u = zeros(nframe, ncol);
for i = 1:m
    u = u + x( (i-1)*nframe+1:(i-1)*nframe+nframe, :);
end
u = u./m;

u = fft(u);
u = abs(u(1:nframe/2,:));

feature = [];
m_avg = 4;
for i = 1:ncol
    th = (sum(u(:,i).^m_avg)/(nframe/2)).^(1/m_avg);
    
    %[~, rmax] = quadmaxlocate(20*log10(x(:,i)'), 20*log10(th));
    [ymax, rmax] = maxlocate(20*log10(u(:,i)'), 20*log10(th));
    feature = [feature median(diff(rmax))];
    
    tt = i * nshift / fs;           % default specgram step is NFFT/2 i.e. 128
    F = rmax * fs / nframe;                  % Convert R from bins to Hz
    plot(tt*ones(length(F),1), F', 'x', 'MarkerSize', 4);              % the tracks follow the specgram peaks
end

figure; plot(feature);

%% sub-band pca

fs = 16000;
[b,a,fx,bx,gd]=gammabank(0.35, fs, '', [375 7000]);
nbanks = size(b,1); 


%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\annotated\CarHorn_02r.wav';
%[x, fs] = audioread(clip);
%if size(x, 2) == 2
%    x = sum(x, 2)./2;
%end
%x = resample(x(1:end,1), 16000, fs);
x = x_true;
fs = 16000;

uERB = [];
for each_band = 1:nbanks
    y = filter(b(each_band,:), a(each_band,:), x);
    uERB = [uERB y];
end
clear y;


% count number of overlapping frames
m = 2;
nframe = 512;
nshift = nframe / m;

x = x(1:floor((length(x) - nframe) / nshift) * nshift + nframe); 
x = buffer(x, nframe, nframe-nshift, 'nodelay');
[~, nobsrvs] = size(x);
clear x;


    
% feature alaysis
%for ibank = 1:nbanks
ibank = 16;    

    % generate sub-band PCA
    Ruu = [];
    
    T0 = 0; T0 = floor(T0 * fs)+1;
    T1 = size(uERB, 1);
    
    t = T0;
    while t + nframe < T1
        rxx = xcorr(uERB(t:t+nframe-1, ibank), 'coeff');
        rxx = rxx(nframe+1:end);
        Ruu = [Ruu rxx];    %Ruu is of size nframe-1 x nobsrvs
        t = t + nshift;
    end
    

    %apply PCA for dim reduction
    mu = mean(Ruu);
    [coeff, score, latent] = princomp(Ruu);
    vardist = 100*latent/sum(latent);

    RuuCtd = bsxfun(@minus, Ruu, mu);
    norm(score * coeff' - RuuCtd)
    
    figure; hold on; plot(vardist); grid on;
    figure; 
    subplot(4,2,1); plot(score(:,1)); grid on;
    subplot(4,2,2); plot(score(:,2)); grid on;
    subplot(4,2,3); plot(score(:,3)); grid on;
    subplot(4,2,4); plot(score(:,4)); grid on;
    subplot(4,2,5); plot(score(:,5)); grid on;
    subplot(4,2,6); plot(score(:,6)); grid on;
    subplot(4,2,7); plot(score(:,7)); grid on;
    subplot(4,2,8); plot(score(:,8)); grid on;
    
    
    
    
%    break;
    
%end



%%

score_true = load('carhorn_score.mat');
score_false = load('noise_score.mat');

pt2t = Ruu.' * score_true.score;
pt2f = Ruu.' * score_false.score;

figure; hold on; grid on;
plot(pt2t(:,1), 'r');
plot(pt2f(:,1), 'b');