close all; 
clear all; 
clc;

% I. populate alarm signals, noise signals, set SNR and mixing 
fs = 8000;

% alarm permutations
alarm = dir('./Alarm/*.wav');
alarm_idx = 1:length(alarm);


% noise permutations
noise = dir('./Noise/*.wav');
noise_idx = 1:length(noise);


% snr permutations
snr = [96 -96];


% read all alarms and noises to mem
alarm_fs = zeros(length(alarm), 1);
noise_fs = zeros(length(noise), 1);
alarm_dat = {};
noise_dat = {};

for i = alarm_idx
    [x, alarm_fs(i)] = audioread(strcat('./Alarm/', alarm(i).name));
    x = resample(x(:,1), fs, alarm_fs(i));
    alarm_dat{i} = x;
    clear x;
end
for i = noise_idx
    [x, noise_fs(i)] = audioread(strcat('./Noise/', noise(i).name));
    x = resample(x(:,1), fs, noise_fs(i));
    noise_dat{i} = x;
    clear x;
end

% mixing for different train/test cases
for i = alarm_idx
    for k = snr
        for j = noise_idx
            
            x = alarm_dat{i};
            y = noise_dat{j};
            nh = randi(length(y));
            y = y(1+mod(nh+0:nh+length(x)-1, length(y)));
            
            xw = x(1:floor((length(x)-512)/256)*256+512);
            xw = buffer(xw, 512, 256, 'nodelay');
            xw = sum(xw.^2);

            yw = y(1:floor((length(x)-512)/256)*256+512);
            yw = buffer(yw, 512, 256, 'nodelay');
            yw = sum(yw.^2);
            
            snrw = yw ./ xw;
            snrw = max(snrw);
            x =  x .* sqrt((snrw * (10.^(k/10))));
            
            x = x + y;
            scal = 0.9 / max(abs(x));
            x = x .* scal;
            
            clip = strcat(num2str(i), '_', num2str(k), '_', num2str(j));
            mixed = strcat('./Mixed/', clip, '.wav');
            audiowrite(mixed, x, fs, 'BitsPerSample', 16);
        
        end
             
    end
end
clear x;
clear y;



%% II. Feature extraction
close all; clear all; clc;

files = dir('./Mixed/*.wav');

% alarm permutations
alarm = dir('./Alarm/*.wav');
alarm_tot_cnt = length(alarm);
alarm_idx = randperm(alarm_tot_cnt);


% noise permutations
noise = dir('./Noise/*.wav');
noise_tot_cnt = length(noise);
noise_idx = randperm(noise_tot_cnt);


% snr permutations
snr = [0 -96];
snridx = randperm(length(snr));
snr = snr(snridx);             


% mixing for different train/test cases
ntrain = floor(noise_tot_cnt * 0.6);
mat_train = [];
mat_test = [];
alarmcnt = 0;
nframe = 512;
for i = alarm_idx
    for k = snr
        
        % training set
        for j = noise_idx(1:ntrain)
            
            clip = strcat('./Mixed/', num2str(i), '_', num2str(k), '_', num2str(j), '.wav');
                      
            [x, fs] = audioread(clip);
            %xcorrfea = AlarmXcorr(x, fs, nframe);
            %xcorrfea = iso7731(x, nframe, nframe/2, fs);
            %xcorrfea = xcorrfea(201:end);
            [feaMat, ~, ~, ~] = FeaExtClsCore(x, fs, 'fileLabel_Train_1000.model', nframe/2, nframe); 
            
            if k < -10
                fea = [-ones(size(feaMat,1), 1)  feaMat];  %feaMat(101:100+length(xcorrfea),:)
            else
                fea = [ones(size(feaMat,1), 1) feaMat];
            end
            mat_train = [mat_train; fea];
        end
        
        % test set
        for j = noise_idx(ntrain+1:end)
            
            clip = strcat('./Mixed/', num2str(i), '_', num2str(k), '_', num2str(j), '.wav');
            
            [x, fs] = audioread(clip);           
            %xcorrfea = AlarmXcorr(x, fs, nframe);
            %xcorrfea = iso7731(x, nframe, nframe/2, fs);
            %xcorrfea = xcorrfea(201:end);
            %[feaMat, ~, ~, ~] = FeaExtClsCore(x, fs, 'fileLabel_Train_1000.model', nframe/2, nframe);
            
            if k < -10
                fea = [-ones(size(feaMat,1), 1)  feaMat];
            else
                fea = [ones(size(feaMat,1), 1)  feaMat];
            end
            mat_test = [mat_test; fea];
        end
    end
    alarmcnt = alarmcnt + 1
end

dlmwrite('./Calib/Train.txt', mat_train, 'delimiter',' ','precision','%.8f');
dlmwrite('./Calib/Test.txt', mat_test, 'delimiter',' ','precision','%.8f');


%% III. Train/Test and Evaluate the results
ml_exe = 'boost.exe real';
ml_train_arg = '-p 1 -l 500 -r 1 -d 4 -i RAW';
ml_test_arg = '-t -p 1 -l 500 -r 1 -d 4 -i RAW';

train_cmd = {ml_exe, ml_train_arg, '-m', './Calib/Model.txt', '-f', './Calib/Train.txt > ./Calib/log.txt'};
test_cmd = {ml_exe, ml_test_arg, '-m', './Calib/Model.txt', '-f', './Calib/Test.txt', '-o', './Calib/Result.txt >> ./Calib/log.txt'};


system(strjoin(train_cmd));
system(strjoin(test_cmd));


%% find false alarms by confidence inspection on samples
ModelFile = './Calib/Model.txt';
FrameShift = 256;
FrameLength = 512;

alarm = 1;
snr = -5;

for noise = 1:10   % total 116
    clip = strcat('./Mixed/',num2str(alarm), '_', num2str(snr), '_', num2str(noise), '.wav');
    [x, fs] = audioread(clip);
    [feaMat,confidence,time,pow] = FeaExtClsCore(x, fs, ModelFile, FrameShift, FrameLength);

    figure;
    subplot(2,1,1);
    spectrogram(x, 512, 256, 512, fs, 'yaxis');
    subplot(2,1,2); grid on;
    plot(time, confidence);
end

%% plot false alarm curve
figure; hold on; grid on;
falsealarm_avg = [];

[loop,acc,recall,precision,falsealarm] = getEvalResult('./Calib/Result.txt');
falsealarm_avg = [falsealarm_avg; falsealarm];
plot(falsealarm);









%% Manual selection of tonality feature 
close all; clear all; clc;

figure;
nframe = 128;

clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise\WGN.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Mixed\6_-96_11.wav'
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise_aside\buccaneer2.wav';
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Audience\Audience_0000000003.wav';
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Speech\Speech_0000000009.wav';
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Speech\Speech_0000000061.wav';
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Speech\Speech_0000000063.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise2\Speech_0000000061.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise\WGN.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise_aside\Train_Station_binaural(0.00-30.00s)_voiceremoved.wav';

[x, fs] = audioread(clip);
x = resample(x(1:end,1), 16000, fs);
fs = 16000;
%fea = AlarmXcorr(x, fs, nframe);
fea = iso7731(x, nframe, nframe/2, fs);

subplot(2,1,1);
plot(fea, '-r'); hold;
mu = mean(fea(21:end));
me = median(fea(21:end));
n = length(fea(21:end));
plot(20+1:20+n, mu*ones(n,1), '--r');
plot(20+1:20+n, me*ones(n,1), '-.r');


clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Mixed\12_0_7.wav'
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Phone\Phone_0000000770.wav'; 
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Horn\Horn_0000000636.wav';   
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Timer\Timer_0000000663.wav'; 
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Sms\Sms_0000001040.wav';    

%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Broadcast\Broadcast_0000000671.wav';%fixed 
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000882.wav';%fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Horn\Horn_0000000260.wav';%trade sensitivity for precision
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Horn\Horn_0000000258.wav';

%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000960.wav';%fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000825.wav';%fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000915.wav';

%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\58016__guitarguy1985__piercer.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\188004__motion-s__police-car-siren.wav'; %failure, short burst
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\159755__conleec__misc-police-siren-in-stop-officer-001.wav'; %voice part, failure
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\159752__conleec__misc-police-siren-int-onboard-002.wav';

[x, fs] = audioread(clip);
x = resample(x(1:end,1), 16000, fs);
x = awgn(x, 96, 'measured');
fs = 16000;
%fea = AlarmXcorr(x, fs, nframe);
fea = iso7731(x, nframe, nframe/2, fs);

subplot(2,1,1);
plot(fea, '-b');
mu = mean(fea(21:end));
me = median(fea(21:end));
n = length(fea(21:end));
plot(20+1:20+n, mu*ones(n,1), '--b');
plot(20+1:20+n, me*ones(n,1), '-.b');
grid;

subplot(2,1,2);
spectrogram(x, 512, 256, 512, fs, 'yaxis'); hold;




%% Horn frequency band statistics
B{1} = [2538 2980 3411 3831];
B{1} = [1517 2168 2549];
B{3} = [500 1000 1500 2000 2500 3000 3500 4000 4500];
B{4} = 375 * (1:12);
B{5} = 415 * (1:11);
B{6} = 420 * (1:11);
B{7} = 500 * (1:9);
B{8} = [2400 2800 3200 3600 4000 4400];
B{9} = [2946 3377 3865 4353 4841];
B{10} = [2600 3000 3400 3800];
B{11} = [2800 3200 3525 3877];

figure; hold on; grid on;
for i = 1:11
    plot(zeros(1, length(B{i})), B{i}, 's', 'MarkerSize', 6);
end





%% Max value Tracking

close all; clear all; clc;


%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\annotated\CarHorn_06r.wav';
clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\annotated\Speech_allison.wav';

[x, fs] = audioread(clip);
if size(x, 2) == 2
   x = sum(x, 2)./2;
end
x = resample(x(1:end,1), 16000, fs);
fs = 16000;

% pre-emphsis to flatten the spectrum
PolyB = [1 -0.99];
PolyA = 1;
x = filter(PolyB, PolyA, x); 



nframe = 512*16;
nshift = 256;

subplot(2,1,1);
spectrogram(x, nframe, nframe - nshift, nframe, fs, 'yaxis');
hold on;



[temp] = XcorrFeature(x, 512*16, 256, fs);
subplot(2,1,2); plot(rms(temp,1)); grid on;

%%


x = x(1:floor((length(x) - nframe) / nshift) * nshift + nframe); 
x = buffer(x, nframe, nframe-nshift, 'nodelay');
[~, ncol] = size(x);
x = fft(x.*repmat(hamming(nframe), 1, ncol));
%x = fft(x);
x = abs(x(1:nframe/2,:));

%slicep = zeros(nframe/2, 1) * NaN;
%slicen = zeros(nframe/2, 32) * NaN;

%pen = zeros(ncol,1);
%spread = zeros(ncol,1);
feature_zeta = zeros(ncol,1);
m_avg = 2;
for i = 1:ncol
    
    th = (sum(x(:,i).^m_avg)/(nframe/2)).^(1/m_avg);
    %[~, rmax] = quadmaxlocate(20*log10(x(:,i)'), 20*log10(th));
    [ymax, rmax] = maxlocate(20*log10(x(:,i)'), 20*log10(th));
    
    %if isnan(ymax)
    %else
    %    [ymax, yidx] = sort(ymax, 2, 'descend');
        %if length(ymax) > M
        %    rmax = rmax(yidx(1:M));
        %end
    %end
    
    %slice = zeros(nframe/2, 1) * NaN;
    %slice(rmax) = 1;    
    %pen(i) = sum(isnan(slice - slicep));
    %slicep = slice;
    
%     slicen(:,1:end-1) = slicen(:,2:end);
%     slicen(:,end) = isnan(slice - slicep);
%     spread(i) = sum((32-sum(slicen')) ~= 0);
    
    
    %calculate experimental feature
    %xm = (th.^2) * ones(nframe/2, 1);
    %xm(rmax) = ymax.^2;    
    %ita = median(diff(rmax)) / length(rmax);
    %feature_zeta(i) = ita * geomean(xm)/mean(xm);
    
    ita = max(diff(rmax)) / numel(rmax);
    [sesf, ratio] = SegEnSpecFlatness(x(:,i));
    feature_zeta(i) = sesf;
    
    tt = (nframe/2 + i * nshift) / fs;           % default specgram step is NFFT/2 i.e. 128
    F = rmax * fs / nframe;                  % Convert R from bins to Hz
    %subplot(2,1,1);
    %plot(tt*ones(length(F),1), F', 'x', 'MarkerSize', 2);              % the tracks follow the specgram peaks
end


tc = 0.02;
alpha = 1-exp(-nshift/fs/tc);
feature_zeta_sm = zeros(size(feature_zeta));
for i = 2:length(feature_zeta)
    feature_zeta_sm(i) = feature_zeta(i-1) * alpha + (1 - alpha) * feature_zeta_sm(i-1);
end
subplot(2,1,2);
plot(feature_zeta_sm); xlim([1 ncol]); grid on;


%% Test of ISO7731 feature
close all; clear all; clc

fs = 16000;
f = 1000;
t = 1;
n = 128;
m = n/2;

nsps = t*fs;
x1 = sin(2*pi*(f/fs)*(0:nsps-1));
x2 = 0.9*sin(2*pi*((f-100)/fs)*(0:nsps-1));
x = (x1 + x2)/2;

ba = [200 1400]; 
ba = ba./fs.*n;
ba(:,1) = ceil(ba(:,1));
ba(:,2) = floor(ba(:,2));
ba

% cut integer number of sliding window from the time series
x = x(1:floor((length(x) - n) / m) * m + n);
x = buffer(x, n, n-m, 'nodelay');

    % zero centering
    u = fft(x(:,30));
    mu = sum(x(:,30));
    u(1) = u(1) - mu;
     
        v = zeros(size(u));
        % v must be complex conjugate
        v(ba(1):ba(2)) = u(ba(1):ba(2));
        v(n+2-ba(2):n+2-ba(1)) = u(n+2-ba(2):n+2-ba(1));
        
        % circular convolution
        rvv = real(fft(v.*conj(v)));
        [~, vmax] = quadmaxlocate(rvv', min(rvv));
        vd = diff(vmax);
        y = std(vd)/mean(vd); 


