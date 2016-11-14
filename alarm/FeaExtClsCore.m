function [feaMat,confidence,time,pow] = FeaExtClsCore(x,fs,modelFile, blk, nfft)
% feature extraction and classification core function
% x:          data
% fs:         sampling rate (typically 16kHz)
% modelFile:  if the model is given, the code not only extracts features but
%             also does classification otherwise only feature extraction is
%             performed.
% blk:        Length of block (frameShift)
% nfft        FFT length, same as frame length
% 
% feaMat:     return the feature matrix (a set of features for every frame)
% confidence: the probability of the frame being speech
% time:       time of each frame in second
% pow:        pow of each input frame
%
% The instant frame feature is included by default. It can be changed by
% setting the variable isFrameFeaIncluded. "frame" indicates instant values
% whereas "win" indicates smoothed value either by averaging over a fixed
% window,e.g. 2 second, or by exponential smoothing.
%
% features calculated using a fixed FFT and overlap
%
% NOTE:
% 1, The Harmonicity feature is normalized by default. A more effective 
% feature is an unnormalized calculation. However, it is energy dependent. 
% Therefore to use it we need leveling AGC.
% 2, the feature "ve" (voice band energy) is also dependent on energy and
% it works the best if we have leveling AGC.
% 3, 5e-7/4 used to calculate "be2" is essentially a voice gating, which is
% energy dependent and works the best if we haev leveling AGC.
% 2 and 3 are used anyway in case they help.

writefile = 1;
if nargin<2
    error('Not enough input arguments!');
end
if nargin == 2
    task = 'extract';
end
if nargin >= 3
    task = 'classify';
    model = loadModel(modelFile);
end

isFrameFeaIncluded = 'y';

if (writefile == 1)
fid = fopen('fea.txt', 'w');
if(fid < 0) 
    'issue with file open' 
end
fidCla = fopen('cla.txt', 'w');
if(fidCla < 0) 
    'issue with file open' 
end
end

%% frame sizes
frameLen   = nfft;
frameShift = blk;
frameOvlp  = frameLen-frameShift;
x = [zeros(frameOvlp,1); x];
nFrame     = floor((length(x)-frameLen)/frameShift)+1;
nFFT       = frameLen;
halfNFFT   = nFFT/2;
blockTime  = frameShift/fs;
time       = single(frameLen/2:frameShift:frameLen/2+(nFrame-1)*frameShift)/fs;

%% features
% The time constant of an exponential moving average is the amount of time
% for the smoothed response of a unit set function to reach 1-1/e, which is
% around 63.2% of the final final. 
% The relationship between this time constant tau and the smoothing factor
% alpha is given by the formula: alpha = 1 - exp(-deltaT/tau), where
% deltat is the sampling time interval.

% for this feature, it is better to have less overlap and longer fft size
specFlux.tc              = 1;  % time constant in second
specFlux.BlockWeight     = 1-exp(-blockTime/specFlux.tc);
specFlux.frameVal        = zeros(nFrame,1, 'single');
specFlux.winVal          = zeros(nFrame,1, 'single');

powSpecSkew.nFramePerWin = round(round(2*fs)/frameShift); % need 2 seconds, 1 is not so discriminative
powSpecSkew.frameVal     = zeros(nFrame,1, 'single');
powSpecSkew.winVal       = zeros(nFrame,1, 'single');

pauseCount.nFramePerWin = round(round(1*fs)/frameShift); % 1 second also works, 2 sec a bit better
pauseCount.framePow     = zeros(nFrame,1, 'single');
pauseCount.winVal       = zeros(nFrame,1, 'single');

zeroCross.tc            = 1;
zeroCross.BlockWeight   = 1-exp(-blockTime/zeroCross.tc);
zeroCross.nFramePerWin  = round(round(2*fs)/frameShift); % 1 second also works, 2 sec a bit better
zeroCross.frameVal      = zeros(nFrame,1);
zeroCross.winVal        = zeros(nFrame,4);               
% skew rate, discriminative
% median-to-mean, discriminative
% average, not discriminative
% std, discriminative

MFCC.tc                 = 0.5;
MFCC.BlockWeight        = 1-exp(-blockTime/MFCC.tc);
MFCC.order              = 8;
MFCC.numBanks           = 19;
MFCC.maxFreq            = fs/2;
MFCC.preEmph            = preEmphVec(halfNFFT);
% MFCC.nFramePerWin       = round(round(2*fs)/frameShift);
MFCC.frameVal           = zeros(nFrame,MFCC.order, 'single');
MFCC.winVal             = zeros(nFrame,2,MFCC.order, 'single');      
% mean
% std

Harmonicity.tc           = 0.3;
Harmonicity.BlockWeight  = 1-single(exp(-blockTime/Harmonicity.tc));
Harmonicity.f            = fs/nFFT*(0:halfNFFT-1);
Harmonicity.frameVal     = zeros(nFrame,2, 'single');
Harmonicity.winVal       = zeros(nFrame,4, 'single');    
% mean strength
% std strength

maxFreq                  = 1250; % 3300 % 5000
SHR                      = myInitSHR(fs,nFFT,75,450,maxFreq,ones(1,maxFreq));
SHR.tc                   = 0.3;
SHR.BlockWeight          = 1-exp(-blockTime/SHR.tc);
% SHR.nFramePerWin         = round(round(2*fs)/frameShift);  
SHR.maxFreq              = maxFreq;
SHR.frameVal             = zeros(nFrame,2, 'single');
% harm
% pitch
SHR.winVal               = zeros(nFrame,4, 'single');               % harm is much more useful than pitch
% harm mean
% harm std
% pitch mean
% pitch std

bandEng.tc               = 0.3;
bandEng.BlockWeight      = 1-exp(-blockTime/bandEng.tc);
load bandFreq2
load BandTilt
bandEng.bandFreq         = [120 200 400 800 1600 3200 8000];
bandEng.bandFreq2        = bandFreq2; 
bandEng.nBands           = length(bandEng.bandFreq)-1;
bandEng.nBands2          = length(bandEng.bandFreq2);
bandEng.voiceBand        = [0 1 1 1 1 0];
bandEng.lfBand           = [0 1 1 1 0 0];
bandEng.hfBand           = [0 0 0 0 1 1];
bandEng.fullBand         = [0 1 1 1 1 1];
bandEng.binFreq          = fs/nFFT*(0:halfNFFT-1);
bandEng.bandMatrix       = zeros(bandEng.nBands,halfNFFT);
bandEng.bandMatrix2      = zeros(bandEng.nBands2,halfNFFT);
bandEng.bandMatrix3       = zeros(bandEng.nBands,bandEng.nBands2);
for b=1:bandEng.nBands
    bandEng.bandMatrix(b,:) = (bandEng.binFreq>=bandEng.bandFreq(b))&(bandEng.binFreq<bandEng.bandFreq(b+1));
end
% calculate banding matrix from 100 Hz up to fs/2, power preserved
[tmp,bandEng.bandFreq2,~,~] = filtbankm(bandEng.nBands2,nFFT,fs,60,fs/2,'ey');
bandEng.bandTilt = interp1(bandFreq2,BandTilt,bandEng.bandFreq2);
clear BandTilt bandFreq2
bandEng.bandMatrix2         = tmp(:,1:halfNFFT);  %discard the Nyquist freq
for b=1:bandEng.nBands
    bandEng.bandMatrix3(b,:) = (bandEng.bandFreq2>=bandEng.bandFreq(b))&(bandEng.bandFreq2<bandEng.bandFreq(b+1));
end
bandEng.talkSens            = 5e-7/4;
bandEng.frameVal            = zeros(nFrame,bandEng.nBands+3);
bandEng.winVal              = zeros(nFrame,2,bandEng.nBands+3);
% band energy: 1 - bandEng.nBands
% voice band total energy: bandEng.nBands+1
% band engergy normalized over max: bandEng.nBands+2, not used yet, need noise estimate  

bandFlux.tc                 = 0.3;  % time constant in second
bandFlux.BlockWeight        = 1-exp(-blockTime/bandFlux.tc);
bandFlux.frameVal           = zeros(nFrame,1);
bandFlux.winVal             = zeros(nFrame,2); 
% mean
% std
bandFlux.FluxSmoothTime     = 0.0315;    % (sec) Smooth flux
bandFlux.FluxFloorUpTime    = 22.4944;   % (sec) Flux floor rise rate
bandFlux.FluxFloorDownTime  = 0.155;     % (sec) Flux floor drop rate
bandFlux.f                  = fs/nFFT*(0:halfNFFT-1);
bandFlux.FluxBands          = find(bandFlux.f<1000 & bandFlux.f>100); % Bands to use (0-1kHz)
bandFlux.FluxThreshold      = 3;          % (dB) Flux Lower cutoff
bandFlux.FluxMax            = 15;         % (dB) Flux Upper cutoff
bandFlux.FluxSmoothAlpha    = 1-exp(-frameShift/fs/bandFlux.FluxSmoothTime);
bandFlux.FluxFloorUpAlpha   = 1-exp(-frameShift/fs/bandFlux.FluxFloorUpTime);
bandFlux.FluxFloorDownAlpha = 1-exp(-frameShift/fs/bandFlux.FluxFloorDownTime);
bandFlux.preBands = ones(length(bandFlux.FluxBands),1);
bandFlux.FluxSmooth = 0;
bandFlux.FluxFloor = 0;
bandFlux.FluxDist = 0;


%% classification results
cls = zeros(nFrame,1);
score = zeros(nFrame,1);
confidence = zeros(nFrame,1);
pow = zeros(nFrame,1);

%% apply windowing and fft
x = single(x(1:(nFrame-1)*frameShift+frameLen));
y = buffer(x,frameLen,frameOvlp,'nodelay');
dataSpec = single(fft(y.*repmat(hamming(nFFT),1,nFrame))/nFFT); % here we use nFFT to scale instead of sqrt(nfft) because we want to make sure the total energy is less than 1 when the signal amplitude is less than 1.
dataSpec = dataSpec(1:halfNFFT,:);
powSpecdB = single(20*log10(abs(dataSpec)+eps));

%% frame-by-frame extraction
% not optimized for Matlab matrix operation but easier to implement on-line
% classifier
for m = 1:nFrame
    
    pow(m) = single(sum(abs(dataSpec(:,m)).^2)*2); %*2 because of halfNfft
    
    % averaged squared L2 norm of weighted flux
    tmp1 = dataSpec(:,m);
    if m==1
        tmp2 = zeros(size(tmp1));
    else
        tmp2 = dataSpec(:,m-1);
    end
    specFlux.frameVal(m) = compSpecflux(abs(tmp1),abs(tmp2),2,0);
    [tmp,~] = compStatistics(specFlux.frameVal(m).^2,specFlux.BlockWeight,specFlux.winVal(max(m-1,1)).^2,0);
    specFlux.winVal(m) = sqrt(tmp);
    
    % skew of spectral power
    k = (1:halfNFFT)';
    powSpecSkew.frameVal(m) = (sum(k.*powSpecdB(:,m))*halfNFFT-sum(powSpecdB(:,m))*sum(k))/(halfNFFT*sum(k.^2)-sum(k)^2);
    tmp = powSpecSkew.frameVal(max(m-powSpecSkew.nFramePerWin+1,1):m);
    powSpecSkew.winVal(m) = log(abs(sum(single(tmp-mean(tmp)).^3)+eps)/length(tmp));
    
    % pause count
    data = y(end-frameShift+1:end,m);
    pauseCount.framePow(m) = single(sum(data.^2)/frameShift);
    tmp = pauseCount.framePow(max(m-pauseCount.nFramePerWin+1,1):m);
    thld = single(0.25*sum(tmp)/length(tmp));
    pauseCount.winVal(m) = single(sum(tmp<=thld))/length(tmp);
    
    % skew of zero cross rate, mean-to-median zero cross rate, average zero
    % cross rate, std zero cross rate
    data = y(end-frameShift+1:end,m);
    zeroCross.frameVal(m) = compZCR(data);
    [zeroCross.winVal(m,3),zeroCross.winVal(m,4)] = compStatistics(zeroCross.frameVal(m),zeroCross.BlockWeight,zeroCross.winVal(max(m-1,1),3),zeroCross.winVal(max(m-1,1),4));
    zcr = zeroCross.frameVal(max(m-zeroCross.nFramePerWin+1,1):m);
    tmp = zcr-mean(zcr);
    zeroCross.winVal(m,1) = sum(tmp.^3)./((sum(tmp.^2)).^1.5+eps); % skew rate'
    zeroCross.winVal(m,2) = median(zcr)/(mean(zcr)+eps); % difficult to do it recursively due to median
    % zeroCross.winVal(m,3) = mean(zcr);
    % zeroCross.winVal(m,4) = std(zcr,1);   
    
    % avr MFCC (size MFCC_ORDER*1))
    % std MFCC (size MFCC_ORDER*1))
    MFCC.frameVal(m,:) = compMFCC(abs(dataSpec(:,m).*MFCC.preEmph ),MFCC.maxFreq,MFCC.numBanks,MFCC.order);
    [MFCC.winVal(m,1,:),MFCC.winVal(m,2,:)] = compStatistics(MFCC.frameVal(m,:).',MFCC.BlockWeight,squeeze(MFCC.winVal(max(m-1,1),1,:)),squeeze(MFCC.winVal(max(m-1,1),2,:)));
    % tmp = MFCC.frameVal(max(m-MFCC.nFramePerWin+1,1):m,:);
    % MFCC.winVal(m,1,:) = mean(tmp,1);
    % MFCC.winVal(m,2,:) = std(tmp,1,1);
    
    % harmonicity strength mean, std
    tmp = dataSpec(:,m);
    [Harmonicity.frameVal(m,1),Harmonicity.frameVal(m,2)] = compHarmonicity(tmp, Harmonicity.f);
    % gm:correction
    Harmonicity.frameVal(m,1) = 0.001*Harmonicity.frameVal(m,1);
    Harmonicity.frameVal(m,2) = 2.5*Harmonicity.frameVal(m,2);
    [Harmonicity.winVal(m,1),Harmonicity.winVal(m,2)] = compStatistics(Harmonicity.frameVal(m,1),Harmonicity.BlockWeight,Harmonicity.winVal(max(m-1,1),1),Harmonicity.winVal(max(m-1,1),2));
    [Harmonicity.winVal(m,3),Harmonicity.winVal(m,4)] = compStatistics(Harmonicity.frameVal(m,2),Harmonicity.BlockWeight,Harmonicity.winVal(max(m-1,1),3),Harmonicity.winVal(max(m-1,1),4));
    % tmp = Harmonicity.frameVal(max(m-Harmonicity.nFramePerWin+1,1):m);
    % Harmonicity.winVal(m,1) = mean(tmp);
    % Harmonicity.winVal(m,2) = std(tmp,1);
    
    % SHR harm and pitch mean and std
    tmp = dataSpec(:,m);
    [harm,pitch,~]=myFuncSHR(abs(tmp),SHR);
    SHR.frameVal(m,1) = 0.1*harm; % a bit normalization so that it is at the level of [0 1] approximately.
    SHR.frameVal(m,2) = 0.001*pitch; % a bit normalization so that it is at the level of [0 1] approximately.
    [SHR.winVal(m,1),SHR.winVal(m,2)] = compStatistics(SHR.frameVal(m,1),SHR.BlockWeight,SHR.winVal(max(m-1,1),1),SHR.winVal(max(m-1,1),2));
    [SHR.winVal(m,3),SHR.winVal(m,4)] = compStatistics(SHR.frameVal(m,2),SHR.BlockWeight,SHR.winVal(max(m-1,1),3),SHR.winVal(max(m-1,1),4));
    % tmp = SHR.frameVal(max(m-SHR.nFramePerWin+1,1):m,1);
    % SHR.winVal(m,1) = mean(tmp);
    % SHR.winVal(m,2) = std(tmp,1);
    % tmp = SHR.frameVal(max(m-SHR.nFramePerWin+1,1):m,2);
    % SHR.winVal(m,3) = mean(tmp);
    % SHR.winVal(m,4) = std(tmp,1);
    
    % band eng mean + voice band eng mean + noise removed band energy sum
%     ve = 1+0.2*log10(bandEng.voiceBand*(bandEng.bandMatrix*(abs(dataSpec(:,m)).^2))+1e-10); % so the range is between [-1,1]
%     be2 = bandEng.bandTilt.'.*(bandEng.bandMatrix2*(abs(dataSpec(:,m)).^2));
%     be = bandEng.bandMatrix3*be2;
%     be2bak = be2;
%     noisePowLevel = 0;  % estimated from noise reduction module
%     be2 = sum(max(0,be2 - 2*noisePowLevel)./(be2 + bandEng.talkSens))/bandEng.nBands2;    
%     % bandEng.frameVal(m,:) = [be/(sum(abs(dataSpec(:,m)).^2+1e-10));ve;be2]; % note here: it is different from be/sum(be+1e-10) because sum(be) = sum((abs(dataSpec(5:end,m)).^2)) (the first four bins are not included in the be.
%     bandEng.frameVal(m,:) = [be/(sum(be2bak+1e-10));ve;be2];
    be = (bandEng.bandTilt).'.*(full(bandEng.bandMatrix2)*(abs(dataSpec(:,m)).^2));
    % ve = 1+0.2*log10(bandEng.voiceBand*bandEng.bandMatrix3*be+1e-10); % so the range is between [-1,1]
    % gm: correction
    ve = 1+0.2*log10(bandEng.voiceBand*bandEng.bandMatrix3*be+1e-10)+0.13; % so the range is between [-1,1]
    hfe = (bandEng.hfBand*bandEng.bandMatrix3*be)/(bandEng.fullBand*bandEng.bandMatrix3*be+1e-10); % a feature for impulsive noise vs. speech/noise
    noisePowLevel = 0;  % estimated from noise reduction module
    be2 = sum(max(0,be - 2*noisePowLevel)./(be + bandEng.talkSens))/bandEng.nBands2; % bandEng.talkSens is actually a gating for voice, which is energy dependent, needs AGC!
    % gm: correction
    be2 = be2*1.05;
    bandEng.frameVal(m,:) = [bandEng.bandMatrix3*be/(sum(be+1e-10));ve;hfe;be2];
    % tmp = bandEng.frameVal(max(m-bandEng.nFramePerWin+1,1):m,:);
    % bandEng.winVal(m,1,:) = mean(tmp,1);
    % bandEng.winVal(m,2,:) = std(tmp,1,1);
    [bandEng.winVal(m,1,:),bandEng.winVal(m,2,:)] = compStatistics(bandEng.frameVal(m,:).',bandEng.BlockWeight,squeeze(bandEng.winVal(max(m-1,1),1,:)),squeeze(bandEng.winVal(max(m-1,1),2,:)));
    
    % feature: selected band flux with smoothing
    bandFlux = compBandFlux(bandFlux, abs(dataSpec(:,m)).^2);
    bandFlux.frameVal(m) = bandFlux.Flux;
    [bandFlux.winVal(m,1),bandFlux.winVal(m,2)] = compStatistics(bandFlux.frameVal(m),bandFlux.BlockWeight,bandFlux.winVal(max(m-1,1),1),bandFlux.winVal(max(m-1,1),2));
    
    % TODO:
    % the first features in myMLCoreCalculate
    
    % classification
    if strcmp(task,'classify')
        if strcmpi(isFrameFeaIncluded,'n')
            fea = [...                
                squeeze(bandEng.winVal(m,1,:)).',...      % 1-9
                squeeze(bandEng.winVal(m,2,:)).',...      % 10-18
                bandFlux.winVal(m,:),...                  % 19-20   
                squeeze(MFCC.winVal(m,1,:)).',...         % 21-28
                squeeze(MFCC.winVal(m,2,:)).',...         % 29-36                
                SHR.winVal(m,:),...                       % 37-40
                Harmonicity.winVal(m,:),...               % 41-44
                specFlux.winVal(m),...                    % 45
                powSpecSkew.winVal(m),...                 % 46
                pauseCount.winVal(m),...                  % 47
                zeroCross.winVal(m,:),...                 % 48-51                
                ];                    
        else
            fea = [...
                bandEng.frameVal(m,:),...                 % 1-9
                squeeze(bandEng.winVal(m,1,:)).',...      % 10-18
                squeeze(bandEng.winVal(m,2,:)).',...      % 19-27
                bandFlux.frameVal(m),...                  % 28
                bandFlux.winVal(m,:),...                  % 29-30
                MFCC.frameVal(m,:),...                    % 31-38
                squeeze(MFCC.winVal(m,1,:)).',...         % 39-46
                squeeze(MFCC.winVal(m,2,:)).',...         % 47-54
                SHR.frameVal(m,:),...                     % 55-56
                SHR.winVal(m,:),...                       % 57-60
                Harmonicity.frameVal(m,:),...             % 61-62
                Harmonicity.winVal(m,:),...               % 63-66
                specFlux.frameVal(m),...                  % 67
                specFlux.winVal(m),...                    % 68      
                powSpecSkew.frameVal(m),...               % 69
                powSpecSkew.winVal(m),...                 % 70
                pauseCount.winVal(m),...                  % 71
                zeroCross.frameVal(m),...                 % 72
                zeroCross.winVal(m,:),...                 % 73-76
                ];                                       
        end
        
        % if you want to be perfectly the same as the c code (boost64bit.exe or boost32bit.exe),
        % you need to truncate the feature values to %.3f
        % dlmwrite('tmp.txt',fea,'delimiter',' ','precision', '%.3f')
        % fea = dlmread('tmp.txt');
        [cls(m),score(m),confidence(m)] = boostClassify(fea,model);
if (writefile == 1)
        fprintf(fid, '%1.8f\n', fea);
        fprintf(fidCla, '%1.8f\n', confidence(m));
end
     end
end
if (writefile == 1)
fclose(fid);
fclose(fidCla);
end
dlmwrite('pow.txt', pow, 'precision', '%1.8f');
dlmwrite('tim.txt', time, 'precision', '%1.8f');

%% output feature vectors
if strcmpi(isFrameFeaIncluded,'n')
    % only win Level    
    feaMat = [...
        squeeze(bandEng.winVal(:,1,:)),...
        squeeze(bandEng.winVal(:,2,:)),...
        bandFlux.winVal,...
        squeeze(MFCC.winVal(:,1,:)),...
        squeeze(MFCC.winVal(:,2,:)),...
        SHR.winVal,...
        Harmonicity.winVal,...
        specFlux.winVal,...
        powSpecSkew.winVal,...
        pauseCount.winVal,...
        zeroCross.winVal,...
        ];
else
    % include frame level (instant features)
    feaMat = [...        
        bandEng.frameVal,...
        squeeze(bandEng.winVal(:,1,:)),...
        squeeze(bandEng.winVal(:,2,:)),...
        bandFlux.frameVal,...
        bandFlux.winVal,...
        MFCC.frameVal,...
        squeeze(MFCC.winVal(:,1,:)),...
        squeeze(MFCC.winVal(:,2,:)),...
        SHR.frameVal,...
        SHR.winVal,...
        Harmonicity.frameVal,...
        Harmonicity.winVal,...
        specFlux.frameVal,...
        specFlux.winVal,...
        powSpecSkew.frameVal,...
        powSpecSkew.winVal,...
        pauseCount.winVal,...
        zeroCross.frameVal,...
        zeroCross.winVal,...
        ];
end

% if strcmpi(task,'extract') % remove the first initial frames to get clean features
%     dropFrame = floor(fs*1/frameShift);
%     feaMat = feaMat(dropFrame+1:end,:);
%     % cls = cls(dropFrame+1:end);
%     % score = score(dropFrame+1:end);
%     % confidence = confidence(dropFrame+1:end);
%     % time = time(dropFrame+1:end,:);
% end

%% plot
% figure;plot([-22:1:-6],histc(powSpecSkew.winVal,[-22:1:-6]),'b');grid on
% figure;plot([0:0.05:1],histc(pauseCount.winVal,[0:0.05:1]),'b');grid on
% figure;plot([-1:0.05:1],histc(zeroCross.winVal(:,1),[-1:0.05:1]),'b');grid on
% figure;plot([0:0.1:1.5],histc(zeroCross.winVal(:,2),[0:0.1:1.5]),'b');grid on
% figure;plot([0:0.05:0.5],histc(zeroCross.winVal(:,3),[0:0.05:0.5]),'b');grid on
% figure;plot([0:0.03:0.3],histc(zeroCross.winVal(:,4),[0:0.03:0.3]),'b');grid on
