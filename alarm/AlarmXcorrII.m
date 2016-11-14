function [features, xpeaks] = AlarmXcorrII(x, fs, nframe)

%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Mixed\9_-3_1.wav'; %fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000960.wav';%fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Alarm\Alarm_0000000825.wav';%fixed
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categoraized\Alarm\Alarm_0000000915.wav';%fixed

%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Speech\Speech_0000000738.wav';
%clip = 'D:\Workspace\Checkout\Alarm\audioFeaExt\data\categorized\Speech\Speech_0000000011.wav';

%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\58016__guitarguy1985__piercer.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\188004__motion-s__police-car-siren.wav'; %failure, short burst
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Alarms_FreeSoundOrg\159755__conleec__misc-police-siren-in-stop-officer-001.wav'; %voice part, failure

%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise\m9_car_6db.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise\WGN.wav';
%clip = 'D:\Workspace\Checkout\Alarm\BenchMark\Noise_aside\city_streets_stereo_44.1k.wav';

%[y, fs] = audioread(clip);
%y = resample(y(1:end,1), 16000, fs);
%fs = 16000;

% bandpass filter for alarm specific frequency region
order = 6;
[b1,a1] = butter(order, [675 1633]/(fs/2), 'bandpass');
[b2,a2] = butter(order, [800 1600]/(fs/2), 'bandpass');
[b3,a3] = butter(order, [800 1750]/(fs/2), 'bandpass');
[b4,a4] = butter(order, [600 1200]/(fs/2), 'bandpass');
[b5,a5] = butter(order, [500 1500]/(fs/2), 'bandpass');

x1 = filter(b1,a1,x);
x2 = filter(b2,a2,x);
x3 = filter(b3,a3,x);
x4 = filter(b4,a4,x);
x5 = filter(b5,a5,x);
x = (x1+x2+x3+x4+x5)/5;



m = 16;
%nframe = 512;
nshift = nframe / m;
x = x(1:floor((length(x) - nframe) / nshift) * nshift + nframe); 
x = buffer(x, nframe, nframe-nshift, 'nodelay');
[~, ncol] = size(x); 

x = fft(x);
x = abs(x(1:nframe/2,:));


xpeak = zeros(ncol,1);
for i = 1:ncol
    [~, rmax] = max( x(:,i) );
    xpeak(i) = rmax;
end    

% exponential filtering the feature
tc = 0.016*4;
alpha = 1-exp(-nshift/fs/tc);
xpeaks = zeros(size(xpeak));
for i = 2:length(xpeak)
    xpeaks(i) = xpeak(i-1) * alpha + (1 - alpha) * xpeaks(i-1);
end


% calculate autocorrelation in window sliding
nxcorr = 768;
xpeaksf = [rand(nxcorr-m/2,1); xpeaks];
xpeaksf = buffer(xpeaksf, nxcorr, nxcorr-m/2, 'nodelay');
feature = zeros(size(xpeaksf,2), 1);
for i = 1:size(xpeaksf,2)
    [rxx, ~] = xcorr(xpeaksf(:,i) - mean(xpeaksf(:,i)), 'coeff');  
    [~, xmax] = quadmaxlocate(rxx(nxcorr:end)', 0.2);
    xmax = [1 xmax];
    feature(i) = std(diff(xmax))/mean(diff(xmax));
    feature(isnan(feature)) = 1;
    feature(feature == 0) = 1;
end

% exponential filtering the feature
tc = 0.016/2;
alpha = 1-exp(-nshift/fs/tc);
features = zeros(size(feature));
for i = 2:length(feature)
    features(i) = feature(i-1) * alpha + (1 - alpha) * features(i-1);
end
