function [SHR]=myInitSHR(Fs,fftlen,minf0,maxf0,maxFrequency,weightVector)
if nargin<3
    minf0=62.5;  % assume 32 ms segment (256 samples for 8K sampling rate), which ensures at least two periods for each segment
    maxf0=450;   % for singing voice, this may be not enough
    maxFrequency = 1250; % a historical value (see Hermes, 1988). frequencies above this limit will be ignored
    weightVector = ones(1,maxFrequency);
else
 %   minf0 = param(1);  % assume 32 ms segment (256 samples for 8K sampling rate), which ensures at least two periods for each segment
 %   maxf0 = param(2);      
 %   maxFrequency = param(3); % a historical value (see Hermes, 1988). frequencies above this limit will be ignored   
 %   weightVector = param(4);
end
%--------------- some thresholds specific to the algorithm -------------------------------
maxlogf=log2(maxf0/2); minlogf=log2(minf0/2); % the search region to compute SHR is as low as 0.5 minf0.
N=floor(maxFrequency/minf0); % maximum number harmonics
N=N-mod(N,2);
% N = N*2;   % tenically, we need to use 2N instead of N. However, ...
%N=N*4; %In fact, in most cases we don't need to multiply N by 4 and get equally good results yet much faster.
%----------------- derive linear and log frequency scale ----------------
frequency=Fs*(1:fftlen/2)/fftlen; % we ignore frequency 0 here since we need to do log transformation later and won't use it anyway.
frequency_idx=find(frequency<=maxFrequency);
frequency=frequency(frequency_idx);
logf=single(log2(frequency));
%clear frequency;
min_bin=logf(end)-logf(end-1); % the minimum distance between two points after interpolation
interp_logf=logf(1):min_bin:logf(end);
i=(2:N);
shift_distance = min(length(interp_logf)-1,round(log2(i)/min_bin));

linear_f = 2.^interp_logf;
% map the weight vector
interp_amplitude_weight = interp_logf;

for i=1:length(interp_amplitude_weight)
    interp_amplitude_weight(i) = weightVector(min(maxFrequency,max(1,round(linear_f(i)))));
end
upperbound=find(interp_logf>=maxlogf); % find out the index of upper bound of search region on the log frequency scale.
upperbound=upperbound(1);% only the first element is useful
lowerbound=find(interp_logf>=minlogf); % find out the index of lower bound of search region on the log frequency scale.
lowerbound=lowerbound(1);

SHR.LOG_SPECTRUM_ON = 1;
SHR.frequency_idx = frequency_idx;
SHR.logf = logf;
SHR.interp_logf = interp_logf;
%SHR.weightVector = ones(fftlen/2,1);
SHR.interp_amplitude_weight = interp_amplitude_weight;
SHR.lowerbound = lowerbound;
SHR.upperbound = upperbound;
SHR.N = N;
SHR.shift_distance = shift_distance;
SHR.numPitch = 1; % default is single pitch
SHR.bins_per_octave = round(1/min_bin);
