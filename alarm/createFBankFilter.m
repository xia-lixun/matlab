function [filters,startFFT,len] = createFBankFilter(maxFreq,maxBin,numBanks)

maxMel = Hz2Mel(maxFreq);
hz_per_bin = maxFreq/maxBin;

% central frequency in mel
centreMel = (0:numBanks+1)*maxMel/(numBanks+1);

% the start/central/end mel of one filter
melStart = centreMel(1:numBanks);
melCentre = centreMel(2:numBanks+1);
melEnd = centreMel(3:numBanks+2);

% slopes, linear in Mel scale
riseSlope = 1./(melCentre-melStart);
fallSlope = 1./(melCentre-melEnd);

% convert to the the scale of fft bin 
startFFT = floor(Mel2Hz(melStart)/hz_per_bin+1)+1;
startFFT(1) = 1;
endFFT = floor(Mel2Hz(melEnd)/hz_per_bin)+1;
endFFT(endFFT>maxBin) = maxBin;
len = endFFT-startFFT+1;

filters = zeros(maxBin,numBanks);

for k = 1:numBanks
    tmp = zeros(len(k),1);
    mel = Hz2Mel((startFFT(k)-1+(0:len(k)-1))*hz_per_bin);
    idx = (mel<=melCentre(k));
    tmp(idx) = riseSlope(k)*(mel(idx)-melStart(k));
    tmp(~idx) = fallSlope(k)*(mel(~idx)-melCentre(k))+1;
    filters(startFFT(k):endFFT(k),k) = tmp;
end

% filters = [];
% for k = 1:numBanks
%     filters{k} = zeros(len(k),1);
%     mel = Hz2Mel((startFFT(k)-1+(0:len(k)-1))*hz_per_bin);
%     idx = (mel<=melCentre(k));
%     filters{k}(idx) = riseSlope(k)*(mel(idx)-melStart(k));
%     filters{k}(~idx) = fallSlope(k)*(mel(~idx)-melCentre(k))+1;
% end
    
function mel = Hz2Mel(f)
mel = 1127*log(1+f/700);

function f = Mel2Hz(mel)
if mel==0
    f = 0;
else
    f = (exp(mel/1127)-1)*700;
end
