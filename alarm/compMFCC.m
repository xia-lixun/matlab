function mfcc = compMFCC(spec,maxFreq,numBanks,MFCC_ORDER)
% spec is magnitude spectra.
% spec in each block/frame is in column

MDCT_SIZE = numBanks*4;
maxBin = size(spec,1);
M = size(spec,2);

% for real implementation, the filters can be made shorter than maxBin
% because the beginning and ending are zeroes. Here we use the same length
% because malab is better at matrix calculation than for loop.
[filters,startFFT,len] = createFBankFilter(maxFreq,maxBin,numBanks);

spec = repmat(spec.',[1,1,numBanks]);
spec = permute(spec,[2 3 1]);
fbank = squeeze(log(sum(spec.*repmat(filters,[1,1,M]))+1e-5));

normFactor = sqrt(2/numBanks);
mdct_coeff = cos(2*pi*((0:numBanks-1)'*2+1)*(1:MFCC_ORDER)/MDCT_SIZE);
if M>1
    mfcc = squeeze(sum(repmat(fbank.',[1,1,MFCC_ORDER]).*permute(repmat(mdct_coeff,[1 1 M]),[3 1 2]),2))*normFactor;
else
    mfcc = sum(repmat(fbank.',1,8).*cos(2*pi*((0:numBanks-1)'*2+1)*(1:MFCC_ORDER)/MDCT_SIZE))*normFactor;
end
mfcc = mfcc.';




% MDCT_SIZE = numBanks*4;
% maxBin = length(spec);
% 
% % for real implementation, the filters can be made shorter than maxBin
% % because the beginning and ending are zeroes. Here we use the same length
% % because malab is better at matrix calculation than for loop.
% [filters,startFFT,len] = createFBankFilter(maxFreq,maxBin,numBanks);
% 
% spec = repmat(spec(:),1,numBanks);
% fbank = log(sum(spec.*filters)+1e-5);
% 
% normFactor = sqrt(2/numBanks);
% mfcc = sum(repmat(fbank.',1,MFCC_ORDER).*cos(2*pi*((0:numBanks-1)'*2+1)*(1:MFCC_ORDER)/MDCT_SIZE))*normFactor;
% 




