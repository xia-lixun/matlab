function []=writeBandEngFiles(nFFT,fs,foldername)

load bandFreq2
load BandTilt

nBands2 = length(bandFreq2);
halfNFFT = nFFT/2;
% calculate banding matrix from 100 Hz up to fs/2, power preserved
[tmp,bandFreq21,~,~] = filtbankm(nBands2, nFFT,fs,60,fs/2,'ey');

bandTilt = interp1(bandFreq2,BandTilt,bandFreq21);
fid = fopen([foldername 'bandEng_bandFreq2.dat'], 'w');
fprintf(fid, '%1.8ff,\n', bandFreq21);
fclose(fid);
fid = fopen([foldername 'bandEng_bandTilt.dat'], 'w');
fprintf(fid, '%1.8ff,\n', bandTilt);
fclose(fid);
clear bandFreq2
clear BandTilt

bandMatrix2         = tmp(:,1:halfNFFT); 
bandMatrix2 = full(bandMatrix2);
fid = fopen([foldername 'bandEng_bandMatrix2.dat'], 'w');
for i=1:nBands2
fprintf(fid, '%1.8ff,\n', bandMatrix2(i,:));
end
fclose(fid);

end