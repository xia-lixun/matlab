clear all
close all

blk = 256;
nfft = 512;
%[x,fs] = audioread('../../Data/VAD/speechmixtrainnoise.wav');
[x48,fs] = audioread('../../Data/VAD/speechmixtrainnoise48.wav');
%[x48,fs] = audioread('../../Data/VAD/speechmixtrainnoise44.wav');

%% if decimating, uncomment the following
filt = dlmread('../../Data/VAD/dec3filt.dat');
nh = length(filt);
x48 = [zeros(nh-1, 1); x48];
rat = 3;
lenin = length(x48);
lenout = lenin/rat;
x = zeros(lenout,1);
for i=1:lenout-nh
    inp = x48((i-1)*rat+1:(i-1)*rat+nh);
    x(i) = sum(inp .* filt);
end
fs = fs/3;
%%
[feaMat,confidence,time] = FeaExtClsCore(x,fs,'../../Data/VAD/fileLabel_Train_1000.model', blk, nfft);
p = [];
vad = zeros(size(confidence));
for k = 1:size(confidence,1)
    [p,vad(k)] = compTransmitVAD(p,confidence(k), fs, blk);
end
dlmwrite('smvad.txt', vad', 'precision', '%1.8f');

% figure;plot(time,vad,'color',[43 129 86]/255);;hold on;grid on
% plot(time,confidence,'r')
% plot([1:length(x)]/fs,x);
% xlabel('Time [sec]'); 
% legend(gca,'signal','vad','confidence of speech');
% axis([0 90 -1 1.5])
