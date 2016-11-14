function compareFeat(fnum)
mat = dlmread('fea.txt');
vst = dlmread('D:\VST\VADVST\vstfea.txt');
% vst = dlmread('D:\\AudioMulch 2.2.4\\vstfea.txt');
len = length(mat);
numfeats = 76;

matdat = mat(fnum:numfeats:end);
vstdat = vst(fnum:numfeats:end);
len2 = min(length(matdat), length(vstdat));
subplot(211), plot(1:len2, matdat(1:len2), 1:len2, vstdat(1:len2)), legend('mat', 'vst');
subplot(212), plot(matdat(1:len2)-vstdat(1:len2));
end