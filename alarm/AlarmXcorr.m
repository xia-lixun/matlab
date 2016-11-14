% Evaluation of alarm detection based on correlation method
function fea_sm = AlarmXcorr(x, fs, nframe)

% II. preparing parameters
%nframe = 512;
nshift = nframe/2;

% bandpass filter for alarm specific frequency region
order = 6;
[b1,a1] = butter(order, [430  1330]/(fs/2), 'bandpass');
[b2,a2] = butter(order, [500  1400]/(fs/2), 'bandpass');
[b3,a3] = butter(order, [500  1700]/(fs/2), 'bandpass');
[b4,a4] = butter(order, [600  1300]/(fs/2), 'bandpass');
[b5,a5] = butter(order, [700  1000]/(fs/2), 'bandpass');

[b6,a6] = butter(order, [700  1600]/(fs/2), 'bandpass');
[b7,a7] = butter(order, [700  1700]/(fs/2), 'bandpass');
[b8,a8] = butter(order, [900  1200]/(fs/2), 'bandpass');
[b9,a9] = butter(order, [1000 1400]/(fs/2), 'bandpass');

[b10,a10] = butter(order, [1100 1500]/(fs/2), 'bandpass');
[b11,a11] = butter(order, [1300 2700]/(fs/2), 'bandpass');
[b12,a12] = butter(order, [1500 1900]/(fs/2), 'bandpass');
[b13,a13] = butter(order, [1600 2000]/(fs/2), 'bandpass');
[b14,a14] = butter(order, [1600 3200]/(fs/2), 'bandpass');
% [b15,a15] = butter(order, [2000 2400]/(fs/2), 'bandpass');
% [b16,a16] = butter(order, [2000 2500]/(fs/2), 'bandpass');
% [b17,a17] = butter(order, [2100 2500]/(fs/2), 'bandpass');
[b15,a15] = butter(order, [2050 2450]/(fs/2), 'bandpass');
[b16,a16] = butter(order, [2100 2700]/(fs/2), 'bandpass'); 
%[b19,a19] = butter(order, [3100 3400]/(fs/2), 'bandpass');


% add rules for horns
[b17,a17] = butter(order, [1875 2400]/(fs/2), 'bandpass');
[b18,a18] = butter(order, [2400 2800]/(fs/2), 'bandpass');
[b19,a19] = butter(order, [2800 3200]/(fs/2), 'bandpass');
[b20,a20] = butter(order, [3200 3600]/(fs/2), 'bandpass');
[b21,a21] = butter(order, [3500 3700]/(fs/2), 'bandpass');
[b22,a22] = butter(order, [3600 4000]/(fs/2), 'bandpass');  % 4000
[b23,a23] = butter(order, [4000 4400]/(fs/2), 'bandpass'); 



y1 = filter(b1,a1,x);
y2 = filter(b2,a2,x);
y3 = filter(b3,a3,x);

y4 = filter(b4,a4,x);
y5 = filter(b5,a5,x);
y6 = filter(b6,a6,x);

y7 = filter(b7,a7,x);
y8 = filter(b8,a8,x);
y9 = filter(b9,a9,x);

y10 = filter(b10,a10,x);
y11 = filter(b11,a11,x);
y12 = filter(b12,a12,x);

% add rules for horns
y13 = filter(b13,a13,x);
y14 = filter(b14,a14,x);
y15 = filter(b15,a15,x);
y16 = filter(b16,a16,x);
y17 = filter(b17,a17,x);
y18 = filter(b18,a18,x);
y19 = filter(b19,a19,x);
y20 = filter(b20,a20,x);
y21 = filter(b21,a21,x);
y22 = filter(b22,a22,x);
y23 = filter(b23,a23,x);

%frame cross correlation
% for each frame shift event, we divide the sliding length in to M portions
% each portion we calculate the cross correlation between the frames
% each frame is 512 samples, for fs=16000 we can detect periodicity of tone
% down to fs/nframe = 31.25Hz. Is that compliant to ISO7731?
fea1 = acorr_fea(y1, nframe, nshift);
fea2 = acorr_fea(y2, nframe, nshift);
fea3 = acorr_fea(y3, nframe, nshift);

fea4 = acorr_fea(y4, nframe, nshift);
fea5 = acorr_fea(y5, nframe, nshift);
fea6 = acorr_fea(y6, nframe, nshift);

fea7 = acorr_fea(y7, nframe, nshift);
fea8 = acorr_fea(y8, nframe, nshift);
fea9 = acorr_fea(y9, nframe, nshift);

fea10 = acorr_fea(y10, nframe, nshift);
fea11 = acorr_fea(y11, nframe, nshift);
fea12 = acorr_fea(y12, nframe, nshift);

% add for horns
fea13 = acorr_fea(y13, nframe, nshift);
fea14 = acorr_fea(y14, nframe, nshift);
fea15 = acorr_fea(y15, nframe, nshift);
fea16 = acorr_fea(y16, nframe, nshift);
fea17 = acorr_fea(y17, nframe, nshift);
fea18 = acorr_fea(y18, nframe, nshift);
fea19 = acorr_fea(y19, nframe, nshift);
fea20 = acorr_fea(y20, nframe, nshift);
fea21 = acorr_fea(y21, nframe, nshift);
fea22 = acorr_fea(y22, nframe, nshift);
fea23 = acorr_fea(y23, nframe, nshift);

fea = (fea1.*fea2.*fea3.*fea4.*fea5.*fea6.*fea7.*fea8.*fea9.*fea10.*fea11.*fea12.*fea13.*fea14.*fea15.*fea16.*fea17.*fea18.*fea19.*fea20.*fea21.*fea22.*fea23).^(1/23);

% exponential filtering the feature
tc = 0.016*18;
alpha = 1-exp(-nshift/fs/tc);
fea_sm = zeros(size(fea));
for i = 2:length(fea)
    fea_sm(i) = fea(i-1) * alpha + (1 - alpha) * fea_sm(i-1);
end
fea_sm = fea_sm';



