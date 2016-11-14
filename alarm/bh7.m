function win = bh7(nfft)

% use 7-term blackman harris
a0=2.712203605850388e-001;
a1=4.334446123274422e-001;
a2=2.180041228929303e-001;
a3=6.578534329560609e-002;
a4=1.076186730534183e-002;
a5=7.700127105808265e-004;
a6=1.368088305992921e-005;

k = (0:nfft-1)/nfft;
k = k';

win = a0; 
win = win - a1*cos(2*pi*k);
win = win + a2*cos(2*2*pi*k);
win = win - a3*cos(3*2*pi*k);
win = win + a4*cos(4*2*pi*k);
win = win - a5*cos(5*2*pi*k);
win = win + a6*cos(6*2*pi*k);

%cpg = sum(win)/nfft;
%y = sum(x.')/m;
%y = (y.') .* win;
%y = abs(fft(y))/cpg;
%semilogy(y(1:nfft/2), 'r');
%y = sum(y.*y) / length(y)

end