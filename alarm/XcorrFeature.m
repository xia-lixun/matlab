function [TempMat] = XcorrFeature(x, nframe, nshift, fs)


uERB = x;
TempMat = [];

% feature alaysis
for ibank = 1:1%nbanks
%ibank = 16;

    % generate sub-band auto-correlation (ACF)
    Temp = [];
    
    T0 = 0; T0 = floor(T0 * fs)+1;
    T1 = size(uERB, 1);
    
    t = T0;
    while t + nframe < T1
        
        %rxx = xcorr(uERB(t:t+nframe-1, ibank), 'coeff');
        %rxx = rxx(nframe:end);
        
        v = fft(uERB(t:t+nframe-1, ibank));
        v(1) = 0;
        rxx = real(ifft(v.*conj(v)));
        rxx = rxx ./ rxx(1);
        
        Temp = [Temp rxx(2:end)];
        
        %Ruu = [Ruu rxx];    %Ruu is of size nframe-1 x nobsrvs
        t = t + nshift;
    end
    
    TempMat = [TempMat; Temp];
end


end