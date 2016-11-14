function [Sesf, MaxRatio] = SegEnSpecFlatness(Spec)


NFFT2 = length(Spec);


%F0 = 250:500:8000-500;
%F1 = F0 + 500;

F0 = 187.5:375:8000-375;
F1 = F0 + 375;
Nseg = numel(F0);

% for fs = 16000, we have 15 segments
% estimate energy of each segment first
SegEng = zeros(Nseg, 1);
SegEsf = zeros(Nseg, 1);

F0 = round(F0/8000 * NFFT2);
F1 = round(F1/8000 * NFFT2);

SegEngMax = 0;
for i = 1:Nseg
    
    SegEng(i) = sum(Spec(F0(i):F1(i)).^2);
    if SegEng(i) > SegEngMax
        SegEngMax = SegEng(i);
    end
    
    XkHat = Spec(F0(i):F1(i)) ./ sum(Spec(F0(i):F1(i)));
    Log2XkHat = log2(XkHat);
    ita = -sum(Log2XkHat .* XkHat) / log2(F1(i) - F0(i) + 1);
    SegEsf(i) = 2.^ita - 1;
    
end

Sesf = sum(SegEsf .* SegEng) / sum(SegEng);
MaxRatio = SegEngMax / sum(SegEng);

end