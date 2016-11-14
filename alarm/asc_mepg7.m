function [asc, ass] = asc_mepg7(pxx)

NFFT2 = length(pxx);

fk = zeros(size(pxx));
for temp = 1:NFFT2;
    fk(temp) = temp;
end
fk = (fk-1)./NFFT2 * 8000;
%fk = log2(fk ./ 1000);

asc = sum(fk .* pxx) / sum(pxx);


ass = sum(((fk-asc).^2) .* pxx) / sum(pxx);
ass = sqrt(ass);

end