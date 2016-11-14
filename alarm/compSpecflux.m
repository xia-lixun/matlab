function flux = compSpecflux(spec,prespec,varargin)
% default:weighted specflux used in ACA
% an alternative calculation given in paper:
% "automated speech/other discrimination for loudness monitoring"
% spec is magnitude spectral, not power spectra or complex spectra

N = size(spec,1);
M = size(spec,2);

if isempty(varargin)
    opt = 1;emph = 1;
else
    opt = varargin{1};emph = varargin{2};
end

if opt == 1
    % magnitude sum
    en = sum(abs(spec));
    preen = sum(abs(prespec));
else
    % power sum
    en = sum(abs(spec).^2);
    preen = sum(abs(prespec).^2);
end
if emph == 1
    preEmph = preEmphVec(N);
else
    preEmph = ones(N,1);
end

weight = (en+preen)/2/N+eps;

if opt == 1
    flux = sum(abs((spec-prespec).*repmat(preEmph,1,M)))./weight/(N*2);    
else
    flux = sqrt(sum(abs((spec-prespec).*repmat(preEmph,1,M)).^2)./weight/(N*2));
end