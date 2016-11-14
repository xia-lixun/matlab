function [p, Flux] = compBandFlux(p, Bands)
% Computes the difference between selected bands from the current and 
% previous frame and the result is smoothed. The FluxFloor tracks the 
% minimum. Returned result is flux in dB. 
%
% INPUTS
%  - Bands     : Band power 
%
% OUTPUTS
%  - p.Flux    : (dB) Spectral flux

flux = sqrt(sum(abs(Bands(p.FluxBands) - p.preBands)));
p.preBands = Bands(p.FluxBands);
p.FluxSmooth = flux*p.FluxSmoothAlpha + (1-p.FluxSmoothAlpha)*p.FluxSmooth;
if flux < p.FluxFloor
   p.FluxFloor = flux*p.FluxFloorDownAlpha + (1-p.FluxFloorDownAlpha)*p.FluxFloor;
else
   p.FluxFloor = flux*p.FluxFloorUpAlpha + (1-p.FluxFloorUpAlpha)*p.FluxFloor;
end

% Don't let flux floor get higher than smoothed flux (can happen due to
% slow movement of floor). 
if (p.FluxFloor > p.FluxSmooth), p.FluxFloor = p.FluxSmooth; end

% Spectral flux is the distance of the smoothed flux from flux floor
p.FluxDist = 20*log10(p.FluxSmooth+eps) - 20*log10(p.FluxFloor+eps);
p.Flux = (min(max(p.FluxThreshold,p.FluxDist),p.FluxMax) - p.FluxThreshold)/(p.FluxMax-p.FluxThreshold);

Flux = p.Flux;

end