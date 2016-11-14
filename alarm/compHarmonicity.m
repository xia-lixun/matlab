function [F0, Strength] = compHarmonicity(Bins, Freq, bnorm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARMONICITY
% Estimates the harmonic strength of an audio frame ranging from
% below 0 for an inharmonic sound (e.g. white noise), to 1 for a
% perfectly harmonic sound.
%
% Input Bins in columns, must be a spectrum using windowing to prevent
% leakage!!
%
% HarmonicityStrength
% Computing the harmonic strength begins by finding potential F0 candidates
% from the peaks in the frequency spectrum (Bins) between p.F0min and
% p.F0max and above the average power of the spectrum. Refined potential F0
% candidates are then found by quadratically interpolating the peaks, from
% which binary masks are formed around the harmonic locations up to
% p.maxHarmonics. The harmonic strength is calculated by summing the
% spectrum at the regions defined by the binary mask and subtracting the
% summed regions outside the mask up. The value is then normalised by
% the total power in the up to frequency defined by p.maxHarmonics. This
% is calculated for each of the peaks and p.Strength is set to the highest
% harmonic strength value, and the corresponding peak set to p.F0. The
% harmonic Strength of the signal ranges from below 0 for an inharmonic
% sound (e.g. white noise), to above 0.5 for harmonic sounds.

% indicate harmonicity should be normalized, default = 'normalized'
if nargin == 2
    bnorm = 1;
end

p.maxHarmonics = 100; % Max # harmonics to use in sum, 20 ?? or maybe determined by maxFreq
p.maxFreq = Freq(end);
p.F0min = 70;         % (Hz) Min F0 for peak search
p.F0max = 310;        % (Hz) Max F0 for peak search
p.minPower = -50;     % (dB) Min power before calc strength
p.halfMaskWidth = 8;  % Mask width at each harmonic is 2*F0/p.halfMaskWidth

lFrames = size(Bins,2);
df = Freq(2)-Freq(1); % freq resolution
[~,F0BinMin] = min(abs(p.F0min-Freq));
[~,F0BinMax] = min(abs(p.F0max-Freq));
BinMax = find(Freq<=p.maxFreq,1,'last');
Bins = abs(Bins);


%% F0 estimation
% PEAK PICKING
% include one point before and afterwards to identify peaks at edges
% please change peaks = find(peakmap(:,m))+F0BinMin-1; to peaks =
% find(peakmap(:,m))+F0BinMin-2; in Line 55
% F0Bins = Bins(max(F0BinMin-1,0):min(F0BinMax+1,lBins),:);
% F0Bins2 = Bins(max(F0BinMin,0):min(F0BinMax,lBins),:);
% peakmap = Maxima(F0Bins,mean(F0Bins2));
F0Bins = Bins(F0BinMin:F0BinMax,:);
peakmap = Maxima(F0Bins,mean(F0Bins));
Strength = zeros(lFrames,1);
F0 = zeros(lFrames,1);

for m = 1:lFrames
    peaks = find(peakmap(:,m))+F0BinMin-1;
    % lowPower = (10*log10(sum(Bins(:,m))) < p.minPower);
    lowPower = 0;
    %% Mask calculations
    if isempty(peaks) || lowPower
        % No peaks or
        Strength(m) = 0;
        F0(m) = 0;
    else
        strength = zeros(size(peaks));
        for k = 1:length(peaks)
            pk = peaks(k);
            
            % Form harmonic binary mask
            lHarmonics = min(p.maxHarmonics,floor(BinMax/pk)); % no. of harmonics
            halfMask = max(0,round(pk/p.halfMaskWidth));
            lM = min(BinMax, round(lHarmonics*pk+halfMask));
            M = zeros(lM,1);
            for n = 1:lHarmonics
                    % FIXED MASK
                    M(max(1,round(n*pk-halfMask)):min(BinMax,round(n*pk+halfMask))) = 1;
            end
            
            % Reward harmonicity, punish inharmonic regions
            % strength(k) = (sum(M(F0BinMin:lM).*Bins(F0BinMin:lM,m)) - sum(~M(F0BinMin:lM).*Bins(F0BinMin:lM,m)))/sum(Bins(F0BinMin:lM,m));
            % strength(k) = (sum((M(F0BinMin:lM).*Bins(F0BinMin:lM,m)).^2) - sum((~M(F0BinMin:lM).*Bins(F0BinMin:lM,m)).^2))/sum(Bins(F0BinMin:lM,m).^2);
            % strength(k) = (sum(M(F0BinMin:lM).*Bins(F0BinMin:lM))) - (sum(~M(F0BinMin:lM).*Bins(F0BinMin:lM)));
            % strength(k) = (sum(M(F0BinMin:lM).*Bins(F0BinMin:lM,m)));
            
            % feature strength(k) = (sum(M(F0BinMin:lM).*Bins(F0BinMin:lM,m))) 
            % is more discriminative for voiced and tonal noise because
            % voiced part usually stands out of noise but this feature is dependent 
            % on total energy and therefore needs AGC! If there is AGC,
            % please use this feature instead of normalized Harmonicity as
            % below.
            if bnorm == 1
                strength(k) = (sum((M(F0BinMin:lM).*Bins(F0BinMin:lM,m)).^2))/sum(Bins(F0BinMin:lM,m).^2);
            else
                strength(k) = sum((M(F0BinMin:lM).*Bins(F0BinMin:lM,m)));
            end
            peaks(k) = pk;
        end
        
        % Get highest harmonic strength
        [Strength(m),ind] = max(strength);
        F0(m) = (peaks(ind)-1)*df + Freq(1);
        
    end
end

%% Helper Functions

function peaks = Maxima(x, thresh)
%function peaks = Maxima(x, [thresh])
% Finds local maxima of vector x for points greater than absolute
% threshhold thresh and outputs binary col vector peaks. thresh default is
% zero.

if (nargin < 2)
    thresh = zeros(1,size(x,2));
end

lx = size(x,1);

% leading edge
lead = (x > [x(1,:); x(1:(lx-1),:)]);

% falling edge - note >= sign to catches plateaus
fall = (x >= [x(2:lx,:); 1+x(lx,:)]);

peaks = (x>repmat(thresh,lx,1)).*lead.*fall;



