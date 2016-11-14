function [harm,pitch,interp_amplitude]=myFuncSHR(amplitude,SHR)
frequency_idx = SHR.frequency_idx;
logf = SHR.logf;
interp_logf = SHR.interp_logf;
interp_amplitude_weight = SHR.interp_amplitude_weight;
lowerbound = SHR.lowerbound;
upperbound = SHR.upperbound;
N = SHR.N;
shift_distance = SHR.shift_distance;
numPitch = SHR.numPitch;
%%%% NOTE: the amplitude input is assume to be FFT_Len/2, i.e. without the DC component. If different, please modify InitSHR2 accordingly
amplitude = amplitude(2:frequency_idx(end)+1);

if (SHR.LOG_SPECTRUM_ON)
amplitude(amplitude<=0)=7.6294e-006;
%amplitude = amplitude ./max(amplitude);
amplitude=single(20*log10(amplitude)); 
end
%amplitude(1:5) = 0;
interp_amplitude=single(interp1(logf,amplitude,interp_logf,'linear'));
%interp_amplitude = interp_amplitude.*interp_amplitude_weight;  % weight the spectrum
interp_amplitude=interp_amplitude-min(interp_amplitude);  % normalize, reduce the impact of low amplitude frequencies

%% spectrum sum and diff
% init shsodd and shseven
if 0
shsodd = interp_amplitude;
shseven = [interp_amplitude(shift_distance(1)+1:end) zeros(1,shift_distance(1))];
for i = 3:2:N
   shsodd = shsodd + [interp_amplitude(shift_distance(i-1)+1:end) zeros(1,shift_distance(i-1))];
   shseven = shseven + [interp_amplitude(shift_distance(i)+1:end) zeros(1,shift_distance(i))];
end
difference=shseven-shsodd;
if 0
figure(7)
subplot(4,1,1),plot(interp_amplitude)
subplot(4,1,2),plot(shsodd)
subplot(4,1,3),plot(shseven)
subplot(4,1,4),plot(difference)
end
interp_amplitude=interp_amplitude-min(interp_amplitude);
shsodd = interp_amplitude;
shseven = [interp_amplitude(shift_distance(1)+1:end) zeros(1,shift_distance(1))];
for i = 3:2:N
   shsodd = shsodd + [interp_amplitude(shift_distance(i-1)+1:end) zeros(1,shift_distance(i-1))];
   shseven = shseven + [interp_amplitude(shift_distance(i)+1:end) zeros(1,shift_distance(i))];
end
difference=shseven-shsodd;
figure(8)
subplot(4,1,1),plot(interp_amplitude)
subplot(4,1,2),plot(shsodd)
subplot(4,1,3),plot(shseven)
subplot(4,1,4),plot(difference)
pause
end


% more efficient implementation, but less readable
difference = [interp_amplitude(shift_distance(1)+1:end) zeros(1,shift_distance(1))] - interp_amplitude;
for i = 3:2:N
    difference = difference + [interp_amplitude(shift_distance(i)+1:end) zeros(1,shift_distance(i))] - [interp_amplitude(shift_distance(i-1)+1:end) zeros(1,shift_distance(i-1))];
end
%difference - difference1,pause
if 0
figure(8)
subplot(2,1,1),plot(interp_amplitude)
subplot(2,1,2),plot(difference)
pause
end

% peak picking process
if 0
%[harm,pitch]=max(difference(lowerbound:upperbound)); % only find two maxima
[harm] = multi_harm(difference,interp_logf,interp_amplitude, numPitch);
harm=max(harm/N,7.6294e-006) % a bit normalization
%HARMnorm=HARMraw/sum(interp_amplitude);
pitch =1;% 2^interp_logf(pitch+lowerbound-1)*2;
else
%%% multi-pitch tracking
[harm,pitch_idx]=multi_pitch(difference(lowerbound:upperbound),numPitch,SHR.bins_per_octave); % only find two maxima
harm=max(harm/N,7.6294e-006); % a bit normalization
%HARMnorm=HARMraw/sum(interp_amplitude);
pitch = 2.^interp_logf(pitch_idx+lowerbound-1)*2;
%plot(difference(lowerbound:upperbound))
%vline(pitch_idx)
%[harm,pitch_idx pitch],pause

end
    
if 0
    [idx;pitch]
figure(3)
subplot(2,1,1),plot(interp_amplitude)
subplot(2,1,2),plot(difference/N), hold, 
plot([idx(1) idx(1)],[min(difference) max(difference)],'r'), 
plot([idx(2) idx(2)],[min(difference) max(difference)],'k'),hold off,pause
end
    
end

function [harm,pitch] = multi_pitch(difference, N,ob)
% seach multiple pitches 
% difference: the input SHR spectrum
% N: number of pitch candidates
pitch = zeros(1,N);
harm = zeros(1,N);
[harm(1), pitch(1)] = max(difference);
    % remove the 1/4f0 peak
    s = max(1,pitch(1)-ob-5);
    %s = 1;
    e = min(pitch(1)-ob+5,length(difference));
difference(s:e) = min(difference);
for n = 2:N    
    
%    plot(difference),pause
    s = max(1,pitch(n-1)-5);
    %s = 1;
    e = min(pitch(n-1)+5,length(difference));
    difference(s:e) = min(difference);
    [harm(n), pitch(n)] = max(difference);
end
%pitch(1)-pitch(2)
end
function [harm] = multi_harm(difference,interp_logf,interp_amplitude, N)
% once the first peak is found at f0, find the second peak at the high
% harmonics, 2f0, 3f0, etc
% difference: the input SHR spectrum
% N: number of pitch candidates
pitch = zeros(1,N);
harm = zeros(1,N);
[harm(1), idx1] = max(difference);
pitch1 = 2^interp_logf(idx1)*2;
[harmt, idx2] = max(difference(idx1+10:end));
idx2=idx2+idx1+10-1;
%pitch2 = 2^interp_logf(idx2+idx1+10-1)*2;
%[pitch1 pitch2 ]
%pitch2/pitch1
tmp = interp_logf(idx1)+log2(2);
idx = find(interp_logf>=tmp);
idx=idx(1);

harm(2)=max(difference(idx-3:idx+3));
if 0
subplot(2,1,1),plot(interp_logf,interp_amplitude)
subplot(2,1,2),plot(interp_logf,difference)
vline([interp_logf(idx1) interp_logf(idx2)]),pause
end
end