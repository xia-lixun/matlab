f_fft_Bark = 1+(13.*atan(0.76.*f_fft./1000) + 3.5.*atan(((f_fft./1000)./7.5).^2)); 

m=1; BarkBinNum=[];
for n = 1:length(f_fft);
    if fix(f_fft_Bark(n)) > m;
        BarkBinNum(m) = n;
        m=m+1;
        BarkBinNum(m) = N_FFT/2+1;
    end;        
end;
%Extension of the absolute Bark values to 26 numbers:
BarkBinNum = [1 BarkBinNum];
Bin_fBark=(BarkBinNum(1:end-1)+BarkBinNum(2:end))*0.5;
fBark=Bin_fBark*Fs/N_FFT;

%Calculation of the relative Bark values (25 numbers):
BarkBinNumDiff = diff(BarkBinNum);

Y_Power_Mono=sum((abs(Y_State).*fft_scaling*2).^2,3); % Level -> Power -> Sum Power
GroundRef = 1.0e-020; Y_Power_Mono(:)=max(Y_Power_Mono(:),GroundRef);
Bm=zeros(length(BarkBinNum)-1, size(Y_State,2));

% Grouping fft bins into critical bands
% Acoustic Calibration

% MiddleEar Transfer Function
MiddleEar;
MiddleEar_P= 10.^(MiddleEar_Bark*0.1); % Linear Power Domain
HPEQ_P=HPEQ;  % Linear Power Domain
for k = 1: size(Y_Power_Mono,2) % loop in frames
    
        Y=Y_Power_Mono(:,k);
        
        % Acoustic Calibeation: MiddleEar, Headphone
        if HPEQ_EN
          Y=Y.*HPEQ_P.*MiddleEar_P';  
        end
    for n = 1:length(BarkBinNum)-1  % loop in barks
         Bm(n,k) =sum(Y(BarkBinNum(n):BarkBinNum(n+1)-1));
    end
end

Bm_Raw=Bm;
% Apply short-term smoothing
 frameRate = Fs/(0.5*N_FFT);
    [att,rel] = calLoudnessTC(1000,frameRate,0);
% % DEC loudness model parameters
% at_DEC = 100;
% rt_DEC  = 300;
% AT_DEC = 1 - exp( -2.2 *1000 * (N_FFT/2) / (at_DEC * Fs   ) );
% RT_DEC = 1 - exp( -2.2 *1000 * (N_FFT/2) / (at_DEC * Fs   ) );

    
% if strcmp(p_ver,'_LL_DEC_Inst') ~=1   
    for k = 2: size(Y_Power_Mono,2) % loop in frames

    %         % Barks in attack phase
    %         at_id= Bm_Raw(:,k-1) < Bm_Raw(:,k);
    %         Bm(at_id,k)=Bm_Raw(at_id,k)*(1-AT_DEC) + Bm(at_id, k-1)*AT_DEC;
    %         % Barks in release phase
    %         Bm(~at_id,k)=Bm_Raw(~at_id,k)*(1-RT_DEC) + Bm(~at_id, k-1)*RT_DEC;

            % Barks in attack phase
            at_id= Bm_Raw(:,k-1) < Bm_Raw(:,k);
            Bm(at_id,k)=Bm_Raw(at_id,k)*(1-att) + Bm(at_id, k-1)*att;
            % Barks in release phase
            Bm(~at_id,k)=Bm_Raw(~at_id,k)*(1-rel) + Bm(~at_id, k-1)*rel;




    end
% end

% Spread Function
S = zeros(length(BarkBinNum)-1);
for k = 1:(length(BarkBinNum)-1);
    for n = 1:(length(BarkBinNum)-1);
        S(k,n) = 15.81 + ((25-10)/2)*((k-n)+0.474)-((25+10)/2)*(sqrt(1+((k-n)+0.474)^2));
    end;
end;
% Linearization
S = (10.^(S./10));

% Apply spread function to critical band.
Cm=S*Bm;
% Cm=Bm;

% Spectral flatness measurement
SFMm = zeros(length(BarkBinNum)-1,size(Y_Power_Mono,2));

for k=1:size(Y_Power_Mono,2)
    
    for n = 1:length(BarkBinNum)-1
    AMm=mean(Y_Power_Mono(BarkBinNum(n):BarkBinNum(n+1)-1,k));
    GMm_dB = 10*sum(log10(Y_Power_Mono(BarkBinNum(n):BarkBinNum(n+1)-1,k)))./BarkBinNumDiff(n);
    SFMm(n,k) = GMm_dB - 10*log10(AMm); 
    
    end
    
end
SFM_dB_Max = -60;
Om = zeros(length(BarkBinNum)-1,size(Y_Power_Mono,2));
Alpham_M=zeros(length(BarkBinNum)-1,size(Y_Power_Mono,2));
for k = 1:size(Y_Power_Mono,2)
    for n = 1:length(BarkBinNum)-1;
        Alpham = min((SFMm(n,k)/SFM_dB_Max),1);
        %Berechnung der Offsets (=O=Geometrische Gewichtung):
        Om(n,k) = Alpham*(14.5+n) + (1-Alpham)*5.5;
        Alpham_M(n,k)=Alpham;
    end;
end;

% Tm = 10*log10(Cm)-Om;
% Test Disable Offset
Tm = 10*log10(Cm);%-Om; 

%  renormalization of the overall side:
% Ce = ones(size(Tm));
% One = ones(size(Tm,1),1);
% for k = 1:size(Y_Power_Mono,2);
%     Ce(:,k) = 10*log10(S * One);
% end
% Tm = Tm - Ce;

%   Calibration, 
%   HPEQ before spread ? 
%   Senstivity before spread as well, if SPL-dependent
%   spread function is applied
P_Ref = HP_Sensitivity;
% Tm_spl = Tm + 90-P_Ref;
Tm_spl = Tm-HP_Sensitivity;

% Absolute threshold of hearing
M = zeros(1,fix(Fs/2));  %M =Level in [dB], the absolute threshold of hearing (=hearing threshold)
for f = 1:Fs/2;
    M(f) = 3.64*((f/1000)^-0.8) - 6.5*exp(-0.6*((f/1000)-3.3)^2) + 10^-3*((f/1000)^4);
end;

%Absolute hearing threshold within the individual frequency groups:
M_Bark = zeros(length(BarkBinNum)-1,1);
for n = 1:length(BarkBinNum)-1;
    if f_fft(BarkBinNum(n+1)) > Fs/2;
        END = floor(Fs/2);
    else
        END = round(f_fft(BarkBinNum(n+1)));
    end;
    %%Using the absolute minimum of the resting threshold within each frequency group:
    %%(After Robinson)
    % M_Bark (n) = min (M (round (f_Fft (BarkBinNum (n))): END));
    % Using the average values of the absolute resting threshold within each frequency group:
    % (After Johnston)
    M_Bark(n) = mean(M(round(f_fft(BarkBinNum(n)))+1:END));
end;


MaskTh_SPLm = zeros(size(Tm_spl));

for k = 1:size(Y_Power_Mono,2)
    
    MaskTh_SPLm(:,k) = max(Tm_spl(:,k),M_Bark);    
end