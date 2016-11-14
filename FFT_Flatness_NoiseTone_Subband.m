close all;
N_FFT=1024;
% hann_window = sqrt(hann(N_FFT));
% hann_window([1:N_FFT/2]+N_FFT/2) = sqrt(   1-hann_window(1:N_FFT/2).^2  );
fft_w=hamming(N_FFT,'periodic');
w_scaling=sum(fft_w);

% Subband Parameters
band_width=16; 
band_i=3:band_width:N_FFT/2-band_width;


f=20:0.4:1000;
flatness=zeros(1,length(band_i));


% smooth parameter for noise
attack   = 1 - exp( -2.2 *1000 * 512 / (1 * 44100   ) );
release = 1 - exp( -2.2 *1000 * 512 / (1000 * 44100) );
noist_T=1.5;

% non-window


% Windowed

fa=0.3*sin(2*pi*1700*(0:1/44100:(N_FFT-1)/44100));
% fa=fa+0.1*sin(2*pi*(f(i)*1.1)*(0:1/44100:(N_FFT-1)/44100));

    NOISE  = zeros(noist_T/0.1,N_FFT); % Tmax/Release
    for in = 1: noist_T/0.1-1  % Tmax/Release

    noise=rand(1,N_FFT);

    bh=fft(noise.*fft_w')/w_scaling;
    NOISE(in,:)=bh;
    end

    for ib=1:N_FFT/2
        for in = 1: noist_T/0.1-1 
             if  abs(NOISE(in,ib))  >NOISE(end,ib)
                 T=attack;
             else
                 T=release;
             end

             NOISE(end,ib)  = (1-T)*NOISE(end,ib)  + T*abs(NOISE(in,ib)) ; 

        end

    end    
    
bh=fft(fa.*fft_w')/w_scaling;

% mix the FFT of noise and tone
bh= bh+ NOISE(end,:);

for i=band_i;
P_bh=abs(bh(i:i+band_width-1)).^2;

gM=sum(10*log10(P_bh))/(band_width);
aM=10*log10(mean(P_bh));

flatness((i+band_width-3)/band_width)=gM-aM;


figure (1), plot(10*log10(P_bh));
ylim([-90, 0])
close 1;

end

figure,
grid on, hold on;
plot(flatness);
% figure,
% plot(10*log10(P_bh));
% 
% 
% fa=0.1*sin(2*pi*f(i)*(0:1/44100:(N_FFT-1)/44100));
% bh=fft(fa.*fft_w');
% Pbh=20*log10(2*abs(bh)/w_scaling);
% i_tone=find(Pbh(1:end/2)==max(Pbh));
% 
% e_tone=2*sum(abs(bh(i_tone-1:i_tone+1)))/w_scaling;
% E_TONE=20*log10(e_tone);
% 
% e_lobe=2*abs(bh([i_tone-2 i_tone+2]))/w_scaling;
% crest=E_TONE-20*log10(e_lobe);
% crest_1(i)=min(crest);
% 
% e_tone=2*sum(abs(bh(i_tone)))/w_scaling;
% E_TONE=20*log10(e_tone);
% 
% e_lobe=2*abs(bh([i_tone-1 i_tone+1]))/w_scaling;
% crest=E_TONE-20*log10(e_lobe);
% crest_2(i)=min(crest);
% end
