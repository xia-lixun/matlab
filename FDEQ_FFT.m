% input: audio_in
% output:   Y_State

N_FFT=1024;
hann_window = sqrt(hann(N_FFT));
hann_window([1:N_FFT/2]+N_FFT/2) = sqrt(   1-hann_window(1:N_FFT/2).^2  );

STFT_step=N_FFT/2;

Y_State=zeros(N_FFT/2+1, fix(length(audio_in)/STFT_step)-2, size(audio_in,2));

% Windowed STFT for music @ reference volume
out_last=zeros(STFT_step,2);
for i=1: length(audio_in)/STFT_step-2       % loop in frames
    ind= [1:STFT_step]' + (i-1)*STFT_step;
 
    for ch=1:size(audio_in,2)
        in_overlap=[audio_in(ind,ch); audio_in(ind+STFT_step,ch)].*hann_window;
        Y=fft(in_overlap,N_FFT);
        Y_State(:,i,ch)=Y(1:end/2+1);

    end
    
end



