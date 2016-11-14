function [p,tx] = compTransmitVAD(p,conf, fs, blk)
% 1, p.VoiceRatio is the smoothed confidence of voice
% 2, p.CountDown imposes a hold time on (p.VoiceRatio>0.5) so that a very
% short pause will be removed
% 3, however, with only p.CountDown, the short burst (sometimes a false alarm) will
% also hold for a long time. We use a long term condition to judge if
% it is within a continuous speech, which requires a removal of short
% pause, or it is just a pure short burst, by using a long-term counter
% p.VoiceCount. If the speech has only lasted for a very short period
% of time, we allow for a fast decay of p.VoiceRatio.

if isempty(p)    
    p.blk = blk;
    p.fs = fs;
    p.blkTime = p.blk/p.fs;
    p.holdTime = 350e-3;
    
    p.VoiceRatio = 0;
    p.VoiceRatioAttTime = 10e-3;
    p.VoiceRatioRelTime = 500e-3;
    p.VoiceAtt = exp(-p.blkTime/p.VoiceRatioAttTime);
    p.VoiceRel = exp(-p.blkTime/p.VoiceRatioRelTime);
    
    p.CountDown = 0;
    p.CountMax = round(p.holdTime/p.blkTime);
    p.dec = 1;
    p.dec_delta = 0.9;
    p.VoiceCount = 0;          % an indication how long the good voice has lasted.
    p.ContinuousVoiceTime = 700e-3;
    p.ContinuousVoice = round(p.ContinuousVoiceTime/p.blkTime);
    p.Decrement = 1;    
end

if (p.VoiceRatio>0.5)
    p.dec = 1;
    p.Decrement = 1;
    p.CountDown = p.CountMax;
    p.VoiceCount = p.VoiceCount + 1;
else
    p.CountDown = max(p.CountDown-p.Decrement,0);
    p.VoiceCount = max(p.VoiceCount-1,0);
    if p.CountDown == 0
        p.VoiceCount = 0;
    end
end

if conf>p.VoiceRatio
    p.VoiceRatio = p.VoiceAtt*p.VoiceRatio + (1-p.VoiceAtt)*conf;
else
    if (p.CountDown == 0 || p.VoiceCount<p.ContinuousVoice)  % p.VoiceCount<p.ContinuousVoice helps to reduce a short busrt vad accuracy in time
        p.dec = p.dec*p.dec_delta;
        p.Decrement = p.Decrement + 1;
    end
    p.VoiceRatio = p.VoiceRel*p.dec*p.VoiceRatio + (1-p.VoiceRel*p.dec)*conf;
end

tx = p.CountDown>0;

%% for datalogging
% persistent a logData
% if (isempty(a))
%     a = 1;
%     logData = [];
% end
% logData.CountDown(a) = p.CountDown;
% logData.VoiceRatio(a) = p.VoiceRatio;
% logData.VoiceCount(a) = p.VoiceCount;
% a = a+1;