function [loop,acc,recall,precision,falsealarm] = getEvalResult(evalFile)

% acc: (positive correct nr + negative correct nr)/nr of all samples
% recall: positive correct nr./nr of all positive samples (sensitivity, or hit rate)
% precision: positive correct nr./nr of all classified positive samples, i.e. (1-false alarm)
% false alarm = nr of false positive/nr of all negative samples
% false alarm is derived according to https://en.wikipedia.org/wiki/Sensitivity_and_specificity

fid = fopen(evalFile);
fgetl(fid);
fgetl(fid);
t = fgetl(fid);
idx = 0;
while (1)
    t = fgetl(fid);
    tmp = sscanf(t,'%d %f %f %f');
    if strcmp(t(1),'-')
        break;
    end
    idx = idx+1;
    loop(idx) = tmp(1);
    acc(idx) = tmp(2);
    recall(idx) = tmp(3);
    precision(idx) = tmp(4);
    accp = acc(idx)/100;
    recallp = recall(idx)/100;
    precisionp = precision(idx)/100;
    falsealarm(idx) = 100*(((1-precisionp)/precisionp)/((1/precisionp+(1-recallp)/recallp-1/accp)/(1/accp-1)+(1-precisionp)/precisionp));
end
fclose(fid);



