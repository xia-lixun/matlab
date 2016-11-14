% fpr: false positive rate, as x-axis
% tpr: true positive rate, as y-axis
% idx: total number of samples, dimension
% opt: optimimal probability threshold
function [fpr, tpr, idx, opt] = getROC2(evalFile)


tstart = tic;
fid = fopen(evalFile);
fgetl(fid);
fgetl(fid);
t = fgetl(fid);

idx = 0;
while (1)
    t = fgetl(fid);
    if strcmp(t(1),'-')
        break;
    end
    idx = idx+1;
end

t = fgetl(fid); % flush the header line
idx = 0;
while(1)
    t = fgetl(fid);
    if ~ischar(t)
        break;
    end
 
    tmp = sscanf(t,'%d %d %f %f');
    idx = idx + 1;
    
    label(idx) = tmp(2);
    prob(idx) = tmp(4);
    
end
fclose(fid);
time = toc(tstart);
disp('file read'); time



% map probability of {-1, 1} to {1}
prob(label==-1) = 1 - prob(label==-1);


% sort the array
tstart = tic;
[probAscend, reorder] = sort(prob);
time = toc(tstart);

label = label(reorder);
label = label > 0;
p= sum(label);
n = numel(label) - p;
disp('data sort'); time

% false positive rate = false_positive / (false_positive + true_negative)
% true positive rate = true_positive / (true_positive + false_negative)

%tp = [cumsum(label, 'reverse') 0];
%fp = [cumsum(not(label), 'reverse') 0];
tstart = tic;
tp = [fliplr(cumsum(label(end:-1:1))) 0]; 
fp = [fliplr(cumsum(not(label(end:-1:1)))) 0];
fpr = fp / n;
tpr = tp / p;
time = toc(tstart);
disp('ROC calculation'); time


% plot the curve and the optimal point
figure; hold on;
plot(fpr, tpr);

[pmax, I] = max((1-fpr).*tpr);
plot(fpr(I), tpr(I), 's');
grid on;
xlabel('False Positive Rate (False Alarm Rate) - x100 %');
ylabel('True Positive Rate (Sensitivity) - x100 %');
title(strcat('ROC of ', evalFile));
opt = probAscend(I);
