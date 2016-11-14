function [newMean, newStd] = compStatistics(Data, w, oldMean, oldStd)
% Statistics online calculation

% Update Mean
newMean = (oldMean + w.*Data)./(1+w);

% Updated Std
newStd = sqrt((oldStd.^2+abs(oldMean-newMean).^2 + w.*abs(Data-newMean).^2) ./ (1+w));
