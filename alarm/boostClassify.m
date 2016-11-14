function [cls,score,confidence] = boostClassify(fea,model)
% used for online classification given a feature vector and a model, and it
% also supports multiple feature vectors input for classification
% rewrite CAudioClassifer::audioClassify in audioClassifier.cpp (audio.sln)
% in Matlab
% boostEval is another version of this function, which also calculates the
% overall accuracy, recall, precision etc.
% input:
% fea: a feature vector 1*N or vectors M*N, where N is the size of feature
% vector
% mode: a struct include the following fields
%       sigA: a scalar
%       sigB: a scalar
%     feaIdx: nLoop*1
%        dir: nLoop*1
%      alpha: nLoop*1
%       thld: nLoop*1
%      range: nLoop*1

nLoop = size(model.feaIdx,1);
nSample = size(fea,1);
strongHypo = zeros(nSample,1, 'single');
for n = 1:nLoop
    % weak Hypo, one feature per time
    weakHypo = (fea(:,model.feaIdx(n))-single(model.thld(n)))*model.dir(n)/(single(model.range(n))+eps);
    weakHypo(weakHypo>1) = 1.0;
    weakHypo(weakHypo<-1) = -1.0;
    % accumulate weak Hypothesis
    strongHypo = strongHypo+sum(weakHypo,2)*single(model.alpha(n));
end

prob = 1./(1+exp(single(model.sigA)*strongHypo+single(model.sigB)));

cls = ((strongHypo>0)-0.5)*2;
score = strongHypo;
% the confidence that the feature set is from the first class
confidence = prob;

% % if you need confidence that the feature set is from the classified class
% confidence(strongHypo<0) = 1-prob(strongHypo<0);