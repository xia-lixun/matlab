function [prob,clsLab,estClsLab,boostScore,acc,recall,precis] = boostEval(modelFile,feaClsFile)
% used to batch benchmark the training set and test set.
% boost model interpretation
% sigmoid function fitting parameters sigA sigB
% loop, index (started from 0), dir, alpha, threshold, range
% rewrite CRealBoost::Eval in boost.sln in Matlab here
trainData = dlmread(feaClsFile); 
clsLab = trainData(:,1);
fea = trainData(:,2:end);
model = loadModel(modelFile);

% results using different number of loops
acc = zeros(size(model.feaIdx,1),1);
recall = zeros(size(acc)); 
precis = zeros(size(acc)); 
nPos = sum(clsLab>0);       % postive samples in the data
nSample = size(clsLab,1);
nLoop = size(model.feaIdx,1);
strongHypo = zeros(nSample,1);
for n = 1:nLoop
    n
    % weak Hypo, one feature per time
    weakHypo = (fea(:,model.feaIdx(n))-model.thld(n))*model.dir(n)/model.range(n);
    weakHypo(weakHypo>1) = 1.0;
    weakHypo(weakHypo<-1) = -1.0;
    % accumulate weak Hypothesis
    strongHypo = strongHypo+sum(weakHypo,2)*model.alpha(n);
    predPos = sum(strongHypo>0); % positive prediction
    pt = sum(strongHypo>0 & clsLab>0); % positive correct
    pf = sum(strongHypo<0 & clsLab>0); % positive incorrect
    nt = sum(strongHypo<0 & clsLab<0); % negative correct
    nf = sum(strongHypo>0 & clsLab<0); % negative incorrect
    acc(n) = 100*(pt+nt)/nSample;
    recall(n) = 100*pt/nPos;
    precis(n) = 100*pt/predPos;
end

% results for the last loop (nLoop) using sigmod function
prob = 1./(1+exp(model.sigA*strongHypo+model.sigB));
prob = abs((prob-0.5)*2); %?? in boost.ext test mode but inconsistent with aca.exe classify mode, where it uses (1-prob)
estClsLab = ((strongHypo>0)-0.5)*2;
boostScore = strongHypo;

