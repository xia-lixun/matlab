%% test the function of boost 
% using a extracted feature vector with class labels(1/-1)
% after the training the result can be output to " -o " file
% the result can also be obtained by rerun the boost using " -t " option.
% the two results will be slightly different because the latter reads the
% data from file which has been truncated to %.3f.
% boost32 and boost64 are the same code compiled at different bit width, so
% the results will be slightly different.
clear all
close all
cmd = '..\boost32bit.exe real -p 1 -l 500 -r 5 -d 4 -i RAW -m speech_versus_others_model.txt -f speech_versus_others_train.txt -o speech_versus_others_train_results.txt';
system(cmd)
cmd = '..\boost32bit.exe real -t -p 1 -l 500 -r 5 -i RAW -m speech_versus_others_model.txt -f speech_versus_others_test.txt -o speech_versus_others_test_results.txt';
system(cmd)

% boost eval, should be idential to training results "*results.txt"
% test the train feature set, should be very high
clear all
close all
modelFile = 'speech_versus_others_model.txt';
feaClsFile = 'speech_versus_others_train.txt';
[prob,clsLab,estClsLab,boostScore,acc,recall,precis] = boostEval(modelFile,feaClsFile);
sum((clsLab-estClsLab)~=0) % number of wrong classification

% test the test feature set, lower than train set scores
clear all
close all
modelFile = 'speech_versus_others_model.txt';
feaClsFile = 'speech_versus_others_test.txt';
[prob,clsLab,estClsLab,boostScore,acc,recall,precis] = boostEval(modelFile,feaClsFile);
sum((clsLab-estClsLab)~=0) % number of wrong classification

% test a wav file assuming the features vectors have been extracted
% acomic.wav is the wav file
% acomic.txt is the classification results from ACA.exe using "classify"
% option. In Matlab, we omit the smoothing in ACA. Therefore, the results
% are raw confidence, slightly different from acomic.txt.
clear all 
close all
% load model
modelFile = 'speech_versus_others_model.txt';
model = loadModel(modelFile);
% load features
feaFile = 'acomic_features.txt';
fea = dlmread(feaFile);
[cls,score,confidence] = boostClassify(fea(:,:),model);

%% self generated example
clear all
close all

rng(0,'twister') % for reproducibility
fea = rand(2000,20);
clsLab = ((sum(fea(:,1:5),2) > 2.5)-0.5)*2;

trainFea = [clsLab,fea];
dlmwrite('example_train.txt',trainFea,'delimiter',' ');
cmd = '..\boost32bit.exe real -p 1 -l 500 -r 5 -d 4 -i RAW -m example_model.txt -f example_train.txt -o example_train_results.txt';
system(cmd)
% cmd = 'boost32bit.exe real -t -p 1 -l 500 -r 5 -i RAW -m example_model.txt -f example_train.txt -o example_train_results.txt';
% system(cmd)

% boost eval, should be idential to training results "*results.txt"
clear all
close all
modelFile = 'example_model.txt';
feaClsFile = 'example_train.txt';
[prob,clsLab,estClsLab,boostScore,acc,recall,precis] = boostEval(modelFile,feaClsFile);
sum((clsLab-estClsLab)~=0)  % number of wrong classification