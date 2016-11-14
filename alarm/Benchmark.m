close all; 
clear all; 
clc;

%% delete all files generated
a = dir('./Mixed/Test/*.wav');
    for i = 1:length(a)
        delete(strcat('./Mixed/Test/', a(i).name));
    end


a = dir('./Mixed/Train/*.wav');
    for i = 1:length(a)
        delete(strcat('./Mixed/Train/', a(i).name));
    end

mkdir('./Mixed');
mkdir('./Mixed/Train');
mkdir('./Mixed/Test');

delete('./Feat/Test/*.txt');
delete('./Feat/Train/*.txt');

delete('./Calib/Test/*.txt');
delete('./Calib/Train/*.txt');
delete('./Calib/*.txt');

%% populate alarm signals, noise signals and set SNR 
fs = 16000;
nfea = 76;

alarm = dir('./Alarm/*.wav');
alarm_tot_cnt = length(alarm);
alarm_train_cnt = floor(alarm_tot_cnt * 0.6);
alarm_test_cnt = alarm_tot_cnt - alarm_train_cnt;

alarm_idx = randperm(alarm_tot_cnt);
alarm_train_idx = alarm_idx(1:alarm_train_cnt);    
alarm_test_idx = alarm_idx(alarm_train_cnt+1:end); 
clear alarm_idx;


noise = dir('./Noise/*.wav');
noise_tot_cnt = length(noise);
noise_train_cnt = floor(noise_tot_cnt * 0.6);
noise_test_cnt = noise_tot_cnt - noise_train_cnt;

noise_idx = randperm(noise_tot_cnt);
noise_train_idx = noise_idx(1:noise_train_cnt);
noise_test_idx = noise_idx(noise_train_cnt+1:end);

              

%% sanity check
% generate the training set of the design matrix
% and collect the features on-th-fly
for i = alarm_train_idx
            [x, fsx] = audioread(strcat('./Alarm/', alarm(i).name));
            x = resample(x(:,1), fs, fsx);
    
            mixed = strcat('./Mixed/Train/a', num2str(i), '.wav');
            audiowrite(mixed, x, fs, 'BitsPerSample', 16);
            
            feacmd = strcat({'VAD.exe '}, mixed, {' ./Feat/Train/a'}, num2str(i), '.txt');
            system(feacmd{1});
            
            feamat = load(strcat('./Feat/Train/a', num2str(i), '.txt'));
            nobser = length(feamat)/nfea;
            feamat = reshape(feamat, nfea, nobser)';
            feamat = [ones(nobser, 1) feamat];
            dlmwrite(strcat('./Calib/Train/a', num2str(i), '.txt'), feamat, 'delimiter',' ','precision','%.8f');
end

for i = alarm_test_idx
            [x, fsx] = audioread(strcat('./Alarm/', alarm(i).name));
            x = resample(x(:,1), fs, fsx);
            
            mixed = strcat('./Mixed/Test/a', num2str(i), '.wav');
            audiowrite(mixed, x, fs, 'BitsPerSample', 16);
            
            feacmd = strcat({'VAD.exe '}, mixed, {' ./Feat/Test/a'}, num2str(i), '.txt');
            system(feacmd{1});
            
            feamat = load(strcat('./Feat/Test/a', num2str(i), '.txt'));
            nobser = length(feamat)/nfea;
            feamat = reshape(feamat, nfea, nobser)';
            feamat = [ones(nobser, 1) feamat];
            dlmwrite(strcat('./Calib/Test/a', num2str(i), '.txt'), feamat,  'delimiter',' ','precision','%.8f');
end

for i = noise_train_idx
            [x, fsx] = audioread(strcat('./Noise/', noise(i).name));
            x = resample(x(:,1), fs, fsx);
            
            mixed = strcat('./Mixed/Train/n', num2str(i), '.wav');
            audiowrite(mixed, x, fs, 'BitsPerSample', 16);
            
            feacmd = strcat({'VAD.exe '}, mixed, {' ./Feat/Train/n'}, num2str(i), '.txt');
            system(feacmd{1});
            
            feamat = load(strcat('./Feat/Train/n', num2str(i), '.txt'));
            nobser = length(feamat)/nfea;
            feamat = reshape(feamat, nfea, nobser)';
            feamat = [-ones(nobser, 1) feamat];
            dlmwrite(strcat('./Calib/Train/n', num2str(i), '.txt'), feamat,  'delimiter',' ','precision','%.8f');
end

for i = noise_test_idx
            [x, fsx] = audioread(strcat('./Noise/', noise(i).name));
            x = resample(x(:,1), fs, fsx);
            
            mixed = strcat('./Mixed/Test/n', num2str(i), '.wav');
            audiowrite(mixed, x, fs, 'BitsPerSample', 16);
            
            feacmd = strcat({'VAD.exe '}, mixed, {' ./Feat/Test/n'}, num2str(i), '.txt');
            system(feacmd{1});
            
            feamat = load(strcat('./Feat/Test/n', num2str(i), '.txt'));
            nobser = length(feamat)/nfea;
            feamat = reshape(feamat, nfea, nobser)';
            feamat = [-ones(nobser, 1) feamat];
            dlmwrite(strcat('./Calib/Test/n', num2str(i), '.txt'), feamat,  'delimiter',' ','precision','%.8f');
end
clear x;
clear y;




%% merge indivisual train and test file in calib folder for two
calibtraintxt = dir('./Calib/Train/*.txt');
calibtesttxt = dir('./Calib/Test/*.txt');

trainmat = [];
testmat = [];
for i = 1:length(calibtraintxt)
    indvmat = dlmread(strcat('./Calib/Train/', calibtraintxt(i).name));
    trainmat = [trainmat; indvmat]; 
end
clear indvmat;
dlmwrite('./Calib/Train.txt', trainmat, 'delimiter', ' ', 'precision','%.8f');
clear trainmat;

for i = 1:length(calibtesttxt)
    indvmat = dlmread(strcat('./Calib/Test/',calibtesttxt(i).name));
    testmat = [testmat; indvmat]; 
end
clear indvmat;
dlmwrite('./Calib/Test.txt', testmat, 'delimiter', ' ', 'precision','%.8f');
clear testmat;


%% train and test
ml_exe = 'boost.exe real';
ml_train_arg = '-p 1 -l 500 -r 5 -d 4 -i RAW';
ml_test_arg = '-t -p 1 -l 500 -r 5 -d 4 -i RAW';

train_cmd = {ml_exe, ml_train_arg, '-m', './Calib/Model.txt', '-f', './Calib/Train.txt'};
test_cmd = {ml_exe, ml_test_arg, '-m', './Calib/Model.txt', '-f', './Calib/Test.txt', '-o', './Calib/Result.txt'};


system(strjoin(train_cmd));
system(strjoin(test_cmd));






