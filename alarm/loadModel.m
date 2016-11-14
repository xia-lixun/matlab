function model = loadModel(modelFile)
% the model file output from boost.ext

modelData = textread(modelFile); % get dimension
nLoop = size(modelData,1)-1;
modelData = textread(modelFile,'%.6f'); % get float but lose dimension
model.sigA = modelData(1);
model.sigB = modelData(2);
modelData = reshape(modelData(3:end),[],nLoop)';
model.feaIdx = round(modelData(:,2))+1;% C indexing starts from 0!
model.dir = modelData(:,3);
model.alpha = modelData(:,4);
model.thld = modelData(:,5);
model.range = modelData(:,6);