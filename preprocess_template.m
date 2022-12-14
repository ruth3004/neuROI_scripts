addpath('../../neuRoi')
%% Clear variables
clear all
%% Step01 Configure experiment and image processing parameters
% Load 
% Experiment parameters
expInfo.name = '20221026_NT0012_f3_13dpf_testNTvsRModors';
expInfo.frameRate = 30;
expInfo.odorList = {'Ala','Ser','TDCA','Ctrl','GCA','TCA','Cad','SA'};
expInfo.nTrial = 4;
expInfo.nPlane = 4;

expSubDir = expInfo.name;

% Raw data
rawDataDir =pwd;


%%
listAllFiles = struct2cell(dir(rawDataDir));
listAllFiles = listAllFiles(1,:);
[indx, tf] = listdlg('PromptString','Select all the files you want to preprocess: ','ListString',listAllFiles,'ListSize',[300,250]);
% res = True;
% while res == True
    rawFileList = listAllFiles(indx);
%     nFiles = length(rawFileList);
%     res = questdlg("")


% Data processing configuration
% Directory for saving processing results
resultDir = fullfile(rawDataDir,'results');

% Directory for saving binned movies

binDir = fullfile(resultDir,'binned');

%% Step02 Initialize NrModel with experiment configuration
myexp = NrModel(rawDataDir,rawFileList,resultDir,...
                expInfo);
%% Step03a (optional) Bin movies
% Bin movie parameters
binParam.shrinkFactors = [1, 1, 2];
binParam.trialOption = struct('process',true,'noSignalWindow',[1 4]);
binParam.depth = 8;
for planeNum=1:myexp.expInfo.nPlane
myexp.binMovieBatch(binParam,binDir,planeNum);
end
%% Step03b (optional) If binning has been done, load binning
%% parameters to experiment
%read from the binConfig file to get the binning parameters
binConfigFileName = 'binConfig.json';
binConfigFilePath = fullfile(binDir,binConfigFileName);
myexp.readBinConfig(binConfigFilePath);
%% Step04 Calculate anatomy maps
% anatomyParam.inFileType = 'raw';
%anatomyParam.trialOption = {'process',true,'noSignalWindow',[1 24]};
anatomyParam.inFileType = 'binned';
anatomyParam.trialOption = [];
for planeNum=1:myexp.expInfo.nPlane
    myexp.calcAnatomyBatch(anatomyParam,planeNum);
end
%% Step04b If anatomy map has been calculated, load anatomy
%% parameters to experiment
anatomyDir = myexp.getDefaultDir('anatomy');
anatomyConfigFileName = 'anatomyConfig.json';
anatomyConfigFilePath = fullfile(anatomyDir,anatomyConfigFileName);
myexp.readAnatomyConfig(anatomyConfigFilePath);
%% Step05 Align trial to template
templateRawName = myexp.rawFileList{1};
% plotFig = false;
% climit = [0 0.5];
for planeNum=1:myexp.expInfo.nPlane
    myexp.alignTrialBatch(templateRawName,[],...
                          'planeNum',planeNum,...
                          'alignOption',{'plotFig',false});
end
%% Save experiment configuration
expFileName = strcat('experimentConfig_',expInfo.name,'.mat');
expFilePath = fullfile(resultDir,expFileName);
save(expFilePath,'myexp')
