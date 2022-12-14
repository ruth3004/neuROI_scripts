
clear all 
close all 

%% Load experiment

expName = '20220426_RM0008_130hpf_fP1_f3';
resultDir = fullfile(pwd,'results');
planeNum = 4;
planeString = NrModel.getPlaneString(planeNum);
traceResultDir = fullfile(resultDir,'time_trace',planeString);

expFilePath = fullfile(resultDir,sprintf('experimentConfig_%s.mat',expName));
foo = load(expFilePath);
myexp = foo.myexp;
disp(myexp.expInfo)
fileNameArray = myexp.rawFileList;

%%% Sort names by odor 
odorList = myexp.expInfo.odorList;
%fileNameArraySorted = shortcut.sortFileNameArray(fileNameArray,'odor',odorList);
fileNameArraySorted = fileNameArray;
nTrials = myexp.expInfo.nTrial;

%%% Load time trace matrices

traceResultArray = struct('timeTraceMat',{},'roiArray',{},...
                          'roiFilePath',{},'rawFilePath',{});
appendix = sprintf('_frame%dtoInfby4',planeNum);
for k=1:length(fileNameArraySorted)
    fileName = fileNameArraySorted{k};
    timeTraceFilePath = shortcut.getTimeTraceFilePath(traceResultDir,fileName,appendix);
    foo = load(timeTraceFilePath);
    traceResultArray(k) = foo.traceResult;
end

% Keep only the ROIs that appear in all trials
[commonRoiTagArray,timeTraceMatList,idxMat] = analysis.findCommonRoi(traceResultArray);

%% OPTIONAL - To delete a trial
% if contains(expName, 'Ruth-20200924-f3') 
%     timeTraceMatList(:,3) = [];
%     
% end
% 
% calcium_trace_matrix = zeros([size(timeTraceMatList{1}) size(timeTraceMatList,2)]);
% savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'ForClaire');
% mkdir(savepath_dFtraces);
% 
% for ii=1:size(timeTraceMatList,2)
%     calcium_trace_matrix(:,:,ii)=cell2mat(timeTraceMatList(ii));
% end
% savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'ForClaire');
% 
% save(fullfile(savepath_dFtraces,'calcium_traces_matrix_unfiltered.mat'),'calcium_trace_matrix')
%% Plot heatmapsd

figure
for tt = 1:length(timeTraceMatList)
    time_traces= timeTraceMatList{tt};
    subplot(nTrials,length(odorList),tt)
    imagesc(time_traces)
   

    %caxis([0 300]) 
end
colorbar
%% Define baseline
%bl_range = 60:90; % in frames
bl_range = 80:110; % in frames
%% Cleaning data based on baselines 
all_baselines = cellfun(@(x) x(:,bl_range),timeTraceMatList(1,:),'UniformOutput',false);
all_baselines_mean       = nanmean(cell2mat(all_baselines),2);

% Detecting outliers and NaN Rois
[~, bl_outliers]=rmoutliers(nanmean(all_baselines_mean,2));
[~, bl_nans]=rmmissing(all_baselines_mean,1);

% merge filters 
bl_filter = or(bl_nans, bl_outliers);

%% Do you want to duplicate baseline? WITH FILTER
duplicate = 1; % Yes: (1), no:(0) 
% duplicate baseline and attach to trimmed time trace
if duplicate == 1
timeTraces_modified = cellfun(@(x,y) cat(2,x(~bl_filter,:),y(~bl_filter,...
    bl_range(1):end)), all_baselines, timeTraceMatList, 'UniformOutput', false);
else
timeTraces_modified = cellfun(@(x) x(~bl_filter,bl_range(1):end), timeTraceMatList, 'UniformOutput', false);
end 

    figure
        for tt = 1:length(timeTraces_modified)
        time_traces= timeTraces_modified{tt};
        subplot(nTrials,length(odorList),tt)
        title(tt)
        imagesc(time_traces)
        %caxis([0 300])
        end
    colorbar
% normalize time trace to baseline mean (dF traces) 
baselines_means = cellfun(@(x) nanmean(x(~bl_filter),2), all_baselines,'UniformOutput', false);

timeTraces_norm = cellfun(@(x,y) (x-y)./y, timeTraces_modified,baselines_means,'UniformOutput', false);

    figure
        for tt = 1:length(timeTraces_norm)
            time_traces= timeTraces_norm{tt};
            subplot(nTrials,length(odorList),tt)
            imagesc(time_traces)
            caxis([0 2])
        end
    colorbar

 without_filter = 0  ; 
%% Do you want to duplicate baseline? WITHOUT FILTER
% duplicate baseline and attach to trimmed time trace    
without_filter = 0;
duplicate = 1; % Yes: (1), no:(0) 
if duplicate == 1
timeTraces_modified = cellfun(@(x,y) cat(2,x,y(:,bl_range(1):end)), all_baselines, timeTraceMatList, 'UniformOutput', false);
else
timeTraces_modified = cellfun(@(x) x(:,bl_range(1):end), timeTraceMatList, 'UniformOutput', false);
end 

    figure
        for tt = 1:length(timeTraces_modified)
        time_traces= timeTraces_modified{tt};
        subplot(6,4,tt)
        imagesc(time_traces)
        %caxis([0 300])
        end
    colorbar
    

baselines_means = cellfun(@(x) nanmean(x,2), all_baselines,'UniformOutput', false);

 % normalize time trace to baseline mean (dF traces) 
timeTraces_norm = cellfun(@(x,y) (x-y)./y, timeTraces_modified,baselines_means,'UniformOutput', false);

    figure
        for tt = 1:length(timeTraces_norm)
            time_traces= timeTraces_norm{tt};
            subplot(6,4,tt)
            imagesc(time_traces)
            caxis([0 10])
        end
    colorbar   
    
    
        
%% saving dF traces
savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'dF_traces_doublebl');
mkdir(savepath_dFtraces);


for ii = 1:length(timeTraces_norm)
    dF_traces = timeTraces_norm{ii};
    if without_filter == 1 
            file2save = append(replace(fileNameArraySorted{ii}, '_.tif','_'),'dF_traces.mat');
    else 
            file2save = append(replace(fileNameArraySorted{ii}, '.tif','_'),'dF_traces.mat');
    end
    save(fullfile(savepath_dFtraces,file2save),'dF_traces')
end 
        

%% saving dF traces as matrix FOR CLAIRE
matrix2save = zeros([size(timeTraces_norm{1}) size(timeTraces_norm,2)]);
savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'ForClaire');
mkdir(savepath_dFtraces);

for ii=1:size(timeTraces_norm,2)
    matrix2save(:,:,ii)=cell2mat(timeTraces_norm(ii));
end

if duplicate ==1 & without_filter ==1
    save(fullfile(savepath_dFtraces,'dF_traces_matrix_repl_baseline_unfiltered.mat'),'matrix2save')
elseif duplicate ==1 & without_filter ==0
    save(fullfile(savepath_dFtraces,'dF_traces_matrix_repl_baseline.mat'),'matrix2save')
else 
    save(fullfile(savepath_dFtraces,'dF_traces_matrix_unfiltered.mat'),'matrix2save')
end
%%

% reshape(cell2mat(timeTraces_norm)
% if duplicate = 1
% save(fullfile(savepath_dFtraces,'dF_traces_repl_baseline.mat'),'timeTraces_norm')

