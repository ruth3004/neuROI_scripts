
clear all 
close all 


%% Load experiment
% resultDir = 'C:\Users\montruth\Documents\Work backup\20210202\f4\results';
resultDir = 'C:\Users\montruth\Documents\Work backup\20200924\f3\results';
% expName = 'Ruth-20210202-f4';
expName = 'Ruth-20200924-f3';
% planeNum = 3;
planeNum = 1;
planeString = NrModel.getPlaneString(planeNum);
traceResultDir = fullfile(resultDir,'time_trace',planeString);

expFilePath = fullfile(resultDir,sprintf('experimentConfig_%s.mat',expName));
foo = load(expFilePath);
myexp = foo.myexp;
disp(myexp.expInfo)
fileNameArray = myexp.rawFileList;

%% Sort names by odor 
odorList = myexp.expInfo.odorList;
fileNameArraySorted = shortcut.sortFileNameArray(fileNameArray,'odor',odorList);
nTrials = myexp.expInfo.nTrial;

%% Load time trace matrices

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



%% Defining baseline
idx = 9;
time_traces = timeTraceMatList{idx};
figure
imagesc(time_traces)

plot(nanmean(time_traces,1))

% bl_range = 80:100;  % in frames
bl_range = 50:80;  % in frames

%% Cleaning data and creating dF traces 
dFtracesMatList;

for tt = 1:length(timeTraceMatList)
time_traces= timeTraceMatList{tt};
baseline = time_traces(:,bl_range); % remove nan columns

%     figure
%     imagesc(baseline)
%     colorbar
%     caxis([0 200])

%%% Filtering noisy neurons and NaN ROIs

[~, bl_outliers]=rmoutliers(nanmean(baseline,2));
[~, bl_nans]=rmmissing(baseline,1);


% merge filters 
bl_filter = or(bl_nans, bl_outliers);

%     figure
%     imagesc(baseline(~bl_filter,:))
%     colorbar
%     caxis([0 200])

%apply filter
bl_filtered = baseline(~bl_filter,:);


%%% Replicating baseline

bl_means = mean(bl_filtered,2);
bl_stds = std(bl_filtered,0,2);

bl_rnd = zeros(size(bl_filtered));

% a = 5;    %mean
% b = 500;  %std
% y = a.*randn(1000,1) + b;

for ii = 1:size(bl_filtered,1)
    bl_rnd(ii,:) = bl_means(ii).*rand(size(bl_filtered,2),1)+bl_stds(ii);
end

% check if mean and std true
% check_mean = bl_means-mean(bl_rnd,2);
% check_std = bl_stds-std(bl_rnd,0,2);
% 
% bl_replica = [bl_rnd, bl_filtered];
% 
% figure
% imagesc(bl_replica)
% colorbar
% caxis([0 200])

% check_mean2 = bl_means-mean(bl_replica,2);
% check_std2 = bl_stds-std(bl_replica,0,2);



%%% Adding new baseline to time trace and cut 

% time_traces_modified1 = time_traces;
% time_traces_modified1 = [bl_rnd,time_traces(~bl_filter,bl_range(1):end)];

time_traces_modified = time_traces;
time_traces_modified = [bl_filtered,time_traces(~bl_filter,bl_range(1):end)];


%%% Normalizing traces 

% time_traces_norm1 = time_traces_modified1-mean(bl_means,2);
% 
% figure
% imagesc(time_traces_norm1)

time_traces_norm = time_traces_modified-mean(bl_means,2);

% figure
% imagesc(time_traces_norm)


%% Save dF traces 
 
savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'dF_traces');

mkdir(savepath_dFtraces);


for ii = 1:length(dFtracesMatList)
    dF_traces = dFtracesMatList{ii};
    file2save = append(replace(fileNameArraySorted{ii}, '.tif','_'),'dF_traces.mat');
    save(fullfile(savepath_dFtraces,file2save),'dF_traces')
end 



end
%%

figure
for ii = 1: size(timeTraceMatList,2)
    bl_mean = zeros(1,size(timeTraceMatList{1},1));
    figure
    ii
    for jj = 1:size(timeTraceMatList{1},1)
        jj
        dF_traces_raw = timeTraceMatList{ii};
        bl_mean(jj)= nanmean(dF_traces_raw(jj,bl),2);
        plot(jj, bl_mean(jj), 'ok')
        hold on 
    end
end
blthr  = 3*nanmean(bl_mean,'all');
line([1,size(dF_traces_raw,1)], [blthr blthr])

%selRois = find(blMean<blthr);

% dF_traces_raw_filtered = dF_traces_raw(selRois);


% rnd_bl = rand


%% Define baseline
baseline = 85:130; 
resp_window = length(baseline):length(baseline)+75;



%% median F0
dFtracesMatList = cell(size(timeTraceMatList));
figure

for ii = 1:length(timeTraceMatList)
    dF_traces_raw = timeTraceMatList{ii};

    % F0 = median(dF_traces_raw(:,baseline),2);
    F0 = dF_traces_raw(:,baseline);
    F0_median = nanmedian(nanmean(F0,1),2);
    
    
        
    
    dF_traces_abs = ((dF_traces_raw(:,baseline(1):end)- F0_median)./F0_median);
    dF_traces = dF_traces_abs;
    %dF_traces = dF_traces_abs/max(dF_traces_abs, [], "all");

    % delete nans
    dF_traces = rmmissing(dF_traces, 1);
    dFtracesMatList{ii} =dF_traces;
    
   

        respMatList = zeros(size(dF_traces,1),1);
        for jj = 1:size(dF_traces,1)
             th = 5*std(dF_traces(jj,1:length(baseline)));
             resp = find(dF_traces(jj,resp_window)>th,1);
             if ~isempty(resp)
                 respMatList(jj) = resp+length(baseline);
 
             end
        end 

    
    
    
    %std too large?
   
    
%     figure; plot(1:size(dF_traces,2), dF_traces(jj,:))
%     hold on 
%     plot(resp, dF_traces(jj,resp), 'r*')
    
    %plt heatmaps
    subplot(nTrials,length(odorList), ii)
    imagesc(dF_traces)
    hold on 
    plot(respMatList, 1:size(respMatList),'r*')
    traceTitle = char(fileNameArraySorted{ii});
    title(traceTitle(27:34))

    grid off
    
    
end

%% filtering non-responsive neurons 






%% save dF_traces

savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'dF_traces');

mkdir(savepath_dFtraces);


for ii = 1:length(dFtracesMatList)
    dF_traces = dFtracesMatList{ii};
    file2save = append(replace(fileNameArraySorted{ii}, '.tif','_'),'dF_traces.mat');
    save(fullfile(savepath_dFtraces,file2save),'dF_traces')
end 




