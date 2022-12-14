%% Load experiment
%Input
expPath = "W:\scratch\gfriedri\montruth\2P_RawData\2022-04-26\f4";
expName = '20220426_RM0008_132hpf_fP1_f4';
duplicate = 1;      % Yes: (1), no:(0) 
filter = 0;         % Yes: (1), no:(0) 
planeArray = [1 2 3 4];

%Extract experiment from NrModel
resultDir = fullfile(expPath,'results');
expFilePath = fullfile(resultDir,sprintf('experimentConfig_%s.mat',expName));
foo = load(expFilePath);
myexp = foo.myexp;
disp(myexp.expInfo)
fileNameArray = myexp.rawFileList;

%Go through all selected planes 
for planeNum = 1: length(planeArray)
    planeString = NrModel.getPlaneString(planeArray(planeNum));
    traceResultDir = fullfile(resultDir,'time_trace',planeString);

    %%% Sort names by odor 
    odorList = myexp.expInfo.odorList;
    fileNameArraySorted = shortcut.sortFileNameArray(fileNameArray,'odor',odorList);
    nTrials = myexp.expInfo.nTrial;

    %%% Load time trace matrices
    traceResultArray = struct('timeTraceMat',{},'roiArray',{},...
                              'roiFilePath',{},'rawFilePath',{});
    appendix = sprintf('_frame%dtoInfby4',planeArray(planeNum));
    for k=1:length(fileNameArraySorted)
        fileName = fileNameArraySorted{k};
        timeTraceFilePath = shortcut.getTimeTraceFilePath(traceResultDir,fileName,appendix);
        foo = load(timeTraceFilePath);
        traceResultArray(k) = foo.traceResult;
    end

    % Keep only the ROIs that appear in all trials
    [commonRoiTagArray,timeTraceMatList,idxMat] = analysis.findCommonRoi(traceResultArray);

    %%% Plot heatmapsd

    figure
    for tt = 1:length(timeTraceMatList)
        time_traces= timeTraceMatList{tt};
        subplot(nTrials,length(odorList),tt)
        imagesc(time_traces)


        %caxis([0 300]) 
    end
    colorbar
    %%% Define baseline
    if planeNum == 1
        prompt = {'Enter baseline start in frames:','Enter baseline stop in frames:'};
        dlgtitle = 'Baseline definition';
        dims = [1 50];
        definput = {'80', '110'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        bl_range =  [str2double(answer{1}):str2double(answer{2})];
    end

    
    %%% Calculating baseline and its mean 
    all_baselines = cellfun(@(x) x(:,bl_range),timeTraceMatList(1,:),'UniformOutput',false);
    all_baselines_mean       = nanmean(cell2mat(all_baselines),2);

    % Detecting outliers and NaN Rois
    [ID_outliers, bl_outliers]=rmoutliers(nanmean(all_baselines_mean,2));
    [ID_nans, bl_nans]=rmmissing(all_baselines_mean,1);

    % merge filters 
    bl_filter = or(bl_nans, bl_outliers);

    % duplicate baseline and attach to trimmed time trace
    if duplicate == 1
        if filter == 1
            timeTraces_modified = cellfun(@(x,y) cat(2,x(~bl_filter,:),y(~bl_filter,...
                bl_range(1):end)), all_baselines, timeTraceMatList, 'UniformOutput', false);
        else 
            timeTraces_modified = cellfun(@(x,y) cat(2,x,y(:,bl_range(1):end)), ...
                all_baselines, timeTraceMatList, 'UniformOutput', false);
        end
    elseif duplicate == 0
        if filter == 1
            timeTraces_modified = cellfun(@(x) x(~bl_filter,bl_range(1):end), timeTraceMatList, 'UniformOutput', false);
        else
        timeTraces_modified = cellfun(@(x) x(:,bl_range(1):end), timeTraceMatList, 'UniformOutput', false);
        end 
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

    %%% normalize time trace to baseline mean (dF traces) 

     if filter == 0 
        baselines_means = cellfun(@(x) nanmean(x,2), all_baselines,'UniformOutput', false);
     elseif filter == 1
        baselines_means = cellfun(@(x) nanmean(x(~bl_filter),2), all_baselines,'UniformOutput', false);
     end
        
    timeTraces_norm = cellfun(@(x,y) (x-y)./y, timeTraces_modified, baselines_means,'UniformOutput', false);

    %%% filter roitags if filer is on 
      if filter == 1 
        roiTags = commonRoiTagArray(~bl_filter);
     elseif filter == 0
         roiTags = commonRoiTagArray;
      end
    
    
        figure
            for tt = 1:length(timeTraces_norm)
                time_traces= timeTraces_norm{tt};
                subplot(nTrials,length(odorList),tt)
                imagesc(time_traces)
                caxis([0 2])
            end
        colorbar

    %%% saving dF traces
    savepath_dFtraces = fullfile(resultDir,'spike_deconvolution',planeString,'dF_traces_doublebl');
    mkdir(savepath_dFtraces);

    for ii = 1:length(timeTraces_norm)
        dF_traces = timeTraces_norm{ii};
        if filter == 1 
                filename = append(replace(fileNameArraySorted{ii}, '.tif','_'),'dF_traces_filtered.mat');
        else 
                filename = append(replace(fileNameArraySorted{ii}, '.tif','_'),'dF_traces.mat');
        end
        save(fullfile(savepath_dFtraces,filename),'dF_traces', 'roiTags')
    end 

end 
