
clear all 
close all 

%% Load experiment
% resultDir = 'C:\Users\montruth\Documents\Work backup\20200924\f3\results\';
% expName = 'Ruth-20200924-f3';
% planeNum = 3;
resultDir = 'W:\scratch\gfriedri\montruth\2P_RawData\2022-04-26\f3\results';
expName = '20220426_RM0008_130hpf_fP1_f3';
planeNum = 2;
nTrialPerOdor = 3;

planeString = NrModel.getPlaneString(planeNum);
predictionDir = fullfile(resultDir,'spike_deconvolution',planeString, 'dF_traces_doublebl');

expFilePath = fullfile(resultDir,sprintf('experimentConfig_%s.mat',expName));
foo = load(expFilePath);
myexp = foo.myexp;
disp(myexp.expInfo)
fileNameArray = myexp.rawFileList;

%% Sort names by odor 
odorList = myexp.expInfo.odorList;
nOdor = length(odorList);
fileNameArraySorted = shortcut.sortFileNameArray(fileNameArray,'odor',odorList);
nTrials = myexp.expInfo.nTrial;


%% Load time trace matrices

spikeRatesArray = cell(size(fileNameArray));
nTrial = length(spikeRatesArray);
suffix = 'full_prediction_';
appendix = '_dF_traces';
for k=1:length(fileNameArraySorted)
    fileName = fileNameArraySorted{k};
    spikeRatesFilePath = fullfile(predictionDir,...
        append(suffix,replace(fileName, '.tif',appendix)));
    foo = load(spikeRatesFilePath);
    spikeRatesArray{k} = foo.spike_rates;
end


%%
%    figure
%         for tt = 1:length(spikeRatesArray)
%             spike_traces= spikeRatesArray{tt};
%             subplot(6,4,tt)
%             imagesc(spike_traces(:,10:250))
%             caxis([0 2])
%         end
%     colorbar

    
zlim = [0 1.5];
nCol = length(odorList)+1;
nRow = nTrialPerOdor;
nSubplot = length(spikeRatesArray);
indMat = reshape(1:nRow*nCol,nCol,nRow).';
cutWindow = 50:300;

figWidth = 1800;
figHeight = 300*nRow;
fig = figure('InnerPosition',[200 500 figWidth figHeight]);
for k=1:nSubplot
    subplot(nRow,nCol,indMat(k))
    imagesc(spikeRatesArray{k})
    % imagesc(timeTraceMatList{k})
    % imagesc(traceResultArray(k).timeTraceMat)
     %ax.Visible = 'off';
    if mod(k,nRow) == 1
        ax = gca;
        odor = shortcut.getOdorFromFileName(fileNameArraySorted{k});
        title(odor);
        set(get(ax,'Title'),'Visible','on');
    end
    caxis(zlim)
    
%     if k<nOdor
%         set(gca,'XTick',[]);
%     else
    ticks=[0:5:cutWindow(2)]*7.5;
    
    for qqq=1:length(ticks)
        xLabels{qqq}=num2str(ticks(qqq)/7.5);
    end    
    set(gca,'XTick',ticks,'XTickLabel',xLabels)
    %xlim([cutWindow(1)*7.5 cutWindow(end)])
end
subplot(nRow,nCol,indMat(nSubplot+1))
caxis(zlim)
colorbar('Location','west')
axis off

%% Define odors to correlate

o1 =1; % File number sorted
o2= 5; % File number sorted



%% Plot correlation across time between two odors
sampling_freq = myexp.expInfo.frameRate/myexp.expInfo.nPlane;
odor_1 = spikeRatesArray{o1};
odor_2 = spikeRatesArray{o2};

figure
 subplot(2,1,1)
for n = 1:size(odor_1,1)
     %lot((1:size(odor_1,2))/sampling_freq, odor_1(n,:)+ n)
     plot((1:size(odor_1,2)), odor_1(n,:)+ n)
     hold on 
%  line( [mt1 mt1],[0 109], 'Color', 'k', 'LineWidth',2) 
%  line( [mt2 mt2],[0 109], 'Color', 'r','LineWidth',2) 
%     
end
axis tight

    ticks=[0:5:cutWindow(2)]*7.5;
    
    for qqq=1:length(ticks)
        xLabels{qqq}=num2str(ticks(qqq)/7.5);
    end    
    set(gca,'XTick',ticks,'XTickLabel',xLabels)


subplot(2,1,2)
for n = 1:size(odor_2,1)
     %plot((1:size(odor_2,2))/sampling_freq, odor_2(n,:)+ n)
     plot((1:size(odor_2,2)), odor_2(n,:)+ n)
     hold on 
%     
%  line( [mt1 mt1],[0 109], 'Color', 'k', 'LineWidth',2) 
%  line( [mt2 mt2],[0 109], 'Color', 'r','LineWidth',2) 
end

    ticks=[0:5:cutWindow(2)]*7.5;
    
    for qqq=1:length(ticks)
        xLabels{qqq}=num2str(ticks(qqq)/7.5);
    end    
    set(gca,'XTick',ticks,'XTickLabel',xLabels)


xlabel('Time in seconds')
axis tight



%% Define window parameters for correlations
window = 4; % in frames
start = 5;    % in seconds
stop = 30;     % in seconds
response = 10*7.5; % in frames

%% Correlation matrix of two odors 

mt= 1:length(odor_1)-window; % moving time window in frames
Mov_Corr_Mat = zeros(length(mt),1);

mt_start = ceil(start*sampling_freq); % in frames
mt_stop = ceil(stop*sampling_freq); % in frames

for ii = 1:length(mt)
       
    time_start = mt(ii);%ceil(mt(ii)*sampling_freq); % in frames

    rates_odor_1 = odor_1(:,time_start:time_start+window-1);
    rates_odor_2 = odor_2(:,time_start:time_start+window-1);

    mean_rates_odor_1 = mean(rates_odor_1,2,'omitnan');
    mean_rates_odor_2 = mean(rates_odor_2,2,'omitnan');

    Mov_Corr_Mat(ii) =  corr2(mean_rates_odor_2,mean_rates_odor_1);
end

figure
plot(1:length(mt), Mov_Corr_Mat)
%view([-90 -90])
% colorbar
% caxis([0 1])

hold on 

 line( [mt_start mt_start],[-0.2 1], 'Color', 'k', 'LineWidth',2) 
 line( [mt_stop mt_stop],[-0.2 1], 'Color', 'r','LineWidth',2) 
 %line( [response response],[-1 1], 'Color', 'm','LineWidth',2, 'LineStyle', '--') 
%ylim([30 60])
% set(gca, 'XTick', 0:sampling_freq/2:length(mt))
    ticks=[0:5:cutWindow(2)]*7.5;
    
    for qqq=1:length(ticks)
        xLabels{qqq}=num2str(ticks(qqq)/7.5);
    end    
    set(gca,'XTick',ticks,'XTickLabel',xLabels)
 

title(['Correlation of ', odorList{(o1-mod(o1,2))/nTrialPerOdor+1},...
    ' and ', odorList{(o2-mod(o2,3))/3+1}, ' at frames ',...
    num2str(mt_start) ' (14.5s) and ' num2str(mt_stop), ' (16.0s)  with window: ',num2str(window)])

xlabel('Time in seconds')
xlim([5 30]*7.5)

%% Correlation matrix of all odors 

mean_rates_all = zeros(size(spikeRatesArray{1},1),length(spikeRatesArray));
time_start = ceil(mt_start); % in frames

for jj = 1:length(spikeRatesArray)
    odor=spikeRatesArray{jj};
    rates_odor = odor(:,time_start:time_start+window-1);
    mean_rates_all(:,jj) = mean(rates_odor,2,'omitnan');
    if jj == 1
        CorrAllMat = mean_rates_all(:,jj);
    else
       CorrAllMat = [CorrAllMat mean_rates_all(:,jj)];
    end
end

A =  corrcoef(CorrAllMat);
  

figure
    imagesc(A)

    colorbar
    caxis([0 1])

    set(gca, 'YTick',2:3:24)
    set(gca, 'YTickLabel',odorList)
    set(gca, 'XTick',2:3:24)
    set(gca, 'XTickLabel',odorList)
xtickangle(gca,45)
     hold on
 for nl = 3.5:3:24
    line([nl nl],[0 24.5],'LineWidth',2, 'Color', 'k')
    line([0 24.5], [nl nl],'LineWidth',2, 'Color', 'k')
 end

 
%% Correlation difference

mean_rates_all = zeros(size(spikeRatesArray{1},1),length(spikeRatesArray));
time_start = ceil(mt_stop); % in frames

   for jj = 1:length(spikeRatesArray)
        odor=spikeRatesArray{jj};
        rates_odor = odor(:,time_start:time_start+window-1);
        mean_rates_all(:,jj) = mean(rates_odor,2,'omitnan');
        if jj == 1
            CorrAllMat = mean_rates_all(:,jj);
        else
           CorrAllMat = [CorrAllMat mean_rates_all(:,jj)];
        end
   end

   [B, p] =  corrcoef(CorrAllMat);
  

figure
f(1) = subplot(1,3,1);
    imagesc(A)  
    colorbar
    caxis([0 1])

     hold on
for li = 1: nOdor-1
    line([0 nOdor*nTrial],[li li]*nTrialPerOdor+0.5,'LineWidth',2, 'Color','k')
    line([li li]*nTrialPerOdor+0.5,[0 nOdor*nTrial],'LineWidth',2, 'Color','k')
end   
 axis square   
 title(['Start correlation at ', sprintf('%.2f', mt_start/sampling_freq),'s']); 

f(2) = subplot(1,3,2);
    imagesc(B)
    colorbar
    caxis([0 1])
    hold on
for li = 1: nOdor-1
    line([0 nOdor*nTrial],[li li]*nTrialPerOdor+0.5,'LineWidth',2, 'Color','k')
    line([li li]*nTrialPerOdor+0.5,[0 nOdor*nTrial],'LineWidth',2, 'Color','k')
end   
     
axis square
 title(['End correlation at ', sprintf('%.2f', mt_stop/sampling_freq),'s']); 
  
     
f(3) = subplot(1,3,3);
    imagesc(B-A)
    colorbar
    
    mycolors = zeros(100,3);
    mycolors(1:50,3) = linspace(1,0,50);
    mycolors(:,2) = zeros(100,1);
    mycolors(51:100,1) = linspace(0,1,50);
    colormap(f(3),mycolors)

for li = 1: nOdor-1
    line([0 nOdor*nTrial],[li li]*nTrialPerOdor+0.5,'LineWidth',2, 'Color','k')
    line([li li]*nTrialPerOdor+0.5,[0 nOdor*nTrial],'LineWidth',2, 'Color','k')
end   
axis square
    
     
set(f, 'YTick',nTrialPerOdor/2:nTrialPerOdor:nTrial, 'YTickLabel',odorList,...
         'XTick',nTrialPerOdor/2:nTrialPerOdor:nTrial,'XTickLabel',odorList,  'FontSize', 14)


 title(['Correlation difference between ', sprintf('%.2f', (mt_stop-mt_start)/sampling_freq),'s (', num2str((mt_stop-mt_start)),'frames).']);      

%% Create matrices: All mean spike rate of rolling window (mean_rates_all) and all correlation matrices along rolling window (CorrAllMat)
windowRow = mt_start:1:mt_stop;
mean_rates_all = zeros(length(windowRow),...            %timesteps
                    size(spikeRatesArray{1},1),...      %roi
                    length(spikeRatesArray));           %trials
                
CorrAllMat = zeros(length(windowRow),...                %timesteps
                    length(spikeRatesArray),...         %correlation matrix
                    length(spikeRatesArray));                  


for win = 1:length(windowRow) % rolling window (timesteps)

    time_start = windowRow(win); % time start in frames

       for jj = 1:length(spikeRatesArray) % along trials 
            odor=spikeRatesArray{jj};
            rates_odor = odor(:,time_start:time_start+window-1);
            mean_rates_all(win,:,jj) = mean(rates_odor,2,'omitnan');
       end
    CorrAllMat(win,:,:) = corrcoef(squeeze(mean_rates_all(win,:,:)));
end
%% Ploting correlations across time against odors (29.05.2021) 


odor1 = 3; 
odor2 = 5;
nt =2;

odor1_list = [odor1*ones(1, 1),(odor1+1)*ones(1, 3),(odor1+2)*ones(1, 3)];%,(odor1+3)*ones(1, 4)]
odor2_list = repmat((odor2-1)*odor2,1,3) + [1 2 3]
corrOd1Od2Mat = zeros(size(CorrAllMat,1),length(odor1_list));


for ii = 1:length(odor1_list)
    for jj= odor2_list
     corrOd1Od2Mat(:,ii) =CorrAllMat(:,odor1_list(ii), jj);
    end
end


figure
%for ii= 1:nTrials
t =  windowRow/sampling_freq; 
corr_mean = nanmean(corrOd1Od2Mat,2)';
corr_std = std(corrOd1Od2Mat,0,2,'omitnan')';

%shaded area of std
curve1 = corr_mean + corr_std;
curve2 = corr_mean - corr_std;
t2 = [t, fliplr(t)];
inBetween = [curve1, fliplr(curve2)];

%plot
%subplot(2, 2,ii)
    fill(t2, inBetween, 'k'); 
    alpha(0.2)
    hold on     
    plot(t,corr_mean, 'Color',my_colors(1,:),'LineWidth',2)
    line([14 19],[0.7 0.7], 'LineWidth', 2,'Color', 'k', 'LineWidth',2)
    
    ylim([-0.1 0.8])
    xlim([5 30])
    xlabel('Time in s')
    ylabel('Correlation coefficient')
    title([odorList{odor1},' vs ',odorList{odor2}])
%end



%%



for kk = 1: length(odor1_list)
subplot(2,3, kk)

    for ii = 1:nTrialPerOdor
        plot((1:size(CorrAllMat,1))+mt_start,CorrAllMat(:,odor1_list(kk)*ii,odor2_list(kk)*ii),'-o','MarkerSize',2)
        hold on 
    end
    title([odorList{odor1_list(kk)},' vs ' odorList{odor2_list(kk)}], 'FontSize',14)
    if kk == length(odor1_list)
        legend('Trial 1','Trial 2','Trial 3','Trial 4')
    end
    ticks=[0:5:cutWindow(2)]*7.5;
    
    for qqq=1:length(ticks)
        xLabels{qqq}=num2str(ticks(qqq)/7.5);
    end    
    set(gca,'XTick',ticks,'XTickLabel',xLabels)
    %axis tight
    xlabel('Time in s')
    ylabel('Correlation coefficient')
end
    

%% Plot correlation rows
counter = 1;
fig = figure;
for win = 1:6:length(windowRow)-3

    time_start = windowRow(win); % in frames
    subplot(2,length(windowRow(win)/2,counter))
    imagesc(squeeze(CorrAllMat(win,:,:)))
    colorbar
    caxis([0 0.8])

     hold on
for li = 1: nOdor-1
    line([0 nOdor*nTrial],[li li]*nTrialPerOdor+0.5,'LineWidth',2, 'Color','k')
    line([li li]*nTrialPerOdor+0.5,[0 nOdor*nTrial],'LineWidth',2, 'Color','k')
end  
     
 axis square   
%  t= char(['Correlation with window ', num2str(window),' at frame ', num2str(win),' (',num2str(win/7.5),'s)']);
 t= char([num2str(windowRow(win)/sampling_freq),' s']);

 title(t); 
 %set(gca, 'YTick',2:3:24, 'YTickLabel',odorList,...
 %        'XTick',2:3:24,'XTickLabel',odorList)
 %   xtickangle(gca,45)
 counter = counter +1 ;   
% saveas(fig, [predictionDir, t,'.png'])
end

%% Average of correlation matrices per trial 


avgmean_rates_all = zeros(length(windowRow),...         %timesteps
                    size(spikeRatesArray{1},1),...      %mean spike rates at that window at that roi
                    length(spikeRatesArray));           %trials

avgCorrAllMat = zeros(length(windowRow),...             %timesteps
                    length(spikeRatesArray)/nTrials,... %correlation matrix
                    length(spikeRatesArray)/nTrials);                  


for win= 1:length(windowRow)

    avgmean_rates_all(win,:,:) = nanmean(permute(reshape(...
                squeeze(mean_rates_all(win,:,:)), ...
                size(squeeze(mean_rates_all(win,:,:)),1),nTrialPerOdor,[]), ...
                [1 3 2]),3);
%     avgmean_rates_all(win,:,:) = nanmean(permute(reshape(...
%                 mean_rates_all(win,:,:), ...
%                 length(mean_rates_all(win,:,:)),nTrialPerOdor,[]), ...
%                 [1 3 2]),3);

    avgCorrAllMat(win,:,:)= corrcoef(squeeze(avgmean_rates_all(win,:,:)));
end
%% Ploting average correlation matrices (odors)
% Averaged correlation matrix across different odors 
CorrAcrossTime = zeros(length(windowRow),1);
Legends = cell(1);
for ii = 1:length(odorList)
odor = ii; 
counter = 1;
my_colors = turbo(8);
% my_colors = brighten(my_colors,-.8)
figure 
 for od = 1:length(odorList)
    
    if odor ~= od 
        for cc = 1:size(avgmean_rates_all,1)         
            CorrAcrossTime(cc)= corr2(avgmean_rates_all(cc,:,odor),...
                        avgmean_rates_all(cc,:,od));
        
        end
        Legends{counter} = [odorList{odor},'-', odorList{od}];
        counter= counter+1;
           plot(windowRow/sampling_freq,CorrAcrossTime,...
        'Color',my_colors(od,:), 'LineWidth', 2)
        hold on 
    end
    
end

legend(Legends)
ylim([0 1])
xlabel('Time in s')
ylabel('Correlation coefficient')
end

%% Ploting average correlation matrices (concentrations)
% Averaged correlation matrix across different concentrations 
CorrAcrossTime = zeros(length(windowRow),1);
Legends = cell(1);
counter = 1;
my_colors = hot(8);
figure 
 for od = 1:2:length(odorList)-1
    
     for cc = 1:size(avgmean_rates_all,1)         
        CorrAcrossTime(cc)= corr2(avgmean_rates_all(cc,:,od),...
                        avgmean_rates_all(cc,:,od+1));
        
        end
        Legends{counter} = [odorList{od},'-', odorList{od+1}];
        counter= counter+1;
        plot(windowRow/sampling_freq,CorrAcrossTime,...
                'Color',my_colors(od,:), 'LineWidth', 2)
        hold on 
    
end

legend(Legends)
ylim([0 1])
xlabel('Time in s')
ylabel('Correlation coefficient')

%% Correlation within same odor (trials)
CorrAcrossTrials = zeros(nTrials,nchoosek(nTrials,2),length(windowRow));
% Legends = cell(1);
% counter = 1;
my_colors = turbo(8);


 for trial = 1:4%1:nTrialPerOdor:length(fileNameArray)
    this_trial = trial;%ceil(trial/nTrialPerOdor) 
    figure 

    Legends = cell(1);

    for cc = 1:size(avgmean_rates_all,1)     
        
        CorrAcrossTrials(this_trial,1,cc)= corr2(squeeze(mean_rates_all(cc,:,trial)),...
            squeeze(mean_rates_all(cc,:,trial+1)));
        CorrAcrossTrials(this_trial,2,cc)= corr2(squeeze(mean_rates_all(cc,:,trial)),...
            squeeze(mean_rates_all(cc,:,trial+2)));
        CorrAcrossTrials(this_trial,3,cc)= corr2(squeeze(mean_rates_all(cc,:,trial)),...
            squeeze(mean_rates_all(cc,:,trial+3)));
        CorrAcrossTrials(this_trial,4,cc)= corr2(squeeze(mean_rates_all(cc,:,trial+1)),...
            squeeze(mean_rates_all(cc,:,trial+2)));
        CorrAcrossTrials(this_trial,5,cc)= corr2(squeeze(mean_rates_all(cc,:,trial+1)),...
            squeeze(mean_rates_all(cc,:,trial+3)));
        CorrAcrossTrials(this_trial,6,cc)= corr2(squeeze(mean_rates_all(cc,:,trial+2)),...
            squeeze(mean_rates_all(cc,:,trial+3)));

             
      end
      Legends{1} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-2'];
      Legends{1} = [odorList{trial},' 1-2'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,1,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2)
      hold on 
      Legends{2} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-3'];
      Legends{2} = [odorList{trial},' 1-3'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,2,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
      
      Legends{3} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-4'];
      Legends{3} = [odorList{trial},' 1-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,3,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', '--')
            
      Legends{4} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 2-3'];
      Legends{4} = [odorList{trial},' 2-3'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,4,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
      
      Legends{5} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 2-4'];
      Legends{5} = [odorList{trial},' 2-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,5,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', '--')
            
      Legends{6} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 3-4'];
      Legends{6} = [odorList{trial},' 3-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossTrials(this_trial,6,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
      
            
                  
      legend(Legends)
      
      
    %ylim([0 1])
    xlabel('Time in s')
    ylabel('Correlation coefficient')
      
    
end

legend(Legends)

%ylim([0 1])
xlabel('Time in s')
ylabel('Correlation coefficient')


%% Mean and std of correlation within same odor (trials)

figure
for ii= 1:nTrials
t =  windowRow/sampling_freq; 
corr_mean = nanmean(squeeze(CorrAcrossTrials(ii,:,:)),1);
corr_std = nanstd(squeeze(CorrAcrossTrials(ii,:,:)),1);

%shaded area of std
curve1 = corr_mean + corr_std;
curve2 = corr_mean - corr_std;
t2 = [t, fliplr(t)];
inBetween = [curve1, fliplr(curve2)];

%plot
subplot(2, 2,ii)
    fill(t2, inBetween, 'k'); 
    alpha(0.2)
    hold on     
    plot(t,corr_mean, 'Color',my_colors(ii,:),'LineWidth',2)
    line([13.9 19],[0.7 0.7], 'LIneWidth', 2,'Color', 'k', 'LineWidth',2)
    
    ylim([-0.1 0.8])
    xlim([5 30])
    xlabel('Time in s')
    ylabel('Correlation coefficient')
    title(odorList{ii})
end


%% Correlation across odors
CorrAcrossOdors = zeros(nTrials,16,length(windowRow));
% Legends = cell(1);
% counter = 1;
my_colors = turbo(8);


 for odor = 1:nTrialPerOdor
    %this_trial = ceil(trial/nTrialPerOdor); 
    figure 

    Legends = cell(1);

    for cc = 1:size(avgmean_rates_all,1)     
        for tr = 1: 16
        CorrAcrossOdors(odor,tr,cc)= corr2(squeeze(mean_rates_all(cc,:,odor)),...
             squeeze(mean_rates_all(cc,:,odor+1)));
%         CorrAcrossOdors(odor,2,cc)= corr2(squeeze(mean_rates_all(cc,:,odor)),...
%             squeeze(mean_rates_all(cc,:,odor+2)));
%         CorrAcrossOdors(odor,3,cc)= corr2(squeeze(mean_rates_all(cc,:,odor)),...
%             squeeze(mean_rates_all(cc,:,odor+3)));
%         CorrAcrossOdors(odor,4,cc)= corr2(squeeze(mean_rates_all(cc,:,odor+1)),...
%             squeeze(mean_rates_all(cc,:,odor+2)));
%         CorrAcrossOdors(odor,5,cc)= corr2(squeeze(mean_rates_all(cc,:,odor+1)),...
%             squeeze(mean_rates_all(cc,:,odor+3)));
%         CorrAcrossOdors(odor,6,cc)= corr2(squeeze(mean_rates_all(cc,:,odor+2)),...
%             squeeze(mean_rates_all(cc,:,odor+3)));
        end
             
      end
      Legends{1} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-2'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,1,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2)
      hold on 
      Legends{2} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-3'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,2,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
      
      Legends{3} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 1-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,3,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', '--')
            
      Legends{4} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 2-3'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,4,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
      
      Legends{5} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 2-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,5,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', '--')
            
      Legends{6} = [odorList{(trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1},' 3-4'];
      plot(windowRow/sampling_freq,squeeze(CorrAcrossOdors(odor,6,:)),...
                'Color',my_colors((trial-mod(trial,nTrialPerOdor))/nTrialPerOdor+1,:), 'LineWidth', 2, 'LineStyle', ':')
                        
      legend(Legends)
      
      
    xlabel('Time in s')
    ylabel('Correlation coefficient')
      
    
end

legend(Legends)

%ylim([0 1])
xlabel('Time in s')
ylabel('Correlation coefficient')


%% Mean and std of correlation across odors

figure
for ii= 1:nTrials
t =  windowRow/sampling_freq; 
corr_mean = nanmean(squeeze(CorrAcrossOdors(ii,:,:)),1);
corr_std = nanstd(squeeze(CorrAcrossOdors(ii,:,:)),1);

% corr_mean = nanmean(squeeze(CorrAcrossTrials(ii,:,:)),1);
% corr_std = nanstd(squeeze(CorrAcrossTrials(ii,:,:)),1);

%shaded area of std
curve1 = corr_mean + corr_std;
curve2 = corr_mean - corr_std;
t2 = [t, fliplr(t)];
inBetween = [curve1, fliplr(curve2)];

%plot
subplot(2, 2,ii)
    fill(t2, inBetween, 'k'); 
    alpha(0.2)
    hold on     
    plot(t,corr_mean, 'Color',my_colors(ii,:),'LineWidth',2)
    line([13.9 19],[0.7 0.7], 'LIneWidth', 2,'Color', 'k', 'LineWidth',2)
    
    ylim([-0.1 0.8])
    xlim([5 30])
    xlabel('Time in s')
    ylabel('Correlation coefficient')
    title(odorList{ii})
end


%% New trial -- 28.04.21 %%%

nOdor = length(odorList);
nTrial = length(spikeRatesArray);
responseWindow = floor([15 17] * 7.5);
responseInd = responseWindow(1):responseWindow(2);
patternArray = cellfun(@(x) mean(x(:,responseInd),2),...
                       spikeRatesArray,'UniformOutput',false);


%%%%
corrMat = zeros(nTrial,nTrial);
sm=0;
for ind=1:nTrial
    for jnd=ind:nTrial
        corrVec = analysis.calcPatternCorrelation(...
            patternArray{ind},...
            patternArray{jnd});
        corrMat(ind,jnd) = corrVec;
    end
end

% plot
fig = figure;
imagesc(helper.makeSymmetricMat(corrMat))
% caxis([0 1])
trialNumVec = linspace(ceil(nTrials/2),nTrial-floor(nTrials/2),nOdor);
xticks(trialNumVec)
xticklabels(odorList)
set(gca,'xaxisLocation','top')
set(gca,'XTick',[])
yticks(trialNumVec)
yticklabels(odorList)
colorbar
set(gca, 'FontSize', 22)
colormap jet
caxis([0 1])


%% Correlation matrix in time 
window = 2; % in frames
start = 7.0;% in seconds
stop = 20.0;% in seconds
sampling_freq = myexp.expInfo.frameRate/myexp.expInfo.nPlane;

o1 = 5; % File number
o2= 6; % File number
odor_1 = spikeRatesArray{o1};
odor_2 = spikeRatesArray{o2};


mt= 1:length(odor_1)-window; % moving time window in frames

mt_start = ceil(start*sampling_freq); % in frames
mt_stop = ceil(stop*sampling_freq); % in frames


 nOdor = length(odorList);
% nTrial = length(spikeRatesArray);
% responseWindow = floor([10 11] * 7.5);
% responseInd = responseWindow(1):responseWindow(2);
% patternArray = cellfun(@(x) mean(x(:,responseInd),2),...
%                        spikeRatesArray,'UniformOutput',false);
 
corrMatTimeSeries = zeros(size(corrMat, 1), size(corrMat,2), mt_stop-mt_start);

for t = 1: length(mt)
    patternArray = cellfun(@(x) mean(x(:,mt(t):mt(t)+2),2),...
                        spikeRatesArray,'UniformOutput',false);
    for ind=1:nTrial
        for jnd=ind:nTrial
            corrVec = analysis.calcPatternCorrelation(...
                patternArray{ind},...
                patternArray{jnd});
            corrMatTimeSeries(ind,jnd,t) = corrVec;
        end
    end
end

%% Plotting 

figure

subplot(od,trial,1)
plot(squeeze(corrMatTimeSeries(o1,o2,:)),'k*-')
