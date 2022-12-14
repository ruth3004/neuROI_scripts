
clear all 
close all 

%% Load experiment
% resultDir = 'C:\Users\montruth\Documents\Work backup\20200924\f3\results\';
% expName = 'Ruth-20200924-f3';
% planeNum = 3;
resultDir = '\\tungsten-nas.fmi.ch\tungsten\scratch\gfriedri\montruth\2P_RawData\2022-01-17\f2\results';
expName = '20220117_RM0012_128hpf_fP16_f2';
planeNum = 2;

planeString = NrModel.getPlaneString(planeNum);
predictionDir = fullfile(resultDir,'spike_deconvolution',planeString, 'dF_traces_doublebl');

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

dFtracesArray = cell(size(fileNameArray));
appendix = '_dF_traces';
for k=1:length(fileNameArraySorted)
    fileName = fileNameArraySorted{k};
    spikeRatesFilePath = fullfile(predictionDir,...
        append(replace(fileName, '.tif',appendix)));
    foo = load(spikeRatesFilePath);
    dFtracesArray{k} = foo.dF_traces;
end




%% Define odors to correlate

o1 = 1; % File number
o2= 2; % File number


% Plot correlation across time between two odors
sampling_freq = myexp.expInfo.frameRate/myexp.expInfo.nPlane;
odor_1 = dFtracesArray{o1};
odor_2 = dFtracesArray{o2};

figure
for n = 1:size(odor_2,1)
     plot((1:size(odor_2,2))/sampling_freq, odor_2(n,:)+ n)
     %plot((1:size(odor_2,2)), odor_2(n,:)+ n)
    hold on 
end

xlabel('Time in frames')
axis tight



% figure
% for n = 1:size(spike_rates,1)
%     plot((1:size(spike_rates,2))/sampling_freq, spike_rates(n,:)+ n)
%     hold on 
% end
% 
% xlabel('Time in s')
% axis tight

%% Define window parameters for correlations
window = 2; % in frames
start = 10.0;    % in seconds
stop = 20.0;     % in seconds
response = 13*sampling_freq; % in frames

%% Correlation matrix of two odors 

mt= 1:length(odor_1)-window; % in frames
Mov_Corr_Mat = zeros(length(mt),1);

mt1 = ceil(start*sampling_freq); % in frames
mt2 = ceil(stop*sampling_freq); % in frames

for ii = 1:length(mt)
       
    time_start = mt(ii);%ceil(mt(ii)*sampling_freq); % in frames

    rates_odor_1 = odor_1(:,time_start:time_start+window-1);
    rates_odor_2 = odor_2(:,time_start:time_start+window-1);

    mean_rates_odor_1 = mean(rates_odor_1,2,'omitnan');
    mean_rates_odor_2 = mean(rates_odor_2,2,'omitnan');

    Mov_Corr_Mat(ii) =  corr2(mean_rates_odor_2,mean_rates_odor_1);
end

figure
imagesc(Mov_Corr_Mat)
view([-90 -90])
colorbar
caxis([0 1])

hold on 

line([0 2], [mt1 mt1], 'Color', 'k', 'LineWidth',2) 
line([0 2], [mt2 mt2], 'Color', 'r','LineWidth',2) 
line([0 2], [response response], 'Color', 'm','LineWidth',2, 'LineStyle', '--') 
%ylim([30 60])
set(gca, 'YTick', 0:sampling_freq*5:length(mt))
y_tick = get(gca, 'YTick');
set(gca, 'YTicklabel',y_tick/sampling_freq)


title(['Correlation of ', odorList{o1},...
    ' and ', odorList{o2}, ' at frames ',...
    num2str(mt1) ' and ' num2str(mt2), ' (window: ',num2str(window),')'])

ylabel('Time in s')

%% Correlation matrix of all odors 

mean_rates_all = zeros(size(dFtracesArray{1},1),length(dFtracesArray));
time_start = ceil(mt1); % in frames

   for jj = 1:length(dFtracesArray)
        odor=dFtracesArray{jj};
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

mean_rates_all = zeros(size(dFtracesArray{1},1),length(dFtracesArray));
time_start = ceil(mt2); % in frames

   for jj = 1:length(dFtracesArray)
        odor=dFtracesArray{jj};
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
     for nl = 3.5:3:24
        line([nl nl],[0 24.5],'LineWidth',2, 'Color', [1 1 1])
        line([0 24.5], [nl nl],'LineWidth',2, 'Color', [1 1 1])
     end
 axis square   
 title(['Correlation at frame ', num2str(mt1),' with window = ', num2str(window),' frames']); 

f(2) = subplot(1,3,2);
    imagesc(B)
    colorbar
    caxis([0 1])
    hold on
     for nl = 3.5:3:24
        line([nl nl],[0 24.5],'LineWidth',2, 'Color', [1 1 1])
        line([0 24.5], [nl nl],'LineWidth',2, 'Color', [1 1 1])
     end
     
axis square
 title(['Correlation at frame', num2str(mt2),' with window = ', num2str(window),' frames']); 
  
     
f(3) = subplot(1,3,3);
    imagesc(B-A)
    colorbar
    
    mycolors = zeros(100,3);
    mycolors(1:50,3) = linspace(1,0,50);
    mycolors(:,2) = zeros(100,1);
    mycolors(51:100,1) = linspace(0,1,50);
    colormap(f(3),mycolors)

     for nl = 3.5:3:24
        line([nl nl],[0 24.5],'LineWidth',2, 'Color',[1 1 1])
        line([0 24.5], [nl nl],'LineWidth',2, 'Color', [1 1 1])
     end
axis square
    
     
set(f, 'YTick',2:3:24, 'YTickLabel',odorList,...
         'XTick',2:3:24,'XTickLabel',odorList)
    xtickangle(f,45)
     
 title(['Correlation difference ', num2str(mt2-mt1),'frames (', num2str((mt2-mt1)/sampling_freq),'s).']);      

%% Correlation rows
windowRow = mt1:1:mt2;
mean_rates_all = zeros(length(windowRow),...               %timesteps
                    size(dFtracesArray{1},1),...      %spike rates
                    length(dFtracesArray));   %odors
                
CorrAllMat = zeros(length(windowRow),...               %timesteps
                    length(dFtracesArray),... %correlation matrix
                    length(dFtracesArray));                  

figure
for win = 1:length(windowRow)

    time_start = windowRow(win); % in frames

       for jj = 1:length(dFtracesArray)
            odor=dFtracesArray{jj};
            rates_odor = odor(:,time_start:time_start+window-1);
            mean_rates_all(win,:,jj) = mean(rates_odor,2,'omitnan');

       end
    CorrAllMat(win,:,:) = corrcoef(squeeze(mean_rates_all(win,:,:)));


    imagesc(squeeze(CorrAllMat(win,:,:)))
    colorbar
    caxis([0 1])

     hold on
     for nl = 3.5:3:24
        line([nl nl],[0 24.5],'LineWidth',2, 'Color', [1 1 1])
        line([0 24.5], [nl nl],'LineWidth',2, 'Color', [1 1 1])
     end
     
 axis square   
 t= char(['Correlation with window ', num2str(window),' at frame ', num2str(win),' (',num2str((win+mt1)/7.5),'s)']);
 title(t); 
 set(gca, 'YTick',2:3:24, 'YTickLabel',odorList,...
         'XTick',2:3:24,'XTickLabel',odorList)
    xtickangle(gca,45)
    
% saveas(fig, [predictionDir, t,'.png'])
drawnow
pause(0.5)
end

%% Average of correlation matrices per trial 


avgmean_rates_all = zeros(length(windowRow),...         %timesteps
                    size(dFtracesArray{1},1),...      %spike rates
                    length(dFtracesArray)/nTrials);   %odors

avgCorrAllMat = zeros(length(windowRow),...             %timesteps
                    length(dFtracesArray)/nTrials,... %correlation matrix
                    length(dFtracesArray)/nTrials);                  


for win= 1:length(windowRow)

    avgmean_rates_all(win,:,:) = nanmean(permute(reshape(...
                squeeze(mean_rates_all(win,:,:)), ...
                size(squeeze(mean_rates_all(win,:,:)),1),3,[]), ...
                [1 3 2]),3);
    avgCorrAllMat(win,:,:)= corrcoef(squeeze(avgmean_rates_all(win,:,:)));


end
%% Ploting average correlation matrices (odors)
% Averaged correlation matrix across different odors 
CorrAcrossTime = zeros(length(windowRow),1);
Legends = cell(1);
figure
for ii = 1:8
    subplot(2,4,ii)
    this_odor = ii; 
    counter = 1;
    my_colors = summer(8);
    % my_colors = brighten(my_colors,-.8)
     
     for od = 1:length(odorList)

        if this_odor ~= od 
            for cc = 1:size(avgmean_rates_all,1)         
                CorrAcrossTime(cc)= corr2(avgmean_rates_all(cc,:,this_odor),...
                            avgmean_rates_all(cc,:,od));

            end
            Legends{counter} = [odorList{this_odor},'-', odorList{od}];
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
CorrAcrossTrials = zeros(nTrials,length(windowRow));
my_colors = summer(8);
figure
trialarray = 1:3:length(fileNameArray);
 for t = 1:length(trialarray)
    trial = trialarray(t);
    subplot(length(odorList)/2,length(odorList)/4,t)
    Legends = cell(1);
    CorrAcrossTrials = zeros(length(windowRow),1);
     for cc = 1:size(avgmean_rates_all,1)     
        
        CorrAcrossTrials(1,cc)= corr2(squeeze(mean_rates_all(cc,:,trial)),...
            squeeze(mean_rates_all(cc,:,trial+1)));
        CorrAcrossTrials(2,cc)= corr2(squeeze(mean_rates_all(cc,:,trial)),...
            squeeze(mean_rates_all(cc,:,trial+2)));
        CorrAcrossTrials(3,cc)= corr2(squeeze(mean_rates_all(cc,:,trial+1)),...
            squeeze(mean_rates_all(cc,:,trial+2)));
       
      end
      Legends{1} = [odorList{(trial-mod(trial,3))/3+1},' 1-2'];
      plot(windowRow/sampling_freq,CorrAcrossTrials(1,:),...
                'Color',my_colors((trial-mod(trial,3))/3+1,:), 'LineWidth', 2)
      hold on 
      Legends{2} = [odorList{(trial-mod(trial,3))/3+1},' 1-3'];
      plot(windowRow/sampling_freq,CorrAcrossTrials(2,:),...
                'Color',my_colors((trial-mod(trial,3))/3+1,:), 'LineWidth', 2, 'LineStyle', ':')
      

      Legends{3} = [odorList{(trial-mod(trial,3))/3+1},' 2-3'];
      plot(windowRow/sampling_freq,CorrAcrossTrials(3,:),...
                'Color',my_colors((trial-mod(trial,3))/3+1,:), 'LineWidth', 2, 'LineStyle', '--')
            
            
      legend(Legends)
      
      
    ylim([0 1])
    xlabel('Time in s')
    ylabel('Correlation coefficient')
      
    
end

legend(Legends)
ylim([0 1])
xlabel('Time in s')
ylabel('Correlation coefficient')





