%1 is motion, 2 is colour
%colour same/different side as motion
% colour before, same and after motion
%
clear all
close all
set(0,'DefaultFigureVisible','on')
t=1;
addpath(genpath('/fs04/jo36/MRH110_EEG_2021/CSDtoolbox/'));
addpath(genpath('/fs04/jo36/MRH110_EEG_2021/eeglab13_6_5b/'));
addpath('function_programs/');
eeglab
onew=1;
no_of_bins = 3;
namestring ={'SameSide','DiffSide'};
namestringCM = {'Motion','Colour'};
addpath('function_programs/');
ssj=[1:10 12:40 42:43 45 47];
%ssj = [1:10 12:38 41:43 45 47];%quarter
%ssj = [2,7,8,9,10,12,14,15,16,17,18,19,21,22,23,25,30:37,42,45,47];
%ssj = [1:20];
%ssj=[6,38,41];%LargeN2
chn = 53;
CSD = 1;
if CSD
    yCPP = [-10,20]; ytickCPP =[-10:5:25];
    yN2pc = [-30,20]; ytickN2pc = [-30:5:20];
    yN2c= [-30,20]; ytickN2c=[-30:5:30];
    yN2i=[-30,20];ytickN2i=[-30:5:30];
    ylim = [-20 20]; ylimBeta=[-2 0.5];ytickBeta=[-2:0.3:0.5];
else
    yCPP = [-4,8]; ytickCPP =[-4:2:8];
    yN2pc = [-4,4]; ytickN2pc = [-4:2:4];
    yN2c= [-4,4]; ytickN2c=[-4:2:4];
    yN2i=[-4,4];ytickN2i=[-4:2:4];
    ylim = [-4 8];ylimBeta=[-0.1 0.1];ytickBeta=[-0.9:1:3];
end
cs = [1:2];
csdiff= [1; 2];
conds=[];
ch_beta = [8 25 8 25];
%combine them
for s = 1:size(ssj,2)
    sT = ssj(s)
    if CSD      
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_CSD_N2cSelect294394' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','subject_folder')
        allBinsAcc =allBins;
       
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_CSD_N2cSelect294394' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','STFT_time','STFT_timer','allBinsB')
        load(['Data/ERPs/group_plots_erp_diff_CSD_DDM_Beta_Bin_Motion' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','STFT_time','STFT_timer','allBinsB')
    else
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_N2cSelect294394' num2str(sT) '_' num2str(chn)],'allBins')
        allBinsAcc= allBins; clear allBins;                
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_N2cSelect294394' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','subject_folder')
        load(['Data/ERPs/group_plots_erp_diff_CSD_DDM_Beta_Bin_Motion' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','STFT_time','STFT_timer','allBinsB')
    end
    
    allBins = allBins;
    % d (third row) is colour (2) and motion(1), c (2nd row) is count and respond
    for bin=1:no_of_bins
        for d=1:length(cs)
            for c=2
                for cong=1:2
                    allstuff_a.RTs{s,c,d,cong,bin} = [allBins.Rts{1,c,d,cong,bin}];
                    allstuff_a.acc{s,c,d,cong,bin} = [allBinsAcc.acc{1,c,d,cong,bin}];
                    sizeRTS(s,d,cong,bin) = length([allstuff_a.RTs{s,c,d,cong,bin}]);
                    CPP_sides(s,d,cong,bin,:)= squeeze(mean(allBins.CPP_side(:,c,d,cong,bin,:),2));
                    CPPr_sides(s,d,cong,bin,:) = squeeze(mean(allBins.CPPr_side(:,c,d,cong,bin,:),2));
                    N2c_sides(s,d,cong,bin,:)=squeeze(mean(allBins.N2c_side(:,c,d,cong,bin,:),2));
                    N2i_sides(s,d,cong,bin,:) =squeeze(mean(allBins.N2i_side(:,c,d,cong,bin,:),2));
                    N2pc_sides(s,d,cong,bin,:)=squeeze(mean(allBins.N2pc_side(:,c,d,cong,bin,:),2));
                    
                    ERP_sides(s,d,cong,bin,:,:) = squeeze(mean(allBins.ERP_side(:,c,d,cong,bin,:,:),2));
                    ERPr_sides(s,d,cong,bin,:,:) = squeeze(mean(allBins.ERPr_side(:,c,d,cong,bin,:,:),2));
                    
                    beta_all(s,d,bin,:,:,:) = squeeze(mean(allBinsB.beta_all(:,c,d,bin,:,:),2));
                    betar_all(s,d,bin,:,:) = squeeze(mean(allBinsB.betar_all(:,c,d,bin,:,:),2));
                    
                    beta_contra_all(s,d,bin,:,:,:) = squeeze(mean(allBinsB.beta_contra(:,c,d,:,:,:),2));
                    betar_contra_all(s,d,bin,:,:) = squeeze(mean(allBinsB.betar_contra(:,c,d,:,:),2));
                    beta_ipsi_all(s,d,bin,:,:,:) = squeeze(mean(allBinsB.beta_ispi(:,c,d,:,:,:),2));
                    betar_ipsi_all(s,d,bin,:,:) = squeeze(mean(allBinsB.betar_ipsi(:,c,d,:,:),2));
                    
                    RT_sides{s,d,cong,bin}= [allBins.Rts{1,c,d,cong,bin}];
                    RT_mean(s,d,cong,bin) = mean([RT_sides{s,d,cong,bin}]);
                    % ERP_sidesO(s,:,:,:,:,:) = squeeze(allBins.ERP_sideO);
                    % ERPr_sidesO(s,:,:,:,:,:) = squeeze(allBins.ERPr_sideO);
                end
            end
        end
    end
end


%% plot reaction time lines
for sj = 1:size(allstuff_a.RTs,1)
    for cong=1:2
        for d=1:2
            for bin=1:no_of_bins
                Rt_sj_times(sj,d,cong,bin) = mean([allstuff_a.RTs{sj,2,d,cong,bin}]);
            end
        end
    end  
end
RT_group_times(:,:,:) = squeeze(mean(Rt_sj_times(:,:,:,:),1));
%% Colors

cyan        = [0.2 0.8 0.8];

mangenta    = [1 0 1];
lightgreen  = [0 1 0];
red         = rgb('DeepPink');
blue        = [0    0.7461    1.0000];
green       = [0    1.0000    0.4961];
black        = [0 0 0 ];
brown       = rgb('Chocolate');
orange      = rgb('Gold');
purple      = rgb('Magenta');
dred         = [ 0.6953    0.1328    0.1328];
dblue        = [0         0    0.5430];
dgreen       = [0    0.3906         0];
dblack       = [0 0 0];
dbrown       = rgb('Maroon');
dorange      = rgb('OrangeRed');
dpurple      =[0.5 0 1];
Cols = [dgreen;dred;dblue;green;red;blue;dpurple;purple];Cols1=Cols; Cols4=Cols;
ColsBeta=[dbrown;dorange;brown;orange];
ColsRT=[red;blue;green;purple];
RTCols=[dblack];0
alphas =[0.4 1 0.2 0.6];
alphaRT=[1 0.7 0.3 0];

close all

figurekk=0;
return
%% Topology Plots
%% ERP scalp topo
if CSD
    maplim = [-15 15];
    mapRlim = [-20 20];
    mapDlim = [-10 10];
else
    maplim = [-1.5 1.5];
    mapRlim = [-2 2];
    mapDlim = [-1 1];
end
%%%%%%%%%%%%%%%%%DA

left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];

plot_mean = zeros([64,1]);

d=1; % motion is 1 and color is 2
for bin=1
    %N2pc
    figurekk=figurekk+1;
    hf1 = figure(figurekk)
    subplot(1,1,1)
    plot_mean(left_hemi) = squeeze(mean(mean(mean(mean(ERP_sides(:,d,:,:,left_hemi,find(t>=300 & t<400))-...
        ERP_sides(:,c,:,:,right_hemi,find(t>=300 & t<400),:),1),3),4),6));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','numbers','plotchans',left_hemi);
    %N2c/N2i
    figurekk=figurekk+1;
    hf2=figure(figurekk)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,d,:,:,:,find(t>=300 & t<400)),1),3),4),6));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','numbers','plotchans',1:64);
    %CPP
    %     hf4=figure(4)
    %     subplot(1,1,1)
    %     plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=400 & t<600),:),1),3),4),6));
    %     topoplot(plot_mean,chanlocs,'maplimits', ...
    %      [min(plot_mean),max(plot_mean)], ...
    %     'electrodes','off','plotchans',1:64);
    %CPPr
    figurekk=figurekk+1;
    hf4=figure(figurekk)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERPr_sides(:,d,:,:,:,find(tr>=-20 & tr<20)),1),3),4),6));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','numbers','plotchans',1:64);
    %BetaWeird
    figurekk=figurekk+1;    
    hf5=figure (figurekk)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(beta_all(:,d,:,:,find(STFT_time>=600 & STFT_time<800)),1),2),3),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(-0.9),max(0)], ...
        'electrodes','numbers','plotchans',1:64);
    %Beta
    figurekk=figurekk+1;    
    hf6=figure (figurekk)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(betar_all(:,d,:,:,find(STFT_timer>=-20 & STFT_timer<20)),1),2),3),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(-0.65),max(0)], ...
        'electrodes','numbers','plotchans',1:64);
    colorbar
end
saveas(hf1,['Figures/N2pc' namestringCM{d} '.png'])
saveas(hf2,['Figures/N2c_c' namestringCM{d} '.png'])
%saveas(hf3,['Figures/N2i_c' namestring{onew} '.png'])
saveas(hf4,['Figures/CPP' namestringCM{d} '.png'])
saveas(hf5,['Figures/Beta' namestringCM{d} '.png'])
saveas(hf6,['Figures/betaR' namestringCM{d} '.png'])

%% ERP
% 3rd is for speed, 2nd is for respond and count
for cong=1:2
    for bin=1:no_of_bins %count only
        dataplottopo(:,:,bin) = squeeze(mean(mean(mean(ERP_sides(:,1,cong,bin,:,:,:),1),2),3));
    end
    figurekk=figurekk+1;    
    figure(figurekk)
    plottopo(dataplottopo,'chanlocs',chanlocs,'limits',[t(1) t(end) ...
        min(min(min(dataplottopo))) max(max(max(dataplottopo)))], ...
        'title',['ERP colourmotion'],'showleg','on','ydir',1)
end
%% ERPr
% 3rd is for colourmotion, 2nd is for respond and count
for cong=1:2
    for bin=1:no_of_bins %count only
        dataplotrtopo(:,:,bin) = squeeze(mean(mean(mean(ERPr_sides(:,1,cong,bin,:,:,:),1),2),3));
    end
    figurekk=figurekk+1;    
    figure(figurekk)
    plottopo(dataplotrtopo,'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
        min(min(min(dataplotrtopo))) max(max(max(dataplotrtopo)))], ...
        'title',['ERP colourmotion'],'showleg','on','ydir',1)
end
%% Beta
% beta , % 3rd is for colourmotion, 2nd is for respond and count
for bin=1:no_of_bins %count only
    dataplottopor(:,:,bin) = squeeze(mean(mean(mean(betar_all(:,1,bin,:,:,:),1),2),3));
end
figurekk=figurekk+1;
figure(figurekk)
plottopo(dataplottopor,'chanlocs',chanlocs,'limits',[STFT_timer(1) STFT_timer(end) ...
    min(min(min(dataplottopor))) max(max(max(dataplottopor)))], ...
    'title',['Beta colourmotion'],'chans',1:64,'showleg','on','ydir',1)
%% RT Behaviour
c=2
clear hs hb
for cong=1:2
    for bin = 1:no_of_bins
        for d=1:2
            temp  =[];
            for s = 1:length(ssj)
                temp = [temp [allstuff_a.RTs{s,c,d,cong,bin}]];
            end
            meanRTs{d,bin} = temp;
            clear temp
        end
    end
end
for cong=1:2
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    hold on
    cs=0;
    
    for bin=1:no_of_bins
        for d=1:2
            cs=cs+1;
            clear de_temp
            de_temp = [allstuff_a.RTs{:,c,d,cong,bin}];
            set(gca,'FontSize',16,'xlim',[300,1600],'xtick',[300:200:1600],'ylim',[0,100],'ytick',[0:20:100]);%,'ylim',[-1.5,0.5]);
            h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
            hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',ColsRT(d,:),'Linewidth',1.5);
            hb(cs) = line([mean(meanRTs{d}) mean(meanRTs{d})],[0 450],'Color',ColsRT(d,:),'Linewidth',1.5,'LineStyle','--');
        end
    end
    ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
    xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
    legend(hs,{'Motion','Colour'},'FontSize',12,'Location','NorthEast','fontweight','bold')
    saveas(hf,['Figures\RT_all_bins' namestring{cong} '.png']);
end
% check within subject
for sj=1:size(allstuff_a.RTs,1)
    a(sj)=mean([allstuff_a.RTs{sj,2,1,1}]);
    b(sj)=mean([allstuff_a.RTs{sj,2,1,2}]);
end
figurekk=figurekk+1;    
hf = figure(figurekk)
plot([a;b],'+-');
ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
%%
% make bar graphs
for d=1:2
    kk=1;
    clear mean_RTs std_RTs
    for cong=1:2
        for bin=1:no_of_bins
            mean_RTs(kk) = mean([allstuff_a.RTs{:,2,d,cong,bin}]);
            std_RTs(kk) = std([allstuff_a.RTs{:,2,d,cong,bin}])/sqrt(length(ssj));
            kk=kk+1;
        end
    end
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    hold on
    cs=0;
    for c=1:length(mean_RTs)
        hs1(c) = bar(c,mean_RTs(c),'Linewidth',1.5);
        set(hs1(c),'FaceColor',ColsRT(max(1,ceil(c/3)),:));
        hs2(c) = errorbar(c,mean_RTs(c),std_RTs(c)'.');
        hs2(c).Color = dblack;
        hs2(c).LineWidth = 2;
    end
    ylabel({'Response Time (ms)'},'FontSize',16,'fontweight','bold')
    axis([0 7 300 1200]);
    xticks([1 2 3 4 5 6]);xticklabels({'SS Fast','SS Med','SS Slow',...
        'DS Fast', 'DS Med','DS Slow'});
    set(gca,'Fontsize',10)
    yticks(300:300:1200);
    saveas(hf,['Figures\RT_bar_graphs' namestringCM{d} '.png']);
end
%% Accuracy Behaviour
c=2
clear hs hb
for cong=1:2
    for bin = 1:no_of_bins
        for d=1:2
            temp  =[];
            for s = 1:length(ssj)
                tempAcc = [allstuff_a.acc{s,2,d,cong,bin}];
                Acc(s,2,d,cong,bin) = sum(tempAcc)/length(tempAcc);
            end
            clear temp tempAcc
        end
    end
end
% make bar graphs
for d=1:2
    kk=1;
    clear mean_RTs std_RTs
    for cong=1:2
        for bin=1:no_of_bins
            meanAcc(kk) = mean(Acc(:,2,d,cong,bin));
            stdAcc(kk) = std(Acc(:,2,d,cong,bin))/sqrt(length(ssj));
            kk=kk+1;
        end
    end
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    hold on
    cs=0;
    for c=1:length(meanAcc)
        hs1(c) = bar(c,meanAcc(c),'Linewidth',1.5);
        set(hs1(c),'FaceColor',ColsRT(max(1,ceil(c/3)),:));
        hs2(c) = errorbar(c,meanAcc(c),stdAcc(c)'.');
        hs2(c).Color = dblack;
        hs2(c).LineWidth = 2;
    end
    ylabel({'Accuracy (%)'},'FontSize',16,'fontweight','bold')
    axis([0 7 0.7 1]);
    xticks([1 2 3 4 5 6]);xticklabels({'SS Fast','SS Med','SS Slow',...
        'DS Fast', 'DS Med','DS Slow'});
    set(gca,'Fontsize',10)
    yticks(0.7:0.1:1);
    saveas(hf,['Figures\Acc_bar_graphs' namestringCM{d} '.png']);
end
%% N2 graphs
%shadedErrorBar
%N2pc
clear hf hss1 hss4
colIndx=1;
for cong=1:2
    kk=1;
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    for d=1
        for bin=1:no_of_bins
            N2pc_c = squeeze(mean(mean(N2pc_sides(:,d,cong,bin,:),3),2));
            meanN2pc = mean(N2pc_c,1);
            stdN2pc = 0*std(N2pc_c,1)/sqrt(size(N2pc_c,1));
            hold on
            hss1(kk) = shadedErrorBar(t,meanN2pc,stdN2pc,'lineprops',{'Color',[Cols(kk,:) alphaRT(bin)],'LineWidth',3,'LineStyle','-'});
            set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],'ylim',yN2pc,'ytick',ytickN2pc);%,'ylim',[-1.5,0.5]);
            hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
                ylim,'Color',[Cols(kk,:) alphaRT(bin)],'LineWidth',1.5,'LineStyle','--');
            ylabel('N2pc Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
            xlabel('Time(ms)','FontName','Arial','FontSize',16','fontweight','bold')
            line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
            %     legend(h,side_tags, ...
            %         'FontSize',16,'Location','NorthWest');
            kk=kk+1;
            colIndx=colIndx+1;            
        end
    end
    legend([hss1(1).mainLine, hss1(2).mainLine, hss1(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
    saveas(hf,['Figures\N2pc' namestringCM{d} '_' namestring{cong} '.png']);
end
%% N2c
% close all
% Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
% alphaRT=[1 1 1 0]; yN2c = [-10 10]; ytickN2c = -10:2:10;ylim = yN2c;
% clear hss3 hss4
% for ssjs=1:43
%     colIndx=1;
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for d=1
%             hold on
%             for bin=1:no_of_bins
%                 N2c_c_grp = squeeze(mean(mean(N2c_sides(ssjs,d,cong,bin,:),2),3))';
%                 meanN2c = mean(N2c_c_grp,1);
%                 stdN2c = 0*ones([1 length(N2c_c_grp)]);%0*std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
%                 hold on
%                 %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
%                 hss3(kk) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',3,'LineStyle','-'});
%                 set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
%                     'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
%                 hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
%                     yN2c,'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',1.5,'LineStyle','--');
%                 ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
%                 line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%                 %     legend(h,side_tags, ...
%                 %         'FontSize',16,'Location','NorthWest');
%                 kk=kk+1;
%                 colIndx=colIndx+1;
%             end
%         end
%         legend([hss3(1).mainLine, hss3(2).mainLine, hss3(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
%         %legend([hss3(1).mainLine, hss3(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartN2c\N2c' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs) '.png']);
%         
%         % finding the timing
%         ts = t(t>0 & t<600);
%         gs = ts(find(meanN2c(t>0 & t<600) == min(meanN2c(t>0 & t<600))))
%     end
% end
%% N2c collapsed across bins
% close all
% Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
% alphaRT=[1 1 1 0]; yN2c = [-10 10]; ytickN2c = -10:2:10;ylim = yN2c;
% clear hss3 hss4
% for ssjs=1:43
%     colIndx=1;
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for d=1
%             hold on
%             N2c_c_grp = squeeze(mean(mean(mean(N2c_sides(ssjs,d,cong,:,:),2),3),4))';
%             meanN2c = mean(N2c_c_grp,1);
%             stdN2c = 0*ones([1 length(N2c_c_grp)]);%0*std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
%             hold on
%             %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
%             hss3(kk) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[Cols(2,:) alphaRT(2)],'LineWidth',3,'LineStyle','-'});
%             set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
%                 'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
%             %hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
%             %    yN2c,'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',1.5,'LineStyle','--');
%             ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
%             xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
%             line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
%             line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%             %     legend(h,side_tags, ...
%             %         'FontSize',16,'Location','NorthWest');
%             kk=kk+1;
%             colIndx=colIndx+1;
%         end
%     legend([hss3(1).mainLine],{'N2c'},'FontSize',12,'Location','NorthWest')
%     %legend([hss3(1).mainLine, hss3(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
%     saveas(hf,['Figures\singlePartmeanN2c\N2c' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs) '.png']);
%     
%     % finding the timing
%     ts = t(t>0 & t<600);
%     gs = ts(find(meanN2c(t>0 & t<600) == min(meanN2c(t>0 & t<600))))
%     end
% end
%%
close all
Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
alphaRT=[1 1 1 0]; yN2c = [-15 4]; ytickN2c = -15:2:4;ylim = yN2c;
clear hss3 hss4
colIndx=1;
for cong=1:2
    kk=1;
    figurekk=figurekk+1;
    hf = figure(figurekk)
    for d=1
        hold on
        for bin=1:no_of_bins
            N2c_c_grp = squeeze(mean(mean(mean(N2c_sides(:,d,cong,bin,:),2),3),4));
            meanN2c = mean(N2c_c_grp,1);
            stdN2c = 0*std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
            hold on
            %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
            hss3(kk) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',3,'LineStyle','-'});
            set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
                'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
            hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
                yN2c,'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',1.5,'LineStyle','--');
            ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
            xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
            line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
            %     legend(h,side_tags, ...
            %         'FontSize',16,'Location','NorthWest');
            kk=kk+1;
            colIndx=colIndx+1;
        end
    end
    legend([hss3(1).mainLine, hss3(2).mainLine, hss3(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
    %legend([hss3(1).mainLine, hss3(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
    saveas(hf,['Figures\N2c' namestringCM{d} '_' namestring{cong} '.png']);
    
    % finding the timing
    ts = t(t>0 & t<600);
    gs = ts(find(meanN2c(t>0 & t<600) == min(meanN2c(t>0 & t<600))))
end
%% N2i
% singleparticipants
% close all
% Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
% alphaRT=[1 1 1 0]; yN2i = [-10 10]; ytickN2i = -10:2:10;ylim = yN2i;
% clear hss3 hss4
% for ssjs=1:43
%     colIndx=1;
%     for cong=1:2
%         kk=1;
%         clear hss3
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for bin = 1:no_of_bins
%             hold on
%             for d=1
%                 N2i_c_grp = squeeze(mean(mean(N2i_sides(ssjs,d,cong,bin,:,:),2),3))';
%                 meanN2i = mean(N2i_c_grp,1);
%                 stdN2i = 0*ones([1 length(N2i_c_grp)]);
%                 hold on
%                 %hss1(kk) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
%                 hss3(kk) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',3,'LineStyle','-'});
%                 set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
%                     'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
%                 hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
%                     yN2i,'Color',[Cols(colIndx,:) alphaRT(d)],'LineWidth',1.5,'LineStyle','--');
%                 ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
%                 line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%                 kk=kk+1;
%                 colIndx=colIndx+1;
%                 %     legend(h,side_tags, ...
%                 %         'FontSize',16,'Location','NorthWest');
%             end
%         end
%         legend([hss3(1).mainLine, hss3(2).mainLine, hss3(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartN2i\N2i' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs) '.png']);
%     end
% end
%% N2i collapsed across bins
% close all
% Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
% alphaRT=[1 1 1 0]; yN2i = [-30 30]; ytickN2i = -30:2:30;ylim = yN2i;
% clear hss3 hss4
% for ssjs=1:43
%     for cong=1:2
%         kk=1;
%         clear hss3
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%             hold on
%             for d=1
%                 N2i_c_grp = squeeze(mean(mean(mean(N2i_sides(ssjs,d,cong,:,:,:),2),3),4))';
%                 meanN2i = mean(N2i_c_grp,1);
%                 stdN2i = 0*ones([1 length(N2i_c_grp)]);
%                 hold on
%                 %hss1(kk) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
%                 hss3(kk) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[Cols(5,:) alphaRT(2)],'LineWidth',3,'LineStyle','-'});
%                 set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
%                     'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
%                 ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
%                 line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
%                 line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%                 kk=kk+1;
%                 colIndx=colIndx+1;
%                 %     legend(h,side_tags, ...
%                 %         'FontSize',16,'Location','NorthWest');
%             end
%         legend([hss3(1).mainLine],{'N2i'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartmeanN2i\N2i' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs) '.png']);
%     end
% end
%%
colIndx=1;
for cong=1:2
    kk=1;
    clear hss3
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    for bin = 1:no_of_bins
        hold on
        for d=1
            N2i_c_grp = squeeze(mean(mean(N2i_sides(:,d,cong,bin,:,:),2),3));
            meanN2i = mean(N2i_c_grp,1);
            stdN2i = 0*std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));
            hold on
            %hss1(kk) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
            hss3(kk) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[Cols(colIndx,:) alphaRT(bin)],'LineWidth',3,'LineStyle','-'});
            set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
                'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
            hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
                yN2i,'Color',[Cols(colIndx,:) alphaRT(d)],'LineWidth',1.5,'LineStyle','--');
            ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
            xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
            line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
            kk=kk+1;
            colIndx=colIndx+1;
            %     legend(h,side_tags, ...
            %         'FontSize',16,'Location','NorthWest');
        end
    end
    legend([hss3(1).mainLine, hss3(2).mainLine, hss3(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
    saveas(hf,['Figures\N2i' namestringCM{d} '_' namestring{cong} '.png']);
end
return
%% CPP single trials
% %% CPP
% close all
% yCPP = [-10 15]; ytickCPP = -10:5:15;ylim = yCPP;
% for ssjs=1:43
%     colIndx=1;
%     clear hss11
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for d=1
%             for bin = 1:no_of_bins
%                 hold on
%                 CPP_c = squeeze(mean(mean(CPP_sides(ssjs,d,cong,bin,:,:),3),2))';
%                 meanCPP = mean(CPP_c,1);
%                 stdCPP = 0*ones([1 length(CPP_c)]);%0*std(CPP_c,1)/sqrt(size(CPP_c,1));
%                 clear CPP_c
%                 % beta_sides = squeeze(mean(mean(beta_side_all(:,c,:),1),5));
%                 hold on
%                 hss1(kk) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(colIndx,:)],'LineWidth',3,'LineStyle','-'});
%                 set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
%                 ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
%                 xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
%                 line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
%                 line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%                 kk=kk+1;
%                 colIndx=colIndx+1;
%                 %     legend(h,side_tags, ...
%                 %         'FontSize',16,'Location','NorthWest');
%             end
%         end
%         legend([hss1(1).mainLine hss1(2).mainLine hss1(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartCPP\CPP'  namestringCM{d} '_' namestring{cong} '_' num2str(ssjs)  '.png']);
%     end
% end
% %% CPP mean single part
% close all
% yCPP = [-10 15]; ytickCPP = -10:5:15;ylim = yCPP;
% for ssjs=1:43
%     clear hss11
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for d=1
%             hold on
%             CPP_c = squeeze(mean(mean(mean(CPP_sides(ssjs,d,cong,:,:,:),3),2),4))';
%             meanCPP = mean(CPP_c,1);
%             stdCPP = 0*ones([1 length(CPP_c)]);%0*std(CPP_c,1)/sqrt(size(CPP_c,1));
%             clear CPP_c
%             % beta_sides = squeeze(mean(mean(beta_side_all(:,c,:),1),5));
%             hold on
%             hss1(kk) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(2,:)],'LineWidth',3,'LineStyle','-'});
%             set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
%             ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
%             xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
%             line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
%             line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%             kk=kk+1;
%             colIndx=colIndx+1;
%             %     legend(h,side_tags, ...
%             %         'FontSize',16,'Location','NorthWest');
%         end
%         legend([hss1(1).mainLine],{'CPP'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartmeanCPP\CPP'  namestringCM{d} '_' namestring{cong} '_' num2str(ssjs)  '.png']);
%     end
% end
%% CPP
colIndx=1;
for cong=1:2
    kk=1;
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    for d=1
        for bin = 1:no_of_bins
            hold on
            CPP_c = squeeze(mean(mean(CPP_sides(:,d,cong,bin,:,:),3),2));
            meanCPP = mean(CPP_c,1);
            stdCPP = 0*std(CPP_c,1)/sqrt(size(CPP_c,1));
            clear CPP_c
            % beta_sides = squeeze(mean(mean(beta_side_all(:,c,:),1),5));
            hold on
            hss1(kk) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(colIndx,:)],'LineWidth',3,'LineStyle','-'});
            set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
            ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
            line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
            kk=kk+1;
            colIndx=colIndx+1;
            %     legend(h,side_tags, ...
            %         'FontSize',16,'Location','NorthWest');
        end
    end
    legend([hss1(1).mainLine hss1(2).mainLine hss1(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
    saveas(hf,['Figures\CPP'  namestringCM{d} '_' namestring{cong} '.png']);
end
% %% single trial CPPr
% close all
% for ssjs=1:43
%     colIndx=1;
%     clear hss11
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         for bin = 1:no_of_bins
%             hold on
%             for d=1
%                 CPPr_c = squeeze(mean(mean(CPPr_sides(ssjs,d,cong,bin,:,:),3),2))';
%                 meanCPPr = mean(CPPr_c,1);
%                 stdCPPr = 0*ones([1 length(CPPr_c)]);%0*std(CPPr_c,1)/sqrt(size(CPPr_c,1));
%                 clear CPPr_c
%                 hold on
%                 hss11(kk) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(colIndx,:)],'LineWidth',3,'LineStyle','-'});
%                 set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
%                     'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
%                 hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
%                     [-10 30],'Color',[RTCols alphaRT(d)],'LineWidth',1.5,'LineStyle','-');
%                 ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
%                 xlabel('from response (ms)','FontName','Arial','FontSize',16)
%                 line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
%                 line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%                 kk=kk+1;
%                 colIndx=colIndx+1;
%                 %     legend(h,side_tags, ...
%                 %         'FontSize',16,'Location','NorthWest');
%             end
%         end
%         legend([hss11(1).mainLine hss11(2).mainLine hss11(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartCPPr\CPPr' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs)  '.png']);
%     end
% end
% %% single trial CPPr mean
% close all
% for ssjs=1:43
%     clear hss11
%     for cong=1:2
%         kk=1;
%         figurekk=figurekk+1;
%         hf = figure(figurekk)
%         hold on
%         for d=1
%             CPPr_c = squeeze(mean(mean(mean(CPPr_sides(ssjs,d,cong,:,:,:),3),2),4))';
%             meanCPPr = mean(CPPr_c,1);
%             stdCPPr = 0*ones([1 length(CPPr_c)]);%0*std(CPPr_c,1)/sqrt(size(CPPr_c,1));
%             clear CPPr_c
%             hold on
%             hss11(kk) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(5,:)],'LineWidth',3,'LineStyle','-'});
%             set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
%                 'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
%             ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
%             xlabel('from response (ms)','FontName','Arial','FontSize',16)
%             line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
%             line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
%             kk=kk+1;
%             colIndx=colIndx+1;
%             %     legend(h,side_tags, ...
%             %         'FontSize',16,'Location','NorthWest');
%         end
%         legend([hss11(1).mainLine],{'CPPr'},'FontSize',12,'Location','NorthWest')
%         saveas(hf,['Figures\singlePartmeanCPPr\CPPr' namestringCM{d} '_' namestring{cong} '_' num2str(ssjs)  '.png']);
%     end
% end
%% CPPr
colIndx=1;
clear hss11
for cong=1:2
    kk=1;
    figurekk=figurekk+1;    
    hf = figure(figurekk)
    for bin = 1:no_of_bins
        hold on
        for d=1
            CPPr_c = squeeze(mean(mean(CPPr_sides(:,d,cong,bin,:,:),3),2));
            meanCPPr = mean(CPPr_c,1);
            stdCPPr = 0*std(CPPr_c,1)/sqrt(size(CPPr_c,1));
            clear CPPr_c
            hold on
            hss11(kk) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(colIndx,:)],'LineWidth',3,'LineStyle','-'});
            set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
                'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
            hss4(kk)=line([mean(mean(RT_group_times(d,cong,bin),1),2) mean(mean(RT_group_times(d,cong,bin),1),2)],...
                [-10 30],'Color',[RTCols alphaRT(d)],'LineWidth',1.5,'LineStyle','-');
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('from response (ms)','FontName','Arial','FontSize',16)
            line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
            kk=kk+1;
            colIndx=colIndx+1;            
            %     legend(h,side_tags, ...
            %         'FontSize',16,'Location','NorthWest');
        end
    end
    legend([hss11(1).mainLine hss11(2).mainLine hss11(3).mainLine],{'Fast','Medium','Slow'},'FontSize',12,'Location','NorthWest')
    saveas(hf,['Figures\CPPr' namestringCM{d} '_' namestring{cong} '.png']);
end
return
%% beta
clear h1 h2 hss4 hf hss3
% contra beta and ipsi beta x distractors resp
figurekk=figurekk+1;    
hf = figure(figurekk)
hold on
d=2;
for bin=1:2
    for c=2
        contra_beta_grp = squeeze(mean(beta_all(:,d,c,8,:,:),2));
        % ipsi_beta_grp = squeeze(mean(beta_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        %   h2((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
        %       {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),
        set(gca,'FontSize',16,'xlim',[-100 1400],'ytick',[-2:0.3:0.5],'xtick',[-100, 0:200:1400],'ylim',[-2 0.5]);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        %  xlabel('Time (ms)','FontName','Arial','FontSize',16)
        % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss3((bin-1)*2+c)=line([mean(RT_group_times(d,cong,bin),2), mean(RT_group_times(d,cong,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]);
        hss4((bin-1)*2+c)=line([mean(RT_group_times(d,cong,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]);
    end
end
legend([h1.mainLine],{'Motion','Colour'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\beta_stim' namestringCM{d} '_' namestring{onew} '.png']);
figurekk=figurekk+1;    
hf = figure(figurekk)
hold on
for bin=1:1
    for c=1:2
        contra_beta_grp = squeeze(mean(betar_all(:,d,c,8,:,:),2));
        % ipsi_beta_grp = squeeze(mean(betar_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        %h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
        %    {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),
        set(gca,'FontSize',16,'xlim',[-600 100],'ytick',[-2:0.3:0.5],'xtick',-600:200:100,'ylim',[-0.8 0.5]);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','-','Color',[RTCols alphaRT(c)]);
    end
end
legend([h1.mainLine],{'Motion','Colour'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\beta_resp' namestring{onew} '.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% N2 and CPP TI
%shadedErrorBar
%N2pc
clear hf hss1
hf = figure(1)
for bin=1:1
    for c=1:2
        N2pc_c = squeeze(mean(mean(N2pc_sidesT(:,c,:,:),4),2));
        meanN2pc = mean(N2pc_c,1);
        stdN2pc = std(N2pc_c,1)/sqrt(size(N2pc_c,1));
        hold on
        (c-1)*2+bin
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanN2pc,stdN2pc,'lineprops',{'Color',[dpurple alphaRT((bin-1)*2+c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],'ylim',yN2pc,'ytick',ytickN2pc);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            ylim,'Color',[dpurple alphaRT((bin-1)*2+c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2pc Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time(ms)','FontName','Arial','FontSize',16','fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine, hss1(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2pcxdistT' namestring{onew} '.png']);
% N2c
clear hss3
hf = figure(3)
for bin = 1:1
    hold on
    for c=1:2
        N2c_c_grp = squeeze(mean(mean(N2c_sidesT(:,c,:,:),4),2));
        meanN2c = mean(N2c_c_grp,1);
        stdN2c = std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[dred alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2c,'Color',[dred alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2cxdistT' namestring{onew} '.png']);
% N2i
clear hss3
hf = figure(2)
for bin = 1:1
    hold on
    for c=1:2
        N2i_c_grp = squeeze(mean(mean(N2i_sidesT(:,c,:,:),4),2));
        meanN2i = mean(N2i_c_grp,1);
        stdN2i = std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[dblue alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2i,'Color',[dblue alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2ixdistT' namestring{onew} '.png']);


%% CPP_CPPr
hf = figure(4)
for bin = 1:1
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPP_sidesT(:,c,:,:),4),2));
        meanCPP = mean(CPP_c,1);
        stdCPP = std(CPP_c,1)/sqrt(size(CPP_c,1));
        clear CPP_c
        % beta_sides = squeeze(mean(mean(beta_side_all(:,c,:),1),5));
        hold on
        (c-1)*2+bin
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine hss1(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPxdistT' namestring{onew} '.png']);
clear hss11
hf = figure(5)
for bin = 1:1
    hold on
    for c=1:2
        CPPr_c = squeeze(mean(mean(CPPr_sidesT(:,c,:,:),4),2));
        meanCPPr = mean(CPPr_c,1);
        stdCPPr = std(CPPr_c,1)/sqrt(size(CPPr_c,1));
        clear CPPr_c
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from response (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'Motion','Colour'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPrxdistT' namestring{onew} '.png']);
%% Pre target Alpha power
alpha_ch = [46 47 48 49 50];
clear hss11
hf = figure(10)
for bin = 1:1
    hold on
    for c=1:2
        alpha_c = squeeze(mean(mean(alpha_sides(:,c,alpha_ch,:),3),2));
        meanalpha = mean(alpha_c,1);
        stdalpha = std(alpha_c,1)/sqrt(size(alpha_c,1));
        clear alpha_c
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0 16],'ytick',[0:2:16]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('time from stimulus (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\alpha_sides' namestring{onew} '.png']);
clear hss11
hf = figure(12)
for bin = 1:1
    hold on
    for c=1:2
        alpha_c = squeeze(mean(mean(alpha_asym_sides(:,c,[19],:),3),2));
        meanalpha = mean(alpha_c,1);
        stdalpha = std(alpha_c,1)/sqrt(size(alpha_c,1));
        clear alpha_c
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0.1 0.1],'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-0.5 0.5],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('time from stimulus (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\alpha_sides_asym' namestring{onew} '.png']);
