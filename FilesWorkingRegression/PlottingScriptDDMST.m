% single-trial mediation
clear all
close all
clc
addpath(genpath('function_programs'));
addpath(genpath('../MediationToolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath(genpath('../BRAVO2-master/'));
%addpath(genpath('../spm5/'));
% define subject
singletrial=1;zscoreTrials=0;
no_of_bins=3;
if singletrial; msoffset=5; msBoffset=8; else;msoffset=0; msBoffset=0; end; % to omit the ERPs from the single trial stats
chn = 53;
CSD = 1;
CPPslopeDA=[];CPPslopeDP=[];CPPtrajDA=[];CPPtrajDP=[];
CPPrslopeDA=[];CPPrslopeDP=[];CPPrtrajDA=[]; CPPrtrajDP =[];
CPPtraj=[];CPPrtraj=[];
condsAll = [];
mediationsAll=[];mediationsAll1=[];mediationsAll2=[];
alphaPower = [];
oldnew={'','old'};
subjIndx = 0;
subjCount=0;
allsubj=[1:10 12:40 42:43 45 47];
ch_CPP = 53; subjIndxStart=1;
Tnames = table();
tempTable=table();
CPPAllss =[]; CPPrAllss = []; RTAllss =[]; stimAllss = []; congAllss=[]; subjectss=[];N2iAllss =[]; N2cAllss=[];responsess=[];RTAllssz=[];
for s=1:size(allsubj,2)
    subjCount=subjCount+1;
    sT = allsubj(s)
    subjIndx = subjCount;
    load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_CSD_N2cSelect296396' num2str(sT) '_' num2str(chn)],'N2cAll','N2iAll','CPPAll','CPPrAll','RTAll','stimAll','congAll','responseAll','mediations','mediationsS','chanlocs','t','tr','subject_folder');
    N2cAllss = [N2cAllss; N2cAll];
    N2iAllss = [N2iAllss; N2iAll];
    CPPAllss = [CPPAllss; CPPAll];
    CPPrAllss = [CPPrAllss; CPPrAll];
    RTAllss = [RTAllss ;RTAll];
    RTAllssz = [RTAllssz; zscore(RTAll)];
    stimAllss = [stimAllss; stimAll];
    congAllss = [congAllss; congAll];
    responsess = [responsess; responseAll];
    sAll = s*ones([1 length(congAll)]);
    subjectss = [subjectss; sAll'];
    subjNo(s) = size(N2cAll,1);
    disp('loaded')
end
%% CPPr 
close all
RTAllsorted =  RTAllss(stimAllss==1 & responsess==1);
CPPrAlltemps = CPPrAllss(stimAllss==1 & responsess==1,:);
[~,indx] = sort(RTAllsorted);
CPPrAllsorted = CPPrAlltemps(indx,:);
plot(tr,CPPrAllsorted(1,:));
size(CPPrAllsorted);
hold all
plot(tr,CPPrAllsorted(100,:));

% plot ERP image
d=1;cong=1;
figure
erpimage(CPPAllss(stimAllss==1 & responsess==1,:)', RTAllss(stimAllss==1 & responsess==1)',t,'tLockedCPP',20,0,'caxis',[0 20],'cbar');%100 moving average window
figure
erpimage(CPPrAllss(stimAllss==1 & responsess==1,:)', RTAllss(stimAllss==1 & responsess==1)',tr,'respLockedCPP',20,0,'caxis',[0 20],'cbar');%100 moving average window
%disp(['discrepancies between beta and actual = ', num2str(sum( mediationsS.Betatrial- mediationsS.Alltrial))])

%% N2c RT z-scored across subjects and within cong/incong
close all
RTAllsortedz =  RTAllssz(stimAllss==1 & responsess==1);
RTAllsorted =  RTAllss(stimAllss==1 & responsess==1);
N2cAlltemps = N2cAllss(stimAllss==1 & responsess==1,:);
CPPAlltemps = CPPAllss(stimAllss==1 & responsess==1,:);
subjAlltemps = subjectss(stimAllss==1 & responsess==1); 
RTAllsortedzSS =  RTAllssz(stimAllss==1 & responsess==1 & congAllss==1);
RTAllsortedSS =  RTAllss(stimAllss==1 & responsess==1 & congAllss==1);
N2cAlltempsSS = N2cAllss(stimAllss==1 & responsess==1 & congAllss==1,:);
CPPAlltempsSS = CPPAllss(stimAllss==1 & responsess==1 & congAllss==1,:,:);
subjAlltempsSS = subjectss(stimAllss==1 & responsess==1 & congAllss==1);
RTAllsortedzDS =  RTAllssz(stimAllss==1 & responsess==1 & congAllss==2);
RTAllsortedDS =  RTAllss(stimAllss==1 & responsess==1 & congAllss==2);
N2cAlltempsDS = N2cAllss(stimAllss==1 & responsess==1 & congAllss==2,:);
CPPAlltempsDS = CPPAllss(stimAllss==1 & responsess==1 & congAllss==2,:);
subjAlltempsDS = subjectss(stimAllss==1 & responsess==1 & congAllss==2);
minN2c = min(N2cAllss(:,451:601)')';
RTAllsortedzN2cRemove=  RTAllssz(stimAllss==1 & responsess==1 & minN2c<0);
RTAllsortedN2cRemove =  RTAllss(stimAllss==1 & responsess==1 & minN2c<0);
N2cAlltempsN2cRemove = N2cAllss(stimAllss==1 & responsess==1 & minN2c<0,:);
CPPAlltempsN2cRemove = CPPAllss(stimAllss==1 & responsess==1 & minN2c<0,:);
subjAlltempsN2cRemove = subjectss(stimAllss==1 & responsess==1 & minN2c<0);
%% N2c images
close all
Cols = [0 0 220; 0 128 255; 153 204 255; 255 128 0; 255 174 80; 255 215 163]/255;
alphaRT=[1 1 1 0]; yN2c = [-7 4]; ytickN2c = -7:2:4;ylim = yN2c;
clear hss3 hss4
colIndx=1;
hf=figure
N2c_c_grp = N2cAlltemps;
meanN2c = mean(N2c_c_grp,1);
stdN2c = 0*std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
rtN2c = mean(RTAllsorted);
hold on
hss3 = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[0, 0.7461, 1.0000],'LineWidth',3,'LineStyle','-'});
set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
    'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
hss4=line(rtN2c,rtN2c,...
    yN2c,'Color',[[0, 0.7461, 1.0000] 1],'LineWidth',1.5,'LineStyle','--');
ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
saveas(hf,['Figures\N2cmean.png']);

%% N2c
close  all
figure
[~,indx] = sort(RTAllsortedz);
N2cAllSorted = movmean(N2cAlltemps(indx,:),1);
RTAllsortedUsed = RTAllsorted(indx,:);
erpimage(N2cAllSorted', RTAllsortedUsed',t,'tLockedN2c',200,0,'caxis',[-5 5],'cbar','nosort');%100 moving average window

figure
[~,indx] = sort(RTAllsortedzSS);
N2cAllSSSorted = movmean(N2cAlltempsSS(indx,:),1);
RTAllsortedSSUsed = RTAllsortedSS(indx,:);
erpimage(N2cAllSSSorted', RTAllsortedSSUsed',t,'tLockedN2cSS',200,0,'caxis',[-5 5],'cbar','nosort');%100 moving average window

figure
[~,indx] = sort(RTAllsortedzDS);
N2cAllDSSorted = movmean(N2cAlltempsDS(indx,:),1);
RTAllsortedDSUsed = RTAllsortedDS(indx,:);
erpimage(N2cAllDSSorted', RTAllsortedDSUsed',t,'tLockedN2cDS',200,0,'caxis',[-5 5],'cbar','nosort');%100 moving average window

%figure
%erpimage(N2cAlltemps25', RTAllsorted25',t,'tLockedN2cDS',2,0,'caxis',[-10 10],'cbar');%100 moving average window
return
%% CPP
figure
[~,indx] = sort(RTAllsortedz);
CPPAllSorted = CPPAlltemps(indx,:);
RTAllsortedUsed = RTAllsorted(indx,:);
erpimage(CPPAllSorted', RTAllsortedUsed',t,'tLockedCPP',100,0,'caxis',[-5 15],'cbar','nosort');%100 moving average window

figure
[~,indx] = sort(RTAllsortedzSS);
CPPAllSSSorted = CPPAlltempsSS(indx,:);
RTAllsortedSSUsed = RTAllsortedSS(indx,:);
erpimage(CPPAllSSSorted', RTAllsortedSSUsed',t,'tLockedCPPSS',100,0,'caxis',[-5 15],'cbar','nosort');%100 moving average window

figure
[~,indx] = sort(RTAllsortedzDS);
CPPAllDSSorted = CPPAlltempsDS(indx,:);
RTAllsortedDSUsed = RTAllsortedDS(indx,:);
erpimage(CPPAllDSSorted', RTAllsortedDSUsed',t,'tLockedCPPDS',100,0,'caxis',[-5 15],'cbar','nosort');%100 moving average window

%% N2c removed
close  all
figure
[~,indx] = sort(RTAllsortedzN2cRemove);
N2cAllSorted = N2cAlltempsN2cRemove(indx,:);
RTAllsortedUsed = RTAllsortedN2cRemove(indx,:);
erpimage(N2cAllSorted', RTAllsortedUsed',t,'tLockedN2c',50,0,'caxis',[-5 5],'cbar','nosort');%100 moving average window
%% plotting the box plot of N2c and CPP for each participant
close all
figure
boxplot(RTAllsorted', subjAlltemps');
title('RT by participant')
xlabel('participant')
ylabel('RT');
figure
boxplot(RTAllsortedSS', subjAlltempsSS');
title('RT by participant Same Side')
xlabel('participant')
ylabel('RT (SS)');
figure
boxplot(RTAllsortedDS', subjAlltempsDS','PlotStyle', 'compact');
title('RT by participant Diff Side')
xlabel('participant')
ylabel('RT (DS)');
figure
hold on
boxplot(RTAllsortedSS', subjAlltempsSS');
boxplot(RTAllsortedDS', subjAlltempsDS','PlotStyle', 'compact');
title('RT by participant')
xlabel('participant')
ylabel('RT (DS and SS together)');
%% find mean of each participant
for s=1:max(subjAlltemps)
    meanRTSS(s) = mean(RTAllsortedSS(subjAlltempsSS ==s)) - mean(RTAllsortedDS(subjAlltempsDS ==s));
    ss(s)=s;
end
figure 
plot(1:43,meanRTSS,'r');
title('RTSS-RTDS by participant')
xlabel('participant')
ylabel('RTSS-RTDS');
%% N2i
close all
RTAllsorted =  RTAllss(stimAllss==1 & responsess==1);
N2iAlltemps = N2iAllss(stimAllss==1 & responsess==1,:);
[~,indx] = sort(RTAllsorted);
N2iAllSorted = N2iAlltemps(indx,:);
plot(t,N2iAllSorted(1,:));
size(N2iAllSorted);
hold all
plot(t,N2iAllSorted(100,:));
close all
% plot ERP image
d=1;cong=1;
figure
erpimage(N2iAllss(stimAllss==2,:)', RTAllss(stimAllss==2)',t,'tLockedCPP',100,0,'caxis',[-1 0],'cbar');%100 moving average window

%% plotting participant RTs boxes

return
%% plot linear regression
%%%%%%%%%%%%%%%%%%%1 is motion, 2 is colour%%%%%%%%%%%%%%%%%%%%%%
% 
% k=0;
% for limits=1
%     clear CPPrslopes RTs t meas b1 
%     k=k+1;
% CPPrslopes = zscore(mediationAll.N2cpeak(mediationAll.run>2  & mediationAll.stim==1 ...
%     & abs(mediationAll.N2cpeak)>(max(mediationAll.N2cpeak)/limits) & abs(mediationAll.CPPrslope)>(max(mediationAll.CPPrslope)/limits)));
% RTs = zscore(mediationAll.CPPrslope(mediationAll.run>2  & mediationAll.stim==1 ...
%     & abs(mediationAll.N2cpeak)>(max(mediationAll.N2cpeak)/limits) & abs(mediationAll.CPPrslope)>(max(mediationAll.CPPrslope)/limits)));
% subjs = (mediationAll.subj_idx(mediationAll.run>2  & mediationAll.stim==1 ...
%     & abs(mediationAll.N2cpeak)>(max(mediationAll.N2cpeak)/limits)  & abs(mediationAll.CPPrslope)>(max(mediationAll.CPPrslope)/limits)));
% CPPrslopesN = CPPrslopes(abs(CPPrslopes)<3 & abs(RTs)<3);
% RTsN = RTs(abs(CPPrslopes)<3 & abs(RTs)<3);
% subjN = subjs(abs(CPPrslopes)<3 & abs(RTs)<3);
% t = table(subjN, RTsN, CPPrslopesN, 'VariableNames',{'subj','dv','iv'});
% meas = table([1 2]','VariableNames',{'Measurements'});
% b1 = fitlme(t,'dv~iv+(1|subj)');
% pvaluesAll(k)=b1.Coefficients.pValue(2);
% limitsAll(k)=limits;
% end
% plot(RTsN, b1.Coefficients.Estimate(2)*RTsN)
% hold on 
% plot(RTsN, CPPrslopesN,'+');
% xlabel('RTs')
% ylabel('CPPrslopes');
% 
% limits = 200;
% % create the data
% for i=1:length(ms)-msoffset
%    mediationAllZs.(ms{i}) = mediationAll.(ms{i})(mediationAll.run>2  & mediationAll.stim==1 ...
%     & abs(mediationAll.N2cpeak)>(max(mediationAll.N2cpeak)/limits) & abs(mediationAll.CPPrslope)>(max(mediationAll.CPPrslope)/limits));
%     
% end
% mediationAllZs.subj_idx = (mediationAll.subj_idx(mediationAll.run>2  & mediationAll.stim==1 ...
%     & abs(mediationAll.N2cpeak)>(max(mediationAll.N2cpeak)/limits)  & abs(mediationAll.CPPrslope)>(max(mediationAll.CPPrslope)/limits)));
% mediationAll = mediationAllZs;
% %% plot linear regression
% %%%%%%%%%%%%%%%%%%%1 is motion, 2 is colour%%%%%%%%%%%%%%%%%%%%%%
clear CPPrslopes RTs
gs=3;
CPPrslopes = zscore(mediationAll.N2cpeak(mediationAll.run>=2  & mediationAll.stim==1));
RTs = zscore(mediationAll.CPPrslope(mediationAll.run>=2 & mediationAll.stim==1));
subjs = mediationAll.subj_idx(mediationAll.run>=2 & mediationAll.stim==1);
CPPrslopesN = CPPrslopes(abs(CPPrslopes)<gs & abs(RTs)<gs);
RTsN = RTs(abs(CPPrslopes)<gs & abs(RTs)<gs);
subjN = subjs(abs(CPPrslopes)<gs & abs(RTs)<gs);
t = table(subjN, RTsN, CPPrslopesN, 'VariableNames',{'subj','dv','iv'});
meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'dv~iv+(1|subj)')
plot(CPPrslopesN, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*CPPrslopesN)
hold on 
plot(CPPrslopesN, RTsN,'+');
xlabel('CPPrslopesN')
ylabel('RTsN');
%% checkwithins
%a = mediationAll.N2cLatency(mediationAll.stim==1 & mediationAll.bin==1);
%b = mediationAll.N2cLatency(mediationAll.stim==1 & mediationAll.bin==2);
%figure
%plot([a-b],'+');
%%
% combine the betas with the mediations
%for i = 4:length(msB)
    
%    mediationAll.(msB{i}) = mediationABeta.(msB{i});
%end

%for i=1:length(ms)-msoffset
%    mediationAll.(ms{i}) = mediationAllERP.(ms{i});
    
%end
%if singletrial
%    mediationAll.Alphapower = alphaPower;
%end

%mediationAll.subjectAll=mediationAllERP.subjectAll;
msA=fieldnames(mediationAll);

%% save the mediationAll to a large matrix
ms= fieldnames(mediationAll);
alldata=[];
kys = 1;
for i=1:length(ms)
    alldata=[alldata mediationAll.(ms{i})];
    for j=1:size(mediationAll.(ms{i}),2)
        if size(mediationAll.(ms{i}),2)>1
            keys{kys-1+j}= [ms{i} '_time_' num2str(j)];
        else
            keys{kys-1+j}= [ms{i}];
        end
    end
    kys = kys+size(mediationAll.(ms{i}),2);
end
%% omit the ERPs from mediationsG
%omitfields = {'CPPrslopeTraj','CPPslopeTraj','ERP','ERPr','erpPeakC','erpPeakI','erpLaterpeak','erprSlope'};
%mediationGs = rmfield(mediationAllG, omitfields);
%omitfields = {'Beta','Betar','BetaI','BetaIr','STFT','STFTr','ERP','ERPr','betarSlope'};
%mediationGBs =  rmfield(mediationAllGB, omitfields);
%% save the file to CSV
patht = '../../Stats/';
pathDDM = '/home/szhou/jo36/MRH110_EEG_2021/DDM_modelling/Data/';
%if singletrial
%    fields = {'CPPslopeTraj','CPPrslopeTraj','Beta','Betar','BetaI','STFT','STFTr','ERP','ERPr','BetaIr'};
%    mediationAll = rmfield(mediationAll,fields);
%end

%% obtain DDM parameters for use in modelling 
% rename the following
% stim=stimulus
% response= hit accuracy
% rt = RT
% subj_idx = subjectAll
mediationDDM = mediationAll;
%oldnames = {'stimulus','hit','RT','subjectAll'}; 
%newnames = {'stim','response','rt','subj_idx'};
%for jj=1:length(oldnames)
%mediationDDM.(newnames{jj}) = mediationDDM.(oldnames{jj});
%end
%mediationDDM = rmfield(mediationDDM,oldnames);
mediationDDM.stim = mediationDDM.stim;
mediationDDM.rt = mediationDDM.RT;
mediationDDM.rt = mediationDDM.rt/1000;
mediationAll.RT = mediationDDM.rt;
%%
if CSD
    struct2csv (mediationAll,[patht 'mediateCSD.csv']);
%    struct2csv (mediationAllA,[patht 'mediateCSDAcc.csv']);
%    struct2csv (mediationGs,[patht 'mediateGs.csv']);
%    struct2csv (mediationGBs,[patht 'mediateGBs.csv']);
    writetable(Tnames,[pathDDM 'subjNames.csv']);
    struct2csv (mediationDDM,[pathDDM 'DDMdataRT.csv']);
%    save(['slopeDataIndividCSD.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
else
    struct2csv (mediationAll,[patht 'mediatenoCSD.csv']);
%    struct2csv (mediationAllA,[patht 'mediatenoCSDAcc.csv']);
%    struct2csv (mediationGs,[patht 'mediatenoGs.csv']);
    writetable(Tnames,[pathDDM 'subjNames.csv']);
    struct2csv (mediationDDM,[pathDDM 'DDMdataRT.csv']);
%    save(['slopeDataIndividnoCSD.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
end



return