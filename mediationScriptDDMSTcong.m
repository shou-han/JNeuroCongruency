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
if singletrial; msoffset=5; msBoffset=8; else; msoffset=0; msBoffset=0; end; % to omit the ERPs from the single trial stats
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
for s=1:size(allsubj,2)
    subjCount=subjCount+1;
    sT = allsubj(s)
            subjIndx = subjCount;
    if CSD
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_CSD_319369' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','chanlocs','t','tr','subject_folder')
        load(['Data/ERPs/group_plots_erp_diff_CSD_DDM_Beta_Bin_Motion' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsB','STFT_time','STFT_timer')
    else
        load(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_319369' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','chanlocs','t','tr','subject_folder')
        %load(['Data/ERPs/group_plots_erp_DDM_bin_dual8Hz_ACC_' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','chanlocs','t','tr','subject_folder')
        load(['Data/ERPs/group_plots_erp_diff_CSD_DDM_Beta_Bin_Motion' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsB','STFT_time','STFT_timer')
    end
    clear ms
    %fields = {'ERP','ERPr','ERPdiff','ERPrdiff'};
    fields = {'ERPdiff','ERPrdiff'};
    %mediationsG = rmfield(mediationsG,fields);
    %mediationsAll = cat(3,mediationsAll,mediationsS);
    runs=[1 2; 3 4];
    if zscoreTrials
        for run=2
            [blah,indx]=find(ismember(mediationsS.run,runs(run,:)) & mediationsS.stim==1); %color or motion
            
            tempa = zscore(mediationsS.CPPrslope(indx));
            tempb= zscore(mediationsS.CPPonset(indx) - nanmean(mediationsS.CPPonset(indx)))/nanstd(mediationsS.CPPonset(indx)); 
            tempc = zscore(mediationsS.N2cpeak(indx));
            tempd = zscore(mediationsS.N2ipeak(indx));
            tempe = zscore(mediationsS.N2cLatency(indx));
            tempf = zscore(mediationsS.N2iLatency(indx));
            tempz = (mediationsS.trial(indx));
            tempBins = mediationsS.bin(indx);
            
%             [blah,indxB]=find(ismember(mediationsSB.run,runs(run,:))); %color or motion
%             tempaa = zscore(mediationsSB.betacrSlope(indxB));
%             tempbb = zscore(mediationsSB.betacAmplitude(indxB));
%             tempzz = (mediationsSB.trial(indxB));
%             tempBinsZ = mediationsSB.bin(indx);
            
            mediationsS.CPPrslope(indx) = tempa;
            mediationsS.CPPonset(indx) = tempb;
            mediationsS.N2cpeak(indx)=tempc;
            mediationsS.N2ipeak(indx)=tempd;
            mediationsS.N2cLatency(indx)=tempe;
            mediationsS.N2iLatency(indx)=tempf;
            mediationsS.Alltrial(indx) = tempz;
            
%             mediationsS.BetaSlope(indxB)=tempaa;
%             mediationsS.BetaAmplitude(indxB)=tempbb;
%             mediationsS.Betatrial(indxB) = tempzz;
            
            clear tempa  tempb tempc tempd tempe tempf tempaa tempbb tempz tempzz
        end
        clear i
    end
    for stim=1:2
        for cong=1:2
            for run=2
                [blah,indx]=find(mediationsS.stim==stim & ismember(mediationsS.run,runs(run,:))); %color or motion
                mediationsS.CPPrslopez(indx) = zscore(mediationsS.CPPrslope(indx));
                mediationsS.N2cpeakz(indx) = zscore(mediationsS.N2cpeak(indx));
                mediationsS.rtz(indx) = zscore(mediationsS.rt(indx));
                mediationsS.cppampz(indx) = zscore(mediationsS.CPPamp(indx));
            end
        end
    end
    for stim=1:2
        for run=2
            for resps =1:2
            [blah,indx]=find(mediationsS.stim==stim  & ismember(mediationsS.run,runs(run,:))); %color or motion
            mediationsS.CPPrslopezC(indx) = zscore(mediationsS.CPPrslope(indx));
            mediationsS.N2cpeakzC(indx) = zscore(mediationsS.N2cpeak(indx));
            mediationsS.rtzC(indx) = zscore(mediationsS.rt(indx));
            mediationsS.cppampzC(indx) = zscore(mediationsS.CPPamp(indx));
            end
        end
    end
    if zscoreTrials
        mediationsS.CPPrslope = zscore(mediationsS.CPPrslope);
         mediationsS.N2cpeak = zscore(mediationsS.N2cpeak);  
          mediationsS.N2cLatency = zscore(mediationsS.N2cLatency); 
    end
    if singletrial;    mediations = mediationsS; else; mediations=mediations;end;
        ms = fieldnames(mediations);
    if s==1
        for i=1:length(ms)
            mediationAll.(ms{i})=[];
                    subjects=[];
        end
    end
   % mediations.CPPrslope = zscore(mediations.CPPrslope);
    for i=1:length(ms)
        clear temp
        mediationAll.(ms{i}) = cat(1, mediationAll.(ms{i}),mediations.(ms{i})');
    end
    for i=1
        temp  = squeeze(mediations.(ms{i}));
        subjects = [subjects allsubj(subjIndx)*ones([1 size(squeeze(temp),2)])];
    end
    tempTable.subjects = subject_folder;
    tempTable.subjNo = allsubj(s);
    Tnames = [Tnames;tempTable];
end
mediationAll.subj_idx=subjects';

disp('loaded')
%disp(['discrepancies between beta and actual = ', num2str(sum( mediationsS.Betatrial- mediationsS.Alltrial))])
%% plot linear mixed effect model 
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

%% congruency data 
clc; close all;
gs=100;
clear CPPrslopes RTs N2cpeaks t
[indx, ~] = find(mediationAll.run>=2  & mediationAll.stim==1 );
N2cpeaks = (mediationAll.N2cpeakz(indx));
CPPrslopes = (mediationAll.CPPrslopez(indx));
CPPamps = (mediationAll.cppampz(indx));
N2cpeaksC = (mediationAll.N2cpeakzC(indx));
CPPrslopesC = (mediationAll.CPPrslopezC(indx));
CPPampsC = (mediationAll.cppampzC(indx));
RTsC = (mediationAll.rt(indx));
RTs = (mediationAll.rt(indx));
subjs = mediationAll.subj_idx(indx);
congs = mediationAll.cong(indx);

N2cpeaks2 = N2cpeaks(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs & abs(zscore(N2cpeaks))<gs);
CPPrslopes2 = CPPrslopes(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs  & abs(zscore(N2cpeaks))<gs);
subj2 = subjs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs  & abs(zscore(N2cpeaks))<gs);
cong2 = congs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs  & abs(zscore(N2cpeaks))<gs);
rt2 = RTs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs  & abs(zscore(N2cpeaks))<gs);


% removed outliers
N2cpeaks = N2cpeaks2; CPPrslopes = CPPrslopes2; subjs = subj2; congs= cong2;RTs = rt2;
clear CPPrslopes2 RTs2 subjs2 cong2 rt2 CPPamps2
t = table(subjs, RTs, CPPrslopes,N2cpeaks, congs,'VariableNames',{'subj','rt','cppslope','n2cpeak', 'cong'});
%% mediation
addpath(genpath('MediationToolbox-master'))
[paths, toplevelstats, firstlevelstats] = mediation(N2cpeaks,RTs,CPPrslopes,'stats');
[pathsCongCPP, toplevelstatsCongCPP, firstlevelstatsCongCPP] = mediation(congs,RTs,CPPrslopes,'stats');
[pathsCongN2c, toplevelstatsCongN2c, firstlevelstatsCongN2c] = mediation(congs,RTs,N2cpeaks,'stats');
toplevelstats
toplevelstatsCongCPP
toplevelstatsCongN2c
return
%% n2cpeak cppslope rt mediation
clc; close all
%regressions
meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'cppslope~n2cpeak')
plot(N2cpeaks, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*N2cpeaks)
hold on
plot(N2cpeaks, CPPrslopes,'+');
xlabel('N2cpeak')
ylabel('CPP slopes');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope')
plot(CPPrslopes, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*CPPrslopes)
hold on
plot(CPPrslopes, RTs,'+');
xlabel('CPPslope')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~n2cpeak')
plot(N2cpeaks, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*N2cpeaks)
hold on
plot(N2cpeaks, RTs,'+');
xlabel('N2cpeak')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope+n2cpeak')


%%
meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'cppslope~n2cpeak+((1+cong)|subj)')
plot(N2cpeaks, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*N2cpeaks)
hold on
plot(N2cpeaks, CPPrslopes,'+');
xlabel('N2cpeak')
ylabel('CPP slopes');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope+(')
plot(CPPrslopes, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*CPPrslopes)
hold on
plot(CPPrslopes, RTs,'+');
xlabel('CPPslope')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~n2cpeak')
plot(N2cpeaks, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*N2cpeaks)
hold on
plot(N2cpeaks, RTs,'+');
xlabel('N2cpeak')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope+n2cpeak')
%% cong CPPrslope RT
meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'n2cpeak~1+cppamps+cong+(1+cong|subj)')

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cong+cppamps+(1|subj)')

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~n2cpeak+cong+cppamps+(1+cong|subj)')

%%

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope+(1+cong|subj)')
plot(CPPrslopes, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*CPPrslopes)
hold on
plot(CPPrslopes, RTs,'+');
xlabel('CPPslope')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~cppslope+cong+(1+cong|subj)')

%% cong N2cpeak RT
meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'n2cpeak~1+(1+cong|subj)')
%%


meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~n2cpeak+(1+cong|subj)')
plot(N2cpeaks, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*N2cpeaks)
hold on
plot(N2cpeaks, RTs,'+');
xlabel('N2cpeak')
ylabel('RT');

meas = table([1 2]','VariableNames',{'Measurements'});
b1 = fitlme(t,'rt~n2cpeak+cong+(1+cong|subj)')

%% congruency two
% %%%%%%%%%%%%%%%%%%%1 is motion, 2 is colour%%%%%%%%%%%%%%%%%%%%%%
clc; close all;
gs=3;
clear CPPrslopes RTs t b1
gs=3;
for cong=1:2
    CPPrslopes = (mediationAll.N2cpeak(mediationAll.run>=2  & mediationAll.cong == cong & mediationAll.stim==1  & mediationAll.response ==1));
    RTs = (mediationAll.CPPrslope(mediationAll.run>=2 & mediationAll.cong == cong & mediationAll.stim==1  & mediationAll.response ==1));
    subjs = mediationAll.subj_idx(mediationAll.run>=2 & mediationAll.cong == cong & mediationAll.stim==1 & mediationAll.response ==1);
    congs = mediationAll.cong(mediationAll.run>=2 & mediationAll.cong == cong & mediationAll.stim==1 & mediationAll.response ==1);
    CPPrslopes2 = CPPrslopes(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs);
    RTs2 = RTs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs);
    subj2 = subjs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs);
    cong2 = congs(abs(zscore(CPPrslopes))<gs & abs(zscore(RTs))<gs);
    % removed outliers
    CPPrslopes = CPPrslopes2; RTs = RTs2; subjs = subj2; congs= cong2;
    clear CPPrslopes2 RTs2 subjs2 cong2
    t = table(subjs, RTs, CPPrslopes,congs, 'VariableNames',{'subj','dv','iv', 'cong'});
    meas = table([1 2]','VariableNames',{'Measurements'});
    b1 = fitlme(t,'dv~iv+(1|subj)')
    plot(CPPrslopes, b1.Coefficients.Estimate(1)+b1.Coefficients.Estimate(2)*CPPrslopes)
    hold on
    plot(CPPrslopes, RTs,'+');
    xlabel('CPPrslopes')
    ylabel('RTs');
end
%% checkwithins


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
    writetable(t,[pathDDM,'statsDataCSD.csv']);
    struct2csv (mediationDDM,[pathDDM 'DDMdataRTCSDbri.csv']);
%    save(['slopeDataIndividCSD.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
else
    struct2csv (mediationAll,[patht 'mediatenoCSD.csv']);
%    struct2csv (mediationAllA,[patht 'mediatenoCSDAcc.csv']);
%    struct2csv (mediationGs,[patht 'mediatenoGs.csv']);
    writetable(Tnames,[pathDDM 'subjNames.csv']);
    writetable(t,[pathDDM,'statsData.csv']);
    struct2csv (mediationDDM,[pathDDM 'DDMdataRT4551200400.csv']);
%    save(['slopeDataIndividnoCSD.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
end



return