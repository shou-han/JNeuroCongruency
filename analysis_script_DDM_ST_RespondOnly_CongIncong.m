% save all the plots in a group and the plot can choose whichever one it wants
% participant 9, choose a different set of electrodes

% don't mind about up/down since there is little behaviour differences,
%suggesting that they have the same level of processing
%count for now; will look at count after respond
% is taken care of.......

% allselect is used to select the trials
%%%%%%%%%%%%%%%%%%%1 is motion, 2 is colour%%%%%%%%%%%%%%%%%%%%%%

% colour same/different side as motion
% colour before, same and after motion
%
function analysis_script_DDM_ST_RespondOnly_CongIncong(single_participants)
%single_participants=5
indexbluedown = {' '};
%check 12 and 21 '017','204',
middlebluedown = {'017','204','011','008','015','007','004','010','022','020','018','007','202','206','208','209','210','211','213','214','216','218','220','222','223','224','225','226','228','229','230'};
disp('start')
oldnew = {'', 'old'};
old=1 ;
chP=53; %naming only
ch_CPP=53;
ch_lr = 15; %45;%ch_N2c(single_participants);% left is all contra using erpN
ch_rl = 20; %51;%ch_N2i(single_participants);%right is all ipsi using erpN
% these two (new ones) didn't have the cap aligned so the CPP was a bit too
% high
consWin = 15;
%% Use Current Source Density transformed erp? 1=yes, 0=no
CSD=1;

addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath('function_programs/');
eeglab

load(['../beg_vals/beg_vals' num2str(single_participants) '.mat'])

chanlocs = readlocs('actiCAP64_2_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
%chanlocs = readlocs('cap64.loc');
TCD_bigdots = {};
path_temp = '../Data/';
skip_step = 10;
%%
side_tags = {'Left','Right'};
%%
duds = [];
%%
if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
end

%% Define channels, having combined Brain Products and Biosemi data
plot_chans = 1:64;
left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];
elec_pairs = [1,32;33,63;34,62;4,31;3,30;37,60;36,61;5,...
    27;38,59;6,28;39,58;7,29;9,26;41,56;8,25;40,57;10,21;...
    42,55;11,22;43,54;12,23;15,20;45,51;14,19;44,52;46,50;...
    47,49;16,18];

tester = zeros(64,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','labels','plotchans',1:64);

tester = zeros(64,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','numbers','plotchans',1:64);

%% Triggers
% ITI,left/right
numch=64;
rtlim=[0.300 1.500];

%%%%%%%%%%%check the channels%%%%%%%%%%%%%%
ch_beta_ci = [8; 25];
ch_contra = [right_hemi;left_hemi];
ch_ipsi = [left_hemi;right_hemi];
BL_erp = [-100,0];
BL_beta = [-100 0];
BL_alpha = [-100 0];
ch_alpha =  [46 47 48 49 50];
% zscore threshold
z_thresh = 3;

%% response time bins parameter and search parameters
no_of_bins = 3;
cs=[1;2];
dPA = [1 ;2];
% plotting parameters
for bin = 1:no_of_bins
    rt_bins_tags{bin} = num2str(bin);
end

%% Groupings
targcodes(1,:,:) =   [101:2:116;102:2:116]'; %color targcodes
targcodes(2,:,:) =   [ 101:108;109:116]';%motion targcodes
%Trial types:
% M Motion, U/D Up/Down, L/R Left/Right, C Colour, R/B Red/Blue, L/R
% Left/Right
% 101 is MULCRL
% 102 is MULCRR
% 103 is MDLCRL
% 104 is MDLCRR
% 105 is MULCBL
% 106 is MULCBR
% 107 is MDLCBL
% 108 is MDLCBR
% 109 is MURCRL
% 110 is MURCRR
% 111 is MDRCRL
% 112 is MDRCRR
% 113 is MURCBL
% 114 is MURCBR
% 115 is MDRCBL
% 116 is MDRCBR
% 117 is REST
targMotionColour(1,:) = [101 103 105 107 110 112 114 116]; %same side
targMotionColour(2,:) = [102 104 106 108 109 111 113 115]; %different sides
%% Start loop
for s=1:length(single_participants)
    sT = allsubj{s};
    %% define the time epochs for beta alpha and ERPs
    fs = 500; %resampled to 500 in the runafew
    
    ts = -0.700*fs:1.500*fs;
    t = ts*1000/fs;
    
    % resp-locked erps
    trs = [-0.5*fs:fs*.05];
    tr = trs*1000/fs;
    
    % generating SSVEP (same frequency)
    %% Load the trials
    clear indx1 indx2
    
    load(['../' path_temp subject_folder{s} '/' allsubj{s} '_alldata'],'erp_LPF_8Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz','erp_LPF_35Hz_CSD',...
        'erp_LPF_8Hz_TI','erp_LPF_8Hz_TI_CSD','erp_LPF_35Hz_TI','erp_LPF_35Hz_TI_CSD',...
        'allRT','allRT_TI','allrespLR','allTrig','allblock_count',...
        'resp_artrej','ET_resp_artrej','all_speed','allCMblock_count','allrun_count');
    allblock_count=(allblock_count-1)/2+1;
    if CSD
        erp_all=double(erp_LPF_35Hz_CSD);
        erp_TI = double(erp_LPF_35Hz_TI_CSD);
    else
        erp_all=double(erp_LPF_35Hz);
        erp_TI = double(erp_LPF_35Hz_TI);
    end
    trig_all = allTrig;
    
    ET_BL_resp_artrej_all = ET_resp_artrej;
    BL_resp_artrej_all = resp_artrej;
    
    allRT_all=allRT;
    % Baseline erp
    
    % find accuracy all the trials which are correct and incorrect
    if ismember(allsubj{s},middlebluedown)
        bluedownmiddle=1;
    else
        bluedownmiddle=0;
    end
    clear middlebluedown
    if bluedownmiddle
        indxC1(1,:) = find(ismember(allTrig,[101:104,109:112]));% Blue
        indxC2(1,:) = find(ismember(allTrig,[105:108 113:116]));%Red
        indxM1(1,:) = find(ismember(allTrig,[103 104 107 108 111 112 115 116]));%Down
        indxM2(1,:) = find(ismember(allTrig,[101 102 105 106 109 110 113 114]));%UP
    else
        indxC2(1,:) = find(ismember(allTrig,[101:104,109:112]));%Red
        indxC1(1,:) = find(ismember(allTrig,[105:108 113:116]));%Blue
        indxM2(1,:) = find(ismember(allTrig,[103 104 107 108 111 112 115 116]));%Up
        indxM1(1,:) = find(ismember(allTrig,[101 102 105 106 109 110 113 114]));%Down
    end
    
    allacc{2} = allrespLR;
    allacc{2}(indxC1(allrespLR(indxC1)==1 |allrespLR(indxC1)==3))=1;
    allacc{2}(indxC1(allrespLR(indxC1)==2 |allrespLR(indxC1)==4))=0;
    allacc{2}(indxC1(allrespLR(indxC1)==5))=nan;
    
    allacc{2}(indxC2(allrespLR(indxC2)==2 | allrespLR(indxC2)==4))=1;
    allacc{2}(indxC2(allrespLR(indxC2)==1 | allrespLR(indxC2)==3))=0;
    allacc{2}(indxC2(allrespLR(indxC2)==5))=nan;
    
    allacc{1} = allrespLR;
    allacc{1}(indxM1(allrespLR(indxM1)==1 | allrespLR(indxM1)==3))=1;
    allacc{1}(indxM1(allrespLR(indxM1)==2 | allrespLR(indxM1)==4))=0;
    allacc{1}(indxM1(allrespLR(indxM1)==5))=nan;
    
    allacc{1}(indxM2(allrespLR(indxM2)==2 | allrespLR(indxM2)==4))=1;
    allacc{1}(indxM2(allrespLR(indxM2)==1 | allrespLR(indxM2)==3))=0;
    allacc{1}(indxM2(allrespLR(indxM2)==5))=nan;
    
    
    %% calculate erps
    %make erps contra ipsi
    for trial=1:size(erp_all,3)
        if allCMblock_count(trial)==1 %motion trial
            if ismember(allTrig(trial),110:116)
                side=2;%right
            else
                side=1;%left
            end
        elseif allCMblock_count(trial)==2 %color trial
            if ismember(allTrig(trial),102:2:116)
                side=2;
            else
                side=1;
            end
        end
        
        erpN(left_hemi,:,trial) = erp_all(ch_contra(side,:),:,trial); %save as new variable since transformation problems
        erpN(right_hemi,:,trial) = erp_all(ch_ipsi(side,:),:,trial);
        erpN(centre_chans,:,trial) = erp_all(centre_chans,:,trial);
        
        erpN_TI(left_hemi,:,trial) = erp_TI(ch_contra(side,:),:,trial); %save as new variable since transformation problems
        erpN_TI(right_hemi,:,trial) = erp_TI(ch_ipsi(side,:),:,trial);
        erpN_TI(centre_chans,:,trial) = erp_TI(centre_chans,:,trial);
        
    end
    clear erp baseline_erp
    erp = erpN;
    erp_TI = erpN_TI;
    disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
    
    
    % Baseline erp
    baseline_erp = mean(erp(:,t>=BL_erp(1) & t<=BL_erp(2),:),2);
    erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp
    
    %Baseline erp TI
    baseline_erpTI = mean(erp_TI(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp_TI = erp_TI-repmat(baseline_erpTI,[1,size(erp_TI,2),1]); % baseline full erp
    
    
    %% make response locked
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    validrlock = zeros(1,length(allRT)); % length of RTs.
    % for each trial, shift the erp signature to make RT = 0;
    for n=1:length(allRT)
        [blah,RTsamp] = min(abs(t*fs/1000-allRT(n)));
        if      RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t)
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            validrlock(n)=1;
        end
    end
    for n=1:length(allRT_TI)
        [blah,RTsamp] = min(abs(t*fs/1000-allRT_TI(n)));
        if      RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t)
            erpr_TI(:,:,n) = erp_TI(:,RTsamp+trs,n);
            validrlockTI(n)=1;
        end
    end
    %% beta times
    STFT_time=[];
    no_of_cycles = 8;
    % % get a broad range, e.g. beta, 20 to 35Hz
    fs_STFT = [20:35]; % stftlen_STFT = 140; % 8 cycles = 280 ms (at mid frequency, i.e. 28 Hz for beta), = 140 samples.
    stftlen_STFT = round((1000/round(median(fs_STFT))*no_of_cycles)/2);
    
    % % get a specific frequency for SSVEP
    %fs_STFT = [25]; % or if you want a particular SSVEP frequency
    %stftlen_STFT = round((1000/fs_STFT*no_of_cycles)/2);
    % for SSVEP frequency make sure it's EXACTLY a particular number of cycles of the frequency.
    % check freq_temp_STFT to make sure SSVEP frequency falls on the range
    
    skip_step = 10;
    cc=1;
    for tt = 1:skip_step:length(ts)-(stftlen_STFT)
        tf = tt:tt+stftlen_STFT-1;
        nfft = length(tf);
        freq_temp_STFT = (0:ceil((nfft+1)/2)-1)*fs/nfft;
        STFT_time(cc) = mean(t(tf));
        cc=cc+1;
    end
    %Response locked STFT time (ms)
    STFT_timer= -500:skip_step*2:100; %anything beyond 500 means that it's not on ts (smallest is -534 before 0)
    %Response locked STFT time in samples
    STFT_timers = -0.5/(skip_step/fs):.100/(skip_step/fs);
    validrlockSB = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT)
        [blah,RTsamp] = min(abs(STFT_time*fs/1000-allRT(n))); % get the sample point of the RT.
        if RTsamp+STFT_timers(1) >0 & RTsamp+STFT_timers(end)<=length(STFT_time) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than last RT point.
            validrlockSB(n)=1;
        end
    end
    
    % remove erps whose N2c is more than 0
%     for trial=1:size(erp,3)
%         erpN2c = squeeze(movmean(erp(ch_lr,:,trial),1));
%         %%%selection:
%         N2cselect = min(erpN2c(t>296 & t<396));
%         if N2cselect<0
%             allselect(trial)=1;
%         else
%             allselect(trial)=1;
%         end
%     end
%     
%     clear N2select    
    %% calculate the conditions
    %bin it
    
    for d = 1:size(cs,1) %1 is attend motion, 2 is attend colour
        ET_BL_resp_artrej = ET_BL_resp_artrej_all;
        BL_resp_artrej = BL_resp_artrej_all;
        allRT = allRT_all;
        
        
        c_ind =  find(ismember(allTrig,squeeze(targcodes(d,:,:))) & allCMblock_count==d );
        disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(c_ind))])
        
        allstuff.RTs{s,d} = allRT(c_ind)*1000/fs;
        
        %prepare for speed bins
        bins_count = zeros([1 no_of_bins]);
        for cong = 1:2 % same or different sides, 1 is same side, 2 is different sides
            
            % calcs the indices of the triggers for counting
            conds{s,1,d,cong,1}  = find(ismember(allTrig,targcodes(d,:,:)) & allCMblock_count==d &...
                ET_BL_resp_artrej & allrun_count<=2 & ismember(allTrig,targMotionColour(cong,:))); %
            conds{s,1,d,cong,2} = conds{s,1,d,cong,1};  %
            
            % calcs the indices of the triggers for response
            condTemp = find(ismember(allTrig,targcodes(d,:,:)) & allCMblock_count==d & ismember(allTrig,targMotionColour(cong,:)) & ...
                allRT>rtlim(1)*fs & allRT<rtlim(2)*fs & validrlock & validrlockSB & ET_BL_resp_artrej & allrun_count>2); % check baseline for subject 3

            %no_of_bins bins
            [RT_temp,indx] = sort(allRT(condTemp));
            gTemp = condTemp(indx);
            rtindx = floor(length(RT_temp)/no_of_bins);
            for bin=1:no_of_bins-1
                binstart= (bin-1)*rtindx+1; %start bin
                binend = (bin-1)*rtindx+rtindx; %end bin
                conds{s,2,d,cong,bin} = gTemp(binstart:binend);
                RTs{s,2,d,cong,bin} = allRT([conds{s,2,d,cong,bin}])*1000/fs;
            end
            conds{s,2,d,cong,no_of_bins} = gTemp(((no_of_bins-1)*rtindx+1):end);
            RTs{s,2,d,cong,no_of_bins} = allRT([conds{s,2,d,cong,no_of_bins}])*1000/fs;
            
            % find conditionsT
            condsT{s,1,d,cong}  = find(ismember(allTrig,targcodes(d,:,:)) & ismember(allTrig,targMotionColour(cong,:)) & allCMblock_count==d & allrun_count<=2 );
            condsT{s,2,d,cong}  = find(ismember(allTrig,targcodes(d,:,:)) & ismember(allTrig,targMotionColour(cong,:)) & allCMblock_count==d & allrun_count>2);
        end
    end
    disp(['Subject ',num2str(s) ': ' ,allsubj{s} ' Condition Valid Trials: ',num2str(length([conds{s,1,:,:,1}])), ...
        ' = ',num2str(round(100*(length([conds{s,1,:,:,1}]))/length([condsT{s,1,:,:}]))),'%'])
    disp(['Subject ',num2str(s) ': ' ,allsubj{s} ' Condition Valid Trials: ',num2str(length([conds{s,2,:,:,:}])), ...
        ' = ',num2str(round(100*(length([conds{s,2,:,:,:}]))/length([condsT{s,2,:,:}]))),'%'])
    clear beta_asym_bins RT_bins RT_bins_all
    %% define search windows
    clear RT_bins RT_bins_all
    allBins.accuracy(1) =round(100*(length([conds{s,1,:,:,:}]))/length([conds{s,1,:,:,:}]));
    % these two (new ones) didn't have the cap aligned so the CPP was a bit too,1,:,:,:}]));
    allBins.accuracy(2) =round(100*(length([conds{s,2,:,:,:}]))/length([conds{s,2,:,:,:}]));
    
    %perform group trials mediation analysis
    j=0;
    %hemiS = [elec_pairs(:,1); centre_chans'; elec_pairs(:,2)];
    %     if CSD
    %         window_searchi = [200 500];
    %         window_searchc = [200 500];
    %         window_searchpc = [200 500];
    %     else
    %         window_searchi = [100 300];
    %         window_searchc = [100 300];
    %         window_searchpc = [100 350];
    %     end
    %% create the data for data analysis and plotting
    for c=2 %respond only
        for d=1:2 %motion and colour
            for cong=1:2
                if d==1 %motion
                    window_searchi = [300 500];
                    window_searchc = [128 560]; %JC: 324-424, S: 128-560
                    window_searchpc = [200 500];
                    window_searchCPP = [-200 -5];
                elseif d==2 %colour
                    window_searchi = [300 500];
                    window_searchc = [128 560];
                    window_searchpc = [150 300];
                    window_searchCPP = [-200 -5];
                end
                for bin=1:no_of_bins
                    j=j+1;
                    kk=1;
                    % ERPS
                    clear coef CPP_smooth CPPr_smooth pretemp N2c_smooth N2i_smooth CPPslope_mean_temp
                    % determined from all the N2cs and N2is
                    pretemp = find(t==0)-1;
                    CPP_mean_temp=[];N2c_mean_temp=[];N2i_mean_temp=[];CPPr_mean_temp=[];RT_temp=[];erp_mean_temp=[];erpr_mean_temp=[];
                    N2pc_mean_temp=[];CPPslope_mean_temp=[];
                    CPP_TI_mean_temp=[];CPPr_TI_mean_temp=[];N2c_TI_mean_temp=[];N2i_TI_mean_temp=[];N2pc_TI_mean_temp=[];
                    if isempty([conds{s,c,:,:}]),break,end
                    CPP_mean_temp = squeeze(mean(erp(ch_CPP,:,[conds{s,c,d,cong,bin}]),1));
                    CPPr_mean_temp = squeeze(mean(erpr(ch_CPP,:,[conds{s,c,d,cong,bin}]),1));
                    if size(CPP_mean_temp,1)==1, CPP_mean_temp=CPP_mean_temp',end
                    
                    CPP_TI_mean_temp = squeeze(mean(erp_TI(ch_CPP,:,[conds{s,c,d,cong,bin}]),1));
                    CPPr_TI_mean_temp = squeeze(mean(erpr_TI(ch_CPP,:,[conds{s,c,d,cong,bin}]),1));
                    
                    N2c_mean_temp = squeeze(erp(ch_lr,:,[conds{s,c,d,cong,bin}]));
                    N2i_mean_temp = squeeze(erp(ch_rl,:,[conds{s,c,d,cong,bin}]));
                    N2pc_mean_temp= squeeze(erp(ch_lr,:,[conds{s,c,d,cong,bin}])-erp(ch_rl,:,[conds{s,c,d,cong,bin}]));
                    
                    N2c_TI_mean_temp = squeeze(erp_TI(ch_lr,:,[conds{s,c,d,cong,bin}]));
                    N2i_TI_mean_temp = squeeze(erp_TI(ch_rl,:,[conds{s,c,d,cong,bin}]));
                    N2pc_TI_mean_temp= squeeze(erp_TI(ch_lr,:,[conds{s,c,d,cong,bin}])-erp_TI(ch_rl,:,[conds{s,c,d,cong,bin}]));
                    
                    RT_temp = allRT([conds{s,c,d,cong,bin}]);
                    erp_mean_temp = squeeze(erp(:,:,[conds{s,c,d,cong,bin}]));
                    erpr_mean_temp = squeeze(erpr(:,:,[conds{s,c,d,cong,bin}]));
                    erp_TI_mean_temp = squeeze(erp_TI(:,:,[conds{s,c,d,cong,bin}]));
                    erpr_TI_mean_temp = squeeze(erpr_TI(:,:,[conds{s,c,d,cong,bin}]));
                    % store the graphs for plotting
                    allBins.subj{s,c,d,cong,bin} = s;
                    allBins.countRespond{s,c,d,cong,bin}=c;
                    allBins.colourMotion{s,c,d,cong,bin} = d;
                    allBins.fastSlow{s,c,d,cong,bin} = bin;
                    allBins.sameDiff{s,c,d,cong,bin} = cong;
                    allBins.Rts{s,c,d,cong,bin} = allRT([conds{s,c,d,cong,bin}])*1000/fs;
                    allBins.acc{s,c,d,cong,bin} = allacc{d}([conds{s,c,d,cong,bin}]);
                    if isempty(erp_mean_temp),break,end
                    allBins.ERP_side(s,c,d,cong,bin,:,:) = squeeze(mean(erp_mean_temp,3));
                    allBins.ERPr_side(s,c,d,cong,bin,:,:) = squeeze(mean(erpr_mean_temp,3));
                    
                    allBins.CPP_side(s,c,d,cong,bin,:) = squeeze(mean(CPP_mean_temp,2));
                    allBins.CPPr_side(s,c,d,cong,bin,:) = squeeze(mean(CPPr_mean_temp,2));
                    
                    allBins.N2c_side(s,c,d,cong,bin,:) = squeeze(mean(N2c_mean_temp,2));
                    allBins.N2i_side(s,c,d,cong,bin,:) = squeeze(mean(N2i_mean_temp,2));
                    allBins.N2pc_side(s,c,d,cong,bin,:) = squeeze(mean(N2pc_mean_temp,2));
                    
                    %Task irrelevant
                    allBins.ERPTI_side(s,c,d,cong,bin,:,:) = squeeze(mean(erp_TI_mean_temp,3));
                    allBins.ERPTIr_side(s,c,d,cong,bin,:,:) = squeeze(mean(erpr_TI_mean_temp,3));
                    
                    allBins.CPPTI_side(s,c,d,cong,bin,:) = squeeze(mean(CPP_TI_mean_temp,2));
                    allBins.CPPTIr_side(s,c,d,cong,bin,:) = squeeze(mean(CPPr_TI_mean_temp,2));
                    
                    allBins.N2cTI_side(s,c,d,cong,bin,:) = squeeze(mean(N2c_TI_mean_temp,2));
                    allBins.N2iTI_side(s,c,d,cong,bin,:) = squeeze(mean(N2i_TI_mean_temp,2));
                    allBins.N2pcTI_side(s,c,d,cong,bin,:) = squeeze(mean(N2pc_TI_mean_temp,2));
                    
                    %
                    RT_side{s,c,d,cong,bin} = allRT([conds{s,c,d,cong,bin}])*1000/fs;
                    RT_st_data{c,d,cong,bin} = RT_side{s,c,d,cong,bin};
                    
                    % statistics on count/respond and motion/colour (1 and 2)
                    mediations.run(j)=c;
                    mediations.stim(j)=d;
                    mediations.RT(j) = mean(RT_temp);
                    mediations.bin(j)=bin;
                    mediations.cong(j)=cong;
                    if size(N2c_mean_temp,1)==1, N2c_mean_temp=N2c_mean_temp', N2i_mean_temp=N2i_mean_temp',N2pc_mean_temp=N2pc_mean_temp',end
                    [mediations.N2c_onset(j)] = obtainOnset_final_subj(-N2c_mean_temp,...
                        window_searchc,2,t,0); %the onset will have the same fields as ERPs
                    [mediations.N2i_onset(j)] = obtainOnset_final_subj(-N2i_mean_temp,...
                        window_searchi,5,t,0); %the onset will have the same fields as ERPs
                    [mediations.N2cpeak(j), mediations.N2cLatency(j)] = obtainPeak(-N2c_mean_temp,...
                        window_searchc,5,t,0); %N2c peak but positive
                    [mediations.N2ipeak(j), mediations.N2iLatency(j)] = obtainPeak(-N2i_mean_temp,...
                        window_searchi,5,t,0); %N2i peak but positive
                    [mediations.N2pcpeak(j), mediations.N2pcLatency(j)] = obtainPeak(-N2pc_mean_temp,...
                        window_searchpc,5,t,0); %N2c peak but positive
                    [mediations.CPPonset(j)] = obtainOnset_final_subj(CPP_mean_temp,...
                        [0,500],15,t,0); %CPP onset
                    [mediations.CPPpeakl(j), CPPPeakTime(j)] = obtainPeak(squeeze(CPP_mean_temp),....
                        [0,1000],5,t,0); %CPP peak
                    if c==2 %respond only
                        CPPsides = allBins.CPPr_side(s,c,d,cong,bin,:);
                        CPPsides = CPPsides(:);
                        [mediations.CPPlevel(j), CPPrTime(j)] = obtainPeak(squeeze(CPPsides),...
                            [-500,100],5,tr,1); %CPPr peak
                        %                     [mediations.CPPrslope(j), mediations.CPPronset(j)] = obtainOnset_slope(squeeze(CPPr_mean_temp), ...
                        %                         [tr(slope_timeframe_index(1)) tr(slope_timeframe_index(2))], tr);
                        [mediations.CPPrslope(j), mediations.CPPronset(j)] = obtainOnset_slope(squeeze(CPPsides), ...
                            window_searchCPP, tr);
                        %timeframesIndices(d,:) = slope_timeframe_index;
                    else
                        mediations.CPPlevel(j)=NaN;mediations.CPPrslope(j)=NaN;mediations.CPPronset(j)=NaN;
                    end
                end
            end
        end
    end
end
%% plot the graphs for the N2c peak
clear h h1
h1 = figure;
c=2;
for d=1
    for bin=1:no_of_bins
        for cong=1:2
            indx = (d-1)*2+bin;
            h(d) = plot(t,squeeze(allBins.N2c_side(s,c,d,cong,bin,:)),'LineWidth',3,'LineStyle','-');hold on
            st = find(t>mediations.N2cLatency(indx),1)-20;
            ed = find(t>mediations.N2cLatency(indx),1)+20;
            r = -mediations.N2cpeak(indx)*ones([1,ed-st+1]); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(t(st:ed), r,'Linewidth',5, 'LineStyle', ':');hold on;
        end
    end
end

set(gca,'FontSize',16,'xlim',[-100,1000]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title([subject_folder{s}, ' N2c'])
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Motion','Colour'},'FontSize',16,'Location','NorthWest');
pause(1)
saveas(h1,['Figures/N2cIndividual/N2c', subject_folder{s} '.jpg'])
%% plot the graphs for each of the slopes
%%%Plot each individual participant's CPPr_slope with time-window varying
%%%per participant
clear h h1
h1 = figure;
c=2;
for d=1
    for bin=1:no_of_bins
        for cong=1:2
            h(d) = plot(tr,squeeze(allBins.CPPr_side(s,c,d,cong,bin,:)),'LineWidth',3,'LineStyle','-');hold on
            st = find(tr>-350,1);ed =find(tr>-10,1);
            coef = polyfit(tr(st:ed),...
                squeeze(allBins.CPPr_side(s,c,d,cong,bin,st:ed))',1);% coef gives 2 coefficients fitting r = slope * x + intercept
            CPP_slope(bin)=coef(1)
            r = coef(1) .* tr(st:ed) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(tr(st:ed), r,'Linewidth',5, 'LineStyle', ':');hold on;
        end
    end
end

set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
xlabel('Time (ms)','FontName','Arial','FontSize',16)
title([subject_folder{s}, ' CPP (resp-locked)'])
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
legend(h,{'Motion','Colour'},'FontSize',16,'Location','NorthWest');
pause(1)
saveas(h1,['Figures/CPPindividual/CPPresp', subject_folder{s} '.jpg'])

%% single trials for DDM
k=0; CPPrAll=[]; RTAll = [];CPPAll=[];stimAll=[];congAll=[];responseAll=[];N2cAll=[]; N2iAll=[]; accAll = [];
clear N2cTemp N2iTemp CPPrTemp tTemp tTempi trTemp    
for c=2
    for d=1:2
        warning off
        if d==1 %motion
            window_searchi = [300 500];
            window_searchc = [319 369]; %JC: 324-424, S: 294-394, 150-600 works
            window_searchpc = [200 500];
            window_searchCPP = [-250 -5];
        elseif d==2 %colour
            window_searchi = [300 500];
            window_searchc = [319 369];
            window_searchpc = [150 300];
            window_searchCPP = [-200 -5];
        end
        window_slope= window_searchCPP;         
        for bin=1:no_of_bins
            for cong=1:2
                if ~isempty(conds{s,c,d,cong,bin})
                    %size(conds_bin{s,c,iti,side,bin},2)
                    for jj=1:size(conds{s,c,d,cong,bin},2)                   
                        k=k+1;
                        kk = conds{s,c,d,cong,bin}(jj);
                        mediationsS.run(k) = allrun_count(kk);
                        mediationsS.block(k) = allblock_count(kk);
                        mediationsS.stim(k) = allCMblock_count(kk); %colourMotion
                        mediationsS.bin(k)=bin;
                        mediationsS.cong(k)=cong;
                        mediationsS.response(k)=allacc{d}(kk);
                        mediationsS.rt(k)=allRT(kk)*1000/fs;
                        mediationsS.trial(k)= kk;
                        %N2c N2i CPP
                        N2cTemp = squeeze(movmean(erp(ch_lr,:,kk),1));
                        N2iTemp = squeeze(movmean(erp(ch_rl,:,kk),1));
                        CPPTemp = squeeze(movmean(mean(erp(ch_CPP,:,kk),1),1));
                        CPPTempOnset = squeeze(mean(erp(ch_CPP,:,kk),1));
                        
                        N2pcTemp=squeeze(movmean(erp(ch_lr,:,kk)-erp(ch_rl,:,kk),100));
                        
                        
                        [blah,RTsamp] = min(abs(t*fs/1000-allRT(kk)));
                        CPPrTemp = CPPTemp(RTsamp+trs);
                        CPPrA = squeeze(mean(erp(ch_CPP,RTsamp+trs,kk),1));
                        
                        [mediationsS.N2cpeak(k)] = mean(N2cTemp(t>window_searchc(1) & t<window_searchc(2)));
                        [mediationsS.N2ipeak(k)] = mean(N2iTemp(t>window_searchi(1) & t<window_searchi(2)));%200 to accomodate for N2c
                        [mediationsS.N2pcpeak(k)] = mean(N2pcTemp(t>window_searchpc(1) & t<window_searchpc(2)));
                        
                        mediationsS.CPPamp(k) = mean(CPPTemp(t>window_searchc(1) & t<window_searchc(2)));
                        %%%
                        %tTemp = t(t>window_searchc(1) & t<window_searchc(2)); mediationsS.N2cLatency(k) = tTemp(indxN2c);
                        %tTemp = t(t>window_searchi(1) & t<window_searchi(2)); mediationsS.N2iLatency(k) = tTemp(indxN2i);
                        
                        %tTemp = t(t>window_searchpc(1) & t<window_searchpc(2)); mediationsS.N2pcLatency(k) = tTemp(indxN2pc);
                        
                        %slopes
                        CPPrslopeTemp = movingslope(CPPrTemp,20,4);
                        coef = fitlm(tr(tr>window_slope(1) & tr<window_slope(2)),CPPrTemp(tr>window_slope(1) & tr<window_slope(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                        mediationsS.CPPrslope(k)=coef.Coefficients.Estimate(2);
                        mediationsS.RT(k) = allRT(kk)*1000/fs;
                        mediationsS.CPPlevel(k) = CPPTemp(t==mediationsS.RT(k));
                        [mediationsS.CPPamplitude(k), indxCPPL] = max(CPPrTemp(tr>-100 & tr<100));
                        tTemp = tr(tr>-100 & tr<100); mediationsS.CPPtime(k)= tTemp(indxCPPL);
                        
                        clear tindx onsetIndx
                        CPPonsetTemp = (CPPTempOnset(t>100 & t< mediationsS.RT(k)));
                        onsetIndx = find((diff(sign(CPPonsetTemp)))>0, 1); tindx = (t(t>100 & t<mediationsS.RT(k)));
                        if isempty(onsetIndx); mediationsS.CPPonset(k)=NaN; else %CPP never go below 1
                            mediationsS.CPPonset(k)= tindx(onsetIndx);
                        end
                        
                        clear slopeMax indx CPPslopeTemp
                        
                        N2iAll = [N2iAll; N2iTemp];
                        N2cAll = [N2cAll; N2cTemp];
                        CPPrAll = [CPPrAll;CPPrTemp];
                        RTAll = [RTAll; mediationsS.RT(k)];
                        stimAll = [stimAll; mediationsS.stim(k)];
                        congAll = [congAll; mediationsS.cong(k)];
                        CPPAll = [CPPAll;CPPTemp];
                        responseAll = [responseAll; mediationsS.response(k)];
                    end
                end
            end
        end
    end
end

kk=0;
for d=1:2
            kk=kk+1;
            CPPplot{kk} = CPPAll(stimAll==d & congAll==cong,:);
            CPPrplot{kk}= CPPrAll(stimAll==d & congAll==cong,:);
            RTplot{kk}= RTAll(stimAll==d & congAll==cong);
end
%% save the data for faster plots
if CSD
    save(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_CSD_319369'  num2str(single_participants) '_' num2str(chP)],'RTs','t','plot_chans','tr','chanlocs',...
        'allBins','allsubj','subject_folder', 'allstuff','mediations','RT_side',...
        'conds','condsT','mediationsS','CPPAll','CPPrAll','RTAll','stimAll','congAll','N2cAll','N2iAll','responseAll');
else
    save(['Data/ERPs/group_plots_erp_DDM_bin_ST_dual8Hz_294394'  num2str(single_participants) '_' num2str(chP)],'RTs','t','plot_chans','tr','chanlocs',...
        'allBins','allsubj','subject_folder', 'allstuff','mediations','RT_side',...
        'conds','condsT','mediationsS','CPPAll','CPPrAll','RTAll','stimAll','congAll','N2cAll','N2iAll','responseAll');
end
%Plot_original