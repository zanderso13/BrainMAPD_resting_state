hyperalignment = 0;

load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/trilevel.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/immune_data.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/meds.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/demographics.mat')
load('/Users/zacharyanderson/Documents/GitHub/BrainMAPD_resting_state/300ROI_atlas_obj.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/scan_day.mat')
load(fullfile('/Users/zacharyanderson/Documents/GitHub/BrainMAPD_resting_state/compiled_data/exclusions_based_on_neuro/outcomes.mat'))
load(fullfile('/Users/zacharyanderson/Documents/GitHub/BrainMAPD_resting_state/compiled_data/exclusions_based_on_neuro/CTIall.mat'))


exclusions = {'10001','10034','10041',...
    '10059','10067','10081','10088','10111','10135',...
    '10141','10196','10282','10309','10336','10438','10439',...
    '20085','20108','20644','21238','21257','21268',...
    '20133','20460','20948','20309','20996','21001',...
    '21523'}; % start immune exclusions
%     '10010','10041','10066','10092','10141',...
%     '10176','10186','10256','10264','10271','10274',...
%     '10315','10322','10328','10332','10341','10434',...
%     '10454','10458','10471','20083','20124','20302',...
%     '20507','20530','20630','20897','20902','20903',...
%     '20915','20934','20949','21223','21325','21539'};

fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/conn_matrices_hyperalignment/*.mat'));

% remove mat files from subject list
for badsub = 1:length(exclusions)
    fnames(contains(fnames,exclusions{badsub}),:)=[];
end
    
nsubs = length(fnames);

for sub = 1:length(fnames)
     id=fnames{sub}(83:87);
%     currmotion = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/FD',strcat('*',id,'*.mat')));
%     load(currmotion{1})
    % motion data
%     if isempty(FD)
%         allFD(sub,:)=NaN;
%     else
%         allFD(sub,:)=mean(FD.framewise_displacement);
%     end
    % symptom data
    if isempty(symp(symp.id==str2double(id),:))
        trilevel(sub,:)=NaN;
    else
        trilevel(sub,:)=symp(symp.id==str2double(id),:);
    end
    % immune data
    if isempty(immune.T1BDicsavg(immune.PID==str2double(id)))
        imm_comp(sub,1)=NaN;
    else
        imm_comp(sub,1)=immune.T1BDicsavg(immune.PID==str2double(id));
    end
    % med data
    if isempty(meds.T1SCcmipsychany(meds.PID==str2double(id)))
        curr_meds(sub,1)=0;
    else
        curr_meds(sub,1)=meds.T1SCcmipsychany(meds.PID==str2double(id));
    end
    % demographic data
    if str2double(id) < 20000
        site(sub,1) = 1;
    else
        site(sub,1) = 0;
    end
    % scan site
    if isempty(dem.sex(dem.PID==str2double(id)))
        sex(sub,1)=0;
        dob(sub,1)=NaN;
    else
        sex(sub,1)=dem.sex(dem.PID==str2double(id));
        dob(sub,1)=dem.dob(dem.PID==str2double(id));
    end

    % scan date
    if isempty(scan_day.T1M1date(scan_day.PID==str2double(id)))
        scanday(sub,1)=0;
    else
        scanday(sub,1)=scan_day.T1M1date(scan_day.PID==str2double(id));
    end

    if isempty(lifestress.SEP_Emaj_sum(lifestress.PID==str2double(id)))
        ctidat(sub,1)=NaN;
    else
        ctidat(sub,:)=[lifestress.SEP_Emaj_sum(lifestress.PID==str2double(id)),...
            lifestress.SEP_Lmaj_sum(lifestress.PID==str2double(id))];
    end

    temp = load(fnames{sub});
    
%     seeds(:,:,sub) = [temp.r(27).dat,temp.r(28).dat,...%mofc l and r
%         temp.r(29).dat,temp.r(30).dat,... % anterior ofc l and r
%         temp.r(31).dat,temp.r(32).dat,... % posterior ofc l and r
%         temp.r(33).dat,temp.r(34).dat]; % lateral ofc l and r];
%     seeds_nina(:,:,sub) = [temp.r(145).dat,temp.r(146).dat,...%subgenual ACC
%         temp.r(147).dat,temp.r(148).dat,... % pregenual ACC
%         temp.r(149).dat,temp.r(150).dat]; % superior ACC
%     for targ = 1:length(temp.r)
%         targets_nina(:,targ,sub) = temp.r(targ).dat;
%     end
%     % remove rois that correspond with seeds
%     targets_nina(:,150,:) = [];targets_nina(:,149,:) = [];
%     targets_nina(:,148,:) = [];targets_nina(:,147,:) = [];
%     targets_nina(:,146,:) = [];targets_nina(:,145,:) = [];
%     seeds_hyp(:,:,sub) = [temp.r(27).all_data,temp.r(28).all_data,... %mofc l and r
%         temp.r(29).all_data,temp.r(30).all_data,... % anterior ofc l and r
%         temp.r(31).all_data,temp.r(32).all_data,... % posterior ofc l and r
%         temp.r(33).all_data,temp.r(34).all_data]; % lateral ofc l and r
%     targets(:,:,sub) = [temp.r(43).dat,temp.r(44).dat,... % amyg l and r
%         temp.r(73).dat,temp.r(74).dat,... % caudate l and r
%         temp.r(75).dat,temp.r(76).dat,... % putamen l and r
%         temp.r(77).dat,temp.r(78).dat,... % pallidum l and r
%         temp.r(151).dat,temp.r(152).dat]; % NAcc l and r
    

end
keyboard
% load('300ROI_atlas_obj.mat')
% 
% SMD = strcmp(atl.labels,'SomatomotorDorsal');
% CO = strcmp(atl.labels,'CinguloOpercular');
% A = strcmp(atl.labels,'Auditory');
% DMN = strcmp(atl.labels,'DefaultMode');
% PM = strcmp(atl.labels,'ParietoMedial');
% V = strcmp(atl.labels,'Visual');
% FP = strcmp(atl.labels,'FrontoParietal');
% S = strcmp(atl.labels,'Salience');
% VA = strcmp(atl.labels,'VentralAttention');
% DA = strcmp(atl.labels,'DorsalAttention');
% MTL = strcmp(atl.labels,'MedialTemporalLobe');
% Rew = strcmp(atl.labels,'Reward');
% SML = strcmp(atl.labels,'SomatomotorLateral');

% input what network you want here. You can also do addition to look at
% several networks at once.
% seed to seed associations: 
% VS - mOFC: 
% idx1 = [246,247]; idx2 = [108,116];
% Amyg - mOFC: 
% idx1 = [244,245]; idx2 = [108,116];
% Amyg - VS: 
% idx1 = [246,247]; idx2 = [244,245];
% Network analysis example for connectivity within reward network: idx1 = Rew; idx2 = Rew;
% idx1 = Rew; idx2 = Rew;

% curr_analysis = all_data(idx1,idx2,:);
% 
% for sub = 1:size(curr_analysis,3)
%     conn(sub,1) = mean2(curr_analysis(:,:,sub));
% end
% 
% dat = [symp,array2table(imm_comp),array2table(conn),array2table(sex),array2table(curr_meds),array2table(site),array2table(allFD)]; 
% 
% fprintf(strcat('How many motion exclusions: ',num2str(sum(dat.allFD>0.3))))
% dat(dat.allFD>0.3,:)=[];
% 
% dat(isnan(dat.imm_comp),:) = [];
% if hyperalignment == 1
%     dat = dat(1:nsubs,:);
% end

%% age at first scan calculator

% mean((scanday-dob)/8760); % 8760 is the number of hours in a year
% std((scanday-dob)/8760);

%% hyperalignment
if hyperalignment == 1
    PID_list = dat.id;
    % [alignedvs1,alignedvs2,alignedamyg1,alignedamyg2,alignedofc1,alignedofc2,transformvs1, transformvs2, transformamyg1, transformamyg2, transformofc1, transformofc2] = BrainMAPD_hyperalign(PID_list);
    
end

%% Final cleaning
% dat(dat.Anhedon(:) > 3,:) = []; dat(dat.Anhedon(:) < -3,:) = [];
% dat(dat.GenDis(:) > 3,:) = []; dat(dat.GenDis(:) < -3,:) = [];
% dat(dat.Fears(:) > 3,:) = []; dat(dat.Fears(:) < -3,:) = [];
% dat(zscore(dat.conn)>3,:) = []; dat(zscore(dat.conn)<-3,:) = [];
% dat(dat.imm_comp>3,:) = []; dat(dat.imm_comp<-3,:) = [];


%% compare conventional analyses with same thing on aligned data

% mdl1a = fitlm(dat,'GenDis ~ conn + curr_meds + sex + site');
% mdl1b = fitlm(dat,'GenDis ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl2a = fitlm(dat,'Anhedon ~ conn + curr_meds + sex + site');
% mdl2b = fitlm(dat,'Anhedon ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl3a = fitlm(dat,'Fears ~ conn + curr_meds + sex + site');
% mdl3b = fitlm(dat,'Fears ~ conn + imm_comp + imm_comp*conn + curr_meds + sex + site');
% mdl4 = fitlm(dat,'conn ~ GenDis + Anhedon + Fears + curr_meds + sex + site');
% plotAdjustedResponse(mdl4,'Anhedon')

if hyperalignment == 1
    
    for sub = 1:length(alignedofc2)
        amyg(:,:,sub) = [zscore(alignedamyg1{sub}),zscore(alignedamyg2{sub})];
        vs(:,:,sub) = [zscore(alignedvs1{sub}),zscore(alignedvs2{sub})];
        ofc(:,:,sub) = [zscore(alignedofc1{sub}),zscore(alignedofc2{sub})];
        ha_conn(:,:,sub) = corr(vs(:,:,sub),ofc(:,:,sub));
    end
    
    dat = [dat,array2table(ha_conn)];

%     ha_mdl1 = fitlm(dat,'ha_conn ~ imm_comp + GenDis + imm_comp*GenDis + curr_meds + sex + site')
%     ha_mdl2 = fitlm(dat,'ha_conn ~ imm_comp + Anhedon + imm_comp*Anhedon + curr_meds + sex + site')
%     ha_mdl3 = fitlm(dat,'ha_conn ~ imm_comp + Fears + imm_comp*Fears + curr_meds + sex + site')
    ha_mdl4 = fitlm(dat,'ha_conn ~ GenDis + Anhedon + Fears + curr_meds + sex + site')
    
end
