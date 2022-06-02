
% I want to write an analysis to hyperalign time course data from each seed
% node in the reward network, then create connectivity profiles

fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/ROI300_conn_hyp/*.mat'));

load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/trilevel.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/immune_data.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/meds.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/demographics.mat')
load('/Users/zacharyanderson/Documents/GitHub/BrainMAPD_resting_state/300ROI_atlas_obj.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/scan_day.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/first_year_project_BrainMAPD_clinical_diagnoses_final.mat')

exclusions = {'10001','10041','10034','10041',...
    '10059','10067','10081','10088','10111','10135',...
    '10141','10196','10282','10309','10336','10438','10439',...
    '20085','20108','20644','21238','21257','21268',...
    '20133','20460','20948','20309','20996','21001',...
    '21523',... % start immune exclusions
    '10010','10041','10066','10092','10141',...
    '10176','10186','10256','10264','10271','10274',...
    '10315','10322','10328','10332','10341','10434',...
    '10454','10458','10471','20083','20124','20302',...
    '20507','20530','20630','20897','20902','20903',...
    '20915','20934','20949','21223','21325','21539',... % need to start some final exclusions where data isn't being extracted the way that it should. Eventually I need to dig into this more. But for the data blitz let's just drop them and move on
    '20235','21661'};

% remove mat files from subject list
for badsub = 1:length(exclusions)
    fnames(contains(fnames,exclusions{badsub}),:)=[];
end
    
nsubs = length(fnames);

for sub = 1:length(fnames)
    id=fnames{sub}(70:74);
    
    % symptom data
    if isempty(trilevel(trilevel.id==str2double(id),:))
        symp(sub,:)=NaN;
    else
        symp(sub,:)=trilevel(trilevel.id==str2double(id),:);
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
    % depression SCID
    if isempty(clinical_info.PID(clinical_info.PID==str2double(id)))
        dep(sub,1)=0;
        anx(sub,1)=0;
        com(sub,1)=0;
    else
        dep(sub,1)=clinical_info.dep_life_any(clinical_info.PID==str2double(id));
        anx(sub,1)=clinical_info.anx_life_any(clinical_info.PID==str2double(id));
        com(sub,1)=clinical_info.comorbid_life_dep_anx(clinical_info.PID==str2double(id));
    end
    % scan date
    if isempty(scan_day.T1M1date(scan_day.PID==str2double(id)))
        scanday(sub,1)=0;
    else
        scanday(sub,1)=scan_day.T1M1date(scan_day.PID==str2double(id));
    end
    
end
anydep = dep + com;
load('300ROI_atlas_obj.mat')

SMD = strcmp(atl.labels,'SomatomotorDorsal');
CO = strcmp(atl.labels,'CinguloOpercular');
A = strcmp(atl.labels,'Auditory');
DMN = strcmp(atl.labels,'DefaultMode');
PM = strcmp(atl.labels,'ParietoMedial');
V = strcmp(atl.labels,'Visual');
FP = strcmp(atl.labels,'FrontoParietal');
S = strcmp(atl.labels,'Salience');
VA = strcmp(atl.labels,'VentralAttention');
DA = strcmp(atl.labels,'DorsalAttention');
MTL = strcmp(atl.labels,'MedialTemporalLobe');
Rew = strcmp(atl.labels,'Reward');
SML = strcmp(atl.labels,'SomatomotorLateral');

% This is going to be different than my other script but I'm leaving this text. 
% I want to know which indices correspond to different parts of the region
% seed to seed associations: 
% VS - mOFC: 
% idx1 = [246,247]; idx2 = [108,116];
% Amyg - mOFC: 
% idx1 = [244,245]; idx2 = [108,116];
% Amyg - VS: 
% idx1 = [246,247]; idx2 = [244,245];
% Network analysis example for connectivity within reward network: idx1 = Rew; idx2 = Rew;
idx1 = Rew;


% This is kinda a weird loop, but it will create a four D matrix that is
% voxel by target region by subject by seed region. The goal will then be
% to separate things by roi. Connectivity profiles by each roi can then be
% hyperaligned across subs using Haxby et al, 2011 methods.
for sub = 1:length(fnames)
    load(fnames{sub});
    curr = r(idx1);
%     for region_check = 1:length(curr)
%         
%     check = sum(isnan(corr(curr(1,2).all_data)));
    for roi = 1:sum(idx1)
        % segment out the data we want to correlate. What I wanted to avoid
        % was correlating voxels from one roi, with the average signal in
        % that same roi. So this hopefully indexes them in a somewhat smart
        % way
        temp1 = curr(roi);
        temp2 = curr; temp2(roi)=[];
        for vox = 1:size(temp1.all_data,2) % create the conn profiles
            for targ = 1:size(temp2,2)
                conn_profile(vox,targ,sub,roi) = corr(temp1.all_data(:,vox),temp2(targ).dat(:));
            end
        end
    end
end

% best to remove exclusions here who don't have immune data
conn_profile(:,:,isnan(imm_comp),:)=[];
dep(isnan(imm_comp))=[]; anydep(isnan(imm_comp))=[];
symp(isnan(imm_comp),:)=[]; sex(isnan(imm_comp))=[];
curr_meds(isnan(imm_comp))=[];
imm_comp(isnan(imm_comp))=[];

% create an index and remove 
[B_imm,I_imm] = sort(imm_comp,'ascend');
[B_dep,I_dep] = sort(anydep,'ascend');
conn_profile=conn_profile(:,:,I_imm,:);

% remake inputs to cell arrays to input into hyperalign function
for sub = 1:size(conn_profile,3)
    prof1{sub} = conn_profile(:,:,sub,1); 
    prof2{sub} = conn_profile(:,:,sub,2); 
    prof3{sub} = conn_profile(:,:,sub,3); 
    prof4{sub} = conn_profile(:,:,sub,4); 
    prof5{sub} = conn_profile(:,:,sub,5); 
    prof6{sub} = conn_profile(:,:,sub,6); 
    prof7{sub} = conn_profile(:,:,sub,7); 
    prof8{sub} = conn_profile(:,:,sub,8);
end

[hyp1,tran1] = hyperalign(prof1{:});
[hyp2,tran2] = hyperalign(prof2{:});
[hyp3,tran3] = hyperalign(prof3{:});
[hyp4,tran4] = hyperalign(prof4{:});
[hyp5,tran5] = hyperalign(prof5{:});
[hyp6,tran6] = hyperalign(prof6{:});
[hyp7,tran7] = hyperalign(prof7{:});
[hyp8,tran8] = hyperalign(prof8{:});

%% perform rsa 

for sub1 = 1:length(prof1)
    for sub2 = 1:length(prof1)
        disim1(sub1,sub2) = 1 - corr2(prof1{sub1},prof1{sub2});
        disim2(sub1,sub2) = 1 - corr2(prof2{sub1},prof2{sub2});
        disim3(sub1,sub2) = 1 - corr2(prof3{sub1},prof3{sub2});
        disim4(sub1,sub2) = 1 - corr2(prof4{sub1},prof4{sub2});
        disim5(sub1,sub2) = 1 - corr2(prof5{sub1},prof5{sub2});
        disim6(sub1,sub2) = 1 - corr2(prof6{sub1},prof6{sub2});
        disim7(sub1,sub2) = 1 - corr2(prof7{sub1},prof7{sub2});
        disim8(sub1,sub2) = 1 - corr2(prof8{sub1},prof8{sub2});

        hdisim1(sub1,sub2) = 1 - corr2(hyp1{sub1},hyp1{sub2});
        hdisim2(sub1,sub2) = 1 - corr2(hyp2{sub1},hyp2{sub2});
        hdisim3(sub1,sub2) = 1 - corr2(hyp3{sub1},hyp3{sub2});
        hdisim4(sub1,sub2) = 1 - corr2(hyp4{sub1},hyp4{sub2});
        hdisim5(sub1,sub2) = 1 - corr2(hyp5{sub1},hyp5{sub2});
        hdisim6(sub1,sub2) = 1 - corr2(hyp6{sub1},hyp6{sub2});
        hdisim7(sub1,sub2) = 1 - corr2(hyp7{sub1},hyp7{sub2});
        hdisim8(sub1,sub2) = 1 - corr2(hyp8{sub1},hyp8{sub2});
    end
end

% create average disimilarity matrix
avg_disim = disim5 + disim6; %+ disim3 + disim4 + disim5 + disim6 + disim7 + disim8;
avg_disim = avg_disim ./ 2;

avg_disim_hyp = hdisim1 + hdisim2; %+ hdisim3 + hdisim4 + hdisim5 + hdisim6 + hdisim7 + hdisim8;
avg_disim_hyp = avg_disim_hyp ./ 2;

figure();heatmap(tril(avg_disim)+triu(avg_disim_hyp))

% take a look at transformation matrices

for sub = 1:length(prof1)
    avg_tran1(sub) = mean2(tran1{sub}.T);
    avg_tran2(sub) = mean2(tran2{sub}.T);
    avg_tran3(sub) = mean2(tran3{sub}.T);
    avg_tran4(sub) = mean2(tran4{sub}.T);
    avg_tran5(sub) = mean2(tran5{sub}.T);
    avg_tran6(sub) = mean2(tran6{sub}.T);
    avg_tran7(sub) = mean2(tran7{sub}.T);
    avg_tran8(sub) = mean2(tran8{sub}.T);
end

avg_tran = [avg_tran1',avg_tran2',avg_tran3',avg_tran4',...
    avg_tran5', avg_tran6', avg_tran7', avg_tran8'];

sorted_avg_tran = avg_tran(I_imm,:);

figure(); heatmap(corr(sorted_avg_tran'))

sex = array2table(sex); sex.Properties.VariableNames = {'sex'};
curr_meds = array2table(curr_meds); curr_meds.Properties.VariableNames = {'meds'};
dat_tran = zscore(mean(avg_tran,2)); dat_tran = array2table(dat_tran); dat_tran.Properties.VariableNames = {'transforms'};
imm_comp = array2table(imm_comp); imm_comp.Properties.VariableNames = {'immune'};

dat = [symp,sex,curr_meds,dat_tran,imm_comp];
%avg_tran_tot = mean(avg_tran,2); 
fitlm(dat,'transforms ~ immune')
datablitzmdl = fitlm(dat,'transforms ~ immune + sex + meds + GenDis + Anhedon + Fears')

