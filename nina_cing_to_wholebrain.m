
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/trilevel.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/immune_data.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/meds.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/demographics.mat')
load('/Users/zacharyanderson/Documents/GitHub/BrainMAPD_resting_state/300ROI_atlas_obj.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/scan_day.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/first_year_project_BrainMAPD_clinical_diagnoses_final.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/childhood_trauma.mat')

exclusions = {'10001','10041','10034','10041',...
    '10059','10067','10081','10088','10111','10135',...
    '10141','10196','10282','10309','10336','10438','10439',...
    '20085','20108','20644','21238','21257','21268',...
    '20133','20460','20948','20309','20996','21001',...
    '21523'}; %,... % start immune exclusions
%     '10010','10041','10066','10092','10141',...
%     '10176','10186','10256','10264','10271','10274',...
%     '10315','10322','10328','10332','10341','10434',...
%     '10454','10458','10471','20083','20124','20302',...
%     '20507','20530','20630','20897','20902','20903',...
%     '20915','20934','20949','21223','21325','21539'};


fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/conn_matrices_Nina/*.mat'));

% remove mat files from subject list
for badsub = 1:length(exclusions)
    fnames(contains(fnames,exclusions{badsub}),:)=[];
end
    
nsubs = length(fnames);

for sub = 1:length(fnames)
    id=fnames{sub}(78:82);
    
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
    
    % trauma
    if isempty(traum.SEP_Emaj_sum(traum.PID==str2double(id))) 
        Emaj(sub,1)=NaN;
        Lmaj(sub,1)=NaN;
    else
        Emaj(sub,1)=traum.SEP_Emaj_sum(traum.PID==str2double(id));
        Lmaj(sub,1)=traum.SEP_Lmaj_sum(traum.PID==str2double(id));
    end
end

load('AAL3_1mm.mat')

%% load in whole brain data
% this comes from a bunch of cingulate to whole brain correlations
% initialize fmri_data object

maskimage = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/mask.nii');
dat_pregen = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/mask.nii'); dat_pregen.dat = [];
dat_subgen = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/mask.nii'); dat_subgen.dat = [];
dat_sup = fmri_data('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/nina/mask.nii'); dat_sup.dat = [];

for sub = 1:length(fnames)
    load(fnames{sub});
    pregen_mat .* maskimage.dat;
    dat_pregen.dat(:,sub) = double(pregen_mat);
    dat_subgen.dat(:,sub) = double(subgen_mat);
    dat_sup.dat(:,sub) = double(sup_mat);
    clear subgen_mat sup_mat pregen_mat
end



