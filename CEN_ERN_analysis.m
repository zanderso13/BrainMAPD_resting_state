
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/trilevel.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/immune_data.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/meds.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/demographics.mat')
load('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/outcomes/scan_day.mat')

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
    '20915','20934','20949','21223','21325','21539',... 
    '20235','21661'};

fnames = filenames(fullfile('/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/conn_matrices_AIB/*ERN*.mat'));

% remove mat files from subject list
for badsub = 1:length(exclusions)
    fnames(contains(fnames,exclusions{badsub}),:)=[];
end
    
nsubs = length(fnames);

for sub = 1:length(fnames)
    id=fnames{sub}(72:76);
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

    load(fnames{sub});
    all_data(:,:,sub) = corr_mat;
end

% calculate average CEN/ERN connectivity

% curr_data = atanh(all_data);
% curr_data(curr_data==Inf) = 1;

% for sub = 1:size(curr_data,3)
%     temp = reshape(tril(curr_data(:,:,sub)),[size(curr_data,1)^2,1]);
%     mean_conn(sub,1) = double(mean(temp));
% end

for sub = 1:size(all_data,3)
    mean_conn(sub,1) = mean2(all_data(:,:,sub));
end

imm_comp = array2table(imm_comp); imm_comp.Properties.VariableNames = {'immune'};
mean_conn = array2table(mean_conn); mean_conn.Properties.VariableNames = {'conn'};
sex = array2table(sex); sex.Properties.VariableNames = {'sex'};
site = array2table(site); site.Properties.VariableNames = {'site'};
curr_meds = array2table(curr_meds); curr_meds.Properties.VariableNames = {'meds'};
dat = [symp, imm_comp, mean_conn, sex, site, curr_meds];

dat_male = dat(dat.sex==0,:);
dat_female = dat(dat.sex==1,:);

mdl = fitlm(dat,'conn ~ immune')
mdl1 = fitlm(dat, 'conn ~ immune + sex + meds')
mdl2 = fitlm(dat, 'conn ~ GenDis + Anhedon + Fears + meds + sex')

