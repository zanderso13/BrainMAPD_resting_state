% note: matlab doesn't like tsv's, so I've renamed all files to have .txt
% extension. This makes it so that the readtable function works
% appropriately.

confounddir = '/Users/zacharyanderson/Documents/ACNlab/ETSI/fmriprep';
outputdir = '/Users/zacharyanderson/Documents/ACNlab/ETSI/sub_confounds';

cd(confounddir)
id_length = 5; % length of sub id. Ex: sub-10001
number_of_images_dropped = 0; % fsl wants the extra confounds file to match time course data after it removes images

fnames = filenames(fullfile('sub-*/ses-1/func/sub-*_task-rest_run-1_desc-confounds_timeseries.txt'));


for sub = 1:length(fnames)
    id = fnames{sub}(1:id_length);
    T = readtable(fnames{sub});

    GSR = T.global_signal;
    WM = T.white_matter;
    CSF = T.csf;
    FD = T.framewise_displacement;
    DVARS = T.dvars;
    
    outliers = table2array(T(:,contains(T.Properties.VariableNames,'motion')));
    transx = table2array(T(:,contains(T.Properties.VariableNames,'trans_x')));
    transy = table2array(T(:,contains(T.Properties.VariableNames,'trans_y')));
    transz = table2array(T(:,contains(T.Properties.VariableNames,'trans_z')));
    rotx = table2array(T(:,contains(T.Properties.VariableNames,'rot_x')));
    roty = table2array(T(:,contains(T.Properties.VariableNames,'rot_y')));
    rotz = table2array(T(:,contains(T.Properties.VariableNames,'rot_z')));
    
    totmotion = [transx,transy,transz,rotx,roty,rotz,transx.^2,transy.^2,transz.^2,rotx.^2,roty.^2,rotz.^2,outliers]; 

    allconfounds = [totmotion,WM,CSF,FD];

    allconfounds(isnan(allconfounds)) = 0;
    allconfounds = allconfounds(number_of_images_dropped+1:size(T,1),:);
    confounds_fname = strcat(id,'_run-1_confounds.txt');
%     motionfname = strcat(id,'motion.txt');
%     GSRfname = strcat(id,'_gsr.txt');
%     WMfname = strcat(id,'_wm.txt');
%     CSFfname = strcat(id,'_csf.txt');
%     FDfname = strcat(id,'_fd.txt');
%     DVARSfname = strcat(id,'_dvars.txt');
    
    if ~exist(fullfile(outputdir,id))
        mkdir(fullfile(outputdir,id))
    end

    writematrix(allconfounds,fullfile(outputdir,id,confounds_fname),'Delimiter','tab')
%     writematrix(GSR,fullfile(outputdir,id,GSRfname))
%     writematrix(WM,fullfile(outputdir,id,WMfname))
%     writematrix(CSF,fullfile(outputdir,id,CSFfname))
%     writematrix(FD,fullfile(outputdir,id,FDfname))
%     writematrix(DVARS,fullfile(outputdir,id,DVARSfname))

end
