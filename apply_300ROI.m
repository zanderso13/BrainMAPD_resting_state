% need to apply an atlas object to a series of subjects. This script has
% some dependencies. Names CanlabCore (freely available on Github) and
% other elements in this repository.

function apply_300ROI(PID)

if nargin==0 % defaults just for testing 
    % Define some 
    PID = "10006";
    
end
% 1. Identify where subject files are held. This will be designed to read
% in residual files generated following the application of SPM12 to resting
% state data. This means that preprocessing and data cleaning have already
% occurred at this point.

datadir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels_hyperalignment';
outdir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/conn_300ROI_hyp';
% 2. Load in residual nii's using functions from CanlabCore.

dat = fmri_data(filenames(fullfile(datadir,strcat('sub-',num2str(PID),'/ses-2/run-1/rest/Res_0*nii'))));

% display some QA

% descriptives(dat)

% 3. Load atlas object, currently saved as a .mat so this is easy

load('300ROI_atlas_obj.mat')

% 4. Go into loop that will create a series of vectors for each of the 300
% regions of interest

r = extract_roi_averages(dat,atl);

% 5. The above command yields a 1x300 'region object'. This is effectively
% a series of cells where average timecourse data from each ROI is stored
% in r(ROI_number).dat

% Pull this data into a giant matrix.

for i = 1:length(r)
    temp_matrix(:,i) = r(i).dat(:,1);
end

% 6. corr2 function will now turn this into a 300x300 correlation matrix

corr_mat = corr(temp_matrix);
curr_fname = fullfile(outdir,strcat(num2str(PID), '_matrix.mat'));
save(curr_fname, 'r')

end
