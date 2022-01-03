% need to apply an atlas object to a series of subjects. This script has
% some dependencies. Names CanlabCore (freely available on Github) and
% other elements in this repository.

% 1. Identify where subject files are held. This will be designed to read
% in residual files generated following the application of SPM12 to resting
% state data. This means that preprocessing and data cleaning have already
% occurred at this point.

datadir = '/Users/zacharyanderson/Documents/ACNlab/BrainMAPD/RS/test_firstlevel_out';

% 2. Load in residual nii's using functions from CanlabCore.

dat = fmri_data(filenames(fullfile(datadir, '*/ses-2/run-1/rest/Res*nii')));

% display some QA

descriptives(dat)

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