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

datadir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels_no_gsr';
outdir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/conn_matrices_AIB';
% 2. Load in residual nii's using functions from CanlabCore.

dat = fmri_data(filenames(fullfile(datadir,strcat('sub-',num2str(PID),'/ses-2/run-1/rest/Res_0*nii'))));

% display some QA

% descriptives(dat)

% 3. Load atlas object, currently saved as a .mat so this is easy

load('CEN_atlas.mat')

% 4. Go into loop that will create a series of vectors for each of the 300
% regions of interest

% r = extract_roi_averages(dat,atl);

% For CEN/ERN I'm creating a slightly different workflow. This works by
% just directly pulling from the data object, voxels that are currently
% tagged as within each network. 
network = dat.dat(logical(atl.dat),:);

% 5. One instance of the above command yields a 1x300 'region object'. This is effectively
% a series of cells where average timecourse data from each ROI is stored
% in r(ROI_number).dat. The other, yields the network variable. This
% reflects the timeseries data of all voxels within the network that can
% now be correlated

% Pull this data into a giant matrix.

% for i = 1:length(r)
%     temp_matrix(:,i) = r(i).dat(:,1);
% end

% 6. corr2 function will now turn this into a 300x300 correlation matrix

corr_mat = corr(network');
curr_fname = fullfile(outdir,strcat(num2str(PID), '_CEN_matrix.mat'));
save(curr_fname, 'corr_mat')

end
