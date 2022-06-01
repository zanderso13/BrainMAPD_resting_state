% need to apply an atlas object to a series of subjects. This script has
% some dependencies. Names CanlabCore (freely available on Github) and
% other elements in this repository.

function cingulate_to_wholebrain(PID)

if nargin==0 % defaults just for testing 
    % Define some 
    PID = "10006";
    
end
% 1. Identify where subject files are held. This will be designed to read
% in residual files generated following the application of SPM12 to resting
% state data. This means that preprocessing and data cleaning have already
% occurred at this point.

datadir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels_no_gsr';
outdir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/conn_matrices_Nina';
% 2. Load in residual nii's using functions from CanlabCore.

dat = fmri_data(filenames(fullfile(datadir,strcat('sub-',num2str(PID),'/ses-2/run-1/rest/Res_0*nii'))));

% display some QA

% descriptives(dat)

% 3. Load atlas object, currently saved as a .mat so this is easy

load('AAL3_1mm.mat')

% 4. Go into loop that will create a series of vectors for each of the 300
% regions of interest

r = extract_roi_averages(dat,atl);

% 5. The above command yields a 1x160 'region object'. Now, we want to pull
% those regions that are specific to the cingulate cortex. These will act
% as our seeds. Then we correlate those seeds with the rest of the brain 

subgen = mean([r(145).dat,r(146).dat],2);
pregen = mean([r(147).dat,r(148).dat],2);
sup = mean([r(149).dat,r(150).dat],2);
for voxel = 1:length(dat.dat)
     subgen_mat(voxel,1) = atanh(corr(dat.dat(voxel,:)', subgen));
     pregen_mat(voxel,1) = atanh(corr(dat.dat(voxel,:)', pregen));
     sup_mat(voxel,1) = atanh(corr(dat.dat(voxel,:)', sup));
end

% 6. save resulting variables for each subject

curr_fname = fullfile(outdir,strcat(num2str(PID), '_ACC_matrix.mat'));
save(curr_fname, 'subgen_mat', 'pregen_mat', 'sup_mat')

end