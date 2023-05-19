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

r = atlas2region(atl);

counter = 1;
for reg = 1:length(r)
    if contains(r(reg).shorttitle,'ACC_') == 1
        ACCmask(counter) = r(reg);
        counter = counter + 1;
    end
end

subgen_ACCmask = ACCmask(1:2);
pregen_ACCmask = ACCmask(3:4);
superior_ACCmask = ACCmask(5:6);

ACCdat = extract_data(subgen_ACCmask,dat);
mean_ACC_signal = (ACCdat(1).dat + ACCdat(2).dat) ./ 2;

% 5. The above command yields a 1x2 'region object' that is currently specific to the subgenual ACC. 
% Create a loop, that correlates and does r to z transformation on every
% brain voxel with average subgen ACC activation.

corr_dat = dat; corr_dat.dat = []; 
for voxel = 1:length(dat.dat)
     corr_dat.dat(voxel,1) = atanh(corr(dat.dat(voxel,:)', mean_ACC_signal));   
end

% 6. save resulting variables for each subject

curr_fname = fullfile(outdir,strcat(num2str(PID), '_ACCconnectivity.nii'));
write(corr_dat, 'fname', curr_fname)

end