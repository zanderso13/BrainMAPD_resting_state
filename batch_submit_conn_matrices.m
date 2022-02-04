% directories
% first is where your stats files will be output to
datadir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels';
repodir = '/home/zaz3744/repo/';

file_list = filenames(fullfile(datadir,strcat('sub-*/ses-',num2str(ses),'/func/ssub*MID*run-0',num2str(run),'*preproc_bold.nii')));
for i = 1:length(file_list)
    sublist{i} = file_list{i}(63:67);
end
sublist = string(sublist);

cd(scriptdir)
keyboard
for sub = 1:length(new_list)

    apply_300ROI(PID)

%          s = ['#!/bin/bash\n\n'...
%       '#SBATCH -A p30954\n'...
%       '#SBATCH -p short\n'...
%       '#SBATCH -t 00:15:00\n'...  
%       '#SBATCH --mem=30G\n\n'...
%       'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); apply_300ROI(' num2str(PID) '); quit"\n\n'];
%   
%       scriptfile = fullfile(scriptdir, 'func_conn_script.sh');
%       fout = fopen(scriptfile, 'w');
%       fprintf(fout, s);
%     
%     
%       !chmod 777 func_conn_script.sh
%       !sbatch func_conn_script.sh
keyboard     
end

