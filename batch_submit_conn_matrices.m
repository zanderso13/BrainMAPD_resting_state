% directories
% first is where your stats files will be output to
preprocdir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/fmriprep';
datadir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels_no_gsr';
repodir = '/home/zaz3744/repo/';
scriptdir = '/projects/b1108/projects/BrainMAPD_preproc_rest_T1_only/first_levels_no_gsr/quest_submission';

file_list = filenames(fullfile(preprocdir,strcat('sub-*/ses-2/func/sub*preproc_bold.nii')));
for i = 1:length(file_list)
    sublist{i} = file_list{i}(70:74);
end
%sublist = string(sublist);

cd(scriptdir)
keyboard
for sub = 1:length(sublist)
    PID = sublist{sub};
    
%    apply_300ROI(PID)

          s = ['#!/bin/bash\n\n'...
       '#SBATCH -A p30954\n'...
       '#SBATCH -p short\n'...
       '#SBATCH -t 00:30:00\n'...  
       '#SBATCH --mem=64G\n\n'...
       'matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath(''' repodir ''')); cingulate_to_wholebrain(' PID '); quit"\n\n'];
   
       scriptfile = fullfile(scriptdir, 'func_conn_script.sh');
       fout = fopen(scriptfile, 'w');
       fprintf(fout, s);
     
     
       !chmod 777 func_conn_script.sh
       !sbatch func_conn_script.sh
        
end

