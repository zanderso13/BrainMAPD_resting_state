load('300ROI_atlas_obj.mat')

% Possible networks
% 'SomatomotorDorsal','CinguloOpercular','Auditory','DefaultMode'
% 'ParietoMedial','Visual','FrontoParietal','Salience','VentralAttention');
% 'DorsalAttention','MedialTemporalLobe','Reward','SomatomotorLateral');

atl_temp = atl;
ind = find(strcmp(atl.labels,'MedialTemporalLobe'));

for i = 1:300
    if sum(i==ind)
        atl_temp.dat(atl_temp.dat == i) = 1;
    else
        atl_temp.dat(atl_temp.dat == i) = 0;
    end
end

montage(atl_temp)

