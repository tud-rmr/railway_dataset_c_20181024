% Import map data
%
%   Other m-files required: 
%       - setDatasetPaths
%       - loadFilenamePatterns
%       - importFromFile
%       - exportToFile
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, processRawData

%   Author: Hanno Winter
%   Date: 17-Nov-2020; Last revision: 17-Nov-2020

%% Init

% Set path ________________________________________________________________

saved_pwd = pwd;
pwd_parts = strfind(pwd,filesep);

if ~exist('dataset_paths_set','var') || (dataset_paths_set == false)
    
    dataset_paths_set = false;
    
    for i = 1:length(pwd_parts)
        
        if exist('setDatasetPaths','file')
            setDatasetPaths();
            dataset_paths_set = true;
            break;
        elseif exist('_fcns','dir')
            cd('_fcns');
        else
            cd('..')
        end % if
        
    end % for i
    
    cd(saved_pwd);
    
end % if

if ~dataset_paths_set
    error('Couldn''t add dataset folders to path!')
end % if

% Load scripts ____________________________________________________________

loadFilenamePatterns

%% Checks

if ~exist('import_map_data_executed','var')
    import_map_data_executed = false;
elseif import_map_data_executed
    fprintf('Import map data script already executed!\n');
    return  
end % if

%% Map Data
                             
ref_track_map_data_paths = { ... 
                                fullfile('ref_track_map.csv') ... 
                              };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map: (X,Y)-Track-Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ref_track_map, ~, load_flag] = importFromFile(ref_track_map_data_paths,'ref_track_map','DataType','map');
if load_flag == 1 % save to .mat for faster access in the future
    writeToMatFile(ref_track_map,'ref_track_map','ref_track_map');
end % for i

%% Finish script

import_map_data_executed = true;
