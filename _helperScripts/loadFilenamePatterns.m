% Filename defintions for this dataset
%
%   Other m-files required: setDatasetPaths
%   MAT-files required: none
%
%   See also: none

%   Author: Hanno Winter
%   Date: 11-Feb-2020; Last revision: 05-Jun-2020

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

%% Checks

if ~exist('loadFilenamePatterns_loaded','var')
    loadFilenamePatterns_loaded = false;
elseif loadFilenamePatterns_loaded
    % fprintf('Already loaded filname patterns!\n');
    return  
end % if

%% Parameters

gnss_thales_parameters_root_string = 'gnss_thales_parameters';
gnss_inatm200stn_parameters_root_string = 'gnss_inatm200stn_parameters';
imu_inatm200stn_parameters_root_string = 'imu_inatm200stn_parameters';

%% Raw Data

gnss_thales_raw_data_root_string = 'gnss_thales_raw_data';
gnss_inatm200stn_raw_data_root_string = 'gnss_inatm200stn_raw_data';
imu_inatm200stn_raw_data_root_string = 'imu_inatm200stn_raw_data';
ref_inatm200stn_raw_data_root_string = 'ref_inatm200stn_raw_data';

%% Processed Data

gnss_thales_processed_data_root_string = 'gnss_thales_processed_data';
gnss_inatm200stn_processed_data_root_string = 'gnss_inatm200stn_processed_data';
imu_inatm200stn_processed_data_root_string = 'imu_inatm200stn_processed_data';
ref_inatm200stn_processed_data_root_string = 'ref_inatm200stn_internal_ekf_data';

%% Reference Data

ref_inatm200stn_reference_data_root_string = ref_inatm200stn_processed_data_root_string;
ref_ekfFusion_inatm200stn_reference_data_root_string = 'ref_ekfFusion_inatm200stn_data';

%% Finish

loadFilenamePatterns_loaded = true;
