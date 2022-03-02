% Import time limited processed data
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
%   Date: 23-Nov-2020; Last revision: 24-Nov-2020

%% Init

% Settings ________________________________________________________________

start_time = '2018-10-24 12:10:00'; % format: uuuu-MM-dd HH:mm:ss
end_time = '2018-10-24 12:20:00'; % format: uuuu-MM-dd HH:mm:ss

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

if ~exist('import_limited_processed_data_executed','var')
    import_limited_processed_data_executed = false;
elseif import_limited_processed_data_executed
    fprintf('Import limited processed data script already executed!\n');
    return  
end % if

%% Processed Data
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01','session01','session02', etc
%
%   If the data has been chunked into several parts it is import to provide
%   the order of the parts.
%   

gnss_thales_processed_data_paths = { ... 
                                     fullfile('gnss_thales_processed_data_session01.csv'); ... 
                                     fullfile('gnss_thales_processed_data_session02.csv'); ... 
                                     fullfile('gnss_thales_processed_data_session03.csv')  ... 
                                   };

gnss_inatm200stn_processed_data_paths = { ... 
                                          fullfile('gnss_inatm200stn_processed_data_session01.csv'); ... 
                                          fullfile('gnss_inatm200stn_processed_data_session02.csv'); ... 
                                          fullfile('gnss_inatm200stn_processed_data_session03.csv') ... 
                                        };
                         
imu_inatm200stn_processed_data_paths = { ... 
                                         fullfile('imu_inatm200stn_processed_data_session01.csv'); ... 
                                         fullfile('imu_inatm200stn_processed_data_session02.csv'); ...
                                         fullfile('imu_inatm200stn_processed_data_session03_part01.csv'); ...
                                         fullfile('imu_inatm200stn_processed_data_session03_part02.csv'); ...
                                         fullfile('imu_inatm200stn_processed_data_session03_part03.csv') ... 
                                       };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNSS: Thales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gnss_thales_limited_processed_data_root_string = ... 
    [ ... 
      gnss_thales_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(gnss_thales_processed_data_paths,gnss_thales_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_thales_limited_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_thales_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_thales_processed_data_root_string,out_data); 
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNSS: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gnss_inatm200stn_limited_processed_data_root_string = ... 
    [ ... 
      gnss_inatm200stn_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(gnss_inatm200stn_processed_data_paths,gnss_inatm200stn_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_inatm200stn_limited_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_inatm200stn_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_inatm200stn_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imu_inatm200stn_limited_processed_data_root_string = ... 
    [ ... 
      imu_inatm200stn_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(imu_inatm200stn_processed_data_paths,imu_inatm200stn_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[imu_inatm200stn_limited_processed_data_root_string,'_session',sprintf('%02i',i)],imu_inatm200stn_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',imu_inatm200stn_processed_data_root_string,out_data);
clear out_data

%% Finish script

import_limited_processed_data_executed = true;
