% Import raw data
%
%   Other m-files required: 
%       - setDatasetPaths
%       - loadFilenamePatterns
%       - importFromFile
%       - exportToFile
%   MAT-files required: none
%
%   See also: loadFilenamePatterns

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 18-Nov-2020

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

if ~exist('import_raw_data_executed','var')
    import_raw_data_executed = false;
elseif import_raw_data_executed
    fprintf('Import raw data script already executed!\n');
    return  
end % if

%% Raw Data
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01' or 'session01'.
% 
%   If the data has been chunked into several parts it is import to provide
%   the order of the parts.
%                     

gnss_thales_raw_data_paths = { ... 
                               fullfile('2018-10-24_Session01','01_raw','Gnss_ThalesUnknown','thales_gnss.csv'); ... 
                               fullfile('2018-10-24_Session02','01_raw','Gnss_ThalesUnknown','thales_gnss.csv'); ... 
                               fullfile('2018-10-24_Session03','01_raw','Gnss_ThalesUnknown','thales_gnss.csv')  ... 
                             };

gnss_inatm200stn_raw_data_paths = { ... 
                                    fullfile('2018-10-24_Session01','01_raw','ImuGnss_iNatM200Stn','inatm200stn_gnss_01.csv'); ... 
                                    fullfile('2018-10-24_Session02','01_raw','ImuGnss_iNatM200Stn','inatm200stn_gnss_02.csv'); ... 
                                    fullfile('2018-10-24_Session03','01_raw','ImuGnss_iNatM200Stn','inatm200stn_gnss_03.csv')  ... 
                                  };
                         
imu_inatm200stn_raw_data_paths = { ... 
                                   fullfile('2018-10-24_Session01','01_raw','ImuGnss_iNatM200Stn','inatm200stn_imu_01.csv'); ... 
                                   fullfile('2018-10-24_Session02','01_raw','ImuGnss_iNatM200Stn','inatm200stn_imu_02.csv'); ... 
                                   fullfile('2018-10-24_Session03','01_raw','ImuGnss_iNatM200Stn','inatm200stn_imu_03.csv')  ... 
                                 };
                         
% ref_inatm200stn_raw_data_paths = { ... 
%                                    fullfile('2018-10-24_Session01','01_raw','ImuGnss_iNatM200Stn','inatm200stn_internal_ekf_01.csv'); ... 
%                                    fullfile('2018-10-24_Session02','01_raw','ImuGnss_iNatM200Stn','inatm200stn_internal_ekf_02.csv'); ... 
%                                    fullfile('2018-10-24_Session03','01_raw','ImuGnss_iNatM200Stn','inatm200stn_internal_ekf_03.csv') ... 
%                                  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thales (GNSS) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(gnss_thales_raw_data_paths,gnss_thales_raw_data_root_string,'DataType','raw');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_thales_raw_data_root_string,'_session',sprintf('%02i',i)],gnss_thales_raw_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_thales_raw_data_root_string,out_data); 
clear out_data
                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iNat M200 STN (GNSS) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(gnss_inatm200stn_raw_data_paths,gnss_inatm200stn_raw_data_root_string,'DataType','raw');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_inatm200stn_raw_data_root_string,'_session',sprintf('%02i',i)],gnss_inatm200stn_raw_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_inatm200stn_raw_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iNat M200 STN (IMU) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(imu_inatm200stn_raw_data_paths,imu_inatm200stn_raw_data_root_string,'DataType','raw');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[imu_inatm200stn_raw_data_root_string,'_session',sprintf('%02i',i)],imu_inatm200stn_raw_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',imu_inatm200stn_raw_data_root_string,out_data); 
clear out_data

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % iNat M200 STN (reference data) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [out_data, num_files, load_flag] = importFromFile(ref_inatm200stn_raw_data_paths,ref_inatm200stn_raw_data_root_string,'DataType','raw');
% for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
%     exportToFile(out_data(i,:),[ref_inatm200stn_raw_data_root_string,'_session',sprintf('%02i',i)],ref_inatm200stn_raw_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
% end % for i
% assignin('base',ref_inatm200stn_raw_data_root_string,out_data); 
% clear out_data

%% Finish script

import_raw_data_executed = true;
