% Process raw data to allow for an unified usage
%
%   Other m-files required: 
%       - setDatasetPaths
%       - loadFilenamePatterns
%       - importRawData
%       - exportToFile
%   MAT-files required: none
%
%   See also: importRawData

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 18-Nov-2020

%% Settings (edit)



% Select data to be processed by session number ___________________________
%   If the session selector array is empty, all sessions are imported which are 
%   available from the raw data.

gnss_thales_session_selector = [];
gnss_inatm200stn_session_selector = [];
imu_inatm200stn_session_selector = [];
ref_inatm200stn_session_selector = [];



% Time settings ___________________________________________________________

% Date of measurement campaign
date_of_measurement = '24-10-2018'; % format: dd-MM-yyyy

% Leap seconds between GPS time and UTC time at date of measurement campaign
leap_seconds_gps_utc = 18; % in sec

%% Initialization

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

% Prepare calculatons _____________________________________________________

loadFilenamePatterns
importParameters
importRawData

date_of_measurement = datetime(date_of_measurement,'InputFormat','dd-MM-yyyy');
date_of_measurement.TimeZone = 'utc';

epoch_date = date_of_measurement;
cal = calendar(epoch_date.Year,epoch_date.Month);
epoch_day = cal((cal(:,1)~=0) & (cal(:,1)<epoch_date.Day));
epoch_day = epoch_day(end);
epoch_date.Day = epoch_day;

imu_inatm200stn_attitude = [imu_inatm200stn_parameters.MountingRoll_deg,imu_inatm200stn_parameters.MountingPitch_deg,imu_inatm200stn_parameters.MountingYaw_deg];

R_x = @(roll_deg) [1 0 0; 0 cosd(roll_deg) sind(roll_deg); 0 -sind(roll_deg) cosd(roll_deg)];
R_y = @(pitch_deg) [cosd(pitch_deg) 0 -sind(pitch_deg); 0 1 0; sind(pitch_deg) 0 cosd(pitch_deg)];
R_z = @(yaw_deg) [cosd(yaw_deg) sind(yaw_deg) 0; -sind(yaw_deg) cosd(yaw_deg) 0; 0 0 1];

rotation_zyx = @(roll_pitch_yaw_attitude_deg) ... 
                 R_x(roll_pitch_yaw_attitude_deg(1)) * ... 
                 R_y(roll_pitch_yaw_attitude_deg(2)) * ... 
                 R_z(roll_pitch_yaw_attitude_deg(3));
            
processed_data_output_varnames = {};
processed_data_output_filenames = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GNSS: Thales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data = eval(gnss_thales_raw_data_root_string);
sessions = find(cellfun(@(cell) ~isempty(cell),raw_data));
input_varname = gnss_thales_raw_data_root_string;
for session_i = sessions(:)'   
    
    if ~isempty(gnss_thales_session_selector) && ~ismember(session_i,gnss_thales_session_selector)
        continue;
    end % if    
    
    raw_data_i = raw_data{session_i,1};
    
    % _____________________________________________________________________
    
    fprintf(['Processing ''',input_varname,'.mat'' raw data of session ',sprintf('%02i',session_i),'...']);

    time_utc = datetime(raw_data_i.TOW,'ConvertFrom','epochtime','Epoch',epoch_date) - seconds(leap_seconds_gps_utc);
    time_utc.TimeZone = 'utc';
    time_utc.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS';
    time_unix = posixtime(time_utc); % Convert to unix time

    data_selector = (yyyymmdd(time_utc) == yyyymmdd(date_of_measurement));
    unique_time_selector = [true;(diff(time_unix) ~= 0)];
    data_selector = data_selector & unique_time_selector;

    v_ground = sqrt(raw_data_i.vn_m_s.^2+raw_data_i.ve_m_s.^2);
    sigma_v_ground = nan(size(v_ground));
    heading = mod(atan2d(raw_data_i.ve_m_s,raw_data_i.vn_m_s),360);
    sigma_heading = nan(size(v_ground));
   
    vn_std = nan(size(v_ground));
    ve_std = nan(size(v_ground));
    vd_std = nan(size(v_ground));
    
    ellip_altitude = raw_data_i.height_m;
    sigma_ellip_altitude = nan(size(ellip_altitude));
    
    lat_std = nan(size(v_ground));
    lon_std = nan(size(v_ground));
    height_std = nan(size(v_ground));
    
    hdop = nan(size(v_ground));
    vdop = nan(size(v_ground));
    pdop = nan(size(v_ground));
    
    std_major = nan(size(v_ground));
    std_minor = nan(size(v_ground));
    alpha = nan(size(v_ground));

    data_tmp{session_i,1} = timetable( ...
                                       seconds(time_unix(data_selector)), ...
                                       raw_data_i.lat_rad(data_selector,:)*180/pi, ...
                                       raw_data_i.lon_rad(data_selector,:)*180/pi, ...
                                       ellip_altitude(data_selector,:), ... 
                                       raw_data_i.undulation_m(data_selector,:), ...
                                       v_ground(data_selector,:), ... 
                                       heading(data_selector,:), ... 
                                       raw_data_i.vn_m_s(data_selector,:), ...
                                       raw_data_i.ve_m_s(data_selector,:), ...
                                       -raw_data_i.vu_m_s(data_selector,:), ...
                                       lat_std(data_selector,:), ...
                                       lon_std(data_selector,:), ...
                                       height_std(data_selector,:), ...  
                                       sigma_v_ground(data_selector,:), ...  
                                       sigma_heading(data_selector,:), ... 
                                       vn_std(data_selector,:), ...
                                       ve_std(data_selector,:), ...
                                       vd_std(data_selector,:), ...
                                       std_major(data_selector,:), ...
                                       std_minor(data_selector,:), ...
                                       alpha(data_selector,:), ...
                                       hdop(data_selector,:), ...
                                       vdop(data_selector,:), ...
                                       pdop(data_selector,:) ...
                                     );
    data_tmp{session_i,1}.Properties.VariableNames = { ...
                                                       'Latitude_deg', ... 
                                                       'Longitude_deg', ... 
                                                       'AltitudeEllipsoid_m', ... 
                                                       'Undulation_m', ... 
                                                       'VelocityGround_ms', ... 
                                                       'Heading_deg', ...                                                                        
                                                       'VelocityNorth_ms', ... 
                                                       'VelocityEast_ms', ... 
                                                       'VelocityDown_ms', ... 
                                                       'LatitudeSigma_m', ... 
                                                       'LongitudeSigma_m', ... 
                                                       'AltitudeEllipsoidSigma_m', ... 
                                                       'VelocityGroundSigma_ms', ... 
                                                       'HeadingSigma_deg', ... 
                                                       'VelocityNorthSigma_ms', ... 
                                                       'VelocityEastSigma_ms', ... 
                                                       'VelocityDownSigma_ms', ... 
                                                       'ErrorEllipseMajor_m', ...
                                                       'ErrorEllipseMinor_m', ...
                                                       'ErrorEllipseOrientation_deg', ...
                                                       'HDOP', ... 
                                                       'VDOP', ... 
                                                       'PDOP' ... 
                                                     };

    % Finish ______________________________________________________________
        
    fprintf('done!\n');
    
end % for session_i


if exist('data_tmp','var')
    assignin('base',gnss_thales_processed_data_root_string,data_tmp);
else
    warning('processRawData: Thales GNSS data: ''data_tmp'' doesn''t exist! Probably the selected sessions are not available in the raw data?')
end % if

clear raw_data sessions input_varname session_i raw_data_i data_tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GNSS: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data = eval(gnss_inatm200stn_raw_data_root_string);
sessions = find(cellfun(@(cell) ~isempty(cell),raw_data));
input_varname = gnss_inatm200stn_raw_data_root_string;
for session_i = sessions(:)'   
    
    if ~isempty(gnss_inatm200stn_session_selector) && ~ismember(session_i,gnss_inatm200stn_session_selector)
        continue;
    end % if    
    
    raw_data_i = raw_data{session_i,1};
    
    % _____________________________________________________________________
    
    fprintf(['Processing ''',input_varname,'.mat'' raw data of session ',sprintf('%02i',session_i),'...']);

    time_utc = datetime(raw_data_i.GPSSeconds,'ConvertFrom','epochtime','Epoch',epoch_date) - seconds(leap_seconds_gps_utc);
    time_utc.TimeZone = 'utc';
    time_utc.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS';
    time_unix = posixtime(time_utc); % Convert to unix time

    data_selector = (yyyymmdd(time_utc) == yyyymmdd(date_of_measurement));
    data_selector = data_selector & ~ismember(raw_data_i.Position_VelocityType,'NONE');
    unique_time_selector = [true;(diff(time_unix) ~= 0)];
    data_selector = data_selector & unique_time_selector;

    v_ground = sqrt(raw_data_i.VelocityNorth_m_s_.^2+raw_data_i.VelocityEast_m_s_.^2);
    heading = mod(atan2d(raw_data_i.VelocityEast_m_s_,raw_data_i.VelocityNorth_m_s_),360);

    % Start: Calculation of standard deviation of 'v_ground' and 'heading' %%%%
    n_times = length(time_unix);
    P_v = zeros(2,2,n_times);
    P_v(1,1,:) = raw_data_i.VelocityNorthStandardDeviation_m_s_.^2;
    P_v(2,2,:) = raw_data_i.VelocityEastStandardDeviation_m_s_.^2;

    % H_v_ground = d(sqrt(v_n^2+v_e^2))/d(v_n,v_e)
    H_v_ground = zeros(1,2,n_times);
    H_v_ground(1,1,:) = raw_data_i.VelocityNorth_m_s_./v_ground;
    H_v_ground(1,2,:) = raw_data_i.VelocityEast_m_s_./v_ground;

    % H_heading = d(atan2(v_e,v_n)/d(v_n,v_e)
    H_heading = zeros(1,2,n_times);
    H_heading(1,1,:) = raw_data_i.VelocityNorth_m_s_./v_ground.^2;
    H_heading(1,2,:) = -raw_data_i.VelocityEast_m_s_./v_ground.^2;

    sigma_v_ground = zeros(size(v_ground)); 
    sigma_heading = zeros(size(heading));
    for time_i = 1:n_times
        sigma_v_ground(time_i) = sqrt(H_v_ground(:,:,time_i)*P_v(:,:,time_i)*H_v_ground(:,:,time_i)');
        sigma_heading(time_i) = sqrt(H_heading(:,:,time_i)*P_v(:,:,time_i)*H_heading(:,:,time_i)');
    end % for time_i
    % End: Calculation of standard deviation of 'v_ground' and 'heading' %%%%%%

    ellip_altitude = raw_data_i.Altitude_m_+raw_data_i.Undulation_m_;
    sigma_ellip_altitude = raw_data_i.AltitudeStandardDeviation_m_;

    hdop = nan(size(v_ground));
    vdop = nan(size(v_ground));
    
    std_major = nan(size(v_ground));
    std_minor = nan(size(v_ground));
    alpha = nan(size(v_ground));

    data_tmp{session_i,1} = timetable( ...
                                       seconds(time_unix(data_selector)), ...
                                       raw_data_i.Latitude_deg_(data_selector,:), ...
                                       raw_data_i.Longitude_deg_(data_selector,:), ...
                                       ellip_altitude(data_selector,:), ... 
                                       raw_data_i.Undulation_m_(data_selector,:), ...
                                       v_ground(data_selector,:), ... 
                                       heading(data_selector,:), ... 
                                       raw_data_i.VelocityNorth_m_s_(data_selector,:), ...
                                       raw_data_i.VelocityEast_m_s_(data_selector,:), ...
                                       raw_data_i.VelocityDown_m_s_(data_selector,:), ...
                                       raw_data_i.LatitudeStandardDeviation_m_(data_selector,:), ...
                                       raw_data_i.LongitudeStandardDeviation_m_(data_selector,:), ...
                                       sigma_ellip_altitude(data_selector,:), ...  
                                       sigma_v_ground(data_selector,:), ...  
                                       sigma_heading(data_selector,:), ... 
                                       raw_data_i.VelocityNorthStandardDeviation_m_s_(data_selector,:), ...
                                       raw_data_i.VelocityEastStandardDeviation_m_s_(data_selector,:), ...
                                       raw_data_i.VelocityDownStandardDeviation_m_s_(data_selector,:), ...
                                       std_major(data_selector,:), ...
                                       std_minor(data_selector,:), ...
                                       alpha(data_selector,:), ...
                                       hdop(data_selector,:), ...
                                       vdop(data_selector,:), ...
                                       raw_data_i.PositionDilutionOfPrecision(data_selector,:), ...
                                       raw_data_i.SatellitesUsed(data_selector,:), ...
                                       raw_data_i.SatellitesTracked(data_selector,:), ...
                                       raw_data_i.DifferentialAge_s_(data_selector,:), ...
                                       raw_data_i.SolutionAge_s_(data_selector,:), ...
                                       raw_data_i.SolutionStatus(data_selector,:), ...
                                       raw_data_i.Position_VelocityType(data_selector,:) ...
                                     );
    data_tmp{session_i,1}.Properties.VariableNames = { ...
                                                       'Latitude_deg', ... 
                                                       'Longitude_deg', ... 
                                                       'AltitudeEllipsoid_m', ... 
                                                       'Undulation_m', ... 
                                                       'VelocityGround_ms', ... 
                                                       'Heading_deg', ...                                                                        
                                                       'VelocityNorth_ms', ... 
                                                       'VelocityEast_ms', ... 
                                                       'VelocityDown_ms', ... 
                                                       'LatitudeSigma_m', ... 
                                                       'LongitudeSigma_m', ... 
                                                       'AltitudeEllipsoidSigma_m', ... 
                                                       'VelocityGroundSigma_ms', ... 
                                                       'HeadingSigma_deg', ... 
                                                       'VelocityNorthSigma_ms', ... 
                                                       'VelocityEastSigma_ms', ... 
                                                       'VelocityDownSigma_ms', ... 
                                                       'ErrorEllipseMajor_m', ...
                                                       'ErrorEllipseMinor_m', ...
                                                       'ErrorEllipseOrientation_deg', ...
                                                       'HDOP', ... 
                                                       'VDOP', ... 
                                                       'PDOP', ... 
                                                       'SatellitesUsed', ... 
                                                       'SatellitesTracked', ... 
                                                       'DifferentialAge_s', ... 
                                                       'SolutionAge_s', ... 
                                                       'SolutionStatus', ... 
                                                       'PositionVelocityType' ...
                                                     };

    % Fill GPS outages with NaN
    gps_outage_selector = ~ismember(raw_data_i.SolutionStatus(data_selector,:),'SOL_COMPUTED');

    data_tmp{session_i,1}.Latitude_deg(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.Longitude_deg(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.AltitudeEllipsoid_m(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.Undulation_m(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityGround_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.Heading_deg(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityNorth_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityEast_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityDown_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.LatitudeSigma_m(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.LongitudeSigma_m(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.AltitudeEllipsoidSigma_m(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityGroundSigma_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.HeadingSigma_deg(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityNorthSigma_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityEastSigma_ms(gps_outage_selector,:) = nan;
    data_tmp{session_i,1}.VelocityDownSigma_ms(gps_outage_selector,:) = nan;


    % Finish ______________________________________________________________
        
    fprintf('done!\n');
    
end % for session_i

if exist('data_tmp','var')
    assignin('base',gnss_inatm200stn_processed_data_root_string,data_tmp);
else
    warning('processRawData: iNat-M200/STN GNSS data: ''data_tmp'' doesn''t exist! Probably the selected sessions are not available in the raw data?')
end % if

clear raw_data sessions input_varname session_i raw_data_i data_tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMU: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data = eval(imu_inatm200stn_raw_data_root_string);
sessions = find(cellfun(@(cell) ~isempty(cell),raw_data));
input_varname = imu_inatm200stn_raw_data_root_string;
for session_i = sessions(:)'   
    
    if ~isempty(imu_inatm200stn_session_selector) && ~ismember(session_i,imu_inatm200stn_session_selector)
        continue;
    end % if    
    
    raw_data_i = raw_data{session_i,1};
    
    % _____________________________________________________________________
    
    fprintf(['Processing ''',input_varname,'.mat'' raw data of session ',sprintf('%02i',session_i),'...']);
    
    time_utc = datetime(raw_data_i.GPSSeconds,'ConvertFrom','epochtime','Epoch',epoch_date) - seconds(leap_seconds_gps_utc);
    time_utc.TimeZone = 'utc';
    time_utc.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS';
    time_unix = posixtime(time_utc); % Convert to unix time

    data_selector = (yyyymmdd(time_utc) == yyyymmdd(date_of_measurement));
    unique_time_selector = [true;(diff(time_unix) ~= 0)];
    data_selector = data_selector & unique_time_selector;

    acc = (rotation_zyx(imu_inatm200stn_attitude) * ... 
            [raw_data_i.AccelerationX_m_s_s_, ... 
             raw_data_i.AccelerationY_m_s_s_, ... 
             raw_data_i.AccelerationZ_m_s_s_]')';
    turn_rates = (rotation_zyx(imu_inatm200stn_attitude) * ... 
                    [raw_data_i.AngularRateX_deg_s_, ... 
                     raw_data_i.AngularRateY_deg_s_, ... 
                     raw_data_i.AngularRateZ_deg_s_]')';

    data_tmp{session_i,1} = timetable( ...
                                       seconds(time_unix(data_selector)), ...
                                       acc(data_selector,1), ...
                                       acc(data_selector,2), ...
                                       acc(data_selector,3), ...
                                       turn_rates(data_selector,1), ...
                                       turn_rates(data_selector,2), ...
                                       turn_rates(data_selector,3) ...  
                                     );
    data_tmp{session_i,1}.Properties.VariableNames = { ...
                                                       'AccX_mss', ...
                                                       'AccY_mss', ...
                                                       'AccZ_mss', ...
                                                       'TurnRateX_degs', ...
                                                       'TurnRateY_degs', ...
                                                       'TurnRateZ_degs' ...
                                                     };

    % Finish ______________________________________________________________
        
    fprintf('done!\n');
    
end % for session_i

if exist('data_tmp','var')
    assignin('base',imu_inatm200stn_processed_data_root_string,data_tmp);
else
    warning('processRawData: iNat-M200/STN IMU data: ''data_tmp'' doesn''t exist! Probably the selected sessions are not available in the raw data?')
end % if

clear raw_data sessions input_varname session_i raw_data_i data_tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference: iNat M200 STN (internal EKF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% raw_data = eval(ref_inatm200stn_raw_data_root_string);
% sessions = find(cellfun(@(cell) ~isempty(cell),raw_data));
% input_varname = ref_inatm200stn_raw_data_root_string;
% for session_i = sessions(:)'   
%     
%     if ~isempty(ref_inatm200stn_session_selector) && ~ismember(session_i,ref_inatm200stn_session_selector)
%         continue;
%     end % if    
%     
%     raw_data_i = raw_data{session_i,1};
%     
%     % _____________________________________________________________________
%     
%     fprintf(['Processing ''',input_varname,'.mat'' raw data of session ',sprintf('%02i',session_i),'...']);
% 
%     time_utc = datetime(raw_data_i.GPSSeconds,'ConvertFrom','epochtime','Epoch',epoch_date) - seconds(leap_seconds_gps_utc);
%     time_utc.TimeZone = 'utc';
%     time_utc.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS';
%     time_unix = posixtime(time_utc); % Convert to unix time
% 
%     data_selector = (yyyymmdd(time_utc) == yyyymmdd(date_of_measurement));
%     data_selector = data_selector & ~cellfun(@isempty,regexp(raw_data_i.SystemStatus,'0x4100*'));
%     unique_time_selector = [true;(diff(time_unix) ~= 0)];
%     data_selector = data_selector & unique_time_selector;
% 
%     v_ground = sqrt(raw_data_i.VelocityNorth_m_s_.^2+raw_data_i.VelocityEast_m_s_.^2);
%     sigma_v_ground = nan(size(v_ground)); 
%     
%     sigma_latitude = nan(size(v_ground));
%     sigma_longitude = nan(size(v_ground));
%     d_vehicle = nan(size(v_ground));
%     sigma_d_vehicle = nan(size(v_ground));
%     d_track = nan(size(v_ground));
%     sigma_d_track = nan(size(v_ground));
%     v_vehicle = nan(size(v_ground));
%     sigma_v_vehicle = nan(size(v_ground));
%     sigma_ellip_altitude = nan(size(v_ground));
%     sigma_v_north = nan(size(v_ground));
%     sigma_v_east = nan(size(v_ground));
%     sigma_v_down = nan(size(v_ground));
% 
%     % Start: Transform attitude to standard vehicle coordinates %%%%%%%%%%%%%%%
%     eulerToQuatXYZ = @(r,p,y) [ ...
%                                 cos(r/2).*cos(p/2).*cos(y/2) - sin(r/2).*sin(p/2).*sin(y/2), ...
%                                 sin(r/2).*cos(p/2).*cos(y/2) + cos(r/2).*sin(p/2).*sin(y/2), ...
%                                 -sin(r/2).*cos(p/2).*sin(y/2) + cos(r/2).*sin(p/2).*cos(y/2), ...
%                                 cos(r/2).*cos(p/2).*sin(y/2) + sin(r/2).*sin(p/2).*cos(y/2)  ... 
%                             ];
%     quatToEulerXYZ = @(q) [ ... 
%                             atan2( -2*(q(3,:).*q(4,:) - q(2,:).*q(1,:)) , q(1,:).^2 - q(2,:).^2 - q(3,:).^2 + q(4,:).^2 ); ...
%                                                                asin( max(min(2*(q(2,:).*q(4,:) + q(3,:).*q(1,:)),1),-1) ); ...
%                             atan2( -2*(q(2,:).*q(3,:) - q(4,:).*q(1,:)) , q(1,:).^2 + q(2,:).^2 - q(3,:).^2 - q(4,:).^2 )  ... 
%                           ];
% 
%     orientation_vec = [raw_data_i.Roll_deg_,raw_data_i.Pitch_deg_,raw_data_i.Yaw_deg_] * pi/180;
%     orientation_q = eulerToQuatXYZ(orientation_vec(:,1),orientation_vec(:,2),orientation_vec(:,3));
%     correction_q = [1/sqrt(2),0,0,-1/sqrt(2)];
%     q_new = [ ... 
%               correction_q(1) -correction_q(2) -correction_q(3) -correction_q(4); ... 
%               correction_q(2)  correction_q(1) -correction_q(4)  correction_q(3); ... 
%               correction_q(3)  correction_q(4)  correction_q(1) -correction_q(2); ... 
%               correction_q(4) -correction_q(3)  correction_q(2)  correction_q(1)  ... 
%             ]*orientation_q';
%     orientation_new = quatToEulerXYZ(q_new)'*180/pi;                            
%     % End: Transform attitude to standard vehicle coordinates %%%%%%%%%%%%%%%%%
%     
%     heading = mod(atan2d(raw_data_i.VelocityEast_m_s_,raw_data_i.VelocityNorth_m_s_),360);
%     sigma_heading = nan(size(v_ground));
%     
%     sigma_roll = nan(size(v_ground));
%     sigma_pitch = nan(size(v_ground));
%     sigma_yaw = nan(size(v_ground));
% 
%     acc = (rotation_zyx(imu_inatm200stn_attitude) * ... 
%             [raw_data_i.AccelerationX_m_s_s_, ... 
%              raw_data_i.AccelerationY_m_s_s_, ... 
%              raw_data_i.AccelerationZ_m_s_s_]')';
%     turn_rates = (rotation_zyx(imu_inatm200stn_attitude) * ... 
%                     [raw_data_i.AngularRateX_deg_s_, ... 
%                      raw_data_i.AngularRateY_deg_s_, ... 
%                      raw_data_i.AngularRateZ_deg_s_]')';
% 
%                  
%     std_major = nan(size(v_ground));
%     std_minor = nan(size(v_ground));
%     alpha = nan(size(v_ground));
%                  
%     data_tmp{session_i,1} = timetable( ...
%                                        seconds(time_unix(data_selector,:)), ...
%                                        raw_data_i.Longitude_deg_(data_selector,:), ...
%                                        raw_data_i.Latitude_deg_(data_selector,:), ...
%                                        raw_data_i.Altitude_m_(data_selector,:), ...
%                                        d_vehicle(data_selector,:), ... 
%                                        d_track(data_selector,:), ... 
%                                        v_vehicle(data_selector,:), ... 
%                                        v_ground(data_selector,:), ... 
%                                        heading(data_selector,:), ... 
%                                        raw_data_i.VelocityNorth_m_s_(data_selector,:), ...
%                                        raw_data_i.VelocityEast_m_s_(data_selector,:), ...
%                                        raw_data_i.VelocityDown_m_s_(data_selector,:), ...
%                                        sigma_latitude(data_selector,:), ...
%                                        sigma_longitude(data_selector,:), ...
%                                        sigma_ellip_altitude(data_selector,:), ...
%                                        sigma_d_vehicle(data_selector,:), ...
%                                        sigma_d_track(data_selector,:), ... 
%                                        sigma_v_vehicle(data_selector,:), ... 
%                                        sigma_v_ground(data_selector,:), ...  
%                                        sigma_heading(data_selector,:), ... 
%                                        sigma_v_north(data_selector,:), ...
%                                        sigma_v_east(data_selector,:), ...
%                                        sigma_v_down(data_selector,:), ...
%                                        std_major(data_selector,:), ...
%                                        std_minor(data_selector,:), ...
%                                        alpha(data_selector,:), ...
%                                        acc(data_selector,1), ...
%                                        acc(data_selector,2), ...
%                                        acc(data_selector,3), ...
%                                        turn_rates(data_selector,1), ...
%                                        turn_rates(data_selector,2), ...
%                                        turn_rates(data_selector,3), ...  
%                                        orientation_new(data_selector,1), ... 
%                                        orientation_new(data_selector,2), ...                                                   
%                                        orientation_new(data_selector,3), ...
%                                        sigma_roll(data_selector,:), ... 
%                                        sigma_pitch(data_selector,:), ... 
%                                        sigma_yaw(data_selector,:) ... 
%                                      );
%     data_tmp{session_i,1}.Properties.VariableNames = { ...
%                                                        'Latitude_deg', ... 
%                                                        'Longitude_deg', ... 
%                                                        'AltitudeEllipsoid_m', ... 
%                                                        'DistanceVehicle_m', ... 
%                                                        'DistanceTrack_m', ... 
%                                                        'VelocityVehicle_ms', ... 
%                                                        'VelocityGround_ms', ... 
%                                                        'Heading_deg', ...                                                                        
%                                                        'VelocityNorth_ms', ... 
%                                                        'VelocityEast_ms', ... 
%                                                        'VelocityDown_ms', ... 
%                                                        'LatitudeSigma_m', ... 
%                                                        'LongitudeSigma_m', ... 
%                                                        'AltitudeEllipsoidSigma_m', ...
%                                                        'DistanceVehicleSigma_m', ... 
%                                                        'DistanceTrackSigma_m', ... 
%                                                        'VelocityVehicleSigma_ms', ... 
%                                                        'VelocityGroundSigma_ms', ... 
%                                                        'HeadingSigma_deg', ... 
%                                                        'VelocityNorthSigma_ms', ... 
%                                                        'VelocityEastSigma_ms', ... 
%                                                        'VelocityDownSigma_ms', ...
%                                                        'ErrorEllipseMajor_m', ...
%                                                        'ErrorEllipseMinor_m', ...
%                                                        'ErrorEllipseOrientation_deg', ...
%                                                        'AccX_mss', ...
%                                                        'AccY_mss', ...
%                                                        'AccZ_mss', ...
%                                                        'TurnRateX_degs', ...
%                                                        'TurnRateY_degs', ...
%                                                        'TurnRateZ_degs', ...
%                                                        'Roll_deg', ...
%                                                        'Pitch_deg', ...
%                                                        'Yaw_deg', ...
%                                                        'RollSigma_deg', ...
%                                                        'PitchSigma_deg', ...
%                                                        'YawSigma_deg' ...
%                                                      };                                                                 
% 
%                                                          
%     % Finish ______________________________________________________________
%         
%     fprintf('done!\n');
%     
% end % for session_i
% 
% if exist('data_tmp','var')
%     assignin('base',ref_inatm200stn_processed_data_root_string,data_tmp);
% else
%     warning('processRawData: iNat-M200/STN REF data: ''data_tmp'' doesn''t exist! Probably the selected sessions are not available in the raw data?')
% end % if
% 
% clear raw_data sessions input_varname session_i raw_data_i data_tmp


%% Synchronize data
    
gnss_thales_data = eval(gnss_thales_processed_data_root_string);
gnss_inatm200stn_data = eval(gnss_inatm200stn_processed_data_root_string);
imu_inatm200stn_data = eval(imu_inatm200stn_processed_data_root_string);
% ref_inatm200stn_data = eval(ref_inatm200stn_processed_data_root_string);

for session_i = 1:size(gnss_thales_data,1)

    fprintf(['Synchronizing data of session ',sprintf('%02i',session_i),'...\n']);

    % Set common start and end time ___________________________________

    % Get intersection time limits
    common_start_time_i = max([ ... 
                                gnss_thales_data{session_i,1}.Time(1), ...  
                                gnss_inatm200stn_data{session_i,1}.Time(1), ... 
                                imu_inatm200stn_data{session_i,1}.Time(1), ... 
                                %ref_inatm200stn_data{session_i,1}.Time(1) ... 
                              ]);
    common_end_time_i = min([ ... 
                              gnss_thales_data{session_i,1}.Time(end), ... 
                              gnss_inatm200stn_data{session_i,1}.Time(end), ... 
                              imu_inatm200stn_data{session_i,1}.Time(end), ... 
                              %ref_inatm200stn_data{session_i,1}.Time(end) ... 
                            ]);

    gnss_thales_data{session_i,1} = limitTt(gnss_thales_data{session_i,1},common_start_time_i,common_end_time_i);
    gnss_inatm200stn_data{session_i,1} = limitTt(gnss_inatm200stn_data{session_i,1},common_start_time_i,common_end_time_i);
    imu_inatm200stn_data{session_i,1} = limitTt(imu_inatm200stn_data{session_i,1},common_start_time_i,common_end_time_i);
    %ref_inatm200stn_data{session_i,1} = limitTt(ref_inatm200stn_data{session_i,1},common_start_time_i,common_end_time_i);

    % End synchronizing data ______________________________________________

    fprintf('done!\n');

end % for session_i  

assignin('base',gnss_thales_processed_data_root_string,gnss_thales_data);
assignin('base',gnss_inatm200stn_processed_data_root_string,gnss_inatm200stn_data);
assignin('base',imu_inatm200stn_processed_data_root_string,imu_inatm200stn_data);
%assignin('base',ref_inatm200stn_processed_data_root_string,ref_inatm200stn_data);

clear gnss_thales_data gnss_inatm200stn_data imu_inatm200stn_data ref_inatm200stn_data

%% Save data

if 1    

    exportToFile(eval(gnss_thales_processed_data_root_string),gnss_thales_processed_data_root_string,gnss_thales_processed_data_root_string,'ChunkSize',300);
    exportToFile(eval(gnss_inatm200stn_processed_data_root_string),gnss_inatm200stn_processed_data_root_string,gnss_inatm200stn_processed_data_root_string,'ChunkSize',300);
    exportToFile(eval(imu_inatm200stn_processed_data_root_string),imu_inatm200stn_processed_data_root_string,imu_inatm200stn_processed_data_root_string,'ChunkSize',300);
    %exportToFile(eval(ref_inatm200stn_processed_data_root_string),ref_inatm200stn_processed_data_root_string,ref_inatm200stn_processed_data_root_string,'ChunkSize',300);

end % if

%% Post-Processing

import_processed_data_executed = true;
