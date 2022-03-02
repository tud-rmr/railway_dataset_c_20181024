% Plot processed output data of simulink EKF Fusion model
%
%   Other m-files required: 
%       * setDatasetPaths
%       * importReferenceData
%       * importProcessedData
%       * getLatLonErrorEllipsePoints
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, importReferenceData, 
%             importProcessedData, processEkfFusionData

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 24-Nov-2020

%% Init

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

%% EKF positioning

if(1)
    session_selector = [1 2 3];
    
    plot_error_ellipses = 1; % enable or disable error-ellipse plots
    ts_error_ellipses_inatm200stn = 1; % distance between error-ellipse plots in sec
    max_major_error_ellipses = 10; % maximum major semi-axis size of error ellipses to plot
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    importReferenceData
    importMapData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        
        min_lon = min(gnss_inatm200stn_processed_data{session_i}.Longitude_deg);
        max_lon = max(gnss_inatm200stn_processed_data{session_i}.Longitude_deg);
        min_lat = min(gnss_inatm200stn_processed_data{session_i}.Latitude_deg);
        max_lat = max(gnss_inatm200stn_processed_data{session_i}.Latitude_deg);        
        ref_track_map_selector = ( (ref_track_map.Longitude_deg > min_lon) ...
                                    & ... 
                                   (ref_track_map.Longitude_deg < max_lon)) ... 
                                 & ...      
                                 ( (ref_track_map.Latitude_deg > min_lat) ... 
                                    & ...
                                   (ref_track_map.Latitude_deg < max_lat));
        plot_ref_track_map = any(ref_track_map_selector);
        
        if plot_error_ellipses
            error_ellipses_inatm200stn_indices = 1:ceil(ts_error_ellipses_inatm200stn/seconds(mean(diff(ref_ekfFusion_inatm200stn_data{session_i}.Time)))):length(ref_ekfFusion_inatm200stn_data{session_i}.Time);
            error_ellipses_inatm200stn_indices = error_ellipses_inatm200stn_indices(ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMajor_m(error_ellipses_inatm200stn_indices)' < max_major_error_ellipses);
            [~,ref_ekfFusion_inatm200stn_error_ellipse] = ... 
                getLatLonErrorEllipsePoints( ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}{error_ellipses_inatm200stn_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMajor_m(error_ellipses_inatm200stn_indices)*2.45, ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseMinor_m(error_ellipses_inatm200stn_indices)*2.45, ... 
                                             ref_ekfFusion_inatm200stn_data{session_i}.ErrorEllipseOrientation_deg(error_ellipses_inatm200stn_indices), ... 
                                             360 ... 
                                           );
        end % if
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Postions (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name);
        clear h_plot    
        h_plot = gobjects(0);
        
        if verLessThan('matlab', '9.5')
                                    
            hold on; grid on;
            
            h_plot(end+1) = plot(gnss_thales_processed_data{session_i}.Longitude_deg,gnss_thales_processed_data{session_i}.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales GNSS Receiver');
            %h_plot(end+1) = plot(ref_inatm200stn_internal_ekf_data{session_i}.Longitude_deg,ref_inatm200stn_internal_ekf_data{session_i}.Latitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
            h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_data{session_i}.Longitude_deg,ref_ekfFusion_inatm200stn_data{session_i}.Latitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
            if plot_ref_track_map
                h_plot(end+1) = plot(ref_track_map.Longitude_deg(ref_track_map_selector),ref_track_map.Latitude_deg(ref_track_map_selector),'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
            end % if
            if plot_error_ellipses
                h_plot(end+1) = plot(ref_ekfFusion_inatm200stn_error_ellipse(2,:),ref_ekfFusion_inatm200stn_error_ellipse(1,:),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%)');
            end % if
            
            % You can try
            % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
            % to create satellite plots with Google Maps API
            if exist('plot_google_map')            
                plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
            end % if            
            
            h_legend = legend(h_plot);
            set(h_legend,'Location','best')
            xlabel('longitude [deg]')
            ylabel('latitude [deg]')
            
        else
            
            h_plot(end+1) = geoplot(gnss_thales_processed_data{session_i}.Latitude_deg,gnss_thales_processed_data{session_i}.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales GNSS Receiver'); hold on;
            %h_plot(end+1) = geoplot(ref_inatm200stn_internal_ekf_data{session_i}.Latitude_deg,ref_inatm200stn_internal_ekf_data{session_i}.Longitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
            h_plot(end+1) = geoplot(ref_ekfFusion_inatm200stn_data{session_i}.Latitude_deg,ref_ekfFusion_inatm200stn_data{session_i}.Longitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');  
            if plot_ref_track_map
                h_plot(end+1) = geoplot(ref_track_map.Latitude_deg(ref_track_map_selector),ref_track_map.Longitude_deg(ref_track_map_selector),'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
            end % if
            if plot_error_ellipses
                h_plot(end+1) = geoplot(ref_ekfFusion_inatm200stn_error_ellipse(1,:),ref_ekfFusion_inatm200stn_error_ellipse(2,:),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF 2.45\sigma-Error-Ellipse (\approx95%)');
            end % if
            h_legend = legend(h_plot);
            set(h_legend,'Location','best')
            geobasemap('topographic'); % 'satellite','topographic','none'
           
        end        
        
    end % for session_i
end % if

%% EKF motion data

if(1)
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importReferenceData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
%         imu_inatm200stn_utc_time = datetime(seconds(imu_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
%         ref_inatm200stn_utc_time = datetime(seconds(ref_inatm200stn_internal_ekf_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_utc_time = datetime(seconds(ref_ekfFusion_inatm200stn_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Motion data (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
%         h_plot(end+1) = errorbar(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.DistanceVehicle_m,ref_ekfFusion_inatm200stn_data{session_i}.DistanceVehicleSigma_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.DistanceTrack_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('d_{track} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
%         h_plot(end+1) = errorbar(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.VelocityVehicle_ms,ref_ekfFusion_inatm200stn_data{session_i}.VelocityVehicleSigma_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.VelocityVehicle_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('v_{vehicle} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
%         h_plot(end+1) = plot(ref_inatm200stn_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.AccX_mss,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.AccX_mss,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('a_{vehicle} [m/s^2]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if

%% Distance and Speed Data

if(1)  
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    importReferenceData        
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        
        % Plot-Calculations _______________________________________________
        
        gnss_inatm200stn_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
%         ref_inatm200stn_utc_time = datetime(seconds(ref_inatm200stn_internal_ekf_data{session_i}.Time),'ConvertFrom','posixtime');
        ref_ekfFusion_utc_time = datetime(seconds(ref_ekfFusion_inatm200stn_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time = datetime(seconds(gnss_thales_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time_selector = isbetween(gnss_thales_utc_time,ref_ekfFusion_utc_time(1),ref_ekfFusion_utc_time(end));
        gnss_thales_utc_time = gnss_thales_utc_time(gnss_thales_utc_time_selector);
        
        % Plot ____________________________________________________________
        figure_name = ['EKF Vehicle Distances and Speeds (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.DistanceVehicle_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('d_{vehicle} [m]')

        clear h_plot    
        h_plot = gobjects(0);   
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.VelocityGround_ms(gnss_thales_utc_time_selector),'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales GNSS Receiver');
        h_plot(end+1) = plot(gnss_inatm200stn_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityGround_ms,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn GNSS Receiver');
%         h_plot(end+1) = plot(ref_inatm200stn_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.VelocityGround_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.VelocityGround_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.Heading_deg(gnss_thales_utc_time_selector),'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales GNSS Receiver');
        h_plot(end+1) = plot(gnss_inatm200stn_utc_time,gnss_inatm200stn_processed_data{session_i}.Heading_deg,'--','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn GNSS Receiver');
%         h_plot(end+1) = plot(ref_inatm200stn_utc_time,ref_inatm200stn_internal_ekf_data{session_i}.Heading_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn (internal EKF)');
        h_plot(end+1) = plot(ref_ekfFusion_utc_time,ref_ekfFusion_inatm200stn_data{session_i}.Heading_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_legend = legend(h_plot);
        set(h_legend,'Location','best')
        xlabel('UTC time')
        ylabel('heading angle [deg]')

        linkaxes([ax1,ax2,ax3],'x');
        
    end % for session_i
end % if
