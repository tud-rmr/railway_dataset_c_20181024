% Plot processed raw data
%
%   Other m-files required: setDatasetPaths, importProcessedData
%   MAT-files required: none
%
%   See also: loadFilenamePatterns, importProcessedData

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 10-Dec-2020

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

%% Reference Track-Map

if(1)
    % Pre-Calculations ____________________________________________________
    importMapData
       
    % Plot Routine ________________________________________________________
    figure_name = 'Reference Track-Map';
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);
    
    if verLessThan('matlab', '9.5')
        
        clear h_plot    
        h_plot = gobjects(0);
        ax1 = subplot(4,1,1); hold on; grid on;
        set(ax1,'OuterPosition',[0,0.62, 1, .40]);    
        h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'-','Color',[0, 0.4470, 0.7410],'LineWidth',1.5,'MarkerSize',10,'DisplayName','Center Line');
        h_plot(end+1) = plot(ref_track_map.Longitude_RightRail_deg,ref_track_map.Latitude_RightRail_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Right Rail');
        h_plot(end+1) = plot(ref_track_map.Longitude_LeftRail_deg,ref_track_map.Latitude_LeftRail_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Left Rail');
        h_plot(end+1) = plot(ref_track_map.Longitude_deg(1,:),ref_track_map.Latitude_deg(1,:),'c.','LineWidth',1.5,'MarkerSize',15,'DisplayName',['Start (d=',num2str(ref_track_map.TrackDistance_m(1)),'m)']);
        h_plot(end+1) = plot(ref_track_map.Longitude_deg(end,:),ref_track_map.Latitude_deg(end,:),'c*','LineWidth',1.5,'MarkerSize',10,'DisplayName',['End (d=',num2str(ref_track_map.TrackDistance_m(end)),'m)']);

        % You can try
        % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
        % to create satellite plots with Google Maps API
        if exist('plot_google_map')
            plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
        end % if 
        h_legend = legend(h_plot);
        set(h_legend,'Location','northwest');
        xlabel('longitude [deg]')
        ylabel('latitude [deg]')
     
    else
        
        clear h_plot    
        h_plot = gobjects(0);
        ax1 = subplot(4,1,1);
        set(ax1,'OuterPosition',[0,0.62, 1, .40]);    
        h_plot(end+1) = geoplot(ref_track_map.Latitude_deg,ref_track_map.Longitude_deg,'-','Color',[0, 0.4470, 0.7410],'LineWidth',1.5,'MarkerSize',10,'DisplayName','Center Line'); hold on;
        h_plot(end+1) = geoplot(ref_track_map.Latitude_RightRail_deg,ref_track_map.Longitude_RightRail_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Right Rail');
        h_plot(end+1) = geoplot(ref_track_map.Latitude_LeftRail_deg,ref_track_map.Longitude_LeftRail_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Left Rail');
        h_plot(end+1) = geoplot(ref_track_map.Latitude_deg(1,:),ref_track_map.Longitude_deg(1,:),'k.','LineWidth',1.5,'MarkerSize',15,'DisplayName',['Start (d=',num2str(ref_track_map.TrackDistance_m(1)),'m)']);
        h_plot(end+1) = geoplot(ref_track_map.Latitude_deg(end,:),ref_track_map.Longitude_deg(end,:),'k*','LineWidth',1.5,'MarkerSize',10,'DisplayName',['End (d=',num2str(ref_track_map.TrackDistance_m(end)),'m)']);
        geobasemap('topographic');
        
        h_legend = legend(h_plot);
        set(h_legend,'Location','northwest');
        
    end % if
    
    clear h_plot    
    h_plot = gobjects(0);       
    ax2 = subplot(4,1,2); hold on; grid on;
    set(ax2,'OuterPosition',[0,0.41, 1, .15]);
    h_plot(end+1) = plot(ref_track_map.TrackDistance_m,ref_track_map.Roll_deg,'LineWidth',1.5);
    xlabel('d [m]')
    ylabel('roll [deg]')
    
    clear h_plot
    h_plot = gobjects(0);
    ax3 = subplot(4,1,3); hold on; grid on;
    set(ax3,'OuterPosition',[0,0.19, 1, .15]);
    h_plot(end+1) = plot(ref_track_map.TrackDistance_m,ref_track_map.Pitch_deg,'LineWidth',1.5);
    xlabel('d [m]')
    ylabel('pitch [deg]')
    
    clear h_plot
    h_plot = gobjects(0);
    ax4 = subplot(4,1,4); hold on; grid on;
    set(ax4,'OuterPosition',[0,0, 1, .15]);
    h_plot(end+1) = plot(ref_track_map.TrackDistance_m,ref_track_map.Heading_deg,'LineWidth',1.5);
    xlabel('d [m]')
    ylabel('heading [deg]')
    
    linkaxes([ax2,ax3,ax4],'x');
end % if

%% GNSS Positions

if(1)        
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    importMapData
    
    % Plot Routine ________________________________________________________    
    for session_i = session_selector(:)'
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time = datetime(seconds(gnss_thales_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
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
        
        % Plot ____________________________________________________________        
        figure_name = ['GNSS Lat-Lon-Positions (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        h_fig = figure('Name',figure_name);
        
        clear h_plot    
        h_plot = gobjects(0);
        
        if verLessThan('matlab', '9.5')
            
            hold on; grid on;
            h_plot(end+1) = plot(gnss_thales_processed_data{session_i}.Longitude_deg,gnss_thales_processed_data{session_i}.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
            h_plot(end+1) = plot(gnss_inatm200stn_processed_data{session_i}.Longitude_deg,gnss_inatm200stn_processed_data{session_i}.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');        
            %h_plot(end+1) = plot(gnss_thales_processed_data{session_i}.Longitude_deg([1 end]),gnss_thales_processed_data{session_i}.Latitude_deg([1 end]),'.','LineWidth',1.5,'MarkerSize',20,'DisplayName','Thales Receiver (End Points)');
            %h_plot(end+1) = plot(gnss_inatm200stn_processed_data{session_i}.Longitude_deg([1 end]),gnss_inatm200stn_processed_data{session_i}.Latitude_deg([1 end]),'.','LineWidth',1.5,'MarkerSize',20,'DisplayName','iNatM200Stn (End Points)');        
            if plot_ref_track_map
                h_plot(end+1) = plot(ref_track_map.Longitude_deg(ref_track_map_selector),ref_track_map.Latitude_deg(ref_track_map_selector),'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
            end % if
            
            % You can try
            % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
            % to create satellite plots with Google Maps API
            if exist('plot_google_map')
                plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
            end % if    
            
            legend(h_plot);
            xlabel('longitude [deg]')
            ylabel('latitude [deg]')
            
            dcm_obj = datacursormode(h_fig);
            set(dcm_obj,'UpdateFcn',{@myCustomDataTipFcn,{gnss_thales_utc_time,gnss_inat_utc_time}})
            
        else
            
            h_plot(end+1) = geoplot(gnss_thales_processed_data{session_i}.Latitude_deg,gnss_thales_processed_data{session_i}.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver'); hold on
            h_plot(end+1) = geoplot(gnss_inatm200stn_processed_data{session_i}.Latitude_deg,gnss_inatm200stn_processed_data{session_i}.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');        
            if plot_ref_track_map
                h_plot(end+1) = geoplot(ref_track_map.Latitude_deg(ref_track_map_selector),ref_track_map.Longitude_deg(ref_track_map_selector),'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
            end % if
            
            geobasemap('topographic');
            legend(h_plot);
            
            dcm_obj = datacursormode(h_fig);
            set(dcm_obj,'UpdateFcn',{@myCustomDataTipFcn,{gnss_thales_utc_time,gnss_inat_utc_time}})
        
        end % if
    
    end % for session_i
end % if

%% GNSS (X,Y,Z)-Plot

if(1)
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
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
                
        % Plot ____________________________________________________________
        figure_name = ['GNSS X-Y-Z-Positions (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold on; grid on;

        clear h_plot    
        h_plot = gobjects(0);         
        
        h_plot(end+1) = plot3(gnss_thales_processed_data{session_i}.Longitude_deg,gnss_thales_processed_data{session_i}.Latitude_deg,gnss_thales_processed_data{session_i}.AltitudeEllipsoid_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot3(gnss_inatm200stn_processed_data{session_i}.Longitude_deg,gnss_inatm200stn_processed_data{session_i}.Latitude_deg,gnss_inatm200stn_processed_data{session_i}.AltitudeEllipsoid_m,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        if plot_ref_track_map
            h_plot(end+1) = plot3(ref_track_map.Longitude_deg(ref_track_map_selector),ref_track_map.Latitude_deg(ref_track_map_selector),ref_track_map.AltitudeEllipsoid_m(ref_track_map_selector),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
        end % if
        
        legend(h_plot);
        xlabel('longitude [deg]')
        ylabel('latitude [deg]')
        zlabel('height_{ellip.} [m]')
        
        view(3);
    end % for session_i
end % if

%% GNSS Altitudes

if(1)
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    importMapData
    
    % Plot Routine ________________________________________________________    
    for session_i = session_selector(:)'        
        % Plot-Calculations _______________________________________________        
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time = datetime(seconds(gnss_thales_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Altitudes from Ellipsoid (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold on; grid on;

        clear h_plot    
        h_plot = gobjects(0);         
        
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.AltitudeEllipsoid_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.AltitudeEllipsoid_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        
        legend(h_plot);
        xlabel('UTC time')
        ylabel('height_{ellip.} [m]')
    
    end % for session_i
end % if

%% GNSS Speeds

if(1)
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData       
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time = datetime(seconds(gnss_thales_processed_data{session_i}.Time),'ConvertFrom','posixtime');
                
        % Plot ____________________________________________________________
        figure_name = ['GNSS NED-Velocities (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;    
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.VelocityNorth_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityNorth_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{north} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.VelocityEast_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityEast_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{east} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.VelocityDown_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityDown_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{down} [m/s]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% Ground Speeds (GNSS)

if(1)   
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'        
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        gnss_thales_utc_time = datetime(seconds(gnss_thales_processed_data{session_i}.Time),'ConvertFrom','posixtime');
                
        % Plot ____________________________________________________________
        figure_name = ['GNSS Speed over Ground (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(2,1,1); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(2,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_thales_utc_time,gnss_thales_processed_data{session_i}.Heading_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Thales Receiver');
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.Heading_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('heading [deg]')

        linkaxes([ax1,ax2],'x');
    end % for session_i
end % if

%% GNSS Position Uncertainties

if(1)    
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        % Plot-Calculations _______________________________________________
        gnss_inat_utc_time = datetime(seconds(gnss_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Position Uncertainties (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.LatitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('\sigma_{latitue} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.LongitudeSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('\sigma_{longitude} [m]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_inat_utc_time,gnss_inatm200stn_processed_data{session_i}.AltitudeEllipsoidSigma_m,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        h_legend = legend(h_plot);
        set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('\sigma_{height} [m]')

        linkaxes([ax1,ax2,ax3],'x');    
    end % for session_i
end % if

%% IMU Accelerations

if(1)    
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________
    for session_i = session_selector(:)'
        % Plot-Calculations _______________________________________________        
        imu_utc_time = datetime(seconds(imu_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Accelerations (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.AccX_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_x [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.AccY_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');    
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_y [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.AccZ_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');   
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_z [m/s^2]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% IMU Turn Rates

if(1)  
    session_selector = [1 2 3];
    
    % Pre-Calculations ____________________________________________________
    importProcessedData
    
    % Plot Routine ________________________________________________________        
    for session_i = session_selector(:)'
        % Plot-Calculations _______________________________________________
        imu_utc_time = datetime(seconds(imu_inatm200stn_processed_data{session_i}.Time),'ConvertFrom','posixtime');
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Turn Rates (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateX_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_x [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateY_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_y [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_utc_time,imu_inatm200stn_processed_data{session_i}.TurnRateZ_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','iNatM200Stn');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_z [deg/s]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% Helper Functions

function txt = myCustomDataTipFcn(pointDataTip,event_obj,utc_time)
% txt = myCustomDataTipFcn(pointDataTip,event_obj,time_utc)
%

data_index = pointDataTip.Cursor.DataIndex;
pos = get(event_obj,'Position');

switch event_obj.Target.DisplayName
    case {'Thales Receiver'}
        time_str = datestr(utc_time{1}(data_index));
    case {'iNatM200Stn'}
        time_str = datestr(utc_time{2}(data_index));
    otherwise
        time_str = 'not available';
end

txt = { ...
        ['X: ',num2str(pos(1))], ...
        ['Y: ',num2str(pos(2))],  ... 
        ['Time: ',time_str] ... 
      };

end 
