function [downstreamSpeed,upstreamSpeed,downDensityRatio,upDensityRatio] = calculate_updownSW(event_start,event_end,Event_number)
    %Calculates the solar wind beam speed by first plotting the distribution in Vx-Vy with magnitude
    %V and angle of phi from distribution function. Then find the max histcount, and fit a normal
    %distribution around the peak for the speed.
    %Plot?
    plotting = 1;
    
    
    window = 1500; %time window in milliseconds, total centered on center time
    specie = 'i';
    peak_angle = 3;
    peak_radius = 3; %6 works for all but Event 13 
    %The SW beam is narrow, therre is also heavier ions thatt appear faster than the ion SW due to
    %the instrument not being able to differentiate
    
    %Take the time interval before and after center time
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    eventStart = datenum(event_start,formatIn);
    eventEnd = datenum(event_end,formatIn);
    
    down_centerTime = eventStart + seconds(3);
    up_centerTime = eventEnd - seconds(3);
    
    [downstreamSpeed,downDensityRatio] = calculate_solarwindBeam(Event_number, down_centerTime,window,peak_angle,peak_radius,specie,plotting);
    [upstreamSpeed,upDensityRatio] = calculate_solarwindBeam(Event_number, up_centerTime,window,peak_angle,peak_radius,specie,plotting);
    
end
function[SolarwindBeamSpeed,densityRatio] = calculate_solarwindBeam(Event_number, centerTime,window,peak_angle,peak_radius,specie,plotting)
    timeLeft = centerTime - milliseconds(window);
    timeRight = centerTime + milliseconds(window);
    
    %Get the Distribution Data
    
        eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
        %eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
        eventDataDirectory = '/Users/andrewvu/data/Event Data/';
        eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
        
    %   if exist(eventDataDirectoryFileName,'file') ~= 2
    %     load(eventDataDirectoryFileName) %#ok<LOAD>
    
    
    
    %
    %         save(eventDataDirectoryFileName,...
    %             'mms1_mec_timedata_raw',...
    %             'mms1_mec_rdata_raw',...
    %             'mms2_mec_timedata_raw',...
    %             'mms2_mec_rdata_raw',...
    %             'mms3_mec_timedata_raw',...
    %             'mms3_mec_rdata_raw',...
    %             'mms4_mec_timedata_raw',...
    %             'mms4_mec_rdata_raw',...
    %             'mms1_fgm_timedata_raw',...
    %             'mms1_fgm_bdata_raw',...
    %             'mms2_fgm_timedata_raw',...
    %             'mms2_fgm_bdata_raw',...
    %             'mms3_fgm_timedata_raw',...
    %             'mms3_fgm_bdata_raw',...
    %             'mms4_fgm_timedata_raw',...
    %             'mms4_fgm_bdata_raw',...
    %             'mms1_fgm_timedata_srvy',...
    %             'mms1_fgm_bdata_srvy',...
    %             'fpi_e_timedata',...
    %             'fpi_e_ndata',...
    %             'fpi_e_vdata',...
    %             'fpi_e_tparadata',...
    %             'fpi_e_tperpdata',...
    %             'fpi_e_edata',...
    %             'fpi_e_espectdata',...
    %             'fpi_e_pressdata',...
    %             'fpi_i_timedata',...
    %             'fpi_i_ndata',...
    %             'fpi_i_vdata',...
    %             'fpi_i_tparadata',...
    %             'fpi_i_tperpdata',...
    %             'fpi_i_edata',...
    %             'fpi_i_espectdata',...
    %             'fpi_i_pressdata',...
    %             'time',...
    %             'phi_vector',...
    %             'theta_vector',...
    %             'energy_vector',...
    %             'dist')
    %   else
    %        load(eventDataDirectoryFileName) %#ok<LOAD>
    %   end
    
    
    
    load(eventDataDirectoryFileName) %#ok<LOAD>
%     [time,phi_vector,theta_vector,energy_vector,dist] = load_dist(timeLeft,timeRight,1,'brst',specie);
    
    %Times
%     formatIn='yyyy-mm-dd HH:MM:SS.FFF';
%     tstart = datenum(timeLeft,formatIn);
%     tend = datenum(timeRight,formatIn);
%     
    
    %Find start and end indices.
    start_index = find(time >= timeLeft, 1);
    end_index = find(time >= timeRight, 1);
    
    %Convert the datetime to date String, and then crop to our event timeframe
    time = time(start_index:end_index);
    phi_vector = phi_vector(start_index:end_index,:);
    energy_vector = energy_vector(start_index:end_index,:);
    dist = dist(:,:,:,start_index:end_index);
    
    %Change Angles
%     phi_vector = mod(phi_vector + 180,360); %ensure between [0 and 360)
%     theta_vector = theta_vector - 90; %range from -90 to 90.
    
    
    
    
    
    
    
    
    
    %Calculate Bins
    r_bins = [0,299792458 ./1000 .* sqrt(1 - 1./(energy_vector(1,:)/939513712.93574440+1).^2)];
    angle_bins = 0:360/32:360;
    
    
    X = r_bins'*cosd(angle_bins);
    Y = r_bins'*sind(angle_bins);
    
    v = [];
    phi = [];
    
    for i=1:length(time)
        
        [phi_indices_perTime,theta_indices_perTime,energy_indices_perTime] = ind2sub(size(dist(:,:,:,i)),find(dist(:,:,:,i)));
        
        for j=1:length(phi_indices_perTime)
            if (theta_indices_perTime(j) >= 7 && theta_indices_perTime(j) <= 10) && (phi_indices_perTime(j) <= 4 || phi_indices_perTime(j) >= 29)
                phi = [phi, phi_vector(i,phi_indices_perTime(j))];
                [vx,vy,vz,v_mag] = getVelocity...
                    (energy_vector(i,energy_indices_perTime(j)),theta_vector(theta_indices_perTime(j)),phi_vector(i,phi_indices_perTime(j)));
                v = [v;[vx,vy,vz,v_mag]];
            end
        end
        
    end
    
    % Calculate Magnitude of velocity components
    C = histcounts2(v(:,4),phi',[r_bins inf],[angle_bins inf]);
    if plotting == 1
        figure
        h = pcolor(X,Y,C) ;
        set(gcf,'color','w');
        set(h, 'EdgeColor', 'none');
        myColorMap = parula(256); myColorMap(1,:) = 1;
        colormap(myColorMap); colorbar
        
        axis equal tight
        A = datestr(timeRight);
        title(strcat('MMS 1 FPI Ion (xy) - ',datestr(timeLeft),'-',A(13:20))); xlabel('Vx'); ylabel('Vy')
        
        plot_name =  strcat('VelocityDistribution_',datestr(timeLeft),'-',A(13:20),'.png');
        print(gcf,'-dpng','-r300',plot_name);
    end
    
    %Remove small counts in bins
    C(C<2) = 0;
    
    %Calculate peak velocity of distribution of solar wind
    [Crow,Ccol] = ind2sub(size(C),find(C==max(C,[],'all')));
    [Crow,I] = max(Crow);
    Ccol = Ccol(I);
    peak_rowIndex_Range = Crow-peak_radius:Crow+peak_radius;
    peak_colIndex_Range = Ccol-peak_angle:Ccol+peak_angle;
    peak_rowIndex_Range(peak_rowIndex_Range > 33) = [];
    %     if Crow > 17
    %         peak_rowIndex_Range(peak_rowIndex_Range < 17) = [];
    %     else
    %         peak_rowIndex_Range(peak_rowIndex_Range > 17) = [];
    %     end
    
    if Ccol > 17
        peak_colIndex_Range(peak_colIndex_Range < 17) = [];
    else
        peak_colIndex_Range(peak_colIndex_Range > 17) = [];
    end
    
    peak_C = C(peak_rowIndex_Range,peak_colIndex_Range);
    peak_Vx_bins = X(peak_rowIndex_Range,peak_colIndex_Range);
    peak_Vy_bins = Y(peak_rowIndex_Range,peak_colIndex_Range);
    peak_Vmag_bins = (peak_Vx_bins.^2 + peak_Vy_bins.^2).^(1/2);
    
    peak_V = [];
    
    for i=1:length(peak_rowIndex_Range)
        for j=1:length(peak_colIndex_Range)
            for k=1:peak_C(i,j)
                peak_V = [peak_V peak_Vmag_bins(i,j)];
            end
        end
    end
    nsw = sum(peak_C,'all');
    nt = sum(C,'all');
    densityRatio = nt/(nt-nsw);
    speedDist = fitdist(peak_V','Normal');
    SolarwindBeamSpeed = speedDist.mu;
end

function [vx,vy,vz,v_mag] = getVelocity(E,theta,phi)
    %Calculate velocity
    restEnergy = 939513712.93574440;
    v = 299792458 /1000 * sqrt(1 - 1/(E/restEnergy+1)^2);
    
    %Calculate Components & Magnitude
    vx = v*cosd(phi)*cosd(theta);
    vy = v*sind(phi)*cosd(theta);
    vz = v*sind(theta);
    
    v_mag = sqrt(vx^2 + vy^2 + vz^2);
end
