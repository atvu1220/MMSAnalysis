function [downstreamSpeed,upstreamSpeed,downDensityRatio,upDensityRatio] = calculate_updownSW3D(event_start,event_end)
    %Calculates the solar wind beam speed by first plotting the distribution in Vx-Vy with magnitude
    %V and angle of phi from distribution function. Then find the max histcount, and fit a normal
    %distribution around the peak for the speed.
    
    window = 1500; %time window in milliseconds, total centered on center time
    specie = 'i';
    peak_angle = 3;
    peak_radius = 6; %6 works for all but Event 13
    
    %Take the time interval before and after center time
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    eventStart = datenum(event_start,formatIn);
    eventEnd = datenum(event_end,formatIn);
    
    down_centerTime = eventStart + seconds(3);
    up_centerTime = eventEnd - seconds(3);
    
    [downstreamSpeed,downDensityRatio] = calculate_solarwindBeam(down_centerTime,window,peak_angle,peak_radius,specie);
    [upstreamSpeed,upDensityRatio] = calculate_solarwindBeam(up_centerTime,window,peak_angle,peak_radius,specie);
    
end
function[SolarwindBeamSpeed,densityRatio] = calculate_solarwindBeam(centerTime,window,peak_angle,peak_radius,specie)
    timeLeft = centerTime - milliseconds(window);
    timeRight = centerTime + milliseconds(window);
    
    %Get the Distribution Data
    [time,phi_vector,theta_vector,energy_vector,dist] = load_dist(timeLeft,timeRight,1,'brst',specie);
    
    %Calculate Bins
    r_bins = [0,299792458 ./1000 .* sqrt(1 - 1./(energy_vector(1,:)/939513712.93574440+1).^2)];
    angle_bins = 0:360/32:360;
    
    
    X = r_bins'*cosd(angle_bins);
    Y = r_bins'*sind(angle_bins);
    
    v = [];
    phi = [];
    nt = 0;
    for i=1:length(time)
        
        [phi_indices_perTime,theta_indices_perTime,energy_indices_perTime] = ind2sub(size(dist(:,:,:,i)),find(dist(:,:,:,i)));
        
        for j=1:length(phi_indices_perTime) %Changed theta for 2 channels instead of 1, 10/2019
            if (theta_indices_perTime(j) >= 7 && theta_indices_perTime(j) <= 10) && (phi_indices_perTime(j) <= 4 || phi_indices_perTime(j) >= 29)
                phi = [phi, phi_vector(i,phi_indices_perTime(j))];
                [vx,vy,vz,v_mag] = getVelocity...
                    (energy_vector(i,energy_indices_perTime(j)),theta_vector(theta_indices_perTime(j)),phi_vector(i,phi_indices_perTime(j)));
                v = [v;[vx,vy,vz,v_mag]];
                %Check common factor for all dist numbers.
            end
        end
        length(phi_indices_perTime)
        nt = nt + length(phi_indices_perTime)
    end
    %max dist, probabably SW beam
    [phi3Dpeak,theta3Dpeak,energy3Dpeak] = ind2sub(size(sum(dist,4)),find(sum(dist,4) == max(sum(dist,4),[],'all')))
    % Calculate Magnitude of velocity components
    C = histcounts2(v(:,4),phi',[r_bins inf],[angle_bins inf]);
    figure
    h = pcolor(X,Y,C) ;
    set(gcf,'color','w');
    set(h, 'EdgeColor', 'none');
    myColorMap = parula(256); myColorMap(1,:) = 1;
    colormap(myColorMap); colorbar
    
    axis equal tight
    A = datestr(timeRight);
    title(strcat('MMS 1 FPI Ion (xy) - ',datestr(timeLeft),'-',A(13:20))); xlabel('Vx'); ylabel('Vy')
    
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
    
    plot_name =  strcat('VelocityDistribution_',datestr(timeLeft),'-',A(13:20),'.png');
    print(gcf,'-dpng','-r300',plot_name);
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
