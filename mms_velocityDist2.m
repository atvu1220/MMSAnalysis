%Calculates Velocity Distributions
%This version uses magnitude of vx and vy on polar distribution plot, instead of v_mag.

close all; clear

time_center = '2018-01-12 01:52:20.000'; %Center time
window = 1500; %time window in milliseconds, total centered on center time
specie = 'i';
peak_angle = 5;
peak_radius = 10;

%Take the time interval before and after center time
formatIn='yyyy-mm-dd HH:MM:SS.FFF'; timeCenter = datenum(time_center,formatIn);

timeLeft = timeCenter - milliseconds(window);
timeRight = timeCenter + milliseconds(window);

%Get the Distribution Data
[time,phi_vector,theta_vector,energy_vector,dist] = load_dist(timeLeft,timeRight,1,'brst',specie);

%Increase bins by adding half bins
% energy_vector2 = circshift([energy_vector(1,:),30000],-1,2);
% energy_vector2 = (energy_vector(1,:) + energy_vector2(1:end-1))./2;
% energy_vector = min(energy_vector,[],'all'):(max(energy_vector,[],'all')-min(energy_vector,[],'all'))/100:max(energy_vector,[],'all');

%Calculate Bins
r_bins = [0,299792458 ./1000 .* sqrt(1 - 1./(energy_vector(1,:)/939513712.93574440+1).^2)];
% r_bins = min(r_bins):(max(r_bins)-min(r_bins))/16:max(r_bins);
angle_bins = 0:360/32:360;
% angle2_bins = 0:180/32:180;

X = r_bins'*cosd(angle_bins);
Y = r_bins'*sind(angle_bins);

% [XX,YY,ZZ] = meshgrid(r_bins,angle_bins,angle2_bins);
% [XX,YY,ZZ] = sphere(32);
v = [];
phi = [];
theta = [];

for i=1:length(time)
    
    [phi_indices_perTime,theta_indices_perTime,energy_indices_perTime] = ind2sub(size(dist(:,:,:,i)),find(dist(:,:,:,i)));
    
    for j=1:length(phi_indices_perTime)
        if (theta_indices_perTime(j) >= 8 && theta_indices_perTime(j) <= 9) && (phi_indices_perTime(j) <= 4 || phi_indices_perTime(j) >= 29)
            phi = [phi, phi_vector(i,phi_indices_perTime(j))];
            % theta = [theta, theta_vector(theta_indices_perTime(j))];
            [vx,vy,vz,v_mag] = getVelocity...
                (energy_vector(i,energy_indices_perTime(j)),theta_vector(theta_indices_perTime(j)),phi_vector(i,phi_indices_perTime(j)));
            v = [v;[vx,vy,vz,v_mag]];
            theta = [theta theta_vector(theta_indices_perTime(j))];
        end
    end
    
end

% Calculate Magnitude of velocity components
vxy = (v(:,1).^2 + v(:,2).^2).^(1/2);

Cphi = histcounts2(v(:,4),phi',[r_bins inf],[angle_bins inf]);
h = pcolor(X,Y,Cphi) ;
set(gcf,'color','w');
set(h, 'EdgeColor', 'none');
myColorMap = parula(256); myColorMap(1,:) = 1;
colormap(myColorMap); colorbar
xlabel({'Vx','[km/s]'})
ylabel({'Vy','[km/s]'})

axis equal tight
A = datestr(timeRight);
title(strcat('MMS 1 FPI Ion (xy) - ',datestr(timeLeft),'-',A(13:20))); xlabel('Vx'); ylabel('Vy')



% %Calculate peak velocity of distribution of solar wind?
[Crow,Ccol] = ind2sub(size(Cphi),find(Cphi==max(Cphi,[],'all')));
peak_rowIndex_Range = Crow-peak_radius:Crow+peak_radius;
peak_colIndex_Range = Ccol-peak_angle:Ccol+peak_angle;
peak_rowIndex_Range(peak_rowIndex_Range > 33) = [];
peak_rowIndex_Range(peak_rowIndex_Range < 17) = [];



peak_C = Cphi(peak_rowIndex_Range,peak_colIndex_Range);
peak_Vx_bins = X(peak_rowIndex_Range,peak_colIndex_Range);
peak_Vy_bins = Y(peak_rowIndex_Range,peak_colIndex_Range);
peak_Vmag_bins = (peak_Vx_bins.^2 + peak_Vy_bins.^2).^(1/2);
peak_phi_bins = 180 + atand(peak_Vy_bins./peak_Vx_bins);


peak_V = [];
peak_phi = [];
for i=1:length(peak_rowIndex_Range)
    for j=1:length(peak_colIndex_Range)
        for k=1:peak_C(i,j)
            peak_V = [peak_V peak_Vmag_bins(i,j)];
            peak_phi = [peak_phi peak_phi_bins(i,j)];
        end
    end
end
speedDist = fitdist(peak_V','Normal');
SolarwindBeamSpeed = speedDist.mu
phiDist = fitdist(peak_phi','Normal');
SolarwindBeamPhi = phiDist.mu


plot_name =  strcat('VelocityDistribution_',datestr(timeLeft),'-',A(13:20),'.png');
print(gcf,'-dpng','-r300',plot_name);
% % % 
% % % %%%%%%%%%%%%%%%%%%%% THeta Calculation
% % % peak_angle = 2;
% % % peak_radius = 6;
% % % angle_bins = -90:180/16:90;
% % % % angle2_bins = 0:180/32:180;
% % % 
% % % X = r_bins'*cosd(angle_bins);
% % % Y = r_bins'*sind(angle_bins);
% % % 
% % % Ctheta = histcounts2(v(:,4),theta',[r_bins inf],[angle_bins inf]);
% % % figure
% % % h2 = pcolor(X,Y,Ctheta) ;
% % % set(gcf,'color','w');
% % % set(h2, 'EdgeColor', 'none');
% % % myColorMap = parula(256); myColorMap(1,:) = 1;
% % % colormap(myColorMap); colorbar
% % % xlabel({'Vy','[km/s]'})
% % % ylabel({'Vz','[km/s]'})
% % % 
% % % axis equal tight
% % % A = datestr(timeRight);
% % % title(strcat('MMS 1 FPI Ion (yz) - ',datestr(timeLeft),'-',A(13:20))); xlabel('Vy'); ylabel('Vz')
% % % 
% % % 
% % % 
% % % % %Calculate peak velocity of distribution of solar wind?
% % % [Crow,Ccol] = ind2sub(size(Ctheta),find(Ctheta==max(Ctheta,[],'all')));
% % % peak_rowIndex_Range = Crow-peak_radius:Crow+peak_radius;
% % % peak_colIndex_Range = Ccol-peak_angle:Ccol+peak_angle;
% % % peak_rowIndex_Range(peak_rowIndex_Range > 33) = [];
% % % peak_rowIndex_Range(peak_rowIndex_Range < 17) = [];
% % % 
% % % 
% % % 
% % % peak_C = Ctheta(peak_rowIndex_Range,peak_colIndex_Range);
% % % peak_Vx_bins = X(peak_rowIndex_Range,peak_colIndex_Range);
% % % peak_Vy_bins = Y(peak_rowIndex_Range,peak_colIndex_Range);
% % % peak_Vmag_bins = (peak_Vx_bins.^2 + peak_Vy_bins.^2).^(1/2);
% % % peak_theta_bins = 180 + atand(peak_Vy_bins./peak_Vx_bins);
% % % 
% % % 
% % % peak_V = [];
% % % peak_theta = [];
% % % for i=1:length(peak_rowIndex_Range)
% % %     for j=1:length(peak_colIndex_Range)
% % %         for k=1:peak_C(i,j)
% % %             peak_V = [peak_V peak_Vmag_bins(i,j)];
% % %             peak_theta = [peak_theta peak_theta_bins(i,j)];
% % %         end
% % %     end
% % % end
% % % 
% % % 
% % % thetaDist = fitdist(peak_theta','Normal');
% % % 
% % % % v_bin_counts = histcounts(peak_V,[r_bins inf]);
% % % % figure
% % % % bar(r_bins,v_bin_counts);
% % % 
% % % SolarwindBeamTheta = thetaDist.mu




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