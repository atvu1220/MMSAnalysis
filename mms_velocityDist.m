%Calculates Velocity Distributions
%Uses V_mag as r magnitude. and theta for theta
close all
clear 
time_center = '2018-01-09 08:35:02.000'; %Center time
window = 2500; %time window in milliseconds, total centered on center time
specie = 'i';

%Take the time interval before and after center time
formatIn='yyyy-mm-dd HH:MM:SS.FFF';
timeCenter = datenum(time_center,formatIn);

timeLeft = timeCenter - milliseconds(window);
timeRight = timeCenter + milliseconds(window);

%Get the Distribution Data
[time,phi_vector,theta_vector,energy_vector,dist] = load_dist(timeLeft,timeRight,1,'brst',specie);

phiOffset = [linspace(0,0,16),linspace(360,360,16)];
phi_midvalues = phi_vector + phiOffset;
phi_midvalues = (phi_midvalues + circshift(phi_midvalues,-1,2))/2 - phiOffset + [linspace(0,0,31) 180];
% theta_midvalues = (theta_vector + circshift(theta_vector,-1,2))/2;

v = zeros(length(time),4);


%Calculate Bins for Pcolor

restEnergy = 939513712.93574440;
v_bins = [0,299792458 ./1000 .* sqrt(1 - 1./(energy_vector(1,:)/restEnergy+1).^2)];
phi_bins = 0:360/32:360;
theta_bins = theta_vector;

X = v_bins'*cosd(phi_bins);
Y = v_bins'*sind(phi_bins);

v = [];
phi = [];
energy_indices = [];
phi_indices = [];
theta_indices = [];
for i=1:length(time)
    
    [phi_indices_perTime,theta_indices_perTime,energy_indices_perTime] = ind2sub(size(dist(:,:,:,i)),find(dist(:,:,:,i)));
    
    for j=1:length(phi_indices_perTime)
        phi = [phi, phi_midvalues(i,phi_indices_perTime(j))];
        
%         [vx,vy,vz,v_mag] = getVelocity...
%             (energy_vector(i,energy_indices_perTime(j)),theta_vector(theta_indices_perTime(j)),phi_vector(i,phi_indices_perTime(j)));
        
%         v = [v;[vx,vy,vz,v_mag]];
    end
    theta_indices = [theta_indices; theta_indices_perTime];
    energy_indices = [energy_indices; energy_indices_perTime];
    %Calculate V directly
    %then bin after loop.
end

velocity_midvalues = (v_bins(1:end-1)+v_bins(2:end))/2;
velocity = velocity_midvalues(energy_indices);
%phi = phi_bins(phi_indices);

C = histcounts2(velocity,phi,[v_bins inf],[phi_bins inf]);
h = pcolor(X,Y,C) ;
set(h, 'EdgeColor', 'none');
myColorMap = parula(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar

axis equal tight
A = datestr(timeRight);
title(strcat('MMS 1 FPI Ion (xy) - ',datestr(timeLeft),'-',A(13:20)))
xlabel('Vx')
ylabel('Vy')



%Calculate peak velocity of distribution of solar wind?
v_bin_counts = histcounts(velocity,[v_bins inf]);
figure
bar(v_bins,v_bin_counts);










% n = 6;
% r = (0:n)'/n;
% theta = pi*(-n:n)/n;
% X = r*cos(theta);
% Y = r*sin(theta);
% C = r*cos(2*theta);
% pcolor(X,Y,C)
% axis equal tight 

% figure
% scatter(v(:,1),v(:,2))
% axis equal tight

function [vx,vy,vz,v_mag] = getVelocity(E,theta,phi)

    %Calculate velocity
    restEnergy = 939513712.93574440;
    v = 299792458 /1000 * sqrt(1 - 1/(E/restEnergy+1)^2);

    %Calculate Components & Magnitude
    vx = v*cosd(phi)*cosd(theta);
    vy = v*sind(phi)*cosd(theta);
    vz = v*sind(theta);
    
    v_mag = sqrt(vx^2 + vy^2 + vz^2);

    %Calculate PitchAngle and determine its bin
%     pitchAngle = acos((B_norm(1)*vx + B_norm(2)*vy + B_norm(3)*vz)/v_mag)/pi*180;
%     pitchAngle_bins(index) = floor(pitchAngle/180*pa_n)+1;
%     
%     %Keep track of f
%     distributionFunction(index) = f;
end