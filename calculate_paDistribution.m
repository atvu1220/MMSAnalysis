function [time,pa_bins_values,f_bins,energy_bins_values] = calculate_paDistribution(event_start,event_end,pa_n,specie,energyChannel)
    %Function file to calculate Pitch Angle Distributions
    
    IntegrationWindow = 4/2; %/Number of points on either side of center time to integrate over
    %Parameters
    %     centerTime = '2015-10-16 13:06:59.985'; %Center time
    %     window = 300*2; %time window in seconds, total centered on center time
    %     pa_n=37; %Number of Bins for Pitch angle
    %     specie = 'i';
    
    
    %% Load Dist and FGM data
    
    %Get the Distribution Data
    [time,phi_vector,theta_vector,energy_vector,dist,start_index,end_index] = load_dist(event_start,event_end,1,'brst',specie);
    
    %Load FGM and average over entire time interval
    [Btimedata,Bdata,~,~] = load_fgm(event_start,event_end,1,'brst','DMPA');
    [~,~,Bstart,Bend]=crop(Btimedata,Bdata,event_start,event_end);
    
    
    %% For the pitch angle distribution for each time step.
    
    pa_bins_values = linspace(0,180,pa_n);
    pa_bins = zeros(pa_n,length(time));
    f_bins = zeros(pa_n,length(time));
    B_norm = Bdata(Bstart:Bend,1:3)./Bdata(Bstart:Bend,4); %B norm for each time step
    energy_bins_values = mean(energy_vector(1:2,:));
    for i=(IntegrationWindow+1):1:length(time)-(IntegrationWindow+1)
        
        for m=-IntegrationWindow:1:IntegrationWindow
            integrationStep = i+m;
            
            
            [phi_indices,theta_indices,energy_indices] = ind2sub(size(dist(:,:,energyChannel,integrationStep)),find(dist(:,:,energyChannel,integrationStep)));
           %[phi_indices,theta_indices,energy_indices] = ind2sub(size(dist(:,:,energyChannel,i)),find(dist(:,:,energyChannel,i)));
            
            pitchAngle_bins = zeros(length(phi_indices),1); %stores all pitch angle bins.
            distributionFunction = zeros(length(phi_indices),1); %store distribution function.
            
            for j=1:length(phi_indices) %each nonzero element per time in all of distribution function
                [distributionFunction, pitchAngle_bins] = fillPitchAngleBins(distributionFunction,pitchAngle_bins,...
                    energy_vector(integrationStep,energy_indices(j)),theta_vector(theta_indices(j)),phi_vector(integrationStep,phi_indices(j)),...
                    dist(phi_indices(j),theta_indices(j),energy_indices(j),integrationStep),...
                    B_norm(i,:),pa_n,j);
            end
            
            
            [f_bins,pa_bins] = sumPitchAngleBins(distributionFunction,pitchAngle_bins,f_bins,pa_bins,i);
            
            
            
        end
        
    end
%     imagesc(time,pa_bins_values',(f_bins))
%     shading interp
    
    
    %pcolor(time,pa_bins_values',(f_bins))
    
    
    
    
    
    % set(gcf,'defaultAxesColorOrder',co)
    % set(gcf,'color','w');
    % shading interp
    % colormap(jet)
    %
    % xlim([time(1) time(end)])
    % ylabel({'Pitch Angle';'[\Theta]'},'FontSize', 14)
    % set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    % title(strcat('MMS1 Pitch Angle Distribution:',datestr(timeLeft)), 'FontSize', 18, 'FontWeight', 'normal')
    % datetick('x','keeplimits')
    
    
    
    %% Calculate the total PA distribution over the entire time interval for all energies
    %     %Integrated within time range, plot pitch angle vs Energy
    %     pa_bins_values = linspace(0,180,pa_n);
    %     pa_bins = zeros(pa_n,size(energy_vector,2));
    %     f_bins = zeros(pa_n,size(energy_vector,2));
    %     B_norm = mean(Bdata(Bstart:Bend,1:3))/mean(Bdata(Bstart:Bend,4));
    %
    %     for i=1:size(energy_vector,2)
    %         for j=1:length(time)
    %             [phi_indices,theta_indices] = find(reshape(dist(:,:,i,j),32,16)~=0);
    %
    %             pitchAngle_bins = zeros(length(theta_indices),1); %stores all pitch angle bins.
    %             distributionFunction = zeros(length(theta_indices),1); %store distribution function.
    %
    %             for k=1:length(theta_indices)
    %                 [distributionFunction, pitchAngle_bins] = fillPitchAngleBins...
    %                     (distributionFunction,pitchAngle_bins,...
    %                     energy_vector(j,i),theta_vector(theta_indices(k)),phi_vector(j,phi_indices(k)),...
    %                     dist(phi_indices(k),theta_indices(k),i,j),...
    %                     B_norm,pa_n,k);
    %             end
    %             [f_bins,pa_bins] = sumPitchAngleBins(distributionFunction,pitchAngle_bins,f_bins,pa_bins,i);
    %         end
    %     end
    %
    %
    %
    %
    %     energy_bins_values = mean(energy_vector);
    %     figure
    %     pcolor(pa_bins_values,log10(energy_bins_values)',log10(f_bins'))
    %     set(gcf,'defaultAxesColorOrder',co)
    %     set(gcf,'color','w')
    %     shading interp
    %     colormap(jet)
    %
    %     color_bar = colorbar('Ticks', [-28, -27, -26, -25, -24,-23,-22,-21,-20],...
    %         'TickLabels', {'10^{-28}', '10^{-27}', '10^{-26}', '10^{-25}', '10^{-24}','10^{-23}','10^{-22}','10^{-21}','10^{-20}'},'FontSize', 10);
    %     ylabel(color_bar,{'f(s^3/cm^6'},'FontSize', 12)
    %     ylabel({'Energy';'[eV]'},'FontSize', 14)
    %     yticklabels(num2str([32;100;320;1000;3200;10000],'%2.f'))
    %     xlabel({'Pitch Angle';'[\Theta]'},'FontSize', 14)
    %     set(gca,'XMinorTick','on','linewidth',1.25)
    %     title(strcat('MMS1 Pitch Angle Distribution:',datestr(timeLeft)), 'FontSize', 18, 'FontWeight', 'normal')
    %
    
    
    
    
    %%
    function [distributionFunction, pitchAngle_bins] = fillPitchAngleBins(distributionFunction,pitchAngle_bins,E,theta,phi,f,B_norm,pa_n,index)
        
        %Calculate velocity
        restEnergy = 939513712.93574440;
        v = 299792458 /1000 * sqrt(1 - 1/(E/restEnergy+1)^2);
        
        %Calculate Components & Magnitude
        vx = v*cos(phi/180 * pi)*cos(theta/180 * pi);
        vy = v*sin(phi/180 * pi)*cos(theta/180 * pi);
        vz = v*sin(theta/180*pi);
        
        v_mag = sqrt(vx^2 + vy^2 + vz^2);
        
        %Calculate PitchAngle and determine its bin
        pitchAngle = acos((B_norm(1)*vx + B_norm(2)*vy + B_norm(3)*vz)/v_mag)/pi*180;
        pitchAngle_bins(index) = floor(pitchAngle/180*pa_n)+1;
        
        %Keep track of f
        distributionFunction(index) = f;
    end
    %%
    function [f_bins,pa_bins] = sumPitchAngleBins(distributionFunction,pitchAngle_bins,f_bins,pa_bins,index)
        %Fill in pitch angle data bins
        %Sort by pitchangle bin
        [pitchAngle_bins,paOrder] = sort(pitchAngle_bins);
        distributionFunction = distributionFunction(paOrder);
        
        %Count the number of events in each bin, then sum and save to bin. Sum up the corresponding f's
        for l=unique(pitchAngle_bins)'
            first = find(sort(pitchAngle_bins)==l,1,'first');
            last = find(sort(pitchAngle_bins)==l,1,'last');
            pa_bins(l,index) = last - first+1;
            f_bins(l,index) = sum(distributionFunction(first:last));
            f_bins(f_bins==0) = NaN;
        end
    end
end