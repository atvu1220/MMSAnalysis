%Plots Pitch Angle Distributions per Energy Level
clear




%Parameters
event_start = '2015-10-16 13:06:40.985';
event_end = '2015-10-16 13:07:20.985';


% event_start = '2015-10-16 13:06:59.685';
% event_end = '2015-10-16 13:07:00.285';

pitch_angle_bins = 180;
specie = 'i';
energyRange = 'all';


event_start = 6;



if size(event_start,2) > 3
else
    %Go by Event number then
    [~, ~,~,~,...
          event_start,event_end,...
          ~, ~,...
          ~, ~,...
          ~, ~,...
          ~, ~] = get_eventTimes(event_start);
end

%Organize plots by Date
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Pitch Angle Distributions'
eventFolderName = strcat(event_start(1:10),'/',event_start(12:13),'-',event_start(15:16),energyRange);
if ~exist(eventFolderName,'dir')
    mkdir(eventFolderName)
end
cd(eventFolderName)


%Plot all or just each 4.
if ~strcmp(energyRange,'all')
    

%% Plotting Calculations
figure1 = figure('Position',[1 1 650 850]);
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_gap=1.10;

num_plots = 9;
plot_order = 1;


%first four energy levels + plotting properties
axes('Position',[0.0 0.0 0.001 0.001])
plot_paDistribution(event_start,event_end,pitch_angle_bins,specie,28:32,num_plots,plot_order)
title(strcat('MMS1 Pitch Angle Distribution:',event_start), 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');

colorbar
Ticks = get(colorbar,'Ticks');
set(colorbar,'TickLabels',(strcat(repmat('10^{',size(Ticks,2)',1),num2str(Ticks'),'}')))

set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order +1;
hold on

%4 energy channels per plot
for i=6:-1:0
    energyChannel = 4*i+1:4*(i+1);
    axes('Position',[0.0 0.0 0.001 0.001])
    plot_paDistribution(event_start,event_end,pitch_angle_bins,specie,energyChannel,num_plots,plot_order)
    
    yticks([0 45 90 135 180])
    
    colorbar
    Ticks = get(colorbar,'Ticks');
    set(colorbar,'TickLabels',(strcat(repmat('10^{',size(Ticks,2)',1),num2str(Ticks'),'}')))
    
    set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
    plot_order = plot_order +1;
end

fileName = strcat('PitchAngleDistribution_',specie,'_',num2str(event_start(15:23)),'-',num2str(event_end(15:23)));
print(gcf,'-depsc2', '-loose', strcat(fileName,'.eps'));
else
    
%All Energies



figure1 = figure('Position',[1 1 800 200]);
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_gap=1.10;

num_plots = 1;
plot_order = 1;


energyChannel = 1:32;
axes('Position',[0.0 0.0 0.001 0.001])
plot_paDistribution(event_start,event_end,pitch_angle_bins,specie,energyChannel,num_plots,plot_order)
title(strcat('MMS1 Pitch Angle Distribution:',event_start), 'FontSize', 18, 'FontWeight', 'normal')

yticks([0 45 90 135 180])

colorbar
Ticks = get(colorbar,'Ticks');
set(colorbar,'TickLabels',(strcat(repmat('10^{',size(Ticks,2)',1),num2str(Ticks'),'}')))

hold on
%Plot Magnetic Field Magnitude Overlay
[fgm_timedata,fgm_bdata,~,~] = load_fgm(event_start,event_end,1,'brst','GSE');
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    
    %Find the start and end limits of the event in the data
    start_index = find(fgm_timedata >= tstart, 1);
    end_index = find(fgm_timedata >= tend, 1);
    
    %scale Bmag
    scale = 180/max(fgm_bdata(:,4));
    
    %crop data
    fgm_timedata = fgm_timedata(start_index:end_index,1);
    fgm_bdata = fgm_bdata(start_index:end_index,:);
    plot(fgm_timedata,scale*fgm_bdata(:,4),'LineWidth',1.25)




%set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
datetick('x','keeplimits')




fileName = strcat('PitchAngleDistribution_',specie,'_',num2str(event_start(15:23)),'-',num2str(event_end(15:23)));
print(gcf,'-depsc2', '-loose', strcat(fileName,'.eps'));
end