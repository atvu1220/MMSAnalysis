%Superimposed Epoch Analysis
clear
close all
tic
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Superimposed Epoch'

figure('Position',[1 1 300 900])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_gap=1.11;
numPlots = 8;
plotOrder = 1;
%Manually Choose Events Method
% Event_Number_vector = [2,8,11,26,76,79,81,82,84,89,90,95];
% numberEvents = length(Event_Number_vector);
% timescale_vector = zeros(size(Event_Number_vector));

% for Event_Number = Event_Number_vector
%
%     %Get the event times
%     [event_start,event_end] = get_event(Event_Number);
%
%
%     %Store the timescales in a vector for all events, for later averaging.
%     timescale_vector(Event_Number==Event_Number_vector) = timescale;
%
% end


%Automate finding of events by Event_Type
Event_Number_vector = 2:1:85;
timescale_vector = zeros(size(Event_Number_vector));
Type = 'Substructures HFA';

parfor i=2:length(Event_Number_vector)
    
    %Get the event times
    [Event_Type,Substructure,~,~,event_start,event_end,~,~,~,~,~,~,~,~] = get_eventTimes(Event_Number_vector(i));
    %     [event_start,event_end] = get_event(Event_Number);
    if Substructure == 1
        Substructure_Type = 'Substructures'
    elseif Substructure == 0
        Substructure_Type = 'No Substructure'
    end
    
    if contains(Type,Event_Type) && xor( contains(Type, Substructure_Type) , ~((contains(Type, 'Substructures') && (contains(Type, 'No Substructure') ) )))
        
        %Calculate the timescale of event in datapoints.
        [timescale] = get_Timescale(event_start,event_end);
        
        %Store the timescales in a vector for all events, for later averaging.
        timescale_vector(i) = timescale;
    end
end

timescale_vector = timescale_vector(timescale_vector~=0);
Event_Number_vector = Event_Number_vector(timescale_vector~=0);
numberEvents = length(Event_Number_vector);


%Calculate the mean time scale and construct reference time scale for all events for use of data
%binning
mean_timescale = round(mean(timescale_vector));
timescale_databins = (-mean_timescale/2:1:+mean_timescale/2)*0.15;


%Construct the structure for all events, as a databin
Databins = struct(...
    'time',timescale_databins,...
    'density',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'magneticfield',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'velocity_x',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'velocity_y',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'velocity_z',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'tpara',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'tperp',zeros(length(timescale_databins),length(Event_Number_vector)),...
    'dynamicpressure',zeros(length(timescale_databins),length(Event_Number_vector)));

%This time get all the data and calculate the average parameters.
for Event_Number = Event_Number_vector
    
    %[event_start,event_end] = get_Event(Event_Number);
    
    [~,~,~,~,event_start,event_end,~,~,~,~,~,~,~,~] = get_eventTimes(Event_Number);
    
    [Event] = Event_Structure(event_start,event_end);
    
    %Calculate dynamic pressure rho*v^2
    Event.dynamicpressure = calculate_dynamic_pressure(Event.density,Event.velocity,'i');
    
    
    %Interpolate this one event's data into the mean timescale and put into databins
    Databins.density(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.density,Databins.time)';
    
    Databins.magneticfield(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.magneticfield,Databins.time)';
    
    Databins.velocity_x(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.velocity(:,1),Databins.time)';
    Databins.velocity_y(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.velocity(:,2),Databins.time)';
    Databins.velocity_z(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.velocity(:,3),Databins.time)';
    
    Databins.tpara(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.tpara,Databins.time)';
    Databins.tperp(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.tperp,Databins.time)';
    
    Databins.dynamicpressure(:,Event_Number==Event_Number_vector) = ...
        interp1(Event.time*mean_timescale/(length(Event.time)),Event.dynamicpressure,Databins.time)';
    
    
end
plot_superimposedEpoch(Databins.time,Databins.magneticfield,{'Magnetic Field';'[nT]'},Type,numberEvents,numPlots,plotOrder)
plot_pos = get(gca,'Position');
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.density,{'Density';'[cm^{-3}]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.velocity_x,{'Vx';'[km/s]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.velocity_y,{'Vy';'[km/s]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.velocity_z,{'Vz';'[km/s]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.tpara,{'Tpara';'[eV]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.tperp,{'Tperp';'[eV]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.dynamicpressure,{'Dynamic Pressure';'[nPa]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

fileName = strcat(Type,'_Superimposed Epoch Summary_',num2str(numberEvents));
print(gcf,'-depsc2', '-loose', strcat(string(fileName),'.eps'));
















toc


%% Find the timescale of the event in the number of data points.
function [timescale] = get_Timescale(event_start,event_end)
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    timescale = (tend-tstart)*86400/0.15 +1;
end
%% Get the data, crop the data, and store in a structure.
function [Event_Struct] = Event_Structure(event_start,event_end)
    %Load FPI_i
    [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tpara,fpi_i_tperp,...
        ~,~] = load_fpi(event_start,event_end,1,'brst','i');
    
    
    
    %Crop the data to the specified time range.
    [~,fpi_i_ndata,~,~] = crop(fpi_i_timedata,fpi_i_ndata,event_start,event_end);
    [~,fpi_i_tpara,~,~] = crop(fpi_i_timedata,fpi_i_tpara,event_start,event_end);
    [~,fpi_i_tperp,~,~] = crop(fpi_i_timedata,fpi_i_tperp,event_start,event_end);
    [fpi_i_timedata,fpi_i_vdata,~,~] = crop(fpi_i_timedata,fpi_i_vdata,event_start,event_end);
    
    %load fgm data & crop
    [fgm_timedata, fgm_bdata,~,~]= load_fgm(event_start,event_end,'1','brst');
    [fgm_timedata,fgm_bdata,~,~] = crop(fgm_timedata,fgm_bdata,event_start,event_end);
    %[~,fgm_bdata]=interpxyz(fgm_timedata,fgm_bdata(:,1:3),fpi_i_timedata); %three components only;
    fgm_bdata = interp1(fgm_timedata,fgm_bdata(:,4),fpi_i_timedata); %magnitude only
    
    %Center the time in the event.
    fpi_i_timedata = 86400*(fpi_i_timedata - median(fpi_i_timedata));
    
    %Construct the structure for the event
    Event_Struct = struct(...
        'start',event_start,...
        'end',event_end,...
        'time',fpi_i_timedata,...
        'magneticfield',fgm_bdata,...
        'density',fpi_i_ndata,...
        'velocity',fpi_i_vdata,...
        'tpara',fpi_i_tpara,...
        'tperp',fpi_i_tperp,...
        'dynamicpressure',[]);
end
%% Calculate the quartiles for Superimposed Epoch Analysis
function [first,mean,third] = calculate_quartiles(parameter)
    mean = quantile(parameter, 0.5,2);
    first = quantile(parameter, 0.25,2);
    third = quantile(parameter, 0.75,2);
end
function [] = plot_superimposedEpoch(time,parameter,parameterStrings,Type,numberEvents,numPlots,plotOrder)
    
    [first, mean, third] = calculate_quartiles(parameter);
    
    subplot(numPlots,1,plotOrder);
    
    
    plot(time,first,'DisplayName','1^{st} Quartile','LineWidth',1,'LineStyle','--','Color','k'); hold on
    plot(time,mean,'DisplayName','Mean','LineWidth',1.5,'LineStyle','-','Color','r');
    plot(time,third,'DisplayName','3^{rd} Quartile','LineWidth',1,'LineStyle','--','Color','k');
    
    
    
    xlim([min(time) max(time)])
    ylabel(parameterStrings)
    %     legend('boxoff')
    %     legend('Location','northwest','FontSize',16,'orientation','horizontal')
    
    
    
    if plotOrder ~= numPlots
        set(gca,'XTickLabel',[], 'XMinorTick','on','YMinorTick','on','linewidth',1.5,'FontSize',8)
    else
        set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.5,'FontSize',8)
        xlabel({'Model Time';'[s]'},'FontSize',12)
    end
    
    if plotOrder == 1
        title(string({'Superimposed Epoch Analysis'; strcat({' '}, num2str(numberEvents), ...
            {' '},Type,'s')}),'FontSize',14,'FontWeight', 'normal')
    end
    
    
end