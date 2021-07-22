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
Event_Number_vector = 2:1:151;
Type = 'HFA';
%Type = 'No Substructure HFA';

%% Core Region 
core_timescale_vector = zeros(size(Event_Number_vector));
leading_timescale_vector = zeros(size(Event_Number_vector));
trailing_timescale_vector = zeros(size(Event_Number_vector));
parfor i=2:length(Event_Number_vector)
    
    %Get the event times
    [Event_Type,Substructure,~,~,~,~,left_OuterEdge, left_InnerEdge,...
          right_InnerEdge, right_OuterEdge,~,~,~,~,~,~,~,~] = get_eventTimes(Event_Number_vector(i));
    %     [event_start,event_end] = get_event(Event_Number);
    if Substructure == 1
        Substructure_Type = 'Substructures'
    elseif Substructure == 0
        Substructure_Type = 'No Substructure'
    end
    
    %if contains(Type,Event_Type) && xor( contains(Type, Substructure_Type) , ~((contains(Type, 'Substructures') && (contains(Type, 'No Substructure') ) )))
    if contains(Event_Type,Type)  
        %Calculate the timescale of event in datapoints.
        [leading_timescale] = get_Timescale(left_OuterEdge,left_InnerEdge);
        [core_timescale] = get_Timescale(left_InnerEdge,right_InnerEdge);
        [trailing_timescale] = get_Timescale(right_InnerEdge,right_OuterEdge);
        
        %Store the timescales in a vector for all events, for later averaging.
        leading_timescale_vector(i) = leading_timescale;
        core_timescale_vector(i) = core_timescale;
        trailing_timescale_vector(i) = trailing_timescale;
    end
end
leading_timescale_vector = leading_timescale_vector(leading_timescale_vector~=0);
core_timescale_vector = core_timescale_vector(core_timescale_vector~=0);
trailing_timescale_vector = trailing_timescale_vector(trailing_timescale_vector~=0);
Event_Number_vector = Event_Number_vector(core_timescale_vector~=0);
numberEvents = length(Event_Number_vector);


%Calculate the mean time scale and construct reference time scale for all events for use of data
%binning
mean_leading_timescale = round(mean(leading_timescale_vector));
mean_core_timescale = round(mean(core_timescale_vector));
mean_trailing_timescale = round(mean(trailing_timescale_vector));

leading_timescale_databins = (-mean_leading_timescale/2:1:+mean_leading_timescale/2)*0.15;
core_timescale_databins = (-mean_core_timescale/2:1:+mean_core_timescale/2)*0.15;
trailing_timescale_databins = (-mean_trailing_timescale/2:1:+mean_trailing_timescale/2)*0.15;
total_databins = length(leading_timescale_databins) + length(core_timescale_databins) + length(trailing_timescale_databins);

all_mean_timescale = mean_leading_timescale + mean_core_timescale + mean_trailing_timescale;
all_timescale_databins = (-(all_mean_timescale+2)/2:1:(all_mean_timescale+2)/2)*0.15;

%Construct the structure for all events, as a databin
Leading_Databins = struct(...
    'time',leading_timescale_databins,...
    'density',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'magneticfield',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'velocity_x',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'velocity_y',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'velocity_z',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'tpara',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'tperp',zeros(length(leading_timescale_databins),length(Event_Number_vector)),...
    'dynamicpressure',zeros(length(leading_timescale_databins),length(Event_Number_vector)));

Core_Databins = struct(...
    'time',core_timescale_databins,...
    'density',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'magneticfield',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'velocity_x',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'velocity_y',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'velocity_z',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'tpara',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'tperp',zeros(length(core_timescale_databins),length(Event_Number_vector)),...
    'dynamicpressure',zeros(length(core_timescale_databins),length(Event_Number_vector)));

Trailing_Databins = struct(...
    'time',trailing_timescale_databins,...
    'density',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'magneticfield',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'velocity_x',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'velocity_y',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'velocity_z',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'tpara',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'tperp',zeros(length(trailing_timescale_databins),length(Event_Number_vector)),...
    'dynamicpressure',zeros(length(trailing_timescale_databins),length(Event_Number_vector)));

%This time get all the data and calculate the average parameters.
for Event_Number = Event_Number_vector
    
    %[event_start,event_end] = get_Event(Event_Number);
    
    [Event_Type,Substructure,~,~,~,~,left_OuterEdge, left_InnerEdge,...
          right_InnerEdge, right_OuterEdge,~,~,~,~,~,~,~,~] = get_eventTimes(Event_Number);
    
    [Leading] = Event_Structure(left_OuterEdge,left_InnerEdge);
    [Core] = Event_Structure(left_InnerEdge,right_InnerEdge);
    [Trailing] = Event_Structure(right_InnerEdge,right_OuterEdge);
    
    %Leading
    %Calculate dynamic pressure rho*v^2
    Leading.dynamicpressure = calculate_dynamic_pressure(Leading.density,Leading.velocity,'i');
    %Interpolate this one event's data into the mean timescale and put into databins
    Leading_Databins.density(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.density,Leading_Databins.time)';
    
    Leading_Databins.magneticfield(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.magneticfield,Leading_Databins.time)';
    
    Leading_Databins.velocity_x(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.velocity(:,1),Leading_Databins.time)';
    Leading_Databins.velocity_y(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.velocity(:,2),Leading_Databins.time)';
    Leading_Databins.velocity_z(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.velocity(:,3),Leading_Databins.time)';
    
    Leading_Databins.tpara(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.tpara,Leading_Databins.time)';
    Leading_Databins.tperp(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.tperp,Leading_Databins.time)';
    
    Leading_Databins.dynamicpressure(:,Event_Number==Event_Number_vector) = ...
        interp1(Leading.time*mean_leading_timescale/(length(Leading.time)),Leading.dynamicpressure,Leading_Databins.time)';
    
    
    %Core
    %Calculate dynamic pressure rho*v^2
    Core.dynamicpressure = calculate_dynamic_pressure(Core.density,Core.velocity,'i');
    %Interpolate this one event's data into the mean timescale and put into databins
    Core_Databins.density(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.density,Core_Databins.time)';
    
    Core_Databins.magneticfield(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.magneticfield,Core_Databins.time)';
    
    Core_Databins.velocity_x(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.velocity(:,1),Core_Databins.time)';
    Core_Databins.velocity_y(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.velocity(:,2),Core_Databins.time)';
    Core_Databins.velocity_z(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.velocity(:,3),Core_Databins.time)';
    
    Core_Databins.tpara(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.tpara,Core_Databins.time)';
    Core_Databins.tperp(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.tperp,Core_Databins.time)';
    
    Core_Databins.dynamicpressure(:,Event_Number==Event_Number_vector) = ...
        interp1(Core.time*mean_core_timescale/(length(Core.time)),Core.dynamicpressure,Core_Databins.time)';
    
    
    %Trailing
    %Calculate dynamic pressure rho*v^2
    Trailing.dynamicpressure = calculate_dynamic_pressure(Trailing.density,Trailing.velocity,'i');
    %Interpolate this one event's data into the mean timescale and put into databins
    Trailing_Databins.density(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.density,Trailing_Databins.time)';
    
    Trailing_Databins.magneticfield(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.magneticfield,Trailing_Databins.time)';
    
    Trailing_Databins.velocity_x(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.velocity(:,1),Trailing_Databins.time)';
    Trailing_Databins.velocity_y(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.velocity(:,2),Trailing_Databins.time)';
    Trailing_Databins.velocity_z(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.velocity(:,3),Trailing_Databins.time)';
    
    Trailing_Databins.tpara(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.tpara,Trailing_Databins.time)';
    Trailing_Databins.tperp(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.tperp,Trailing_Databins.time)';
    
    Trailing_Databins.dynamicpressure(:,Event_Number==Event_Number_vector) = ...
        interp1(Trailing.time*mean_trailing_timescale/(length(Trailing.time)),Trailing.dynamicpressure,Trailing_Databins.time)';
    
    
end


%Construct the structure for all events, as a databin
Databins = struct(...
    'time',all_timescale_databins,...
    'density',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'magneticfield',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'velocity_x',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'velocity_y',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'velocity_z',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'tpara',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'tperp',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)),...
    'dynamicpressure',zeros(length(all_timescale_databins)+2,length(Event_Number_vector)));


% Databins.density = [Leading_Databins.density;Core_Databins.density;Trailing_Databins.density];
% Databins.magneticfield = [Leading_Databins.magneticfield;Core_Databins.magneticfield;Trailing_Databins.magneticfield];

fn = fieldnames(Databins); %Store the name of each of the fluxes
for k=2:numel(fn) %loop through the fluxes and perform the arithmetic
    Databins.(fn{k}) = [Leading_Databins.(fn{k}) ; Core_Databins.(fn{k}) ; Trailing_Databins.(fn{k})];
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

plot_superimposedEpoch(Databins.time,Databins.dynamicpressure,{'Dynamic Pressure';'[nPa]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.tpara,{'Tpara';'[eV]'},Type,numberEvents,numPlots,plotOrder)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plotOrder-1), plot_pos(3), plot_pos(4)]);
plotOrder=plotOrder+1;

plot_superimposedEpoch(Databins.time,Databins.tperp,{'Tperp';'[eV]'},Type,numberEvents,numPlots,plotOrder)
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
    
    [~,fgm_bdata] = interpxyz(fgm_timedata,fgm_bdata(:,4),fpi_i_timedata);  %magnitude only
    
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