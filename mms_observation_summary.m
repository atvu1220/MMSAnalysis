%Revised MMS Observation Summary
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear
close all
figure('Position',[1 1 650 900])
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/CompleteSummaryPlots/'
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_gap=1.25;
%cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
%figure('PaperPositionMode', 'auto')
probe_num = '1';
num_plots = 8;
data_type = 'brst';


Event_number = 29
if Event_number == 0
    trange = {'2017-01-14 06:03:45.000','2017-01-14 06:05:20.000'};
    event_start = cell2mat(trange(1));
    event_end = cell2mat(trange(2));
else
    [~, ~,~,~,...
        event_start,event_end,...
        ~, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~,...
        ~, ~] = get_eventTimes(Event_number);
    eventDataFileName = strcat('MMS1_Data_EventNumber_',num2str(Event_number),'.mat');
    %Fucking icloud keeps deleting files every 30 minutes of downloading them
    %eventDataDirectory = '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/Event Data/';
    eventDataDirectory = '/Users/andrewvu/data/Event Data/';
    
    eventDataDirectoryFileName = strcat(eventDataDirectory,eventDataFileName);
    load(eventDataDirectoryFileName)
end
formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(event_start,formatIn);
tend = datenum(event_end,formatIn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%Position Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load mec data
% [mec_timedata,mec_r_gsedata] = load_mec(event_start,probe_num,'srvy');

% %Find the start and end limits of the event in the data
% start_index_r = find(mec_timedata > tstart, 1);
% end_index_r = find(mec_timedata > tend, 1);
%
% %Convert the datetime to date String, and then crop to our event timeframe
% mec_timedata = mec_timedata(start_index_r:end_index_r,1);
% mec_r_gsedata = mec_r_gsedata(start_index_r:end_index_r,:)/6371.2;

%     [~,mec_r_gsedata,~,~] = crop(mec_timedata,mec_r_gsedata,event_start,event_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Load MEC Data
% [mec_timedata_raw, mec_rdata_raw] = load_mec(event_start,1,'srvy');
% [~,mec_r_gsedata,~,~] = crop(mec_timedata_raw,mec_rdata_raw,event_start,event_end);
% 
% 
% %load FGM data
% [fgm_timedata, fgm_bdata, ~, ~] = load_fgm(event_start,event_end,1,'brst');
% 
% %[fgm_timedata_srvy, fgm_bdata_srvy] = load_fgm(event_start,event_end,1,'srvy'); %For Sliding Window
% 
% %Load FPI_e
% [fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
%     fpi_e_edata,fpi_e_espectdata,fpi_e_pressdata] = load_fpi(event_start,event_end,1,'brst','e');
% %Load FPI_i
% [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
%     fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(event_start,event_end,1,'brst','i');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_order = 1;
%Plot Magnetic Fields
plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,num_plots,plot_order)
title('MMS1 Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;


%Plot FPI number density
plot_fpi_number(event_start,event_end,fpi_i_timedata,fpi_i_ndata,num_plots,plot_order,'i')
hold on
plot_fpi_number(event_start,event_end,fpi_e_timedata,fpi_e_ndata,num_plots,plot_order,'e')
hold off
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

%     annotation('textbox',get(gca,'Position'),...
%         'String',{strcat('\sigma=',num2str(coreDensity_STD));strcat('cv=',num2str(coreDensity_CV))},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 10);

%Plot FPI velocity data
plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;


%Plot FPI dynamic pressure
plot_total_pressure(event_start,event_end,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_pressdata,num_plots,plot_order,'i')
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

[shock_normal,r_sc] = calculate_bowshocknormal(event_start,mms1_mec_rdata_raw,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata);

% plot_entropy(event_start,event_end,fpi_i_timedata,fpi_i_ndata,fpi_i_tparadata,fpi_i_tperpdata,'i',num_plots,plot_order); hold on
% %plot_entropy(date_start,date_end,fpi_i_timedata,fpi_e_ndata,fpi_e_tparadata,fpi_e_tperpdata,'e',num_plots,plot_order); hold off
% set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% plot_order = plot_order+1;

%Plot shock angles
% plot_shockangle(event_start,event_end,fgm_timedata,fgm_bdata,num_plots,plot_order,shock_normal)
% set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% plot_order = plot_order+1;

%Plot -VxB
% plot_vcrossb(event_start,event_end,fgm_timedata,fgm_bdata,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order,'n_cs')
% set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% plot_order = plot_order+1;
%leading edge negative edot n_td means towards the discontinuity and
%positive edot n_td means towards the discontinuity

%Plot Temperatures
plot_fpi_temp(event_start,event_end,fpi_i_timedata,fpi_i_tparadata,fpi_i_tperpdata,num_plots,plot_order,'i')
hold on
plot_fpi_temp(event_start,event_end,fpi_e_timedata,fpi_e_tparadata,fpi_e_tperpdata,num_plots,plot_order,'e')
hold off
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

%Energy Spectrum for Ions
plot_fpi_energyspect(event_start,event_end,fpi_i_timedata,fpi_i_edata,fpi_i_espectdata,num_plots,plot_order)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

%Energy Spectrum for Electrons
plot_fpi_energyspect(event_start,event_end,fpi_e_timedata,fpi_e_edata,fpi_e_espectdata,num_plots,plot_order)

set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
datetick('keeplimits')
plot_order = plot_order+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


annotation('textbox',[plot_pos(1)-plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    'String',{datestr(tstart,'YYYY mmm dd'),'X-GSE (Re):', 'Y-GSE (Re):', 'Z-GSE (Re):'},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);
%
r_cols = length(mms1_mec_rdata_raw);
%n_cs = tdnormal(event_start,event_end,fgm_timedata,fgm_bdata,'event'); %(-X) GSE for current sheet normal
annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    'String',{'',num2str(round(mms1_mec_rdata_raw(1,:)/6371.2*10)/10,'%2g\n')},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);

%     annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
%         'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']'),strcat('n_{bs}: [',num2str(shock_normal,'%.4f '),']'),...
%         strcat('Shear Angle (Deg):',num2str(shear_angle,'%2.1f ')),strcat('Duration (s):',num2str(duration,'%2.2f '))},...
%         'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);


%
% annotation('textbox',[plot_pos(1)+2*(plot_pos(3)/r_cols)+0.025, plot_pos(2)-plot_pos(4)*7.5, plot_pos(3), plot_pos(4)],...
%     'String',{'',num2str(round(mec_r_gsedata(1,:)*10)/10,'%2g\n')},...
%     'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);


if Event_number == 0
    plot_name =  strcat('mms',probe_num,'_ObservationSummary_',event_start(1:19),'_',event_end(1:19),'.png');
else
    plot_name =  strcat('mms',probe_num,'_ObservationSummary_','Event_',num2str(Event_number),'.png');
end
%plot_name =  strcat(num2str(Event_number),'_Event_Number');
%print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
print(gcf,'-dpng','-r900', '-loose', plot_name);
%movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
