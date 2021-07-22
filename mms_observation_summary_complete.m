%Matlab CDF Plotter for MMS
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear
close all
figure('Position',[1 1 650 850])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_gap=1.25;
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
%figure('PaperPositionMode', 'auto');
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 9;
data_type = 'brst';
% date = '2018/03/18';




event_start = '2018-03-01 01:05:40.000';

event_start = '2018-03-01 01:03:44.000';
event_end = '2018-03-01 01:06:22.000';

event_start = '2018-03-01 01:03:45.000';
event_end = '2018-03-01 01:04:20.000';
% 
% date_start = '2018-03-01 01:04:30.000';
% date_end = '2018-03-01 01:05:05.000';

event_start = '2018-03-01 01:05:40.000';
event_end = '2018-03-01 01:06:22.000';

event_start = '2018-01-09 08:34:28.000';
event_end = '2018-01-09 08:34:59.000';

% date_start ='2018-03-01 01:18:03.000';
% date_end = '2018-03-01 01:20:33.000';

% date_start = '2018-04-27 19:47:14.000';
% date_end = '2018-04-27 19:48:22.000';



% date_start = '2015-12-28 05:26:50.000';
% date_end = '2015-12-28 05:27:40.000';

event_start = '2018-03-01 01:03:45.000';
event_end = '2018-03-01 01:04:20.000';

date_start = '2018-03-01 01:03:54.600';
date_end = '2018-03-01 01:03:55.100';

event_start    = '2017-1-15 14:24:50.000';
event_end      = '2017-12-15 14:26:10.000';

event_start    = '2017-11-19 11:31:24.000';
event_end      = '2017-11-19 11:31:35.000';

event_start    = '2018-01-12 01:50:21.000';
event_end      = '2018-01-12 01:52:21.000';


event_start    = '2017-12-21 08:08:10.000';
event_end      = '2017-12-21 08:09:25.000';

event_start = strcat('2017-11-04 05:15:45.0');
event_end = strcat('2017-11-04 05:17:00.0');

event_start = strcat('2017-11-04 05:15:45.0');
event_end = strcat('2017-11-04 05:19:30.0');

formatIn='yyyy-mm-dd HH:MM:SS.FFF';

tstart = datenum(event_start,formatIn);
tend = datenum(event_end,formatIn);



%%%%%%%%%%%%%%%%%%%%%%%%%%%Position Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load mec data
%[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,probe_num,'srvy');
% [mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,'srvy');
% [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
% [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
% [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');

%[mms1_mec_timedata_raw,mms1_mec_rdata_raw,~,~] = crop(mms1_mec_timedata_raw,mms1_mec_rdata_raw,event_start,event_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load fgm data
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,data_type);
% [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,data_type);
% [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,data_type);
% [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,data_type);

%Mec data has fewer data points than FGM Epheremis
[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,'srvy');
% [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
% [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
% [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');

%Load FPI_e
[fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
    fpi_e_edata,fpi_e_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'e');
%Load FPI_i
[fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
    fpi_i_edata,fpi_i_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'i');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_order = 1;
%Plot Magnetic Fields
plot_fgm_magnetic(event_start,event_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1:3),num_plots,plot_order)
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

%Plot FPI velocity data
plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

% % % % % %Plot FPI dynamic pressure
% % % % % plot_fpi_dynamic_pressure(event_start,event_end,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,num_plots,plot_order,'i')
% % % % % hold on
% % % % % plot_fpi_dynamic_pressure(event_start,event_end,fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,num_plots,plot_order,'e')
% % % % % hold off
% % % % % set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% % % % % plot_order = plot_order+1;
% % % % % 
% % % % % %Plot -VxB
% % % % % plot_vcrossb(event_start,event_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,fpi_i_timedata,fpi_i_vdata,num_plots,plot_order,'n_cs')
% % % % % set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% % % % % plot_order = plot_order+1;
% % % % % %leading edge negative edot n_td means towards the discontinuity and
% % % % % %positive edot n_td means towards the discontinuity

shock_normal = calculate_bowshocknormal(event_start,mms1_mec_rdata_raw,fpi_i_timedata,fpi_i_ndata,fpi_i_vdata)

%Plot shock angles
plot_shockangle(event_start,event_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,num_plots,plot_order,shock_normal)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;


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
n_cs = tdnormal(event_start,event_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,'event');
annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    'String',{'',num2str(round(mms1_mec_rdata_raw(1,:)/6371.2*10)/10,'%2g\n')},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);

annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05*4 plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)],...
    'String',{strcat('n_{cs}: [',num2str(n_cs,'%.4f '),']'),strcat('n_{bs}: [',num2str(shock_normal,'%.4f '),']')},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);

%
% annotation('textbox',[plot_pos(1)+2*(plot_pos(3)/r_cols)+0.025, plot_pos(2)-plot_pos(4)*7.5, plot_pos(3), plot_pos(4)],...
%     'String',{'',num2str(round(mec_r_gsedata(1,:)*10)/10,'%2g\n')},...
%     'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);


plot_name =  strcat('mms',probe_num,'_ObservationSummary_',event_start(1:19),'_',event_end(1:19),'.pdf');
print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
%movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
% 
% 
% plot_complete_boundary_analysis(event_start,event_end,date_start,date_end,...
% mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
% mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
% mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
% mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
% mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
% mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
% mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
% mms4_mec_timedata_raw,mms4_mec_rdata_raw)
% 
