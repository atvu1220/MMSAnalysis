%Matlab CDF Plotter for MMS
%Andrew Vu 9/19/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear
close all

plot_gap=1;
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Evolution'
%figure('PaperPositionMode', 'auto');
mms_directory = '/Users/andrewvu/data/mms/';

data_type = 'brst';
% date = '2018/03/18';



trange = {'2019-02-16 05:40:45.000','2019-02-16 05:42:00.000'};
event_start = cell2mat(trange(1));
event_end = cell2mat(trange(2));



formatIn='yyyy-mm-dd HH:MM:SS.FFF';

tstart = datenum(event_start,formatIn);
tend = datenum(event_end,formatIn);

probe_order = [3,4,1,2];


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
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'srvy');
[mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'srvy');
[mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'srvy');
[mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'srvy');

%Mec data has fewer data points than FGM Epheremis
[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,'srvy');
[mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
[mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
[mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');

%Load FPI_e
% [mms1_fpi_e_timedata,mms1_fpi_e_ndata,mms1_fpi_e_vdata,mms1_fpi_e_tparadata,mms1_fpi_e_tperpdata,...
%     mms1_fpi_e_edata,mms1_fpi_e_espectdata] = load_fpi(event_start,event_end,1,data_type,'e');
%Load FPI_i
[mms1_fpi_i_timedata,mms1_fpi_i_ndata,~,~,~,...
    ~,~] = load_fpi(event_start,event_end,1,data_type,'e');
[mms2_fpi_i_timedata,mms2_fpi_i_ndata,~,~,~,...
    ~,~] = load_fpi(event_start,event_end,2,data_type,'e');
[mms3_fpi_i_timedata,mms3_fpi_i_ndata,~,~,~,...
    ~,~] = load_fpi(event_start,event_end,3,data_type,'e');
[mms4_fpi_i_timedata,mms4_fpi_i_ndata,~,~,~,...
    ~,~] = load_fpi(event_start,event_end,4,data_type,'i');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolate B to FPI timescale
[~,mms1_fgm_bdata_interp] = interpxyz(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,mms1_fpi_i_timedata);
[~,mms2_fgm_bdata_interp] = interpxyz(mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,mms2_fpi_i_timedata);
[~,mms3_fgm_bdata_interp] = interpxyz(mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,mms3_fpi_i_timedata);
[~,mms4_fgm_bdata_interp] = interpxyz(mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,mms4_fpi_i_timedata);
%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[1 1 650 850])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_order = 1;
num_plots = 5;
%Plot Magnetic Fields
plot_fgm_magnetic(event_start,event_end,mms2_fpi_i_timedata,mms2_fgm_bdata_interp(:,1:3),num_plots,plot_order)
ylabel({'MMS2';'B';'[nT]'},'fontsize',14)
title('MMS FGM Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');
bLimits = ylim;
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fgm_magnetic(event_start,event_end,mms1_fpi_i_timedata,mms1_fgm_bdata_interp(:,1:3),num_plots,plot_order)
ylabel({'MMS1';'B';'[nT]'},'fontsize',14)
ylim(bLimits)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fgm_magnetic(event_start,event_end,mms4_fpi_i_timedata,mms4_fgm_bdata_interp(:,1:3),num_plots,plot_order)
ylabel({'MMS4';'B';'[nT]'},'fontsize',14)
ylim(bLimits)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fgm_magnetic(event_start,event_end,mms3_fpi_i_timedata,mms3_fgm_bdata_interp(:,1:3),num_plots,plot_order)
ylabel({'MMS3';'B';'[nT]'},'fontsize',14)
ylim(bLimits)
timeLimits = xlim;
datetick
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
xlim([timeLimits])
plot_order = plot_order+1;


plot_name =  strcat('mms','_fourSpacecraftFGM_',event_start(1:19),'_',event_end(1:19),'.pdf');
print(gcf,'-dpng','-r300', '-loose', plot_name);







% %Plot FPI number density
% plot_fpi_number(event_start,event_end,mms1_fpi_i_timedata,log(mms1_fpi_i_ndata),num_plots,plot_order,'i')
% hold on
% plot_fpi_number(event_start,event_end,mms2_fpi_i_timedata,log(mms2_fpi_i_ndata),num_plots,plot_order,'i')
% plot_fpi_number(event_start,event_end,mms3_fpi_i_timedata,log(mms3_fpi_i_ndata),num_plots,plot_order,'i')
% plot_fpi_number(event_start,event_end,mms4_fpi_i_timedata,log(mms4_fpi_i_ndata),num_plots,plot_order,'i')
% hold off
% set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
% plot_order = plot_order+1;



figure('Position',[1 1 650 850])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_order = 1;
num_plots = 5;

plot_fpi_number(event_start,event_end,mms2_fpi_i_timedata,(mms2_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'MMS2';'n_i';'[1/cc]'},'fontsize',14)
title('MMS FPI Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');
nLimits = ylim;
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fpi_number(event_start,event_end,mms1_fpi_i_timedata,(mms1_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'MMS1';'n_i';'[1/cc]'},'fontsize',14)
ylim(nLimits)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fpi_number(event_start,event_end,mms4_fpi_i_timedata,(mms4_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'MMS4';'n_i';'[1/cc]'},'fontsize',14)
ylim(nLimits)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;

plot_fpi_number(event_start,event_end,mms3_fpi_i_timedata,(mms3_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'MMS3';'n_i';'[1/cc]'},'fontsize',14)
ylim(nLimits)
timeLimits = xlim;
datetick
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
xlim([timeLimits])
plot_order = plot_order+1;


plot_name =  strcat('mms','_fourSpacecraftFPI_',event_start(1:19),'_',event_end(1:19),'.pdf');
print(gcf,'-dpng','-r300', '-loose', plot_name);








figure('Position',[1 1 850 850])
% co = [0 0 1;
%     1 0 0;
%     0 0 0];
% set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_order = 1;
num_plots = 5;

% left_color = [.5 .5 0];
% right_color = [0 .5 .5];


%Plot Magnetic Fields
% colororder({'b','m'})
yyaxis left 
plot_fgm_magnetic(event_start,event_end,mms2_fpi_i_timedata,mms2_fgm_bdata_interp(:,4),num_plots,plot_order)
bLimits = ylim;
bLimits = [0 20];
ylabel({'MMS2';'|B|';'[nT]'},'fontsize',14)
timeLimits = xlim;
ylim(bLimits)

yyaxis right
plot_fpi_number(event_start,event_end,mms2_fpi_i_timedata,(mms2_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'n';'[1/cc]'},'fontsize',14)
nLimits = ylim;
nLimits = [0 15];
ylim(nLimits);
xlim(timeLimits)
legend off

title('MMS FGM-FPI Observatory Summary', 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
ax=gca;
ax.YAxis(1).Color = [0    0.4470    0.7410];
ax.YAxis(2).Color = [0.8500    0.3250    0.0980];

plot_order = plot_order+1;


yyaxis left
plot_fgm_magnetic(event_start,event_end,mms1_fpi_i_timedata,mms1_fgm_bdata_interp(:,4),num_plots,plot_order)
ylabel({'MMS1';'|B|';'[nT]'},'fontsize',14)
xlim(timeLimits)
ylim(bLimits)

yyaxis right
plot_fpi_number(event_start,event_end,mms1_fpi_i_timedata,(mms1_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'n';'[1/cc]'},'fontsize',14)
xlim(timeLimits)
ylim(nLimits)
legend off

ax=gca;
ax.YAxis(1).Color = [0    0.4470    0.7410];
ax.YAxis(2).Color = [0.8500    0.3250    0.0980];
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;


yyaxis left
plot_fgm_magnetic(event_start,event_end,mms4_fpi_i_timedata,mms4_fgm_bdata_interp(:,4),num_plots,plot_order)
ylabel({'MMS4';'|B|';'[nT]'},'fontsize',14)
xlim(timeLimits)
ylim(bLimits)

yyaxis right
% plot_fpi_number(event_start,event_end,mms4_fpi_i_timedata,median(mms4_fpi_i_ndata)+filter_fpi(mms4_fpi_i_ndata,'i',25,3),num_plots,plot_order,'i')
plot_fpi_number(event_start,event_end,mms4_fpi_i_timedata,mms4_fpi_i_ndata,num_plots,plot_order,'i')
ylabel({'n';'[1/cc]'},'fontsize',14)
xlim(timeLimits)
ylim(nLimits)
legend off

ax=gca;
ax.YAxis(1).Color = [0    0.4470    0.7410];
ax.YAxis(2).Color = [0.8500    0.3250    0.0980];
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+1;



yyaxis left
plot_fgm_magnetic(event_start,event_end,mms3_fpi_i_timedata,mms3_fgm_bdata_interp(:,4),num_plots,plot_order)
ylabel({'MMS3';'|B|';'[nT]'},'fontsize',14)
xlim(timeLimits)
ylim(bLimits)

yyaxis right
plot_fpi_number(event_start,event_end,mms3_fpi_i_timedata,(mms3_fpi_i_ndata),num_plots,plot_order,'i')
ylabel({'n';'[1/cc]'},'fontsize',14)
xlim(timeLimits)
ylim(nLimits)
legend off


datetick
ax=gca;
ax.YAxis(1).Color = [0    0.4470    0.7410];
ax.YAxis(2).Color = [0.8500    0.3250    0.0980];
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);
xlim(timeLimits)
yyaxis left
xlim(timeLimits)
plot_order = plot_order+1;
set(gcf,'position',[1 1 900 850])

plot_name =  strcat('mms','_fourSpacecraftFGM-FPI_',event_start(1:19),'_',event_end(1:19),'.pdf');
print(gcf,'-dpng','-r300', '-loose', plot_name);

