%Plots all 4 spacecrafts FGM data and FPI velocties and densities
clear; close all
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];


%Event Parameters
% date_start =  '2018-01-29 03:36:05.000';
% date_end = '2018-01-29 03:36:25.000';

% date_start =  '2017-12-29 19:11:17.000'
% date_end = '2017-12-29 19:11:40.000'
Event_number = 6;

[Event_Type, Substructure,threshold_std,~,...
    date_start,date_end,...
    ~, left_InnerEdge,...
    right_InnerEdge, ~,...
    leading_leftmost_date, leading_rightmost_date,...
    trailing_leftmost_date, trailing_rightmost_date,...
    ~, ~,...
    ~, ~] = get_eventTimes(Event_number);



%MVAB Analysis Time Range

%[MVAB_start,MVAB_end] = getTimeRange(cursor_info.Position(1),cursor_info.Position(2),'datestr')

MVAB_start = '2018-01-29 03:36:14.792';
MVAB_end = '2018-01-29 03:36:16.792';

MVAB_start = '2018-01-29 03:36:14.318';
MVAB_end = '2018-01-29 03:36:16.731';

%
% MVAB_start ='2018-01-29 03:36:15.257';
% MVAB_end ='2018-01-29 03:36:16.418';




cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Multiplots'
eventFolderName = strcat(date_start(1:10),'/',date_start(12:13),'-',date_start(15:16));
eventFolderName = strcat('Event_Number_',num2str(Event_number));
if ~exist(eventFolderName,'dir')
    mkdir(eventFolderName)
end
cd(eventFolderName)


LMN_coordinates = 0;
%% Load Data from each spacecraft
%%%%MMS1
%load FGM data
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(date_start,date_end,1,'brst');
%Load FPI_i
[mms1_fpi_i_timedata_raw,mms1_fpi_i_ndata_raw,mms1_fpi_i_vdata_raw,~,~,~,~] = load_fpi(date_start,date_end,1,'brst','i');

%%%%MMS2
%load FGM data
[mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(date_start,date_end,2,'brst');
%Load FPI_i
[mms2_fpi_i_timedata_raw,mms2_fpi_i_ndata_raw,mms2_fpi_i_vdata_raw,~,~,~,~] = load_fpi(date_start,date_end,2,'brst','i');

%%%%MMS3
%load FGM data
[mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(date_start,date_end,3,'brst');
%Load FPI_i
[mms3_fpi_i_timedata_raw,mms3_fpi_i_ndata_raw,mms3_fpi_i_vdata_raw,~,~,~,~] = load_fpi(date_start,date_end,3,'brst','i');


%%%%MMS4
%load FGM data
[mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(date_start,date_end,4,'brst');
%Load FPI_i
[mms4_fpi_i_timedata_raw,mms4_fpi_i_ndata_raw,mms4_fpi_i_vdata_raw,~,~,~,~] = load_fpi(date_start,date_end,4,'brst','i');

%% Transform into LMN Coordinates
if LMN_coordinates == 1
    %Calculate LMN Coordinates with MVAB on MMS1
    [MVAB_fgm_timedata,MVAB_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,MVAB_start,MVAB_end);
    [~,eig,MVAB_transform] = mvab(MVAB_fgm_bdata(:,1:3))
    
    %Transformation Time
    mms1_fgm_bdata_mag = mms1_fgm_bdata_raw(:,4);
    mms1_fgm_bdata_raw = mms1_fgm_bdata_raw(:,1:3)*MVAB_transform;
    mms1_fpi_i_vdata_raw = mms1_fpi_i_vdata_raw*MVAB_transform;
    
    mms2_fgm_bdata_mag = mms2_fgm_bdata_raw(:,4);
    mms2_fgm_bdata_raw = mms2_fgm_bdata_raw(:,1:3)*MVAB_transform;
    mms2_fpi_i_vdata_raw = mms2_fpi_i_vdata_raw*MVAB_transform;
    
    mms3_fgm_bdata_mag = mms3_fgm_bdata_raw(:,4);
    mms3_fgm_bdata_raw = mms3_fgm_bdata_raw(:,1:3)*MVAB_transform;
    mms3_fpi_i_vdata_raw = mms3_fpi_i_vdata_raw*MVAB_transform;
    
    mms4_fgm_bdata_mag = mms4_fgm_bdata_raw(:,4);
    mms4_fgm_bdata_raw = mms4_fgm_bdata_raw(:,1:3)*MVAB_transform;
    mms4_fpi_i_vdata_raw = mms4_fpi_i_vdata_raw*MVAB_transform;
    
    Blabel = {'B_l';'B_m';'B_n';'B'};
    Vlabel = {'V_l';'V_m';'V_n'};
else
    mms1_fgm_bdata_mag = mms1_fgm_bdata_raw(:,4);
    mms2_fgm_bdata_mag = mms2_fgm_bdata_raw(:,4);
    mms3_fgm_bdata_mag = mms3_fgm_bdata_raw(:,4);
    mms4_fgm_bdata_mag = mms4_fgm_bdata_raw(:,4);
    
    Blabel =  {'B_x';'B_y';'B_z';'B'};
    Vlabel =  {'V_x';'V_y';'V_z'};
end
%% Plotting B and V separately, but all SC on same plot

%Magnetic Field
plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1),4,1); datetick('keeplimits'); hold on
plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,1),4,1); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,1),4,1); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,1),4,1); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
title('MMS FGM Magnetic Field', 'FontSize',14)
ylabel({string(Blabel(1));'[nT]'})

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,2),4,2); datetick('keeplimits'); hold on
plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,2),4,2); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,2),4,2); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,2),4,2); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Blabel(2));'[nT]'})

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,3),4,3); datetick('keeplimits'); hold on
plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,3),4,3); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,3),4,3); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,3),4,3); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Blabel(3));'[nT]'})

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_mag,4,4); datetick('keeplimits'); hold on
plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_mag,4,4); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_mag,4,4); datetick('keeplimits')
plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_mag,4,4); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Blabel(4));'[nT]'})


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','_B_all_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)



%Ion Density
figure
plot_fpi_number(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_ndata_raw,1,1,'i'); datetick('keeplimits'); hold on
plot_fpi_number(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_ndata_raw,1,1,'i'); datetick('keeplimits')
plot_fpi_number(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_ndata_raw,1,1,'i'); datetick('keeplimits')
plot_fpi_number(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_ndata_raw,1,1,'i'); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
title('MMS FPI Number Density', 'FontSize',14)


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','_n_all_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)



%Ion Velocity
figure
plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits');hold on
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Vlabel(1));'[km/s]'})
title('MMS FPI Ion Bulk Flow', 'FontSize',14)

plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits'); hold on
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Vlabel(2));'[km/s]'})

plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits');hold on
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits')
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits')
legend({'MMS1', 'MMS2', 'MMS3','MMS4'},'FontSize',10)
ylabel({string(Vlabel(3));'[km/s]'})


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','_V_all_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS1 B and Bmag Plot
figure
%Find the lowest and highest points in Bn, then plot a line
%Drawing lines
% line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')

[~,minIndex] = min(mms1_fgm_bdata_raw(:,3));
[~,maxIndex] = max(mms1_fgm_bdata_raw(:,3));


plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_mag,4,1); datetick('keeplimits');
line([(mms1_fgm_timedata_raw(minIndex)) (mms1_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms1_fgm_timedata_raw(maxIndex)) (mms1_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(4));'[nT]'})
title('MMS1: Magnetic Field','FontSize',14)
legend off

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1),4,2); datetick('keeplimits');
line([(mms1_fgm_timedata_raw(minIndex)) (mms1_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms1_fgm_timedata_raw(maxIndex)) (mms1_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(1));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,2),4,3); datetick('keeplimits');
line([(mms1_fgm_timedata_raw(minIndex)) (mms1_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms1_fgm_timedata_raw(maxIndex)) (mms1_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(2));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,3),4,4); datetick('keeplimits');
line([(mms1_fgm_timedata_raw(minIndex)) (mms1_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms1_fgm_timedata_raw(maxIndex)) (mms1_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(3));'[nT]'})
legend off



set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','1','_B_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS2 B and Bmag Plot
figure
%Find the lowest and highest points in Bn, then plot a line
%Drawing lines
% line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')

[~,minIndex] = min(mms2_fgm_bdata_raw(:,3));
[~,maxIndex] = max(mms2_fgm_bdata_raw(:,3));


plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_mag,4,1); datetick('keeplimits');
line([(mms2_fgm_timedata_raw(minIndex)) (mms2_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms2_fgm_timedata_raw(maxIndex)) (mms2_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(4));'[nT]'})
title('MMS2: Magnetic Field','FontSize',14)
legend off

plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,1),4,2); datetick('keeplimits');
line([(mms2_fgm_timedata_raw(minIndex)) (mms2_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms2_fgm_timedata_raw(maxIndex)) (mms2_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(1));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,2),4,3); datetick('keeplimits');
line([(mms2_fgm_timedata_raw(minIndex)) (mms2_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms2_fgm_timedata_raw(maxIndex)) (mms2_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(2));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,3),4,4); datetick('keeplimits');
line([(mms2_fgm_timedata_raw(minIndex)) (mms2_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms2_fgm_timedata_raw(maxIndex)) (mms2_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(3));'[nT]'})
legend off



set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','2','_B_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS3 B and Bmag Plot
figure
%Find the lowest and highest points in Bn, then plot a line
%Drawing lines
% line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')

[~,minIndex] = min(mms3_fgm_bdata_raw(:,3));
[~,maxIndex] = max(mms3_fgm_bdata_raw(:,3));


plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_mag,4,1); datetick('keeplimits');
line([(mms3_fgm_timedata_raw(minIndex)) (mms3_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms3_fgm_timedata_raw(maxIndex)) (mms3_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(4));'[nT]'})
title('MMS3: Magnetic Field','FontSize',14)
legend off

plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,1),4,2); datetick('keeplimits');
line([(mms3_fgm_timedata_raw(minIndex)) (mms3_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms3_fgm_timedata_raw(maxIndex)) (mms3_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(1));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,2),4,3); datetick('keeplimits');
line([(mms3_fgm_timedata_raw(minIndex)) (mms3_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms3_fgm_timedata_raw(maxIndex)) (mms3_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(2));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,3),4,4); datetick('keeplimits');
line([(mms3_fgm_timedata_raw(minIndex)) (mms3_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms3_fgm_timedata_raw(maxIndex)) (mms3_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(3));'[nT]'})
legend off



set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','3','_B_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS4 B and Bmag Plot
figure
%Find the lowest and highest points in Bn, then plot a line
%Drawing lines
% line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')

[~,minIndex] = min(mms4_fgm_bdata_raw(:,3));
[~,maxIndex] = max(mms4_fgm_bdata_raw(:,3));


plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_mag,4,1); datetick('keeplimits');
line([(mms4_fgm_timedata_raw(minIndex)) (mms4_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms4_fgm_timedata_raw(maxIndex)) (mms4_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(4));'[nT]'})
title('MMS4: Magnetic Field','FontSize',14)
legend off

plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,1),4,2); datetick('keeplimits');
line([(mms4_fgm_timedata_raw(minIndex)) (mms4_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms4_fgm_timedata_raw(maxIndex)) (mms4_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(1));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,2),4,3); datetick('keeplimits');
line([(mms4_fgm_timedata_raw(minIndex)) (mms4_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms4_fgm_timedata_raw(maxIndex)) (mms4_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(2));'[nT]'})
legend off

plot_fgm_magnetic(date_start,date_end,mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,3),4,4); datetick('keeplimits');
line([(mms4_fgm_timedata_raw(minIndex)) (mms4_fgm_timedata_raw(minIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms4_fgm_timedata_raw(maxIndex)) (mms4_fgm_timedata_raw(maxIndex))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
ylabel({string(Blabel(3));'[nT]'})
legend off



set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','4','_B_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS1 Interpolation for V&B plot

%Interpolate data with fpi instrument, FGM down to FPI
mms1_fgm_bdata = zeros(length(mms1_fpi_i_timedata_raw),3);
for i=1:3
    mms1_fgm_bdata(:,i)=interp1(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,i),mms1_fpi_i_timedata_raw,'pchip');
end

figure

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fgm_bdata(:,1),3,1); datetick('keeplimits');
ylabel(string(Blabel(1)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits');
ylabel(string(Vlabel(1)))
title('MMS1: Magnetic Field and Velocity','FontSize',14)
legend off

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fgm_bdata(:,2),3,2); datetick('keeplimits');
ylabel(string(Blabel(2)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits');
ylabel(string(Vlabel(2)))
legend off

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fgm_bdata(:,3),3,3); datetick('keeplimits');
ylabel(string(Blabel(3)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms1_fpi_i_timedata_raw,mms1_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits');
ylabel(string(Vlabel(3)))
legend off


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','1','_BV_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS2 Interpolation for V&B plot
%Interpolate data with fpi instrument, FGM down to FPI
mms2_fgm_bdata = zeros(length(mms2_fpi_i_timedata_raw),3);

for i=1:3
    mms2_fgm_bdata(:,i)=interp1(mms2_fgm_timedata_raw,mms2_fgm_bdata_raw(:,i),mms2_fpi_i_timedata_raw,'pchip');
end

figure

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fgm_bdata(:,1),3,1); datetick('keeplimits');
ylabel(string(Blabel(1)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits');
ylabel(string(Vlabel(1)))
title('MMS2: Magnetic Field and Velocity','FontSize',14)
legend off


yyaxis left
plot_fgm_magnetic(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fgm_bdata(:,2),3,2); datetick('keeplimits');
ylabel(string(Blabel(2)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits');
ylabel(string(Vlabel(2)))
legend off

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fgm_bdata(:,3),3,3); datetick('keeplimits');
ylabel(string(Blabel(3)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms2_fpi_i_timedata_raw,mms2_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits');
ylabel(string(Vlabel(3)))
legend off


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','2','_BV_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS3 Interpolation for V&B plot

%Interpolate data with fpi instrument, FGM down to FPI
mms3_fgm_bdata = zeros(length(mms3_fpi_i_timedata_raw),3);

for i=1:3
    mms3_fgm_bdata(:,i)=interp1(mms3_fgm_timedata_raw,mms3_fgm_bdata_raw(:,i),mms3_fpi_i_timedata_raw,'pchip');
end

figure

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fgm_bdata(:,1),3,1); datetick('keeplimits');
ylabel(string(Blabel(1)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits');
ylabel(string(Vlabel(1)))
title('MMS3: Magnetic Field and Velocity','FontSize',14)
legend off


yyaxis left
plot_fgm_magnetic(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fgm_bdata(:,2),3,2); datetick('keeplimits');
ylabel(string(Blabel(2)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits');
ylabel(string(Vlabel(2)))
legend off

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fgm_bdata(:,3),3,3); datetick('keeplimits');
ylabel(string(Blabel(3)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms3_fpi_i_timedata_raw,mms3_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits');
ylabel(string(Vlabel(3)))
legend off


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','3','_BV_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)

%% MMS4 Interpolation for V&B plot

%Interpolate data with fpi instrument, FGM down to FPI
mms4_fgm_bdata = zeros(length(mms4_fpi_i_timedata_raw),3);

for i=1:3
    mms4_fgm_bdata(:,i)=interp1(mms4_fgm_timedata_raw,mms4_fgm_bdata_raw(:,i),mms4_fpi_i_timedata_raw,'pchip');
end

figure

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fgm_bdata(:,1),3,1); datetick('keeplimits');
ylabel(string(Blabel(1)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,1),3,1); datetick('keeplimits');
ylabel(string(Vlabel(1)))
title('MMS4: Magnetic Field and Velocity','FontSize',14)
legend off


yyaxis left
plot_fgm_magnetic(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fgm_bdata(:,2),3,2); datetick('keeplimits');
ylabel(string(Blabel(2)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,2),3,2); datetick('keeplimits');
ylabel(string(Vlabel(2)))
legend off

yyaxis left
plot_fgm_magnetic(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fgm_bdata(:,3),3,3); datetick('keeplimits');
ylabel(string(Blabel(3)))
yyaxis right
plot_fpi_bulkv(date_start,date_end,mms4_fpi_i_timedata_raw,mms4_fpi_i_vdata_raw(:,3),3,3); datetick('keeplimits');
ylabel(string(Vlabel(3)))
legend off


set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w')
plot_name =  strcat('mms','4','_BV_',date_start(1:19),'_',date_end(1:19),'.eps');
print(gcf,'-depsc2', '-loose', plot_name)
